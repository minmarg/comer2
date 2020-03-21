/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/psl.h"
#include "libpro/srcpro/MOptions.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/CuDeviceMemory.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_final.cuh"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_com.h"
#include "CuBatchSS_score_probs.cuh"
#include "CuBatchSS_dynplm.cuh"
#include "CuBatchSS.cuh"

#define CUBATCHSS_SCOREPROBS_DYNPRLM

// -------------------------------------------------------------------------
// constructor
//
CuBatchSS::CuBatchSS()
{
    MYMSG( "CuBatchSS::CuBatchSS", 5 );
}

// // -------------------------------------------------------------------------
// // default constructor
// //
// CuBatchSS::CuBatchSS()
// {
//     throw MYRUNTIME_ERROR("CuBatchSS::CuBatchSS: "
//                 "Default initialization is not allowed.");
// }

// -------------------------------------------------------------------------
// destructor
//
CuBatchSS::~CuBatchSS()
{
    MYMSG( "CuBatchSS::~CuBatchSS", 5 );
}





// -------------------------------------------------------------------------
// CalculateAlnStatisticsDevice: calculate alignment statistics for multiple 
// profile-profile pairs on device
//
void CuBatchSS::CalculateAlnStatisticsDevice(
    cudaStream_t streamproc,
    bool Xuninf,
    size_t ndb1pros,
    size_t ndbCpros,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    //
    int ssemodel,
    float qyeno,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
    //
    CUBSM_TYPE* scores,//in:substitution scores
    CUBDP_TYPE* tmpdpbuffer,//in:temporary information of passed profiles
    float* tmpss2datbuffers,//[out:profile-specific]
    float* dp2alndatbuffers,//[out:profile-specific]
    //summary of profiles passed to phase 2
    unsigned int* attrpassed )
{
    MYMSG( "CuBatchSS::CalculateAlnStatisticsDevice", 4 );
    CalculateStatisticsDeviceHelper(
        streamproc,
        Xuninf,
        ndb1pros, ndbCpros,
        ndb1prosOmtd, ndbCprosOmtd,
        nqyposs, ndb1poss, ndbCposs, dbxpad,
        querposoffset, bdb1posoffset, bdbCposoffset,
        ssemodel,
        qyeno,
        searchspace,
        reflambda, refK,
        expgappedlambda, expgappedK,
        scores,//in:substitution scores
        tmpdpbuffer,//in:temporary information of passed profiles
        tmpss2datbuffers,//[out:temporary profile-specific probabilities]
        dp2alndatbuffers,//[out:profile-specific]
        attrpassed
    );
}

// -------------------------------------------------------------------------
// CalculateStatisticsDeviceHelper: calculate score probabilities, 
// statistical parameters, and alignment significance for multiple 
// significantly similar profile-profile pairs representing part of query 
// and database profiles on device;
// streamproc, CUDA stream for computations;
// Xuninf, flag indicating whether masked positions are uninformative;
// ndb1pros, number of profiles over the db1 positions passed;
// ndbCpros, number of profiles over the complementary dbC positions;
// ndb1prosOmtd, number of files missed up to the beginning of db1 profiles;
// ndbCprosOmtd, number of files missed up to the beginning of dbC profiles;
// nqyposs, number of query positions to process;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// dbxpad, padding along the x axis when using 2D thread blocks;
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
// ssemodel, number of the model for statistical significance estimation;
// qyeno, query ENO;
// searchspace, effective search space for the query under processing;
// reflambda, reference lambda value for ungapped alignments;
// refK, reference K value for ungapped alignments;
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
// scores, the address of the device memory region of calculated scores;
// tmpdpbuffer, device memory allocated for temporary bottom scores 
// which now temporarily contains the information of profiles passed to 
// phase 2;
// tmpss2datbuffers, device memory allocated for temporary probability arrays;
// dp2alndatbuffers, device memory for keeping profile-profile alignment
//  statistics;
// attrpassed, summary of profiles that passed to phase 2;
//
void CuBatchSS::CalculateStatisticsDeviceHelper(
    cudaStream_t streamproc,
    bool Xuninf,
    size_t ndb1pros,
    size_t /*ndbCpros*/,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    //
    int ssemodel,
    float qyeno,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
    //
    CUBSM_TYPE* scores,//in:substitution scores
    CUBDP_TYPE* tmpdpbuffer,//in:temporary information of passed profiles
    float* tmpss2datbuffers,//[out:profile-specific]
    float* dp2alndatbuffers,//[out:profile-specific]
    //summary of profiles passed to phase 2
    unsigned int* attrpassed )
{
    MYMSG( "CuBatchSS::CalculateStatisticsDeviceHelper", 4 );
    const mystring preamb = "CuBatchSS::CalculateStatisticsDeviceHelper: ";
    myruntime_error mre;

#ifdef __DEBUG__
    size_t ndbxposs = ndb1poss + ndbCposs;
    size_t szsspace = nqyposs * ndbxposs;//search space
    if( szsspace < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of search space.");
#endif

    size_t pass2ndbpros = attrpassed[CuDeviceMemory::dgvNPassedPros];
    size_t pass2ndbxposs = attrpassed[CuDeviceMemory::dgvNPosits];
    size_t pass2dbmaxprolen = attrpassed[CuDeviceMemory::dgvMaxProLen];

    if( pass2ndbpros < 1 || pass2ndbxposs < 1 || pass2dbmaxprolen < 1 )
        return;

    try {
#ifdef CUBATCHSS_SCOREPROBS_DYNPRLM
        //execution configuration
        dim3 nthrds(1,1,1);
        dim3 nblcks((unsigned int)pass2ndbpros,1,1);

        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf, "%sKernelCalcStatisticsDynPlm execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, "grid size= (%u,%u,%u) block size= (%u,%u,%u);%s"
                "# query poss= %zu; db poss= %zu db pros= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, 
                NL, nqyposs, pass2ndbxposs, pass2ndbpros );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        CalcStatisticsDynPlm<<<nblcks,nthrds,0,streamproc>>>( 
            Xuninf,
            (uint)ndb1pros,
            (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            ssemodel,
            qyeno,
            searchspace,
            reflambda, refK,
            expgappedlambda, expgappedK,
            scores,
            tmpdpbuffer, 
            tmpss2datbuffers,
            dp2alndatbuffers
        );
#else
        //execution configuration
        dim3 nthrds(CUSS_2DCACHE_DIM,CUSS_2DCACHE_DIM,1);
        dim3 nblcks((pass2dbmaxprolen+nthrds.x-1)/nthrds.x,
                    (nqyposs+nthrds.y-1)/nthrds.y,
                    pass2ndbpros);

        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb + "KernelCalcScoreProbs execution configuration: ";
            sprintf(msgbuf, "grid size= (%u,%u,%u) block size= (%u,%u,%u); "
                "# query poss= %lu; db poss= %lu db pros= %lu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, 
                nqyposs, pass2ndbxposs, pass2ndbpros );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        CalcScoreProbs<<<nblcks,nthrds,0,streamproc>>>( 
            (uint)ndb1pros,
            (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            scores,
            tmpdpbuffer, 
            tmpss2datbuffers
        );
#endif
        MYCUDACHECKLAST;

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}





// TEST METHODS ------------------------------------------------------------
// -------------------------------------------------------------------------
//
