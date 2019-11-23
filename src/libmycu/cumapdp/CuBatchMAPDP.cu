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
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP.cuh"
#include "libmycu/cudp/CuBatchDP_final.cuh"
#include "CuBatchMAPDP_com.h"
#include "CuBatchMAPDP_prnt.cuh"
#include "CuBatchMAPDP_fwd.cuh"
#include "CuBatchMAPDP_bwd.cuh"
#include "CuBatchMAPDP_map_btck.cuh"
#include "CuBatchMAPDP_map_final.cuh"
#include "CuBatchMAPDP.cuh"

#define CUBATCHMAPDP_PARENT_DYNPRLM

// -------------------------------------------------------------------------
// constructor
//
CuBatchMAPDP::CuBatchMAPDP()
{
    MYMSG( "CuBatchMAPDP::CuBatchMAPDP", 5 );
}

// // -------------------------------------------------------------------------
// // default constructor
// //
// CuBatchMAPDP::CuBatchMAPDP()
// {
//     throw MYRUNTIME_ERROR("CuBatchMAPDP::CuBatchMAPDP: "
//                 "Default initialization is not allowed.");
// }

// -------------------------------------------------------------------------
// destructor
//
CuBatchMAPDP::~CuBatchMAPDP()
{
    MYMSG( "CuBatchMAPDP::~CuBatchMAPDP", 5 );
}





// -------------------------------------------------------------------------
// PerformCompleteMAPDynProgDevice: perform complete MAP dynamic 
// programming on device
//
void CuBatchMAPDP::PerformCompleteMAPDynProgDevice(
    cudaStream_t streamproc,
    float logevthld,
    bool modscoresinuse,
    size_t ndb1pros,
    size_t ndbCpros,
    size_t querprosOmtd,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t dbpro1len,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    CUBSM_TYPE* scores,
    CUBSM_TYPE* modscores,
    float* tmpdpdiagbuffers, 
    float* tmpdpbotbuffer,
//float* tmpblockprobscales, 
//float* tmpprobscales,
    unsigned int* maxcoordsbuf,
    char* btckdata,
    float* tmpfwdprobs, 
//float* tmpbwdprobs,
    float* tmpss2datbuffers,//[tmp]
    float* dp2alndatbuffers,//[in]
    char* outalns,//[out]
    unsigned int* attrpassed,
    size_t dbxpad2,
    size_t cbdbalnlen2 )
{
    MYMSG( "CuBatchMAPDP::PerformCompleteMAPDynProgDevice", 4 );
    CalculateFwdBwdProbsMAPDynProgDevice(
        streamproc,
        logevthld,
        modscoresinuse,
        ndb1pros, ndbCpros,
        querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
        dbpro1len,
        nqyposs, ndb1poss, ndbCposs, dbxpad,
        querposoffset, bdb1posoffset, bdbCposoffset,
        scores, modscores, 
        tmpdpdiagbuffers, tmpdpbotbuffer,
//tmpblockprobscales, tmpprobscales,
        tmpfwdprobs,
//tmpbwdprobs,
        tmpss2datbuffers,
        dp2alndatbuffers,
        attrpassed,
        dbxpad2,
        maxcoordsbuf,
        btckdata
    );
    FinalizeMAPDynProgDevice(
        streamproc,
        logevthld,
        ndb1pros,
        querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
        nqyposs, ndb1poss, ndbCposs, dbxpad,
        querposoffset, bdb1posoffset, bdbCposoffset,
        scores,//[in]
        tmpdpdiagbuffers,//[in]
        dp2alndatbuffers,//[in/out]
        attrpassed,
        dbxpad2,
        cbdbalnlen2,
        maxcoordsbuf,//[in]
        btckdata,//[in]
        outalns//[out]
    );
}

// -------------------------------------------------------------------------
// CalculateFwdProbsMAPDynProgDevice: perform dynamic programming of forward 
// probability calculations for multiple profile-profile pairs representing 
// part of query and database profiles on device;
// streamproc, CUDA stream for computations;
// modscoresinuse, indication of using the modular scores;
// ndb1pros, number of profiles over the db1 positions passed;
// ndbCpros, number of profiles over the complementary dbC positions;
// querprosOmtd, number of query profiles up to this query profile;
// ndb1prosOmtd, number of files missed up to the beginning of db1 profiles;
// ndbCprosOmtd, number of files missed up to the beginning of dbC profiles;
// dbpro1len, the largest profile length over all profiles currently being 
// processed;
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
// scores, the address of the device memory region of calculated scores;
// modscores, the address of the device memory region of modular scores;
// tmpdpdiagbuffers, device memory allocated for temporary diagonal scores;
// tmpdpbotbuffer, device memory allocated for temporary bottom DP scores;
// tmpblockprobscales, allocation for temporary diagonal block-specific 
// probability scale factors;
// tmpprobscales, temporary final probability scale factors;
// tmpfwdprobs, device memory allocated for forward probabilities;
// tmpbwdprobs, device memory allocated for backward probabilities;
// tmpss2datbuffers, profile-specific (phase-2) device memory 
// allocated for SS estimation, here used as temporary buffer for total 
// probability for each pair of profiles;
// dp2alndatbuffers, input data of profile-profile alignment statistics;
// attrpassed, summary of profiles that passed to (this) phase 2;
// dbxpad2, padding for the phase-2 db profile data
//
void CuBatchMAPDP::CalculateFwdBwdProbsMAPDynProgDevice(
    cudaStream_t streamproc,
    float logevthld,
    bool modscoresinuse,
    size_t ndb1pros,
    size_t /*ndbCpros*/,
    size_t querprosOmtd,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t dbpro1len,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    CUBSM_TYPE* scores,
    CUBSM_TYPE* modscores,
    float* tmpdpdiagbuffers, 
    float* tmpdpbotbuffer,
//float* /*tmpblockprobscales*/, 
//float* /*tmpprobscales*/,
    float* tmpfwdprobs,
//float* tmpbwdprobs,
    float* tmpss2datbuffers,//[tmp]
    float* dp2alndatbuffers,//[in]
    unsigned int* attrpassed,
    size_t dbxpad2,
    unsigned int* maxcoordsbuf,
    char* btckdata )
{
    MYMSG( "CuBatchMAPDP::CalculateFwdBwdProbsMAPDynProgDevice", 4 );
    const mystring preamb = "CuBatchMAPDP::CalculateFwdBwdProbsMAPDynProgDevice: ";
    myruntime_error mre;

    const float alnextreg = MOptions::GetMINPP();
    const float alnextgapreg = alnextreg * 0.5f;

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

    //dimensions of a block diagonal
    const uint dimblkdiag = CUMAPDP_2DCACHE_DIM_D;

    //maximum number of elements in a diagonal of single elements
    //const uint maxdiagelems = SLC_MIN(nqyposs, dbpro1len);
    //maximum number of block diagonal elements;
    //NOTE: we consider oblique blocks in their diagonals as opposed to 
    // rectangular-block diagonals;
    //const uint maxblkdiagelems = 
    //    CuBatchMAPDP::GetMaxRegBlockDiagonalElems( maxdiagelems, nqyposs, dimblkdiag );
    //----
    //NOTE: use block diagonals, where blocks share a common point 
    // (corner, instead of edge) with a neighbour in a diagonal;
    const uint maxblkdiagelems = CuBatchDP::GetMaxBlockDiagonalElems( 
            dbpro1len, nqyposs, dimblkdiag, CUMAPDP_2DCACHE_DIM_X );

    //execution configuration
    dim3 nthrds(dimblkdiag,1,1);
    dim3 nblcks(maxblkdiagelems,(uint)pass2ndbpros,1);
    //number of regular DIAGONAL block diagonal series, each of given dimensions;
    //rect coords (x,y) are related to diagonal number d by d=x+y-1;
    //then, this number d is divided by the width of the block diagonal 
    // (which is the block length)
    //uint nblkdiags = ((dbpro1len+nqyposs-1)+CUMAPDP_2DCACHE_DIM_X-1)/CUMAPDP_2DCACHE_DIM_X;
    //REVISION: due to the positioning of the first block, the first 1-position diagonal of the 
    // first diagonal block is out of bounds: remove -1
    uint nblkdiags = (uint)(((dbpro1len+nqyposs)+CUMAPDP_2DCACHE_DIM_X-1)/CUMAPDP_2DCACHE_DIM_X);
    //----
    //NOTE: use block DIAGONALS, where blocks share a COMMON POINT 
    // (corner, instead of edge) with a neighbour in a diagonal;
    //the number of series of such block diagonals equals 
    // #regular block diagonals (above) + {(l-1)/w}, 
    // where l is query length (y coord.), w is block width (dimension), and
    // {} denotes floor rounding; {(l-1)/w} is #interior divisions;
    nblkdiags += (uint)((nqyposs-1)/dimblkdiag);

    //the number of the last diagonal starting at x=-dimblkdiag
    uint nintdivs = (uint)((nqyposs-1)/dimblkdiag);
    float nsepds = (float)CUMAPDP_2DCACHE_DIM_D/(float)CUMAPDP_2DCACHE_DIM_X + 1.0f;
    //each division separates a number of diagonals (nsepds); -1 for zero-based indices;
    uint lastdiag = (uint)(nsepds * (float)nintdivs) + 1 - 1;


    try {
        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf, "%sKernelFwdBwdProbsMAPDynProg execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, 
                "grid size= (%u,%u,%u) block size= (%u,%u,%u);%s# query poss= %zu; db poss= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, NL, nqyposs, pass2ndbxposs );
            strbuf += msgbuf;
            sprintf(msgbuf, " db pros= %zu", pass2ndbpros );
            strbuf += msgbuf;
            sprintf(msgbuf, "; # block diags= %u (dim, %ux%u)", 
                    nblkdiags, dimblkdiag, CUMAPDP_2DCACHE_DIM_X );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

#ifdef CUBATCHMAPDP_PARENT_DYNPRLM
        ExecMAPDP_Parent_DynPlm<<<1,1,0,streamproc>>>(
            alnextreg,
            alnextgapreg,
            //
            nblkdiags,
            nblcks, nthrds,
            //
            lastdiag,
            logevthld,
            modscoresinuse,
            (uint)ndb1pros,
            (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            (uint)(pass2ndbxposs + dbxpad2),
            (uint)dbxpad2,
            scores,
            modscores,
            dp2alndatbuffers,
            tmpdpdiagbuffers,
            tmpdpbotbuffer,
            tmpfwdprobs,
//tmpbwdprobs,
            tmpss2datbuffers,
            maxcoordsbuf,
            btckdata
        );
        MYCUDACHECKLAST;
#else
        //launch blocks along block diagonals so that d=x+y-1
        for( uint d = 0; d < nblkdiags; d++) {
            ExecMAPDP_Fwd_Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                d,
                lastdiag,
                logevthld,
                modscoresinuse,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                (uint)pass2ndbxposs + dbxpad2,
                (uint)dbxpad2,
                scores,//[in]
                modscores,//[in]
                dp2alndatbuffers,//[in]
                tmpdpdiagbuffers, 
                tmpdpbotbuffer,
                tmpfwdprobs
            );
            MYCUDACHECKLAST;
        }
        //launch blocks along block diagonals so that d=x+y-1
        for( int d = nblkdiags-1; 0 <= d; d--) {
            ExecMAPDP_Bwd_Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                d,
                lastdiag,
                logevthld,
                modscoresinuse,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                (uint)pass2ndbxposs + dbxpad2,
                (uint)dbxpad2,
                scores,//[in]
                modscores,//[in]
                dp2alndatbuffers,//[in]
                tmpdpdiagbuffers, 
                tmpdpbotbuffer,
                tmpfwdprobs,
//tmpbwdprobs,
                tmpss2datbuffers
            );
            MYCUDACHECKLAST;
        }
        //launch blocks along block diagonals so that d=x+y-1
        for( uint d = 0; d < nblkdiags; d++) {
            ExecMAPDP_MAP_Btck_Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                alnextreg,
                alnextgapreg,
                d,
                lastdiag,
                logevthld,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, //(uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                //(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                (uint)pass2ndbxposs + dbxpad2,
                (uint)dbxpad2,
                tmpfwdprobs,//[in]
                dp2alndatbuffers,//[in]
                tmpdpdiagbuffers, 
                tmpdpbotbuffer,
                tmpss2datbuffers,
                maxcoordsbuf,
                btckdata
            );
            MYCUDACHECKLAST;
        }
#endif

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// FinalizeMAPDynProgDevice: finalize MAP dynamic programming and draw 
// alignments in parallel on device;
// streamproc, CUDA stream for computations;
// ...
//
void CuBatchMAPDP::FinalizeMAPDynProgDevice(
    cudaStream_t streamproc,
    float logevthld,
    size_t ndb1pros,
    size_t querprosOmtd,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    CUBSM_TYPE* scores,//[in]
    float* tmpdpdiagbuffers,//[in]
    float* dp2alndatbuffers,//[in/out]
    unsigned int* attrpassed,
    size_t dbxpad2,
    size_t cbdbalnlen2,
    unsigned int* maxcoordsbuf,//[in]
    char* btckdata,//[in]
    char* outalns )//[out]
{
    MYMSG( "CuBatchMAPDP::FinalizeMAPDynProgDevice", 4 );
    const mystring preamb = "CuBatchMAPDP::FinalizeMAPDynProgDevice: ";
    myruntime_error mre;

    //TODO: add the option of showing SS states in the output
    const bool printsss = MOptions::GetSSSWGT() > 0.0f;

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

    //set db length for alignments, including padding
    size_t dbalnlen2 = pass2ndbxposs + pass2ndbpros * (nqyposs+1);
    dbalnlen2 = ALIGN_UP( dbalnlen2, CUL2CLINESIZE );

    if(cbdbalnlen2 != dbalnlen2)
        throw MYRUNTIME_ERROR( preamb +
            "Inconsistent total alignment length in phase 2.");

    //execution configuration
    dim3 nthrds(CUMAPDP_2DCACHE_DIM_D,1,1);
    dim3 nblcks((uint)pass2ndbpros,1,1);


    try {
        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf, 
                "%sKernelFinalizeMAPDynProg execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, 
                "grid size= (%u,%u,%u) block size= (%u,%u,%u);%s# db pros= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, NL, pass2ndbpros );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        FinalizeMAPDP_MAP<<<nblcks,nthrds,0,streamproc>>>(
            printsss,
            logevthld,
            (uint)ndb1pros,
            (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            (uint)(pass2ndbxposs + dbxpad2),//dblen2
            (uint)dbxpad2,
            (uint)dbalnlen2,
            scores,
            tmpdpdiagbuffers, 
            maxcoordsbuf,
            btckdata,
            dp2alndatbuffers,
            outalns
        );
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
