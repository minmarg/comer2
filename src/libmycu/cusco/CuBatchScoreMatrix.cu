/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/psl.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/TCTXVECT.h"
#include "libpro/srcsco/ScoresAttr.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/SerializedScores.cuh"
#include "libmycu/cupro/SerializedScoresCtor.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresCtor.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "CuBatchScoreMatrix_Init.cuh"
#include "CuBatchScoreMatrix_ssss.cuh"
#include "CuBatchScoreMatrix_cvs2s.cuh"
#include "CuBatchScoreMatrix_hdp1s.cuh"
#include "CuBatchScoreMatrix_Init_hdp1s.cuh"
#include "CuBatchScoreMatrix_ssss_cvs2s.cuh"
#include "CuBatchScoreMatrix_ssss_cvs2s_hdp1s.cuh"
#include "CuBatchScoreMatrix.cuh"

// -------------------------------------------------------------------------
// constructor
//
CuBatchScoreMatrix::CuBatchScoreMatrix(
        Configuration config[],
        int precscale )
:
    AbstractScoreMatrix<CUBSM_TYPE>( 
        AbstractScoreMatrix<CUBSM_TYPE>::CuBatchLSO, 
        config, 
        AbstractScoreMatrix<CUBSM_TYPE>::StatisticsGiven, 
        precscale )
{
    MYMSG( "CuBatchScoreMatrix::CuBatchScoreMatrix", 5 );
}

// -------------------------------------------------------------------------
// Init: data initialization
//
void CuBatchScoreMatrix::Init( int /*querylen*/, int /*sbjctlen*/ )
{
    MYMSG( "CuBatchScoreMatrix::Init", 5 );
    AbstractScoreMatrix<CUBSM_TYPE>::Init( 0, 0 );
}

// -------------------------------------------------------------------------
// default constructor
//
CuBatchScoreMatrix::CuBatchScoreMatrix()
:
    AbstractScoreMatrix<CUBSM_TYPE>()
{
    throw MYRUNTIME_ERROR("CuBatchScoreMatrix::CuBatchScoreMatrix: "
                "Default initialization is prohibited.");
}

// -------------------------------------------------------------------------
// destructor
//
CuBatchScoreMatrix::~CuBatchScoreMatrix()
{
    MYMSG( "CuBatchScoreMatrix::~CuBatchScoreMatrix", 5 );
}





// -------------------------------------------------------------------------
// ComputeScoreMatrix: compute score matrix
//
void CuBatchScoreMatrix::ComputeScoreMatrix( bool /*final*/ )
{
    MYMSG( "CuBatchScoreMatrix::ComputeScoreMatrix", 4 );
}

// -------------------------------------------------------------------------
// ComputeScoreMatrixDevice: compute score matrices for part of query and 
// database profiles on device;
// streamproc, CUDA stream for computations;
// sssscores, serialized pairwise secondary structure state scores;
// ssssattr, attributes of serialized SS scores;
// cvs2scores, serialized context vector scores;
// cvs2sattr, attributes of serialized context vector scores;
// hdp1sTexo, texture object encapsulating serialized HDP1 scores;
// hdpsattr, attributes of HDP1 scores;
// hdp1scoresinuse, flag indicating wether HDP1 scores are in use;
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
// outscores, the address of the device memory region dedicated to scores;
// outmodscores, the address of the device memory region dedicated to 
// modular scores;
//
void CuBatchScoreMatrix::ComputeScoreMatrixDevice(
    cudaStream_t streamproc,
    CUBSM_TYPE* sssscores,
    SerializedScoresAttr ssssattr,
    CUBSM_TYPE* cvs2scores,
    SerializedCVS2ScoresAttr cvs2sattr,
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr hdpsattr,
    bool hdp1scoresinuse,
    size_t nqyposs,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t querposoffset,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* outmodscores )
{
    MYMSG( "CuBatchScoreMatrix::ComputeScoreMatrixDevice", 4 );
    const mystring preamb = "CuBatchScoreMatrix::ComputeScoreMatrixDevice: ";
    myruntime_error mre;

    const float ssswgt = MOptions::GetSSSWGT();
    const float cvswgt = MOptions::GetCVSWGT();
    float hdp1swgt = MOptions::GetADJWGT(); 
    if( ssswgt > 0.0f )
        hdp1swgt *= ADJWGTDG;

    size_t ndbxposs = ndb1poss + ndbCposs;
    size_t szsspace = nqyposs * ndbxposs;//search space

#ifdef __DEBUG__
    if( szsspace < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of search space.");
#endif

    //threads per block; should be multiple of 64 to allow TS to 
    // avoid register memory bank conflicts
//     const unsigned int nthsperblock = 256;
//     //blocks per grid
//     const unsigned int nbkspergrid = (szsspace+nthsperblock-1) / nthsperblock;
//     dim3 nblcks(nbkspergrid,1,1);
//     dim3 nthrds(nthsperblock,1,1);

    //configuration for CalcSMInitSMEM
//     dim3 nthrds(SMINIT_2DCACHE_DIM,SMINIT_2DCACHE_DIM,1);
//     dim3 nblcks(((ndbxposs)+nthrds.x-1)/nthrds.x,(nqyposs+nthrds.y-1)/nthrds.y,1);

    //configuration for CalcSMInitSMEMUnroll2
//     dim3 nthrds(SMINIT_2DCACHE_DIM,SMINIT_2DCACHE_DIM,1);
//     dim3 nblcks(((ndbxposs)+nthrds.x-1)/nthrds.x,(nqyposs+2*nthrds.y-1)/(2*nthrds.y),1);

    //configuration for CalcSMInitSMEMUnroll2x
    dim3 nthrds(SMINIT_2DCACHE_DIM,SMINIT_2DCACHE_DIM,1);
    dim3 nblcks(((unsigned int)(ndbxposs)+2*nthrds.x-1)/(2*nthrds.x),
				((unsigned int)nqyposs+nthrds.y-1)/nthrds.y,1);

    try {
        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sKernelComputeScoreMatrix execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, 
                "grid size= (%u,%u,%u) block size= (%u,%u,%u);%sspace size= %zu; # query poss= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, NL, szsspace, nqyposs );
            strbuf += msgbuf;
            if( ndb1poss ) {
                sprintf(msgbuf, "; # db cached poss= %zu", ndb1poss );
                strbuf += msgbuf;
            }
            if( ndbCposs ) {
                sprintf(msgbuf, "; # db new poss= %zu", ndbCposs );
                strbuf += msgbuf;
            }
            sprintf(msgbuf, " (padding, %zu)", dbxpad );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

//         //{{TEST kernels:
//         nblcks = dim3((nblcks.x+31)/32,(nblcks.y+31)/32,1);
//         nblcks = dim3(100/*4900*/,1,1);
//         for( uint i = 0; i < nqyposs; i+=32 ) {
//             for( uint j = 0; j < 469; j++ ) {
//                 KernelStreamedIterativeTest<<<nblcks,nthrds,0,streamproc>>>( 
//                     i,
//                     (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
//                     querposoffset, bdb1posoffset, bdbCposoffset,
//                     outscores );
//                 MYCUDACHECKLAST;
//             }
//         }
// //         CalcSMInit<<<nblcks,nthrds>>>( 
// //             (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, 
// //             querposoffset, bdb1posoffset, bdbCposoffset );
// //         CalcSMInitSMEM<<<nblcks,nthrds>>>(
// //             (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, 
// //             querposoffset, bdb1posoffset, bdbCposoffset );
//         //}}

        if( hdp1scoresinuse && hdp1swgt > 0.0f ) {
            //{{Initial plus HDP scores
            CalcSM_Init_HDP1S_SMEMUnroll2x<<<nblcks,nthrds,0,streamproc>>>(
                hdp1sTexo,
                hdpsattr,
                hdp1swgt,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
				(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                //NOTE: these pointers must be aligned
                outscores,
                NULL/*outmodscores*/
            );
            MYCUDACHECKLAST;
            //}}
        } else {
            //{{Initial scores
            CalcSM_Init_SMEMUnroll2x<<<nblcks,nthrds,0,streamproc>>>(
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
				(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                //NOTE: these pointers must be aligned
                outscores,
                NULL/*outmodscores*/
            );
            MYCUDACHECKLAST;
            //}}
//             //{{Scoring HDP cluster membership match
//             CalcSM_HDP1S_SMEMUnroll2x<<<nblcks,nthrds,0,streamproc>>>(
//                 hdp1sTexo,
//                 hdpsattr,
//                 MOptions::GetADJWGT() * ADJWGTDG,
//                 (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
//                 (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
//                 //NOTE: these pointers must be aligned
//                 outscores,
//                 NULL/*outmodscores*/
//             );
//             MYCUDACHECKLAST;
//             //}}
        }

//         if( GetHDP1ScoresInUse() && hdp1swgt > 0.0f && ssswgt > 0.0f && cvswgt > 0.0f ) {
//             //{{Scoring secondary structure, context vector, and HDP cluster membership match
//             CalcSM_SSSS_CVS2S_HDP1S_SMEMUnroll2x
//             <<<nblcks,nthrds,ssssattr.szalloc_+cvs2sattr.szalloc_,streamproc>>>(
//                 sssscores,
//                 ssssattr,
//                 ssswgt,
//                 cvs2scores,
//                 cvs2sattr,
//                 cvswgt,
//                 hdp1sTexo,
//                 hdpsattr,
//                 hdp1swgt,
//                 (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
//                 (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
//                 //NOTE: these pointers must be aligned
//                 outscores,
//                 outmodscores
//             );
//             MYCUDACHECKLAST;
//             //}}
//         } else

        if( ssswgt > 0.0f && cvswgt > 0.0f ) {
            //{{Scoring secondary structure and context vector match
            CalcSM_SSSS_CVS2S_SMEMUnroll2x
            <<<nblcks,nthrds,ssssattr.szalloc_+cvs2sattr.szalloc_,streamproc>>>(
                sssscores,
                ssssattr,
                ssswgt,
                cvs2scores,
                cvs2sattr,
                cvswgt,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
				(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                //NOTE: these pointers must be aligned
                outscores,
                outmodscores
            );
            MYCUDACHECKLAST;
            //}}
        } else if( ssswgt > 0.0f ) {
            //{{Scoring secondary structure match
            CalcSM_SSSS_SMEMUnroll2x<<<nblcks,nthrds,ssssattr.szalloc_,streamproc>>>(
                sssscores,
                ssssattr,
                ssswgt,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
				(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                //NOTE: these pointers must be aligned
                outscores,
                outmodscores
            );
            MYCUDACHECKLAST;
            //}}
        } else if( cvswgt > 0.0f ) {
            //{{Scoring context vector match
            CalcSM_CVS2S_SMEMUnroll2x<<<nblcks,nthrds,cvs2sattr.szalloc_,streamproc>>>(
                cvs2scores,
                cvs2sattr,
                cvswgt,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
				(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                //NOTE: these pointers must be aligned
                outscores,
                outmodscores
            );
            MYCUDACHECKLAST;
            //}}
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}





// TEST METHODS ------------------------------------------------------------
// -------------------------------------------------------------------------
// TESTPrintProProScores1: Print calculated profile-profile scores for 
// testing;
//
void CuBatchScoreMatrix::TESTPrintProProScores1(
    char** cachedbdb1pmbeg, 
    char** cachedbdb1pmend, 
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend,
    char** bdbCpmbeg,
    char** bdbCpmend,
    unsigned int dbxpad,
    size_t szscores,
    CUBSM_TYPE* outscores )
{
    MYMSG( "CuBatchScoreMatrix::TESTPrintProProScores1", 3 );
    const mystring preamb = "CuBatchScoreMatrix::TESTPrintProProScores1: ";
    myruntime_error mre;

    //NOTE: SERIAL NUMBER OF THE PROFILE TO PRINT
    int prosernr = 1286;//3347;//2;//2097;//120;//4899;//2097;

    //const size_t nflds = PMBatchProData::GetNoFields();

    size_t ndbcachedposs = PMBatchProData::GetNoPositsFromTo( cachedbdb1pmbeg, cachedbdb1pmend );
    size_t ndbcachedpros = PMBatchProData::GetNoProsFromTo( cachedbdb1pmbeg, cachedbdb1pmend );
    size_t ndb1missedpros = ndbcachedpros;
    size_t ndb1missedposs = ndbcachedposs;
    size_t ndb1pros = 0UL;
    size_t ndbCpros = 0UL;

    const size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmbeg, querypmend );
    size_t ndb1poss = 0UL;
    size_t ndbCposs = 0UL;

    if( bdb1pmbeg && bdb1pmend ) {
        ndb1missedpros = PMBatchProData::GetNoProsFromTo( cachedbdb1pmbeg, bdb1pmbeg );
        ndb1missedposs = PMBatchProData::GetNoPositsFromTo( cachedbdb1pmbeg, bdb1pmbeg );
        ndb1pros = PMBatchProData::GetNoProsFromTo( bdb1pmbeg, bdb1pmend );
        ndb1poss = PMBatchProData::GetNoPositsFromTo( bdb1pmbeg, bdb1pmend );
    }
    if( bdbCpmbeg && bdbCpmend ) {
        ndbCpros = PMBatchProData::GetNoProsFromTo( bdbCpmbeg, bdbCpmend );
        ndbCposs = PMBatchProData::GetNoPositsFromTo( bdbCpmbeg, bdbCpmend );
    }

    if( szscores < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of scores.");

    char* tmpbuf = NULL;

    fprintf( stderr,"%s szscores = %zu\n", preamb.c_str(), szscores );

    tmpbuf = (char*)malloc(szscores * sizeof(char));
    if( !tmpbuf )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

    //char*  pm1dat[nflds];
	//char* pm1dat[pmv2DTotFlds];

    size_t col1ndx = 0;//column index of the identified db profile in the score matrix
    size_t pro1len = 0;//profile length
    int i = -1, j;

    try {

        if( ndb1missedpros < (size_t)prosernr + 1 ) {
            prosernr -= (int)ndb1missedpros;
            if( ndb1pros < (size_t)prosernr + 1 ) {
                prosernr -= (int)ndb1pros;
                if(!( ndbCpros < (size_t)prosernr + 1 ) && bdbCpmbeg ) {
                    col1ndx = ((LNTYPE*)(bdbCpmbeg[pps2DDist]))[i=prosernr] + ndb1poss;
                    pro1len = (size_t)((INTYPE*)(bdbCpmbeg[pps2DLen]))[prosernr];
                }
            }
            else if( bdb1pmbeg ) {
                col1ndx = ((LNTYPE*)(bdb1pmbeg[pps2DDist]))[i=prosernr] - ndb1missedposs;
                pro1len = (size_t)((INTYPE*)(bdb1pmbeg[pps2DLen]))[prosernr];
            }
            fprintf(stderr,"ndb1missedpros= %zu ndb1pros= %zu ndbCpros %zu "
                    "prosernr= %d  col1ndx= %zu pro1len= %zu\n",
                    ndb1missedpros,ndb1pros,ndbCpros,prosernr,col1ndx,pro1len);
        }

//         if( bdb1pmbeg && bdb1pmend ) {
//             memcpy( pm1dat, bdb1pmbeg, nflds * sizeof(void*));
//             for( i = ndb1missedpros; i < prosernr && pm1dat[pps2DLen] < bdb1pmend[pps2DLen]; i++ ) {
//                 pro1len = (size_t)*(INTYPE*)pm1dat[pps2DLen];
//                 col1ndx += pro1len;
//                 PMBatchProData::PMDataNextPro( pm1dat );
//             }
//         }
//         if(( bdbCpmbeg && bdbCpmend ) && 
//            (!bdb1pmend || bdb1pmend[pps2DLen] <= pm1dat[pps2DLen] )) {
//             memcpy( pm1dat, bdbCpmbeg, nflds * sizeof(void*));
//             for( ; i < prosernr && pm1dat[pps2DLen] < bdbCpmend[pps2DLen]; i++ ) {
//                 pro1len = (size_t)*(INTYPE*)pm1dat[pps2DLen];
//                 col1ndx += pro1len;
//                 PMBatchProData::PMDataNextPro( pm1dat );
//             }
//         }

        do {
            if( prosernr != i ) {
                warning((preamb + "Profile with a serial number not found.").c_str());
                break;
            }

//             pro1len = (size_t)*(INTYPE*)pm1dat[pps2DLen];

            MYCUDACHECK( cudaMemcpy( tmpbuf, outscores, szscores, cudaMemcpyDeviceToHost ));
            MYCUDACHECKLAST;

            fprintf( stderr, "## Scores (QxS= %zux%zu):%s##%s", nqyposs, pro1len, NL, NL );
            for( i = 0; i < (int)nqyposs; i++ ) {
                for( j = 0; j < (int)pro1len; j++ ) {
                    fprintf( stderr, " %10.6f", 
                        ((CUBSM_TYPE*)tmpbuf)[i * (ndb1poss + ndbCposs + dbxpad) + col1ndx + j]);
                }
                fprintf( stderr, "%s", NL );
                //for(size_t _k=0;_k<100000UL;_k++)clock();
            }
            fprintf( stderr, "%s", NL );

        } while(0);

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( tmpbuf ) {
        free(tmpbuf);
        tmpbuf = NULL;
    }

    if( mre.isset())
        throw mre;
}
