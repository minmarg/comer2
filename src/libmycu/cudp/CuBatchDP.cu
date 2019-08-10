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
#include "liblib/msg.h"
#include "libpro/srcpro/MOptions.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"
#include "CuBatchDP_init.cuh"
#include "CuBatchDP_init_block2.cuh"
#include "CuBatchDP_init_corr.cuh"
#include "CuBatchDP_init_btck.cuh"
#include "CuBatchDP_final.cuh"
#include "CuBatchDP.cuh"

// #define TEST_CUBATCHDP_INIT_BLOCK2

// -------------------------------------------------------------------------
// constructor
//
CuBatchDP::CuBatchDP()
{
    MYMSG( "CuBatchDP::CuBatchDP", 5 );
}

// // -------------------------------------------------------------------------
// // default constructor
// //
// CuBatchDP::CuBatchDP()
// {
//     throw MYRUNTIME_ERROR("CuBatchDP::CuBatchDP: "
//                 "Default initialization is not allowed.");
// }

// -------------------------------------------------------------------------
// destructor
//
CuBatchDP::~CuBatchDP()
{
    MYMSG( "CuBatchDP::~CuBatchDP", 5 );
}





// -------------------------------------------------------------------------
// PerformDynProgDevice: perform complete dynamic programming on device
//
void CuBatchDP::PerformCompleteDynProgDevice(
    cudaStream_t streamproc,
    CUBSM_TYPE scorethld,
    bool calcbcktrcmatrix,
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
    CUBDP_TYPE* tmpdpdiagbuffers, 
    CUBDP_TYPE* tmpdpbotbuffer,
    unsigned int* maxcoordsbuf,
    char* btckdata,
    unsigned int* attrpassed )
{
    MYMSG( "CuBatchDP::PerformCompleteDynProgDevice", 4 );
    PerformDynProgDevice(
        streamproc,
        calcbcktrcmatrix,
        ndb1pros, ndbCpros,
        querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
        dbpro1len,
        nqyposs, ndb1poss, ndbCposs, dbxpad,
        querposoffset, bdb1posoffset, bdbCposoffset,
        scores,
        tmpdpdiagbuffers, 
        tmpdpbotbuffer,
        maxcoordsbuf,
        btckdata
    );
    FinalizeDynProgDevice(
        streamproc,
        scorethld,
        ndb1pros, ndbCpros,
        ndb1prosOmtd, ndbCprosOmtd,
        dbpro1len,
        nqyposs, ndb1poss, ndbCposs, dbxpad,
        querposoffset, bdb1posoffset, bdbCposoffset,
        scores,
        tmpdpdiagbuffers, 
        tmpdpbotbuffer,
        maxcoordsbuf,
        btckdata,
        //[out:]
        attrpassed
    );
}

// -------------------------------------------------------------------------
// PerformDynProgDevice: perform dynamic programming for multiple 
// profile-profile pairs representing part of query and database profiles on 
// device;
// streamproc, CUDA stream for computations;
// calcbcktrcmatrix, indication of calculating back-tracing matrix;
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
// tmpdpdiagbuffers, device memory allocated for temporary diagonal scores;
// tmpdpbotbuffer, device memory allocated for temporary bottom DP scores;
// maxscoordsbuf, device memory allocated for the coordinates of maximum 
// alignment scores;
// btckdata, device memory allocated for backtracking information;
//
void CuBatchDP::PerformDynProgDevice(
    cudaStream_t streamproc,
    bool /*calcbcktrcmatrix*/,
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
    CUBDP_TYPE* tmpdpdiagbuffers, 
    CUBDP_TYPE* tmpdpbotbuffer,
    unsigned int* maxcoordsbuf,
    char* btckdata )
{
    MYMSG( "CuBatchDP::PerformDynProgDevice", 4 );
    const mystring preamb = "CuBatchDP::PerformDynProgDevice: ";
    myruntime_error mre;

    size_t ndbxposs = ndb1poss + ndbCposs;

#ifdef __DEBUG__
    size_t szsspace = nqyposs * ndbxposs;//search space
    if( szsspace < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of search space.");
#endif

    //dimensions of a block diagonal
    const uint dimblkdiag = CUDP_2DCACHE_DIM_D;

    //maximum number of elements in a diagonal of single elements
    //const uint maxdiagelems = SLC_MIN(nqyposs, dbpro1len);
    //maximum number of block diagonal elements;
    //NOTE: we consider oblique blocks in their diagonals as opposed to 
    // rectangular-block diagonals;
    //const uint maxblkdiagelems = 
    //    CuBatchDP::GetMaxRegBlockDiagonalElems( maxdiagelems, nqyposs, dimblkdiag );
    //----
    //NOTE: use block diagonals, where blocks share a common point 
    // (corner, instead of edge) with a neighbour in a diagonal;
    const uint maxblkdiagelems = CuBatchDP::GetMaxBlockDiagonalElems( 
            dbpro1len/*ndbxposs*/, nqyposs, dimblkdiag, CUDP_2DCACHE_DIM_X );

    //execution configuration
    dim3 nthrds(dimblkdiag,1,1);
    dim3 nblcks(maxblkdiagelems,(uint)(ndb1pros+ndbCpros),1);
    //number of regular DIAGONAL block diagonal series, each of given dimensions;
    //rect coords (x,y) are related to diagonal number d by d=x+y-1;
    //then, this number d is divided by the width of the block diagonal 
    // (which is the block length)
    //uint nblkdiags = ((dbpro1len+nqyposs-1)+CUDP_2DCACHE_DIM_X-1)/CUDP_2DCACHE_DIM_X;
    //REVISION: due to the positioning of the first block, the first 1-position diagonal of the 
    // first diagonal block is out of bounds: remove -1
    uint nblkdiags = (uint)(((dbpro1len+nqyposs)+CUDP_2DCACHE_DIM_X-1)/CUDP_2DCACHE_DIM_X);
    //----
    //NOTE: use block DIAGONALS, where blocks share a COMMON POINT 
    // (corner, instead of edge) with a neighbour in a diagonal;
    //the number of series of such block diagonals equals 
    // #regular block diagonals (above) + {(l-1)/w}, 
    // where l is query length (y coord.), w is block width (dimension), and
    // {} denotes floor rounding; {(l-1)/w} is #interior divisions;
    nblkdiags += (uint)(nqyposs-1)/dimblkdiag;

    //the number of the last diagonal starting at x=-dimblkdiag
    uint nintdivs = (uint)(nqyposs-1)/dimblkdiag;
    float nsepds = (float)CUDP_2DCACHE_DIM_D/(float)CUDP_2DCACHE_DIM_X + 1.0f;
    //each division separates a number of diagonals (nsepds); -1 for zero-based indices;
    uint lastdiag = (uint)(nsepds * (float)nintdivs) + 1 - 1;


#ifdef TEST_CUBATCHDP_INIT_BLOCK2
    //RECTANGULAR block diagonals
    //max number of rectangular blocks in one diagonal
    *const_cast<uint*>(&maxblkdiagelems) = CuBatchDP::GetMaxBlockDiagonalElemsBlock2( 
            dbpro1len, nqyposs, dimblkdiag, CUDP_2DCACHE_DIM_X );
    //execution configuration
    nthrds = dim3(dimblkdiag,1,1);
    nblcks = dim3(maxblkdiagelems,ndb1pros+ndbCpros,1);
    //# rectangular block diagonals
    nblkdiags = (dbpro1len+CUDP_2DCACHE_DIM_X-1)/CUDP_2DCACHE_DIM_X;
    nblkdiags += (nqyposs+CUDP_2DCACHE_DIM_D-1)/CUDP_2DCACHE_DIM_D - 1;
    //-1 for zero-based indices
    lastdiag = (nqyposs+CUDP_2DCACHE_DIM_D-1)/CUDP_2DCACHE_DIM_D - 1;
#endif


    try {
        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sKernelPerformDynProg execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, 
                "grid size= (%u,%u,%u) block size= (%u,%u,%u);%s# query poss= %zu; db poss= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, NL, nqyposs, ndbxposs );
            strbuf += msgbuf;
            sprintf(msgbuf, " (padding, %zu) db pros= %zu", dbxpad, ndb1pros+ndbCpros );
            strbuf += msgbuf;
            sprintf(msgbuf, "; # block diags= %u (dim, %ux%u)", 
                    nblkdiags, dimblkdiag, CUDP_2DCACHE_DIM_X );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        //launch blocks along block diagonals so that d=x+y-1
        for( uint d = 0; d < nblkdiags; d++) {
#ifdef TEST_CUBATCHDP_INIT_BLOCK2
            ExecDPBlock2Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                d,
                lastdiag,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                scores,
                tmpdpdiagbuffers, 
                tmpdpbotbuffer
            );
#elif defined(CUDP_CALC_CORRSCORES_INLINE)
            ExecDP_Corr_Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                d,
                lastdiag,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                scores,
                tmpdpdiagbuffers, 
                tmpdpbotbuffer
            );
#else
//             ExecDPUnroll32x<<<nblcks,nthrds,0,streamproc>>>( 
//                 d,
//                 lastdiag,
//                 (uint)ndb1pros,
//                 (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
//                 (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
//                 (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
//                 scores,
//                 tmpdpdiagbuffers, 
//                 tmpdpbotbuffer
//             );
            ExecDP_Btck_Unroll32x<<<nblcks,nthrds,0,streamproc>>>( 
                d,
                lastdiag,
                (uint)ndb1pros,
                (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
                (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
                (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
                scores,
                tmpdpdiagbuffers, 
                tmpdpbotbuffer,
                maxcoordsbuf,
                btckdata
            );
#endif
            MYCUDACHECKLAST;
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// FinalizeDynProgDevice: finalize dynamic programming by calculating 
// correlation match scores and reducing maximum alignment scores for 
// multiple profile-profile pairs representing part of query and database 
// profiles on device;
// streamproc, CUDA stream for computations;
// ndb1pros, number of profiles over the db1 positions passed;
// ndbCpros, number of profiles over the complementary dbC positions;
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
// tmpdpdiagbuffers, device memory allocated for temporary diagonal scores;
// tmpdpbotbuffer, device memory allocated for temporary bottom DP scores, 
//  it will serve as a temporary output buffer when finalizing DP;
// maxscoordsbuf, device memory allocated for the coordinates of maximum 
// alignment scores;
// btckdata, device memory allocated for backtracking information;
//
void CuBatchDP::FinalizeDynProgDevice(
    cudaStream_t streamproc,
    CUBSM_TYPE scorethld,
    size_t ndb1pros,
    size_t ndbCpros,
    size_t ndb1prosOmtd,
    size_t ndbCprosOmtd,
    size_t /*dbpro1len*/,
    size_t /*nqyposs*/,
    size_t ndb1poss,
    size_t ndbCposs,
    size_t dbxpad,
    size_t /*querposoffset*/,
    size_t bdb1posoffset,
    size_t bdbCposoffset,
    CUBSM_TYPE* scores,
    CUBDP_TYPE* tmpdpdiagbuffers,
    CUBDP_TYPE* tmpdpbotbuffer,
    unsigned int* maxcoordsbuf,
    char* btckdata,
    //[out:]
    unsigned int* /*attrpassed*/ )
{
    MYMSG( "CuBatchDP::FinalizeDynProgDevice", 4 );
    const mystring preamb = "CuBatchDP::FinalizeDynProgDevice: ";
    myruntime_error mre;

    //size_t ndbxposs = ndb1poss + ndbCposs;

//     //dimensions of a block diagonal
//     const uint dimblkdiag = CUDP_2DCACHE_DIM_D;

    //execution configuration
    dim3 nthrds(1,1,1);//one thread
    dim3 nblcks(1,1,1);//one block

    try {
        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sKernelFinalizeDynProg execution configuration: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, 
                "grid size= (%u,%u,%u) block size= (%u,%u,%u);%s# db pros= %zu",
                nblcks.x, nblcks.y, nblcks.z, nthrds.x, nthrds.y, nthrds.z, NL, ndb1pros+ndbCpros );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        //launch parent grid for finalizing DP results
        FinalizeDP<<<nblcks,nthrds,0,streamproc>>>( 
            scorethld,
            (uint)ndb1pros, (uint)ndbCpros,
            (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)bdb1posoffset, (uint)bdbCposoffset,
            scores,
            tmpdpdiagbuffers, 
            tmpdpbotbuffer,
            maxcoordsbuf,
            btckdata
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
