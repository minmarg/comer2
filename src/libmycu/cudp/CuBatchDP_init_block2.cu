/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"
#include "CuBatchDP_init.cuh"
#include "CuBatchDP_init_block2.cuh"

// =========================================================================
// NOTE: THIS KERNEL IS FOR PERFORMANCE EVALUATION ONLY AND IS NOT 
// THOROUGHLY TESTED
// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming;
// NOTE: Version for thread block of one warp!
// blkdiagnum, block diagonal serial number;
// lastydiagnum, last block diagonal serial number along y axis 
// (starting at x=-CUDP_2DCACHE_DIM);
// ndb1pros, number of profiles in the first profile data buffer db1;
// querprosOmtd, number of query profiles up to this query profile;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// nqyposs, number of query positions to process;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// dbxpad, number of padded positions for memory alignment;
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
//

// DP processing layout:
// +-----====----+-----===+--+-
// |  |  |  |  | |  |  |  |  |
// |  |  |  |  | |  |  |  |  |
// |  |  |  |  | |  |  |  |  |
// +--=======----+--======+--+-
// |  |  |  |  | |  |  |  |  |
// +--====-------+--====--+--+-
// (double line indicates current parallel processing)

// -------------------------------------------------------------------------
// ExecDPBlock2Unroll32x: device code for executing dynamic programming 
// using shared memory and 32-fold unrolling along the diagonal of dimension 
// CUDP_2DCACHE_DIM;
// this version of DP uses rectangular thread blocks and the threads of a 
// block perform twice as long as those of a diagonal block 
// (keeping in mind warps execute simultaneously);
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers for last calculated diagonal scores;
// tmpdpbotbuffer, temporary buffers for last calculated bottom scores;
// 
__global__ void ExecDPBlock2Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer )
{
    /*__shared__ */LNTYPE dbprodstCache;//distance in positions to db profile blockIdx.y
    /*__shared__ */INTYPE dbprolenCache;//length of profile blockIdx.y
    __shared__ FPTYPE 
            qrtrpCache[nTDPUsedProTrn][CUDP_2DCACHE_DIM_D],//cache for query transition probabilities
            dbtrpCache[nTDPUsedProTrn][CUDP_2DCACHE_DIM_X];//cache for db transition probabilities
    //cache of scores: use padding of +1 to eliminate bank conflicts when accessing the cache
    __shared__ CUBSM_TYPE
            scoreCache[CUDP_2DCACHE_DIM_D][CUDP_2DCACHE_DIM_X+CUDPIDS];
    __shared__ CUBDP_TYPE
            vertiCache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_D+1],//cache for state scores at the boundary of block
            diag1Cache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_D+1],//cache for state scores of the 1st diagonal
            diag2Cache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_D+1],//cache for state scores of the last diagonal
            bottmCache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_X+1];//cache for state scores of the bottom of the diagonals
    //__shared__ CUBDP_TYPE
    //        maxscCache[CUDP_2DCACHE_DIM_D];//cache for maximum scores of the last processed diagonal
    CUBDP_TYPE maxscCache;
    //
    // blockIdx.y is the profile serial number;
    uint dbfldsndx;
    uint pronr;
    //NOTE: protection against overflow ensured on the host side
    if( blockIdx.y < ndb1pros ) { pronr = blockIdx.y + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = blockIdx.y - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( /*threadIdx.y == 0 && */threadIdx.x == 0 ) {
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[pronr];
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[pronr];
    }

    // blockIdx.x is block serial number s within diagonal blkdiagnum;
    // (x,y) is the bottom-left corner (x,y) coordinates for profile blockIdx.y
    int x, y;
    if( blkdiagnum <= lastydiagnum) {
        //x=bs; y=dw+w-sw -1 (-1, zero-based indices);
        //(b, block's length; w, block's width)
        x = blockIdx.x * CUDP_2DCACHE_DIM_X;
        y = (blkdiagnum + 1 - blockIdx.x) * CUDP_2DCACHE_DIM_D - 1;
    } else {
        //x=(d-d_l)w+ws; y=dw+w-sw -1;
        x = (blockIdx.x + (blkdiagnum-lastydiagnum)) * CUDP_2DCACHE_DIM_X;
        y = (lastydiagnum + 1 - blockIdx.x) * CUDP_2DCACHE_DIM_D - 1;
    }

    //{{use registers efficiently
    dbprodstCache = __shfl_sync(0xffffffff, dbprodstCache, 0);
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}
    //...or use shared memory; bank conflicts arise when accessed
//     __syncthreads();


    //number of iterations for this block to perform;
    //uint makes the compiler to reduce #registers
    uint ilim = GetMaxNoIterationsBlock2( x, y, nqyposs, dbprolenCache );

    if( y < 0 || (int)nqyposs <= (int)(y+1 - blockDim.x) || 
        (int)dbprolenCache <= x /*+ CUDP_2DCACHE_DIM_DpX */ ||
        ilim < 1 )
        //block does not participate: out of profile boundaries
        return;

    int qpos = y - threadIdx.x;//going upwards

    if( 0 <= qpos && qpos < (int)nqyposs )
        // (add querprosOmtd vectors of beginning transitions):
        DPLocCacheTrnProbs( qrtrpCache, 0, querposoffset + qpos + querprosOmtd );
    else
        DPLocInitTrnProbs( qrtrpCache );

    //x is now the position this thread will process
    x += threadIdx.x;

    //db profile position corresponding to the diagonal block's bottom-left 
    // corner in the buffers dc_pm2dvfields_:
    int dbpos = x + dbprodstCache;//going right
    int dblen = ndb1poss + ndbCposs + dbxpad;

#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( 0 <= x && x < (int)dbprolenCache )
            // (add pronr vectors of beginning transitions):
            DPBlock2LocCacheTrnProbs( dbtrpCache, dbfldsndx, dbpos + pronr );
        else
            DPBlock2LocInitTrnProbs( dbtrpCache );
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif


    //dbpos is now the x position of the diagonal block's bottom-left 
    // corner in the score matrix plus the offset determined by thread id:
    dbpos = (blockIdx.y < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;

    //cache scores over the thread block
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        #pragma unroll 4
        for( int i = 0; i < CUDP_2DCACHE_DIM_D; i++ ) {
            //going upwards
            scoreCache[i][threadIdx.x] = CUBDP_Q(0);
            if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x && x < (int)dbprolenCache ) {
                //starting position of line i of the block in the score matrix: 
                qpos = (y-i) * dblen;
                scoreCache[i][threadIdx.x] = scores[qpos + dbpos];
            }
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif

    //cache the vertical from the previous block;
    //the structure of tmpdpdiagbuffers is position-specific 
    //(1D, layed down on the x-axis)
    if( 0 <= x-CUDP_2DCACHE_DIM_D && x < (int)dbprolenCache ) {
        DPLocCacheBuffer( vertiCache, tmpdpdiagbuffers, dbpos-CUDP_2DCACHE_DIM_D, 0, dblen );
        //cache the buffer of maximum scores here
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
        //maxscCache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbpos - CUDP_2DCACHE_DIM_D];
        maxscCache = tmpdpdiagbuffers[qpos + dbpos - CUDP_2DCACHE_DIM_D];
    }
    else {
        DPLocInitCache(vertiCache);
        //maxscCache[threadIdx.x] = CUBDP_Q(0);
        maxscCache = CUBDP_Q(0);
    }
    DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag1Cache, 1/*shift*/);
    DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag2Cache, 1/*shift*/);

    //cache the bottom of the above blocks;
    //the structure of tmpdpbotbuffer is position-specific (1D, along x-axis)
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( CUDP_2DCACHE_DIM_D <= y && 
            0 <= x && x < (int)dbprolenCache ) {
            DPLocCacheBuffer<CUDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer, dbpos, 0, dblen, 1/*shift*/);
        }
        else {
            DPLocInitCache<CUDP_2DCACHE_DIM_X>(bottmCache, 1/*shift*/);
        }
        if( threadIdx.x == 0 ) {
            DPLocInitCache<CUDP_2DCACHE_DIM_X>(bottmCache, 0/*shift*/);
            if( CUDP_2DCACHE_DIM_D <= y && 0 <= x-1 )
                DPLocCacheBuffer<CUDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer, dbpos-1, 0, dblen, 0/*shift*/);
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif


    __syncthreads();

    CUBDP_TYPE (*pdiag1)[CUDP_2DCACHE_DIM_D+1] = diag1Cache;
    CUBDP_TYPE (*pdiag2)[CUDP_2DCACHE_DIM_D+1] = diag2Cache;

    //start calculations for this position with 32x unrolling
    //NOTE: sync inside: do not branch;
    for( int i = 0; i < ilim/*CUDP_2DCACHE_DIM_X*/; i++ ) {
        CUBDP_TYPE val1, val2;
        int dbtrndx = i-(CUDP_2DCACHE_DIM_D-1-threadIdx.x);
        dbtrndx = dbtrndx<=0? 0: (CUDP_2DCACHE_DIM_X<=dbtrndx? CUDP_2DCACHE_DIM_X-1: dbtrndx);
        //
        if( threadIdx.x == blockDim.x-1 && i < CUDP_2DCACHE_DIM_X ) {
            pdiag2[dpdsssStateMM][threadIdx.x+1] = bottmCache[dpdsssStateMM][i];
            pdiag2[dpdsssStateMI][threadIdx.x+1] = bottmCache[dpdsssStateMI][i];
            pdiag2[dpdsssStateIM][threadIdx.x+1] = bottmCache[dpdsssStateIM][i];
            pdiag2[dpdsssStateDN][threadIdx.x+1] = bottmCache[dpdsssStateDN][i];
            pdiag2[dpdsssStateND][threadIdx.x+1] = bottmCache[dpdsssStateND][i];
            //
            pdiag1[dpdsssStateMM][threadIdx.x+1] = bottmCache[dpdsssStateMM][i+1];
            pdiag1[dpdsssStateMI][threadIdx.x+1] = bottmCache[dpdsssStateMI][i+1];
            pdiag1[dpdsssStateIM][threadIdx.x+1] = bottmCache[dpdsssStateIM][i+1];
            pdiag1[dpdsssStateDN][threadIdx.x+1] = bottmCache[dpdsssStateDN][i+1];
            pdiag1[dpdsssStateND][threadIdx.x+1] = bottmCache[dpdsssStateND][i+1];
        }
        if( threadIdx.x == CUDP_2DCACHE_DIM_D-i-1 ) {
            if( i ) {
                pdiag2[dpdsssStateMM][threadIdx.x+1] = vertiCache[dpdsssStateMM][i-1];
                pdiag2[dpdsssStateMI][threadIdx.x+1] = vertiCache[dpdsssStateMI][i-1];
                pdiag2[dpdsssStateIM][threadIdx.x+1] = vertiCache[dpdsssStateIM][i-1];
                pdiag2[dpdsssStateDN][threadIdx.x+1] = vertiCache[dpdsssStateDN][i-1];
                pdiag2[dpdsssStateND][threadIdx.x+1] = vertiCache[dpdsssStateND][i-1];
            }
            pdiag1[dpdsssStateMM][threadIdx.x] = vertiCache[dpdsssStateMM][i];
            pdiag1[dpdsssStateMI][threadIdx.x] = vertiCache[dpdsssStateMI][i];
            pdiag1[dpdsssStateIM][threadIdx.x] = vertiCache[dpdsssStateIM][i];
            pdiag1[dpdsssStateDN][threadIdx.x] = vertiCache[dpdsssStateDN][i];
            pdiag1[dpdsssStateND][threadIdx.x] = vertiCache[dpdsssStateND][i];
        }
        __syncthreads();

        //MM state update (diagonal direction)
        val1 = pdiag2[dpdsssStateMI][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIMp][dbtrndx];
        val2 = pdiag2[dpdsssStateIM][threadIdx.x+1] + 
            qrtrpCache[dptrIMp][threadIdx.x] + dbtrpCache[dptrMMp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateDN][threadIdx.x+1] + 
            qrtrpCache[dptrDMp][threadIdx.x] + dbtrpCache[dptrMMp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateND][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrDMp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateMM][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMMp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        //accessing scores is without bank conflicts because of padding
        val1 = myhdmax/*fmaxf*/( val1 + scoreCache[threadIdx.x][dbtrndx], CUBDP_Q(0));

        __syncthreads();//last point pdiag2 is used for reading in this iteration
        pdiag2[dpdsssStateMM][threadIdx.x] = val1;//WRITE
        //maxscCache[threadIdx.x] = myhdmax/*fmaxf*/( maxscCache[threadIdx.x], val1 );
        maxscCache = myhdmax/*fmaxf*/( maxscCache, val1 );

        //MI state update (up direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMIc][dbtrndx];
        val2 = pdiag1[dpdsssStateMI][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIIc][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateMI][threadIdx.x] = val1;//WRITE

        //IM state update (left direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMIc][threadIdx.x] + dbtrpCache[dptrMMp][dbtrndx];
        val2 = pdiag1[dpdsssStateIM][threadIdx.x] + 
            qrtrpCache[dptrIIc][threadIdx.x] + dbtrpCache[dptrMMp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateIM][threadIdx.x] = val1;//WRITE

        //DN state update (up)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + qrtrpCache[dptrMDp][threadIdx.x];
        val2 = pdiag1[dpdsssStateDN][threadIdx.x+1] + qrtrpCache[dptrDDp][threadIdx.x];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateDN][threadIdx.x] = val1;//WRITE

        //ND state update (left)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + dbtrpCache[dptrMDp][dbtrndx];
        val2 = pdiag1[dpdsssStateND][threadIdx.x] + dbtrpCache[dptrDDp][dbtrndx];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateND][threadIdx.x] = val1;//WRITE

        if( threadIdx.x == 0 && CUDP_2DCACHE_DIM_D <= i+1 ) {
            //WRITE
            //this position is not used by other threads in the current iteration
            bottmCache[dpdsssStateMM][i+1-CUDP_2DCACHE_DIM_D] = pdiag2[dpdsssStateMM][threadIdx.x];
            bottmCache[dpdsssStateMI][i+1-CUDP_2DCACHE_DIM_D] = pdiag2[dpdsssStateMI][threadIdx.x];
            bottmCache[dpdsssStateIM][i+1-CUDP_2DCACHE_DIM_D] = pdiag2[dpdsssStateIM][threadIdx.x];
            bottmCache[dpdsssStateDN][i+1-CUDP_2DCACHE_DIM_D] = pdiag2[dpdsssStateDN][threadIdx.x];
            bottmCache[dpdsssStateND][i+1-CUDP_2DCACHE_DIM_D] = pdiag2[dpdsssStateND][threadIdx.x];
        }
        if( CUDP_2DCACHE_DIM_D-1-threadIdx.x == i+1-CUDP_2DCACHE_DIM_X ) {
            //WRITE to the right boundary (vertical) of the block 
            vertiCache[dpdsssStateMM][i+1-CUDP_2DCACHE_DIM_X] = pdiag2[dpdsssStateMM][threadIdx.x];
            vertiCache[dpdsssStateMI][i+1-CUDP_2DCACHE_DIM_X] = pdiag2[dpdsssStateMI][threadIdx.x];
            vertiCache[dpdsssStateIM][i+1-CUDP_2DCACHE_DIM_X] = pdiag2[dpdsssStateIM][threadIdx.x];
            vertiCache[dpdsssStateDN][i+1-CUDP_2DCACHE_DIM_X] = pdiag2[dpdsssStateDN][threadIdx.x];
            vertiCache[dpdsssStateND][i+1-CUDP_2DCACHE_DIM_X] = pdiag2[dpdsssStateND][threadIdx.x];
        }

        myhdswap( pdiag1, pdiag2 );

        __syncthreads();
    }


    //write the result of calculations for following blocks;
    //write two diagonals;
    if( 0 <= x && x < (int)dbprolenCache ) {
        DPLocWriteBuffer<CUDP_2DCACHE_DIM_D>( vertiCache, tmpdpdiagbuffers, dbpos, 0, dblen );
        //write the buffer of maximum scores
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
        //tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxscCache[threadIdx.x];
        if( CUDP_2DCACHE_DIM_D <= y ) {
            if( tmpdpdiagbuffers[qpos + dbpos] < maxscCache )
                tmpdpdiagbuffers[qpos + dbpos] = maxscCache;
        } else
            tmpdpdiagbuffers[qpos + dbpos] = maxscCache;
    }

    //write the bottom of the diagonal block;
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( 0 <= x && x < (int)dbprolenCache )
            DPLocWriteBuffer<CUDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer, dbpos, 0, dblen );
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif
}
