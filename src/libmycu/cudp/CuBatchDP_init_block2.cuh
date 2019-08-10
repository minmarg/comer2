/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_init_block2_h__
#define __CuBatchDP_init_block2_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"
#include "CuBatchDP_init.cuh"

//device functions for executing dynamic programming

__global__ void ExecDPBlock2Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer 
);


// =========================================================================
// -------------------------------------------------------------------------
// GetMaqxNoIterations: get the maximum number of iterations for the 
// block to perform 
// 
__device__ __forceinline__
int GetMaxNoIterationsBlock2(int x, int y, int quelen, int dbprolen)
{
    x = dbprolen - x;
    if( x <= 0 )
        return 0;
    y -= quelen;
    if( CUDP_2DCACHE_DIM_D <= y )
        return 0;
    y = ( 0 < y )? CUDP_2DCACHE_DIM_D - y: CUDP_2DCACHE_DIM_D;
    if( CUDP_2DCACHE_DIM_X < x )
        x = CUDP_2DCACHE_DIM_X;
    return x + y - 1;
}

// -------------------------------------------------------------------------
// DPBlock2LocInitTrnProbs: initialize transition probabilities
//
__device__ __forceinline__
void DPBlock2LocInitTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_X])
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         smtrpCache[i][threadIdx.x] = CUBDP_Q(-32768);
//     }
    smtrpCache[dptrMMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMIc][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMDp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIIc][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDDp][threadIdx.x] = DPLocGetLog0();
}

// -------------------------------------------------------------------------
// DPBlock2LocCacheTrnProbs: cache transition probabilities to SMEM at 
// position pos
//
__device__ __forceinline__
void DPBlock2LocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_X],
    int fldsndx,
    int pos )
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         if( i == dptrMIc || i == dptrIIc )
//             smtrpCache[i][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i+1, pos );
//         else
//             smtrpCache[i][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i, pos );
//     }
    smtrpCache[dptrMMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMMp, pos );
    smtrpCache[dptrMIc][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMIc, pos+1 );
    smtrpCache[dptrMDp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMDp, pos );
    smtrpCache[dptrIMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrIMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrIIc][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrIIc, pos+1, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrDMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDDp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrDDp, pos, CUDP_TRNPRB_LOGFCT );
}

// -------------------------------------------------------------------------
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
// DPLocInitCache: initialize cache of a DP buffer
// shft, shift of inner-most index for writing values
//
__device__ __forceinline__
void DPLocInitCache( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X+1],
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++ )
        smdpbufCache[i][threadIdx.x+shft] = CUBDP_Q(0);
}

// DPLocCacheBuffer: cache a DP buffers;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to read data from;
// x, x position at which the buffer has to be read;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index for writing values
__device__ __forceinline__
void DPLocCacheBuffer( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X+1],
    const CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        smdpbufCache[i][threadIdx.x+shft] = gmemtmpbuffer[y + x];
    }
}
// DPLocWriteBuffer: write a cahced DP buffer back to GMEM;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to write data to;
// x, x position at which the data has to be written;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index of the cache
__device__ __forceinline__
void DPLocWriteBuffer( 
    const CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X+1],
    CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        gmemtmpbuffer[y + x] = smdpbufCache[i][threadIdx.x+shft];
    }
}

#endif
// -------------------------------------------------------------------------

#endif//__CuBatchDP_init_block2_h__
