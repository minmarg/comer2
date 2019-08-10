/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_init_corr_h__
#define __CuBatchDP_init_corr_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"

//device functions for executing dynamic programming

__global__ void ExecDP_Corr_Unroll32x(
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
// DPLocCacheCorrBuffer: cache a DP buffer of scores for calculating 
// correlation scores;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to read data from;
// x, x position at which the buffer has to be read;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index for writing values
//
__device__ __forceinline__
void DPLocCacheCorrBuffer( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_D],
    const CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++, y += stride ) {
        smdpbufCache[i][threadIdx.x+shft] = gmemtmpbuffer[y + x];
    }
}

// DPLocInitCache: initialize cache of a DP buffer
// shft, shift of inner-most index for writing values
//
__device__ __forceinline__
void DPLocInitCorrCache( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_D],
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++ )
        smdpbufCache[i][threadIdx.x+shft] = CUBDP_Q(0);
}

// -------------------------------------------------------------------------
// DPLocWriteCorrBuffer: write a buffer associated with correlation scores 
// back to GMEM;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to write data to;
// x, x position at which the data has to be written;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index of the cache
//
__device__ __forceinline__
void DPLocWriteCorrBuffer( 
    const CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_D],
    CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++, y += stride ) {
        gmemtmpbuffer[y + x] = smdpbufCache[i][threadIdx.x+shft];
    }
}

// -------------------------------------------------------------------------
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
__device__ __forceinline__
void DPLocCacheCorrBuffer( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X],
    const CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++, y += stride ) {
        smdpbufCache[i][threadIdx.x] = gmemtmpbuffer[y + x];
    }
}
__device__ __forceinline__
void DPLocInitCorrCache( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X] )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++ )
        smdpbufCache[i][threadIdx.x] = CUBDP_Q(0);
}
__device__ __forceinline__
void DPLocWriteCorrBuffer( 
    const CUBDP_TYPE (* __restrict__ smdpbufCache)[CUDP_2DCACHE_DIM_X],
    CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride )
{
    #pragma unroll
    for( int i = 0; i < CUDP_CORR_NSCORES; i++, y += stride ) {
        gmemtmpbuffer[y + x] = smdpbufCache[i][threadIdx.x];
    }
}
#endif//!defined(CUDP_2DCACHE_DIM_DequalsX)


#endif//__CuBatchDP_init_corr_h__
