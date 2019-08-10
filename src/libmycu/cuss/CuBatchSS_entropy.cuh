/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_entropy_h__
#define __CuBatchSS_entropy_h__

#include "extsp/psl.h"
#include "CuBatchSS_com.h"

#define CUSS_RELENT_FAILVAL -1.0f

// =========================================================================
// -------------------------------------------------------------------------
// CalcEntropyProfile: calculate the entropy which is
//    _   lambda s(k)
//   \  e            p(s(k)) s(k) lambda
//   /_
//    k
//
__device__ __forceinline__
void CalcEntropyProfile(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore,
    float lambda )
{
    if( lambda <= 0.0f ) {
        if( threadIdx.y == 0 && threadIdx.x == 0 )
            result[0] = CUSS_RELENT_FAILVAL;
        __syncthreads();
        return;
    }

    float ent = 0.0f;
    float sc = leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN;

    if( minscore <= sc && sc <= maxscore ) {
        float y = lambda * sc;//lambda s(k)
        if( SLC_LOG_SP_MIN < y && y < SLC_LOG_SP_MAX-3.0f )
            ent = __expf(y) * y * scaledprobsCache[threadIdx.y][threadIdx.x];
        //assume an infinitesimally small probability otherwise
    }

    //warp reduce the terms of partial equations (each for threadIdx.y)
    ent += __shfl_down_sync(0xffffffff, ent, 16);
    ent += __shfl_down_sync(0xffffffff, ent, 8);
    ent += __shfl_down_sync(0xffffffff, ent, 4);
    ent += __shfl_down_sync(0xffffffff, ent, 2);
    ent += __shfl_down_sync(0xffffffff, ent, 1);
    //write the partial sum to SMEM
    if( threadIdx.x == 0 )
        tmpsumCache[threadIdx.y] = ent;

    __syncthreads();

    //reduce the partial sums
    if( threadIdx.y == 0 ) {
        ent = 0.0f;
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
            ent = tmpsumCache[threadIdx.x];
        ent += __shfl_down_sync(0xffffffff, ent, 16);
        ent += __shfl_down_sync(0xffffffff, ent, 8);
        ent += __shfl_down_sync(0xffffffff, ent, 4);
        ent += __shfl_down_sync(0xffffffff, ent, 2);
        ent += __shfl_down_sync(0xffffffff, ent, 1);
        if( threadIdx.x == 0 )
            result[0] = ent;
    }

    __syncthreads();
}

// -------------------------------------------------------------------------
// version for unscaled scores and their probabilities
__device__ __forceinline__
void CalcEntropyProfile(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscaledprobsCache[CUSS_N_DIFF_SCORES],
    float minscore,
    float maxscore,
    float lambda )
{
    if( lambda <= 0.0f ) {
        if( threadIdx.y == 0 && threadIdx.x == 0 )
            result[0] = CUSS_RELENT_FAILVAL;
        __syncthreads();
        return;
    }

    float ent = 0.0f;
    float sc = (int)threadIdx.x + CUSS_SCORE_MIN;

    if( threadIdx.y == 0 && threadIdx.x < CUSS_N_DIFF_SCORES && 
        minscore <= sc && sc <= maxscore ) 
    {
        float y = lambda * sc;//lambda s(k)
        if( SLC_LOG_SP_MIN < y && y < SLC_LOG_SP_MAX-3.0f )
            ent = __expf(y) * y * unscaledprobsCache[threadIdx.x];
        //assume an infinitesimally small probability otherwise
    }

    if( threadIdx.y == 0 ) {
        ent += __shfl_down_sync(0xffffffff, ent, 16);
        ent += __shfl_down_sync(0xffffffff, ent, 8);
        ent += __shfl_down_sync(0xffffffff, ent, 4);
        ent += __shfl_down_sync(0xffffffff, ent, 2);
        ent += __shfl_down_sync(0xffffffff, ent, 1);
        if( threadIdx.x == 0 )
            result[0] = ent;
    }

    __syncthreads();
}

#endif//__CuBatchSS_entropy_h__
