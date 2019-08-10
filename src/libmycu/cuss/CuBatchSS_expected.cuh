/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_expected_h__
#define __CuBatchSS_expected_h__

#include "libmycu/cucom/cucommon.h"
#include "CuBatchSS_com.h"

// =========================================================================
// -------------------------------------------------------------------------
// CalcExpectedScoreProfile: calculate the expected score of scores 
// written in a structured SMEM block;
// result, SMEM address to write the result at;
// tmpCache, array for temporary values;
// scaledprobsCache, score probabilities;
// leftndx, starting index of the score for threads 
//  threadIdx.y[0]...threadIdx.y[blockIdx.x-1]
//
__device__ __forceinline__
void CalcExpectedScoreProfile(
    float* __restrict__ result,
    float* __restrict__ tmpCache,
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx )
{
    float expsc = (leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN) *
            scaledprobsCache[threadIdx.y][threadIdx.x];

    //warp reduce the partial terms (each for threadIdx.y)
    expsc += __shfl_down_sync(0xffffffff, expsc, 16);
    expsc += __shfl_down_sync(0xffffffff, expsc, 8);
    expsc += __shfl_down_sync(0xffffffff, expsc, 4);
    expsc += __shfl_down_sync(0xffffffff, expsc, 2);
    expsc += __shfl_down_sync(0xffffffff, expsc, 1);
    //write the partial reductions to SMEM
    if( threadIdx.x == 0 )
        tmpCache[threadIdx.y] = expsc;

    __syncthreads();

    //reduce the partial reductions
    if( threadIdx.y == 0 ) {
        expsc = 0.0f;
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
            expsc = tmpCache[threadIdx.x];
        expsc += __shfl_down_sync(0xffffffff, expsc, 16);
        expsc += __shfl_down_sync(0xffffffff, expsc, 8);
        expsc += __shfl_down_sync(0xffffffff, expsc, 4);
        expsc += __shfl_down_sync(0xffffffff, expsc, 2);
        expsc += __shfl_down_sync(0xffffffff, expsc, 1);
        if( threadIdx.x == 0 )
            result[0] = expsc;
    }

    __syncthreads();
}

#endif//__CuBatchSS_expected_h__
