/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_minmax_h__
#define __CuBatchSS_minmax_h__

#include "libmycu/cucom/cucommon.h"
#include "CuBatchSS_com.h"

// =========================================================================
// -------------------------------------------------------------------------
// GetScoreMinMaxProfile: get the minimum and maximum of scores written in a 
// structured SMEM block;
// result_min, SMEM address where to write the minimum value;
// result_max, SMEM address where to write the maximum value;
// tmpCache, array for temporary values;
// scaledprobsCache, score probabilities;
// leftndx, starting index of the score for threads 
//  threadIdx.y[0]...threadIdx.y[blockIdx.x-1]
//
__device__ __forceinline__
void GetScoreMinMaxProfile(
    float* __restrict__ result_min,
    float* __restrict__ result_max,
    float tmpCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx )
{
    float min = leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN;
    min = scaledprobsCache[threadIdx.y][threadIdx.x] ? min: 0.0f;
    float max = min;

    //warp reduce the partial terms (each for threadIdx.y)
    min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 16));
    min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 8));
    min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 4));
    min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 2));
    min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 1));
    //same for max
    max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 16));
    max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 8));
    max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 4));
    max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 2));
    max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 1));
    //write the partial reductions to SMEM
    if( threadIdx.x == 0 ) {
        tmpCache[threadIdx.y] = min;
        tmpCache[threadIdx.y+CUSS_N_DIFF_SCORES] = max;
    }
    __syncthreads();

    //reduce the partial reductions
    if( threadIdx.y == 0 ) {
        if( threadIdx.x < CUSS_N_DIFF_SCORES ) {
            min = tmpCache[threadIdx.x];
            max = tmpCache[threadIdx.x+CUSS_N_DIFF_SCORES];
        }
        //min
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 16));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 8));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 4));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 2));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 1));
        //max
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 16));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 8));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 4));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 2));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 1));
        if( threadIdx.x == 0 ) {
            result_min[0] = min;
            result_max[0] = max;
        }
    }
    __syncthreads();
}

// -------------------------------------------------------------------------
// GetScoreMinMaxProfile: get the minimum and maximum of scores written in a 
// structured SMEM block; (version for unscaled scores);
// result_min, SMEM address where to write the minimum value;
// result_max, SMEM address where to write the maximum value;
// tmpCache, array for temporary values;
// unscaledprobsCache, probabilities of unscaled scores;
//
__device__ __forceinline__
void GetScoreMinMaxProfile(
    float* __restrict__ result_min,
    float* __restrict__ result_max,
    float tmpCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscaledprobsCache[CUSS_N_DIFF_SCORES] )
{
    float min = (int)threadIdx.x + CUSS_SCORE_MIN;
    min = (threadIdx.y == 0 && threadIdx.x < CUSS_N_DIFF_SCORES && 
            unscaledprobsCache[threadIdx.x]) ? min: 0.0f;
    float max = min;

    if( threadIdx.y == 0 ) {
        //min
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 16));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 8));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 4));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 2));
        min = myhdmin( min, __shfl_down_sync(0xffffffff, min, 1));
        //max
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 16));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 8));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 4));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 2));
        max = myhdmax( max, __shfl_down_sync(0xffffffff, max, 1));
        if( threadIdx.x == 0 ) {
            result_min[0] = min;
            result_max[0] = max;
        }
    }
    __syncthreads();
}

#endif//__CuBatchSS_minmax_h__
