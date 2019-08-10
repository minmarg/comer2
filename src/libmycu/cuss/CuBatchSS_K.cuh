/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_K_h__
#define __CuBatchSS_K_h__

#include "extsp/psl.h"
#include "CuBatchSS_com.h"

#define CUSS_K_ACCRCY 1.0e-4f
#define CUSS_K_MAXIT 32

#define CUSS_K_FAILVAL -1.0f

// limits for the minimum and maximum alignment score when 
// calculating the probabilities of alignments of different 
// lengths;
// NOTE: alignments with scores exceeding these limits are highly 
// unlikely, and their probabilities are most often 
// infinitesimally small, even in cases when lambda is small
#define CUSS_K_ALIGNMENT_SCORE_LIMITS
#define CUSS_K_MIN_ALIGNMENT_SCORE -32
#define CUSS_K_MAX_ALIGNMENT_SCORE  32

#define CUSS_K_SETRESULT(pntr,val) \
    if(threadIdx.y==0 && threadIdx.x==0) pntr[0] = val; \
    __syncthreads();


__device__ __forceinline__
float CalcKSigma(
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscldprobsCache[CUSS_N_DIFF_SCORES],
    float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    const int minunscldscore,
    const int maxunscldscore,
    const int unscldscorerange,
    int leftndx,
    int rghtndx,
    float lambda );

__device__ __forceinline__
float CalcProbOfAlnLengthj(
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscldprobsCache[CUSS_N_DIFF_SCORES],
    float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    const int minunscldscore,
    const int maxunscldscore,
    const int unscldscorerange,
    int leftndx,
    int rghtndx,
    int minscore,
    int maxscore,
    float lambda );

// =========================================================================
// -------------------------------------------------------------------------
// CalcKProfile: calculate the statistical parameter K 
// (Karlin et al. PNAS 87, 1990) given score probabilities in 
// scaledprobsCache, the expected score of the scores, expected, lambda, and 
// relative entropy, relent;
// result represents the SMEM address of the resulting value;
// tmpsumCache, array for temporary sum values;
// leftndx and rghtndx, starting and end index of the score for threads 
//  threadIdx.y[0]...threadIdx.y[blockIdx.x-1];
// (NOTE: This implementation takes about 3x the calculation of lambda, 
// which is satisfactory since it requires more than twice the number of 
// iterations [CalcKSigma])
//
// K is related to the EVD location parameter mu by
//   mu = ln (Kmn) / lambda;
// K is calculated by the formula
//
//       gcd lambda exp( -2 sigma )
//   K = --------------------------
//       H (1 - exp( -gcd lambda ))
//
// where gcd is the greatest common divisor of the scores, H is relative
// entropy, and sigma is a measure related to the probability to obtain an 
// arbitrary alignment of any length >=1 given the scoring scheme;
//
//            _   1 (  _           i lambda     _          )
//   sigma = \   -- ( \   P[j](i) e          + \   P[j](i) )
//           /_   j ( /_                       /_          )
//           j>0      i<0                     i>=0
//
// P[j](i) is a probability to obtain an alignment of score i and length j
//              _
//   P[j](i) = \  P[1](k) P[j-1](i-k)
//             /_
//              k
//
// Three special cases: if high score = gcd and low score = -gcd,
//
//   K = [P[1](gcd) - P[1](-gcd)]^2  / p[1](-gcd)
//
// if high score = gcd and low score <> -gcd,
//
//   K = H (1 - exp(-gcd lambda )) / (gcd lambda)
//
// if high score <> gcd and low score = -gcd,
//                                                          _
//   K = lambda (1 - exp(-gcd lambda )) / (gcd H) squared( \  i P[1](i) )
//                                                         /_
//                                                        i=l:u
//
__device__ __forceinline__
void CalcKProfile(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscldprobsCache[CUSS_N_DIFF_SCORES],
    float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    const int minunscldscore,
    const int maxunscldscore,
    const int unscldscorerange,
    int leftndx,
    int rghtndx,
    float expected,
    float lambda,
    float relent )
{
    if( 0.0f <= expected || lambda <= 0.0f || relent <= 0.0f ) {
        CUSS_K_SETRESULT( result, CUSS_K_FAILVAL);
        return;
    }

    if( minunscldscore == -1.0f && maxunscldscore == 1.0f ) {
        //return [p(-1)-p(1)]^2 / p(-1)
        float minprb = unscldprobsCache[-1-CUSS_SCORE_MIN];
        float maxprb = unscldprobsCache[1-CUSS_SCORE_MIN];
        CUSS_K_SETRESULT( result, __fdividef( SQUARE(minprb - maxprb), minprb ));
        return;
    }
    else if( minunscldscore == -1.0f || maxunscldscore == 1.0f ) {
        if( maxunscldscore != 1.0f ) {
            CUSS_K_SETRESULT( result,
                __fdividef(lambda,relent) * (1.0f - __expf(-lambda)) * SQUARE(expected));
            return;
        }
        else {
            CUSS_K_SETRESULT( result,
                __fdividef(relent,lambda) * (1.0f - __expf(-lambda)));
            return;
        }
    }

    float sigma = CalcKSigma(
        tmpsumCache,
        unscldprobsCache,
        scaledprobsCache,
        minunscldscore,
        maxunscldscore,
        unscldscorerange,
        leftndx,
        rghtndx,
        lambda
    );

    CUSS_K_SETRESULT( result,
        __fdividef( lambda*__expf(-2.0f*sigma), relent*(1.0f - __expf(-lambda))) );
}



// -------------------------------------------------------------------------
// CalcKSigma: calculate sigma required for finding K
//
__device__ __forceinline__
float CalcKSigma(
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscldprobsCache[CUSS_N_DIFF_SCORES],
    float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    const int minunscldscore,
    const int maxunscldscore,
    const int unscldscorerange,
    int leftndx,
    int rghtndx,
    float lambda )
{
    float sigma = 0.0f;//related to the probability of alignment of any length
    float Pj = 1.0f;//probability to obtain arbitrary alignment of length j

    int minscore = minunscldscore;
    int maxscore = maxunscldscore;

    //compute sigma associated with the probability of alignment of any length;
    // scaledprobsCache is reused in each iteration to get probabilities of
    // alignment of length one smaller
    for( int j = 1; j <= CUSS_K_MAXIT && CUSS_K_ACCRCY < Pj; j++ )
    {
        Pj = CalcProbOfAlnLengthj(
            tmpsumCache,
            unscldprobsCache,
            scaledprobsCache,
            minunscldscore,
            maxunscldscore,
            unscldscorerange,
            leftndx,
            rghtndx,
            minscore,
            maxscore,
            lambda
        );

        sigma += __fdividef(Pj, j);

// if(threadIdx.y==0&&threadIdx.x==0)
// printf(" *** minscore= %d && maxscore= %d; j= %d sigma= %.4f Pj= %.4f\n",
// minscore,maxscore,j,sigma,Pj);

        minscore += minunscldscore;
        maxscore += maxunscldscore;
        __syncthreads();
    }

    return sigma;
}

// -------------------------------------------------------------------------
// CalcProbOfAlnLengthj: calculate the probability to obtain arbitrary 
// alignment of length j in parallel
// NOTE: the range of unscaled scores unscldscorerange must NOT exceed the 
// warp size or scale factor CUSS_SCALE_FACTOR;
//
__device__ __forceinline__
float CalcProbOfAlnLengthj(
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float unscldprobsCache[CUSS_N_DIFF_SCORES],
    float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    const int minunscldscore,
    const int maxunscldscore,
    const int unscldscorerange,
    int leftndx,
    int /*rghtndx*/,
    int minscore,
    int maxscore,
    float lambda )
{
    //pji is P[j](i), the probability of alignment of length j and score i;
    //j is the iteration number and i is the score the thread will process;
    float pji = 0.0f;
    float pj = 0.0f;//probability of alignment of length j
    int sc = leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN;

    if( minscore <= sc && sc <= maxscore
#ifdef CUSS_K_ALIGNMENT_SCORE_LIMITS
    &&  CUSS_K_MIN_ALIGNMENT_SCORE <= sc && sc <= CUSS_K_MAX_ALIGNMENT_SCORE
#endif
    ) {
        float etols = (sc < 0)? lambda * (float)sc: 0.0f;
        if( SLC_LOG_SP_MIN < etols )
        {
            int minunscld = (maxscore - sc < unscldscorerange)? 
                            maxunscldscore - (maxscore - sc): minunscldscore;
            int maxunscld = (sc - minscore < unscldscorerange)? 
                            minunscldscore + (sc - minscore): maxunscldscore;
            //P[j](i) = SUM P[j-1](i-k) P[1](k)
            for( int i = minunscld; i <= maxunscld; i++ ) {
                int y = threadIdx.y;
                int x = threadIdx.x - i;
                if( (int)blockDim.x <= x ) { 
                    x -= blockDim.x;
                    if( y + 1 < (int)blockDim.y ) y++; else continue;
                }
                else if( x < 0 ) { 
                    x += blockDim.x;
                    //NOTE: for y==1 the result will not be correct because
                    //scaledprobsCache[y==0] places scores starting from the beginning, but
                    //this range of scores will not be used due to small probabilities;
                    //this should be revised in case of the scale factor down-scaled 4x
                    if( 0 < y ) y--; else continue;
                }
                pji += scaledprobsCache[y][x] * unscldprobsCache[i-CUSS_SCORE_MIN];

// printf(" ^^^ minscore= %d && maxscore= %d; tid.y= %u tid.x= %u sc= %d minunscld= %d maxunscld= %d "
// "i= %d y= %d x= %d pji= %.4f\n",
// minscore,maxscore,threadIdx.y,threadIdx.x,sc,minunscld,maxunscld,i,y,x,pji);
            }

            pj = pji * __expf(etols);

// printf(" ^^^ minscore= %d && maxscore= %d; tid.y= %u tid.x= %u sc= %d minunscld= %d maxunscld= %d "
// "etols= %.4f pji= %.4f pj= %.4f\n",
// minscore,maxscore,threadIdx.y,threadIdx.x,sc,minunscld,maxunscld,etols,pji,pj);
        }
    }

    scaledprobsCache[threadIdx.y][threadIdx.x] = pji;

    //warp reduce the terms of partial equations (each for threadIdx.y)
    pj += __shfl_down_sync(0xffffffff, pj, 16);
    pj += __shfl_down_sync(0xffffffff, pj, 8);
    pj += __shfl_down_sync(0xffffffff, pj, 4);
    pj += __shfl_down_sync(0xffffffff, pj, 2);
    pj += __shfl_down_sync(0xffffffff, pj, 1);
    //write the partial sum to SMEM
    if( threadIdx.x == 0 )
        tmpsumCache[threadIdx.y] = pj;

    __syncthreads();

    //reduce the partial sums
    if( threadIdx.y == 0 ) {
        pj = 0.0f;
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
            pj = tmpsumCache[threadIdx.x];
        pj += __shfl_down_sync(0xffffffff, pj, 16);
        pj += __shfl_down_sync(0xffffffff, pj, 8);
        pj += __shfl_down_sync(0xffffffff, pj, 4);
        pj += __shfl_down_sync(0xffffffff, pj, 2);
        pj += __shfl_down_sync(0xffffffff, pj, 1);
        if( threadIdx.x == 0 )
            tmpsumCache[0] = pj;
    }

    __syncthreads();

    return tmpsumCache[0];
}

#endif//__CuBatchSS_K_h__
