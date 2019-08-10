/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_lambda_h__
#define __CuBatchSS_lambda_h__

#include "extsp/psl.h"
#include "CuBatchSS_com.h"

#define CUSS_LAMBDA_LOWERBND 1.0e-5f
#define CUSS_LAMBDA_UPPERBND 10.0f
#define CUSS_LAMBDA_ACCRCY 1.0e-4f //1.0e-5f
#define CUSS_LAMBDA_MAXIT 32

#define CUSS_LAMBDA_FAILVAL -1.0f

//define to solve the lambda equation represented in a polynomial form;
//leads to a smaller number of iterations but also to reduced 
// instruction throughput per iteration
// #define CUSS_LAMBDA_REPRSNT_POLY
#define CUSS_LAMBDA_POLY_LOWERBND 0.99999f
#define CUSS_LAMBDA_POLY_UPPERBND 4.55e-5f

#define CUSS_LAMBDA_LEADING_THREAD_BEG if(threadIdx.y==0 && threadIdx.x==0) {
#define CUSS_LAMBDA_LEADING_THREAD_END } __syncthreads();


#ifdef CUSS_LAMBDA_REPRSNT_POLY
#define CUSS_LMB_X1 CUSS_LAMBDA_POLY_LOWERBND
#define CUSS_LMB_X2 CUSS_LAMBDA_POLY_UPPERBND
#else
#define CUSS_LMB_X1 CUSS_LAMBDA_LOWERBND
#define CUSS_LMB_X2 CUSS_LAMBDA_UPPERBND
#endif


__device__ __forceinline__
void FindLmbdRootByNRandBisection(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore );

__device__ __forceinline__
void GetFdFlmbdequation(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore,
    float* f, 
    float* df );

__device__ __forceinline__
void GetFlmbdequation(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore,
    float* f );

// =========================================================================
// -------------------------------------------------------------------------
// CalcLambdaProfile: calculate the statistical parameter lambda given 
// score probabilities in scaledprobsCache and the expected score of the 
// scores, expected;
// result, SMEM address to locate the lambda root on exit;
// tmpsumCache, array for temporary sum values;
// leftndx, starting index of the score for threads 
//  threadIdx.y[0]...threadIdx.y[blockIdx.x-1]
//
__device__ __forceinline__
void CalcLambdaProfile(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float expected,
    float minscore,
    float maxscore )
{
    if( 0.0f <= expected || maxscore <= 0.0f ) {
        CUSS_LAMBDA_LEADING_THREAD_BEG
            *result = CUSS_LAMBDA_FAILVAL;
        CUSS_LAMBDA_LEADING_THREAD_END
        return;
    }

    FindLmbdRootByNRandBisection(
        result,
        tmpsumCache,
        scaledprobsCache,
        leftndx,
        minscore,
        maxscore
    );

#ifdef CUSS_LAMBDA_REPRSNT_POLY
    CUSS_LAMBDA_LEADING_THREAD_BEG
        if( result[0] <= 0.0f )
            result[0] = CUSS_LAMBDA_FAILVAL;
        else
            result[0] = -__logf(result[0]);
    CUSS_LAMBDA_LEADING_THREAD_END
#endif
}



// -------------------------------------------------------------------------
// FindLmbdRootByNRandBisection: find the root of a function bracketed 
//  between x1 and x2 using a combination of Newton-Raphson and bisection 
//  methods.
//  A Bisection step is taken whenever Newton-Raphson would take the solution
//  out of bounds, or whenever it is not reducing the size of the brackets
//  rapidly enough.
//  The root will be refined until its accuracy is reached within ±xacc or
//  maximum number of iteration, maxit, has been passed. 
//  The function value and the first derivative of the function are 
//  executed in parallel.
// Computation logic is as found in
//  Numerical Recipes in C by W.H.Press et al.
//
__device__ __forceinline__
void FindLmbdRootByNRandBisection(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore )
{
    __shared__ int flagSM;

    float   f, df, dx, dxold, fl;
    float   xh, xl, rts;
    int n;

    //tmpsumCache[0] is used to pass the argument (lambda) to all the 
    //threads of the block
    CUSS_LAMBDA_LEADING_THREAD_BEG
        result[0] = CUSS_LMB_X1;
    CUSS_LAMBDA_LEADING_THREAD_END

    GetFlmbdequation( result, tmpsumCache,scaledprobsCache,leftndx,minscore,maxscore,  &fl );

    CUSS_LAMBDA_LEADING_THREAD_BEG
        result[0] = CUSS_LMB_X2;
    CUSS_LAMBDA_LEADING_THREAD_END

    GetFlmbdequation( result, tmpsumCache,scaledprobsCache,leftndx,minscore,maxscore,  &f );

    CUSS_LAMBDA_LEADING_THREAD_BEG
        flagSM = ( fl > 0.0f && f > 0.0f ) || ( fl < 0.0f && f < 0.0f );
        if( flagSM )
            result[0] = CUSS_LAMBDA_FAILVAL;
    CUSS_LAMBDA_LEADING_THREAD_END

    if( flagSM )
        return;

    CUSS_LAMBDA_LEADING_THREAD_BEG
        if( fl < 0.0f ) {           //Orient the search so that f (xl) < 0.
            xl = CUSS_LMB_X1;
            xh = CUSS_LMB_X2;
        } else {
            xh = CUSS_LMB_X1;
            xl = CUSS_LMB_X2;
        }
        rts = 0.5f * ( CUSS_LMB_X1 + CUSS_LMB_X2 ); //Initialize the guess for root,
        dxold = fabsf( CUSS_LMB_X2 - CUSS_LMB_X1 ); //the "stepsize before last,"
        dx = dxold;                 //and the last step.
        result[0] = rts;
    CUSS_LAMBDA_LEADING_THREAD_END

    GetFdFlmbdequation( result, tmpsumCache,scaledprobsCache,leftndx,minscore,maxscore,  &f, &df );

    for( n = 0; n < CUSS_LAMBDA_MAXIT; n++ )
    {
        CUSS_LAMBDA_LEADING_THREAD_BEG                          //Bisect if Newton out of range,
            if( ((( rts - xh ) * df - f ) * (( rts - xl ) * df - f ) > 0.0f ) ||
                ( fabs( 2.0f * f ) > fabsf( dxold * df ))) {    //or not decreasing fast enough.
                dxold = dx;
                dx = 0.5f * ( xh - xl );
                rts = xl + dx;
                if( xl == rts ) flagSM = 1;     //Change in root is negligible.
            } else {                            //Newton step acceptable. Take it.
                dxold = dx;
                dx = __fdividef(f, df);
                fl/*temp*/ = rts;
                rts -= dx;
                if( fl/*temp*/ == rts ) flagSM = 1;
            }
            if( fabsf( dx ) < CUSS_LAMBDA_ACCRCY )
                flagSM = 1;                     //Convergence criterion.
            result[0] = rts;
        CUSS_LAMBDA_LEADING_THREAD_END

        if( flagSM )
            break;

        GetFdFlmbdequation( result, tmpsumCache,scaledprobsCache,leftndx,minscore,maxscore,  &f, &df );

        //The one new function evaluation per iteration maintain the bracket on the root.
        if( f < 0.0f )
            xl = rts;
        else
            xh = rts;
    }

    //NOTE: all threads update n
    if( n == CUSS_LAMBDA_MAXIT ) {
        CUSS_LAMBDA_LEADING_THREAD_BEG
            result[0] = CUSS_LAMBDA_FAILVAL;
        CUSS_LAMBDA_LEADING_THREAD_END
    }

    return;
}

// -------------------------------------------------------------------------
// GetFdFlmbdequation: calculate the equation and its derivative
//    _   lambda s(k)
//   \  e            p(s(k)) - 1
//   /_
//    k
// the equation is divided by exp(max s(k) lambda) for numerical stability;
// // // Optionally, the equation can be transformed into the polynomial form:
// // //     exp(max s(k) lambda) * poly(exp(-lambda)),
// // // where poly is a polynomial of exp(-lambda) and has exactly two zeros:
// // // 1 and another one in interval (0,1);
// // // x is then supposed to be exp(-lambda)
//
__device__ __forceinline__
void GetFdFlmbdequation(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore,
    float* f, 
    float* df )
{
    float x = result[0];

    float lf = 0.0f;
    float ldf = 0.0f;
    float sc = leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN;

    if( minscore <= sc && sc <= maxscore ) {
#ifdef CUSS_LAMBDA_REPRSNT_POLY
        sc = maxscore - sc;//division by the max score for numerical stability
        float y = __powf(x, sc-1.0f) * //x^sc=exp[-lambda(max s(k)-s(k))]
                scaledprobsCache[threadIdx.y][threadIdx.x];
        lf = y * x;
        ldf = y * sc;
#else
        sc -= maxscore;//subtract max for numerical stability
        float y = x * sc;//lambda(s(k)-max s(k))
        if( SLC_LOG_SP_MIN < y ) {
            lf = __expf(y) * scaledprobsCache[threadIdx.y][threadIdx.x];
            ldf = lf * sc;
        }
#endif
    }

    //warp reduce the terms of partial equations (each for threadIdx.y)
    lf += __shfl_down_sync(0xffffffff, lf, 16);
    lf += __shfl_down_sync(0xffffffff, lf, 8);
    lf += __shfl_down_sync(0xffffffff, lf, 4);
    lf += __shfl_down_sync(0xffffffff, lf, 2);
    lf += __shfl_down_sync(0xffffffff, lf, 1);
    //same for ldf
    ldf += __shfl_down_sync(0xffffffff, ldf, 16);
    ldf += __shfl_down_sync(0xffffffff, ldf, 8);
    ldf += __shfl_down_sync(0xffffffff, ldf, 4);
    ldf += __shfl_down_sync(0xffffffff, ldf, 2);
    ldf += __shfl_down_sync(0xffffffff, ldf, 1);
    //write the partial sum to SMEM
    if( threadIdx.x == 0 ) {
        tmpsumCache[threadIdx.y] = lf;
        tmpsumCache[threadIdx.y+CUSS_N_DIFF_SCORES] = ldf;
    }
    __syncthreads();

    //reduce the partial sums
    if( threadIdx.y == 0 ) {
        lf = ldf = 0.0f;
        if( threadIdx.x < CUSS_N_DIFF_SCORES ) {
            lf = tmpsumCache[threadIdx.x];
            ldf = tmpsumCache[threadIdx.x+CUSS_N_DIFF_SCORES];
        }
        //lf
        lf += __shfl_down_sync(0xffffffff, lf, 16);
        lf += __shfl_down_sync(0xffffffff, lf, 8);
        lf += __shfl_down_sync(0xffffffff, lf, 4);
        lf += __shfl_down_sync(0xffffffff, lf, 2);
        lf += __shfl_down_sync(0xffffffff, lf, 1);
        //same for ldf
        ldf += __shfl_down_sync(0xffffffff, ldf, 16);
        ldf += __shfl_down_sync(0xffffffff, ldf, 8);
        ldf += __shfl_down_sync(0xffffffff, ldf, 4);
        ldf += __shfl_down_sync(0xffffffff, ldf, 2);
        ldf += __shfl_down_sync(0xffffffff, ldf, 1);
        //subtract the right-hand side
#ifdef CUSS_LAMBDA_REPRSNT_POLY
        float y = __powf(x, maxscore-1.0f);//x^max s(k)=exp[-lambda max s(k)]
        lf -= y * x;
        ldf -= y * maxscore;
#else
        x *= -maxscore;//-lambda max s(k)
        if( SLC_LOG_SP_MIN < x ) {
            x = __expf(x);
            lf -= x;
            ldf += x * maxscore;
        }
#endif
    }
    //save results
    *df = ldf;
    *f = lf;
}

// -------------------------------------------------------------------------
// GetFlmbdequation: same as GetFdFlmbdequation but the derivative is not 
// calculated
//
__device__ __forceinline__
void GetFlmbdequation(
    float* __restrict__ result,
    float tmpsumCache[TIMES2(CUSS_N_DIFF_SCORES)],
    const float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1],
    int leftndx,
    float minscore,
    float maxscore,
    float* f )
{
    float x = result[0];

    float lf = 0.0f;
    float sc = leftndx + (int)threadIdx.x + CUSS_SCALED_SCORE_MIN;

    if( minscore <= sc && sc <= maxscore ) {
#ifdef CUSS_LAMBDA_REPRSNT_POLY
        sc = maxscore - sc;//division by the max score for numerical stability
        //x^sc=exp[-lambda(max s(k)-s(k))]
        lf = __powf(x, sc) * scaledprobsCache[threadIdx.y][threadIdx.x];
#else
        sc -= maxscore;//subtract max for numerical stability
        float y = x * sc;//lambda(s(k)-max s(k))
        if( SLC_LOG_SP_MIN < y )
            lf = __expf(y) * scaledprobsCache[threadIdx.y][threadIdx.x];
#endif
    }

    //warp reduce the terms of partial equations (each for threadIdx.y)
    lf += __shfl_down_sync(0xffffffff, lf, 16);
    lf += __shfl_down_sync(0xffffffff, lf, 8);
    lf += __shfl_down_sync(0xffffffff, lf, 4);
    lf += __shfl_down_sync(0xffffffff, lf, 2);
    lf += __shfl_down_sync(0xffffffff, lf, 1);
    //write the partial sum to SMEM
    if( threadIdx.x == 0 )
        tmpsumCache[threadIdx.y] = lf;

    __syncthreads();

    //reduce the partial sums
    if( threadIdx.y == 0 ) {
        lf = 0.0f;
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
            lf = tmpsumCache[threadIdx.x];
        //lf
        lf += __shfl_down_sync(0xffffffff, lf, 16);
        lf += __shfl_down_sync(0xffffffff, lf, 8);
        lf += __shfl_down_sync(0xffffffff, lf, 4);
        lf += __shfl_down_sync(0xffffffff, lf, 2);
        lf += __shfl_down_sync(0xffffffff, lf, 1);
        //subtract the right-hand side
#ifdef CUSS_LAMBDA_REPRSNT_POLY
        lf -= __powf(x, maxscore);//exp[-lambda max s(k)]
#else
        x *= -maxscore;//-lambda max s(k)
        if( SLC_LOG_SP_MIN < x )
            lf -= __expf(x);
#endif
    }
    //save the result
    *f = lf;
}

#endif//__CuBatchSS_lambda_h__
