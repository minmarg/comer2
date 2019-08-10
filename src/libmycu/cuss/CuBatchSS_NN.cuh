/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_NN_h__
#define __CuBatchSS_NN_h__

#include "extsp/psl.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchSS_com.h"

#define CUSS_NN_MU_OUTPUTFACTOR 30.0f
#define CUSS_NN_SCALE_OUTPUTFACTOR 5.0f

enum Tcuss_NN_Ouputs{
    cuss_NN_E_MU,
    cuss_NN_L_MU,
    cuss_NN_L_SCALE,
    nTcuss_NN_Ouputs
};

// =========================================================================
// -------------------------------------------------------------------------
// PredictParameterbyNNconfEProfile: predict a statistical parameter 
// using an ANN model of one layer and four hidden units with the hyperbolic 
// tangent activation function;
// assuming four input variables under configuration E (prediction based on 
// profile attributes);
// output, SMEM address the result to be written to;
// tmpbuf, SMEM buffer for temporary data;
// inputs, input values;
// weights, ANN weights;
//
__device__ __forceinline__
void PredictParameterbyNNconfEProfile(
    float* __restrict__ output,
    float* __restrict__ tmpbuf,
    const float* __restrict__ inputs,
    const float* __restrict__ weights )
{
    if( threadIdx.x < CUSS_NHIDDEN * CUSS_E_NINPUT )
        tmpbuf[threadIdx.x] = 
            //inputs[threadIdx.x % CUSS_E_NINPUT] * w; when CUSS_E_NINPUT is 2^n:
            inputs[threadIdx.x & (CUSS_E_NINPUT-1)] * 
            //weights[threadIdx.x + 1 + threadIdx.x/CUSS_E_NINPUT]: 1st in the row is the bias
            weights[threadIdx.x+1+(threadIdx.x>>CUSS_LOG_E_NINPUT)];

    __syncthreads();

    if( threadIdx.x < CUSS_NHIDDEN ) {
        tmpbuf[threadIdx.x] = 
            weights[threadIdx.x*(CUSS_E_NINPUT+1)] + //bias term
            tmpbuf[threadIdx.x*CUSS_E_NINPUT] + tmpbuf[threadIdx.x*CUSS_E_NINPUT+1] +
            tmpbuf[threadIdx.x*CUSS_E_NINPUT+2] + tmpbuf[threadIdx.x*CUSS_E_NINPUT+3];
        tmpbuf[threadIdx.x] = 
            //activation function: tanh(x*.5) (.5, steepness):
            (__fdividef(2.0f, 1.0f + __expf(-tmpbuf[threadIdx.x])) - 1.0f) *
            //weight of the next layer (output); the 1st in the row is the bias:
            weights[CUSS_NHIDDEN*(CUSS_E_NINPUT+1) + threadIdx.x+1];
    }

    __syncthreads();

    if( threadIdx.x == 0 ) {
        output[0] = 
            weights[CUSS_NHIDDEN*(CUSS_E_NINPUT+1)] + //bias term
            tmpbuf[threadIdx.x] + tmpbuf[threadIdx.x+1] +
            tmpbuf[threadIdx.x+2] + tmpbuf[threadIdx.x+3];
        output[0] = 
            //activation function: tanh(x*.5) (.5, steepness):
            (__fdividef(2.0f, 1.0f + __expf(-output[0])) - 1.0f) *
            CUSS_NN_MU_OUTPUTFACTOR;
    }

    __syncthreads();
}

// -------------------------------------------------------------------------
// PredictParametersby2NNsConfLProfile: predict statistical parameters 
// by two ANN models of one layer and four hidden units with the hyperbolic 
// tangent activation function SIMULTANEOUSLY;
// assuming four input variables (with one dummy) under configuration L 
// (prediction based on profile compositional similarity);
// output, SMEM address the TWO predictions to be written to;
// tmpbuf, SMEM buffer for temporary data;
// inputs, input values;
// weights, weights of TWO ANNs;
//
__device__ __forceinline__
void PredictParametersby2NNsConfLProfile(
    float* __restrict__ output,
    float* __restrict__ tmpbuf,
    const float* __restrict__ inputs,
    const float* __restrict__ weights )
{
    int offset = (threadIdx.x < CUSS_NHIDDEN*CUSS_L_NINPUT)? 0: CUSS_NHIDDEN+1;

    //the whole warp participates:
    //if( threadIdx.x < 2 * CUSS_NHIDDEN * CUSS_L_NINPUT )
        tmpbuf[threadIdx.x] = 
            //inputs[threadIdx.x % CUSS_L_NINPUT] * w; when CUSS_L_NINPUT is 2^n:
            inputs[threadIdx.x & (CUSS_L_NINPUT-1)] * 
            //weights[threadIdx.x + 1 + threadIdx.x/CUSS_L_NINPUT]: 1st in the row is the bias
            weights[threadIdx.x+1+(threadIdx.x>>CUSS_LOG_L_NINPUT) + offset];

    __syncthreads();

    if( threadIdx.x < 2 * CUSS_NHIDDEN ) {
        offset = (threadIdx.x < CUSS_NHIDDEN)? 0: CUSS_NHIDDEN+1;
        tmpbuf[threadIdx.x] = 
            weights[threadIdx.x*(CUSS_L_NINPUT+1) + offset] + //bias term
            tmpbuf[threadIdx.x*CUSS_L_NINPUT] + tmpbuf[threadIdx.x*CUSS_L_NINPUT+1] +
            tmpbuf[threadIdx.x*CUSS_L_NINPUT+2] + tmpbuf[threadIdx.x*CUSS_L_NINPUT+3];
        tmpbuf[threadIdx.x] = 
            //activation function: tanh(x*.5) (.5, steepness):
            (__fdividef(2.0f, 1.0f + __expf(-tmpbuf[threadIdx.x])) - 1.0f) *
            //weight of the next layer (output); the 1st in the row is the bias:
            weights[((threadIdx.x>>CUSS_LOG_NHIDDEN)+1) * CUSS_NHIDDEN * (CUSS_E_NINPUT+1) +
                    offset + 
                    (threadIdx.x & (CUSS_NHIDDEN-1)) + 1];
    }

    __syncthreads();

    if( threadIdx.x < 2 ) {
        offset = (threadIdx.x == 0)? 0: CUSS_NHIDDEN+1;
        float factor = (threadIdx.x == 0)? CUSS_NN_MU_OUTPUTFACTOR: CUSS_NN_SCALE_OUTPUTFACTOR;
        output[threadIdx.x] = 
            weights[(threadIdx.x+1) * CUSS_NHIDDEN * (CUSS_E_NINPUT+1) + offset] + //bias term
            tmpbuf[CUSS_NHIDDEN*threadIdx.x] + tmpbuf[CUSS_NHIDDEN*threadIdx.x+1] +
            tmpbuf[CUSS_NHIDDEN*threadIdx.x+2] + tmpbuf[CUSS_NHIDDEN*threadIdx.x+3];
        output[threadIdx.x] = 
            //activation function: tanh(x*.5) (.5, steepness):
            (__fdividef(2.0f, 1.0f + __expf(-output[threadIdx.x])) - 1.0f) *
            factor;
    }

    __syncthreads();
}

// -------------------------------------------------------------------------
// FormatInputConfigurationE: format input for a NN for predicting a 
// statistical parameter based on profile attributes (configuration E)
// inputs, SMEM array of input values to the NN;
// 
__device__ __forceinline__
void FormatInputConfigurationE(
    uint nqyposs, uint db1prolen,
    float qyeno, float db1proeno,
    float inputs[4] )
{
    if( threadIdx.x == 0 ) {
        if( qyeno < db1proeno ) {
            //exchange: NN trained asymmetrically
            myhdswap( qyeno, db1proeno );
            myhdswap( nqyposs, db1prolen );
        }
        else if( qyeno == db1proeno && nqyposs < db1prolen )
            myhdswap( nqyposs, db1prolen );

        inputs[0] = 0.05f * qyeno;
        inputs[1] = __fdividef(__logf((float)nqyposs), SLC_LN1K);
        inputs[2] = 0.05f * db1proeno;
        inputs[3] = __fdividef(__logf((float)db1prolen), SLC_LN1K);
    }

    __syncthreads();
}

// -------------------------------------------------------------------------
// FormatInputConfigurationL: format input for a NN for predicting 
// statistical parameters based on profile compositional similarity 
// (configuration L);
// inputs, SMEM array of input values to the NN;
// 
__device__ __forceinline__
void FormatInputConfigurationL(
    uint nqyposs, uint db1prolen,
    float lambda,
    float inputs[4] )
{
    if( threadIdx.x == 0 ) {
        lambda = __fdividef( rintf(lambda*10.0f), 10.0f);
        if( lambda <= 0.0f )
            lambda = 1.0f;
        if( 10.0f < lambda )
            lambda = 10.0f;
        if( nqyposs < db1prolen )
            //exchange: NN trained asymmetrically
            myhdswap( nqyposs, db1prolen );

        inputs[0] = 0.2f * lambda;
        inputs[1] = __fdividef(__logf((float)nqyposs), SLC_LN1K);
        inputs[2] = __fdividef(__logf((float)db1prolen), SLC_LN1K);
        inputs[3] = 0.0f;//dummy last input
    }

    __syncthreads();
}

#endif//__CuBatchSS_NN_h__
