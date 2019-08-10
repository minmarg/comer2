/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_Init_h__
#define __CuBatchScoreMatrix_Init_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchScoreMatrix_com.h"

//constant shift for the initial profile-profile substitution scoring function
#define CONSTINITSHIFT  0.002f //0.0003f

//device variable for testing
extern __device__ unsigned long long d_testvar;

//device functions for computing initial scores

__global__ void CalcSM_Init_SMEMUnroll2x(
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* outscores, 
    CUBSM_TYPE* outmodscores );

__global__ void CalcSMInitSMEMUnroll2(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores, 
    CUBSM_TYPE* outmodscores );

__global__ void CalcSMInitSMEM(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* outmodscores );

__global__ void CalcSMInit(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t queroffset, size_t bdb1offset, size_t bdbCoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* outmodscores );

//test kernel
__global__ void KernelStreamedIterativeTest(
    uint serqposnr,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores );
__global__ void KernelParentStreamedIterativeTest(
    uint serqposnr,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores );

// -------------------------------------------------------------------------
// CalcInitScoreSMEMUnroll2x: calculate two initial scores according to the 
// SMEMUnroll2x approach
//
__device__ __forceinline__
void CalcInitScoreSMEMUnroll2x( 
    const FPTYPE* __restrict__ qrproCache,
    const FPTYPE (* __restrict__ dbproCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ qrposCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbposCache)[SMINIT_2DCACHE_DIM],
    bool col2cond,
    float& score1, float& score2 )
{
    score1 = 0.0f, score2 = 0.0f;
    FPTYPE f1, f2;

    #pragma unroll
    for( int i = 0; i < pmv2DNoElems; i++ ) {
        //exchanged access for qrposCache
        f1 = qrposCache[i][threadIdx.y];
        f2 = qrproCache[i];
        score1 += __fdividef( 
            ( f1 * dbposCache[i][threadIdx.x]), 
            ( 0.5f * (f2 + dbproCache[i][threadIdx.x]))
        );
        if( col2cond ) {
            score2 += __fdividef( 
                ( f1 * dbposCache[pmv2DNoElems+i][threadIdx.x]), 
                ( 0.5f * (f2 + dbproCache[pmv2DNoElems+i][threadIdx.x]))
            );
        }
    }

    MYASSERT(score1>0.0f,"Invalid score1.");
    score1 = __logf(score1) - CONSTINITSHIFT;

    if( col2cond ) {
        MYASSERT(score2>0.0f,"Invalid score2.");
        score2 = __logf(score2) - CONSTINITSHIFT;
    }
}

#endif//__CuBatchScoreMatrix_Init_h__
