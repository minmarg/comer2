/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_cvs2s_h__
#define __CuBatchScoreMatrix_cvs2s_h__

#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchScoreMatrix_com.h"

// CVS2S_SCORESE1S1Step2:fast reference to a cvs2s map!
#define CVS2S_SCORESE1S1Step2

//constant shift for the context vector scoring function
#define CONSTCVSSHIFT  0.5f

//device functions for computing CVS scores

__global__ void CalcSM_CVS2S_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr attr,
    float cvswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    FPTYPE(*roundfunc)(FPTYPE),
    CUBSM_TYPE* __restrict__ outscores, 
    CUBSM_TYPE* __restrict__ outmodscores );

// -------------------------------------------------------------------------
// rounding functors dependent on a execution configuration
struct FtorRoundMAPALN1 {
    __device__ __forceinline__
    FPTYPE operator()(FPTYPE arg) {return arg;}
};
struct FtorRoundMAPALN0 {
    __device__ __forceinline__
    FPTYPE operator()(FPTYPE arg) {return (FPTYPE)(int)arg;}
};

// -------------------------------------------------------------------------
// CalcCVS2ScoreSMEMUnroll2x: calculate two cvs2s scores according to the 
// SMEMUnroll2x approach
//
template <typename Func>
__device__ __forceinline__
void CalcCVS2ScoreSMEMUnroll2x( 
    const SerializedCVS2ScoresAttr& attr,
    const float cvswgt,
    const CUBSM_TYPE* __restrict__ cvs2sCache,
#if !defined(CVS2S_SCORESE1S1Step2)
    const FPTYPE qrenoCache,
    const FPTYPE (* __restrict__ dbenoCache)[SMINIT_2DCACHE_DIM],
#endif
    const FPTYPE* __restrict__ qrcvpriorCache,
    const FPTYPE (* __restrict__ dbcvpriorCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE* __restrict__ qrcvnm2Cache,
    const FPTYPE (* __restrict__ dbcvnm2Cache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ qrcvCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbcvCache)[SMINIT_2DCACHE_DIM],
    Func roundfunc,
    bool col2cond,
    float& score1, float& score2 )
{
    FPTYPE f1;
    score1 = 0.0f, score2 = 0.0f;
#if !defined(CVS2S_SCORESE1S1Step2)
    SerializedCVS2Scores<CUBSM_TYPE> cvs2s(
        const_cast<CUBSM_TYPE*>(cvs2sCache), attr.szalloc_,
        attr.nenos_, attr.ntbls_
    );
#endif
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoCVEls; i++ ) {
        f1 = qrcvCache[i][threadIdx.y];
        score1 += f1 * dbcvCache[i][threadIdx.x];
        if( col2cond )
            score2 += f1 * dbcvCache[pmv2DNoCVEls+i][threadIdx.x];
    }

    score1 += attr.CVS_loKAPPA0_;//add off-diagonal entry of C_2 to the dot product
    //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
    f1 = qrcvnm2Cache[threadIdx.y];
    score1 = (f1 + 1.0f + attr.CVS_loKAPPA0_) *
             (dbcvnm2Cache[0][threadIdx.x] + 1.0f + attr.CVS_loKAPPA0_) - 
             score1 * score1;
    MYASSERT( score1 > 0.0f, "Det 1.");
    //log odds score of normal vectors
    //NOTE: test (wrt sens and aq) and REMOVE typecast to int: CHANGED!
    score1 = attr.CVS_CTERM_ + attr.CVS_PowerNU0_ * __logf(score1) - 
             roundfunc(qrcvpriorCache[threadIdx.y]) - roundfunc(dbcvpriorCache[0][threadIdx.x]);
    //translate score
    score1 = 
#ifdef CVS2S_SCORESE1S1Step2
        //SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2( cvs2sCache, score1, attr.card_, attr.shft_ );
        SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2Boundary( 
            cvs2sCache, score1, attr.card_, attr.shft_,
            attr.key1first_, attr.value1first_, attr.key1last_, attr.value1last_);
#else
        //cvs2s.GetScore( score1, qrenoCache, dbenoCache[0][threadIdx.x] );
        cvs2s.GetScoreE1S1Step2( score1, attr.card_, attr.shft_, 0.f, 0.f );
#endif

    //if( score1 )
    score1 *= cvswgt;

    if( col2cond ) {
        score2 += attr.CVS_loKAPPA0_;//add off-diagonal entry of C_2 to the dot product
        //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
        score2 = (f1 + 1.0f + attr.CVS_loKAPPA0_) *
                (dbcvnm2Cache[1][threadIdx.x] + 1.0f + attr.CVS_loKAPPA0_) - 
                score2 * score2;
        MYASSERT( score2 > 0.0f, "Det 2.");
        //log odds score of normal vectors
        //NOTE: test (wrt sens and aq) and REMOVE typecast to int: CHANGED!
        score2 = attr.CVS_CTERM_ + attr.CVS_PowerNU0_ * __logf(score2) - 
                roundfunc(qrcvpriorCache[threadIdx.y]) - roundfunc(dbcvpriorCache[1][threadIdx.x]);
        //translate score
        score2 = 
#ifdef CVS2S_SCORESE1S1Step2
            //SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2( cvs2sCache, score2, attr.card_, attr.shft_ );
            SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2Boundary(
                cvs2sCache, score2, attr.card_, attr.shft_,
                attr.key1first_, attr.value1first_, attr.key1last_, attr.value1last_);
#else
            //cvs2s.GetScore( score2, qrenoCache, dbenoCache[1][threadIdx.x] );
            cvs2s.GetScoreE1S1Step2( score2, attr.card_, attr.shft_, 0.f, 0.f );
#endif
        //if( score2 )
        score2 *= cvswgt;
    }
}

// -------------------------------------------------------------------------
// CalcCVS2ScoreSMEMUnroll2x_1: calculate one cvs2s score according to the 
// SMEMUnroll2x approach
//
template <typename Func, int OFF, int IND>
__device__ __forceinline__
void CalcCVS2ScoreSMEMUnroll2x_1( 
    const SerializedCVS2ScoresAttr& attr,
    const float cvswgt,
    const CUBSM_TYPE* __restrict__ cvs2sCache,
    const FPTYPE* __restrict__ qrcvpriorCache,
    const FPTYPE (* __restrict__ dbcvpriorCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE* __restrict__ qrcvnm2Cache,
    const FPTYPE (* __restrict__ dbcvnm2Cache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ qrcvCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbcvCache)[SMINIT_2DCACHE_DIM],
    Func roundfunc,
    float& score1 )
{
    FPTYPE f1;
    score1 = 0.0f;
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoCVEls; i++ ) {
        f1 = qrcvCache[i][threadIdx.y];
        score1 += f1 * dbcvCache[OFF+i][threadIdx.x];
    }

    score1 += attr.CVS_loKAPPA0_;//add off-diagonal entry of C_2 to the dot product
    //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
    f1 = qrcvnm2Cache[threadIdx.y];
    score1 = (f1 + 1.0f + attr.CVS_loKAPPA0_) *
             (dbcvnm2Cache[IND][threadIdx.x] + 1.0f + attr.CVS_loKAPPA0_) - 
             score1 * score1;
    MYASSERT( score1 > 0.0f, "Det 1.");
    //log odds score of normal vectors
    //NOTE: test (wrt sens and aq) and REMOVE typecast to int: CHANGED!
    score1 = attr.CVS_CTERM_ + attr.CVS_PowerNU0_ * __logf(score1) - 
             roundfunc(qrcvpriorCache[threadIdx.y]) - roundfunc(dbcvpriorCache[IND][threadIdx.x]);
    //translate score
    score1 = 
        //SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2( cvs2sCache, score1, attr.card_, attr.shft_ );
        SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2Boundary( 
            cvs2sCache, score1, attr.card_, attr.shft_,
            attr.key1first_, attr.value1first_, attr.key1last_, attr.value1last_);

    //if( score1 )
    score1 *= cvswgt;
}

#include "CuBatchScoreMatrix_cvs2s.cu"

#endif//__CuBatchScoreMatrix_cvs2s_h__
