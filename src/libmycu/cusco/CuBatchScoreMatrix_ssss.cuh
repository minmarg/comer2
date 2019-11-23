/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ssss_h__
#define __CuBatchScoreMatrix_ssss_h__

#include "libmycu/cupro/SerializedScoresSM.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchScoreMatrix_com.h"

// SSSS_SCORESP1E1:fast reference to a SSS score table!
#define SSSS_SCORESP1E1

//constant shift for the secondary structure scoring function
#define CONSTSSSSSHIFT  0.05f

// GetSSSSTableIndex: function to get a score table index from SS 
// states and probabilities;
// sss, secondary structure state value,
// prb, probability of a secondary structure state
//
__host__ __device__ __forceinline__
int GetSSSSTableIndex( CHTYPE sss, FPTYPE prb )
{
    //NOTE: uncomment and test (wrt sens and aq) the line below! CHANGED!
    return (int)sss * 10 + int(prb * 9.999f);
    //return (int)sss * 10 + int(prb * 10.0f);
}

//device functions for computing SSS scores

__global__ void CalcSM_SSSS_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr attr,
    float ssswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores, 
    CUBSM_TYPE* __restrict__ outmodscores );

// -------------------------------------------------------------------------
// CalcSSSScoreSMEMUnroll2x: calculate two sss scores according to the 
// SMEMUnroll2x approach
//
__device__ __forceinline__
void CalcSSSScoreSMEMUnroll2x( 
    const SerializedScoresAttr& attr,
    const float ssswgt,
    const CUBSM_TYPE* __restrict__ ssssCache,
#if !defined(SSSS_SCORESP1E1)
    const FPTYPE qrenoCache,
    const FPTYPE (* __restrict__ dbenoCache)[SMINIT_2DCACHE_DIM],
#endif
    const FPTYPE (* __restrict__ qrposCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbposCache)[SMINIT_2DCACHE_DIM],
    bool col2cond,
    float& score1, float& score2 )
{
    int ndxtabrow, ndxtabcol;
    FPTYPE f1, f2;
    score1 = 0.0f, score2 = 0.0f;
#if !defined(SSSS_SCORESP1E1)
    SerializedScoresSM<CUBSM_TYPE> ssss(
        const_cast<CUBSM_TYPE*>(ssssCache), attr.szalloc_,
        attr.nplvs_, attr.nenos_, attr.card_,
        attr.ntbls_, attr.nelems_
    );
#endif
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoSSSps; i++ ) {
        //exchanged access for qrposCache
        f1 = qrposCache[i][threadIdx.y];
        ndxtabrow = GetSSSSTableIndex( i/*s1*/, f1 );
        #pragma unroll
        for( int j = 0; j < pmv2DNoSSSps; j++ ) {
            f2 = dbposCache[j][threadIdx.x];
            ndxtabcol = GetSSSSTableIndex( j/*dbsssCache[0][threadIdx.x]*/, f2 );
            score1 += f1 * f2 * 
#ifdef SSSS_SCORESP1E1
                SerializedScoresSM<CUBSM_TYPE>::GetScoreP1E1( ssssCache, ndxtabrow, ndxtabcol );
#else
                //ssss.GetScoreP1( ndxtabrow, ndxtabcol, 0.f, 0.f, qrenoCache, dbenoCache[0][threadIdx.x] );
                ssss.GetScoreP1E1( ndxtabrow, ndxtabcol, 0.f, 0.f, 0.f, 0.f );
#endif
            //
            if( col2cond ) {
                f2 = dbposCache[pmv2DNoSSSps+j][threadIdx.x];
                ndxtabcol = GetSSSSTableIndex( j/*dbsssCache[1][threadIdx.x]*/, f2 );
                score2 += f1 * f2 * 
#ifdef SSSS_SCORESP1E1
                    SerializedScoresSM<CUBSM_TYPE>::GetScoreP1E1( ssssCache, ndxtabrow, ndxtabcol );
#else
                    //ssss.GetScoreP1( ndxtabrow, ndxtabcol, 0.f, 0.f, qrenoCache, dbenoCache[1][threadIdx.x] );
                    ssss.GetScoreP1E1( ndxtabrow, ndxtabcol, 0.f, 0.f, 0.f, 0.f );
#endif
            }
        }
    }

    f1 = qrposCache[pmv2DNoSSSps-1][threadIdx.y];
    if( -99.f < score1 ) {
        //{{NOTE: test (wrt sens and aq) and REMOVE the block below; CHECKED!
        f2 = dbposCache[pmv2DNoSSSps-1][threadIdx.x];
        score1 += (f1 + 0.1f) * (f2 + 0.1f) * 0.1111111f;//(/9)
        //}}
        score1 -= CONSTSSSSSHIFT;
        score1 *= ssswgt;
    }
//     //NOTE: should not be the case
//     else {
//         score1 = 0.0f;
//     }

    if( col2cond && -99.f < score2 ) {
        //{{NOTE: test (wrt sens and aq) and REMOVE the block below; CHECKED!
        f2 = dbposCache[pmv2DNoSSSps+pmv2DNoSSSps-1][threadIdx.x];
        score2 += (f1 + 0.1f) * (f2 + 0.1f) * 0.1111111f;//(/9)
        //}}
        score2 -= CONSTSSSSSHIFT;
        score2 *= ssswgt;
    }
//     //NOTE: the condition score2 < -99.f should not occur
//     else {
//         score2 = 0.0f;
//     }
}

// -------------------------------------------------------------------------
// CalcSSSScoreSMEMUnroll2x_1: calculate one sss score according to the 
// SMEMUnroll2x approach
//
template <int OFF>
__device__ __forceinline__
void CalcSSSScoreSMEMUnroll2x_1( 
    const float ssswgt,
    const CUBSM_TYPE* __restrict__ ssssCache,
    const FPTYPE (* __restrict__ qrposCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbposCache)[SMINIT_2DCACHE_DIM],
    float& score1 )
{
    int ndxtabrow, ndxtabcol;
    FPTYPE f1, f2;
    score1 = 0.0f;
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoSSSps; i++ ) {
        //exchanged access for qrposCache
        f1 = qrposCache[i][threadIdx.y];
        ndxtabrow = GetSSSSTableIndex( i/*s1*/, f1 );
        #pragma unroll
        for( int j = 0; j < pmv2DNoSSSps; j++ ) {
            f2 = dbposCache[OFF+j][threadIdx.x];
            ndxtabcol = GetSSSSTableIndex( j, f2 );
            score1 += f1 * f2 * 
                SerializedScoresSM<CUBSM_TYPE>::GetScoreP1E1( ssssCache, ndxtabrow, ndxtabcol );
        }
    }

    f1 = qrposCache[pmv2DNoSSSps-1][threadIdx.y];
    if( -99.f < score1 ) {
        //{{NOTE: test (wrt sens and aq) and REMOVE the block below; CHECKED!
        f2 = dbposCache[OFF+pmv2DNoSSSps-1][threadIdx.x];
        score1 += (f1 + 0.1f) * (f2 + 0.1f) * 0.1111111f;//(/9)
        //}}
        score1 -= CONSTSSSSSHIFT;
        score1 *= ssswgt;
    }
}

#endif//__CuBatchScoreMatrix_ssss_h__
