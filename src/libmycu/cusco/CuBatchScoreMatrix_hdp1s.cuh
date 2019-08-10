/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_hdp1s_h__
#define __CuBatchScoreMatrix_hdp1s_h__

#include "extsp/psl.h"
#include "libmycu/cupro/SerializedScoresTM.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchScoreMatrix_com.h"

#define ADJWGTDG 0.45f

#define HDP1_ZZICOV 0.47f
#define HDP1_ZWICOV 0.47f
#define HDP1_XXICOV 0.37f
#define HDP1_XYICOV 1.17f

// HDP1S_SCORESP1E2:fast reference to scores!
#define HDP1S_SCORESP1E2


// GetFactorWrtENOAndHDPPrb: get score factor wrt to the effective 
// number of observations and HDP cluster membership probability of two 
// positions from two profiles;
// ensqry, enssub, effective numbers of observations;
// prbqry, prbsub, HDP cluster membership probabilities;
//
__device__ __forceinline__
float GetFactorWrtENOAndHDPPrb( 
    FPTYPE ensqry, FPTYPE enssub, FPTYPE prbqry, FPTYPE prbsub )
{
//     ensqry -= 23.0f; ensqry *= ensqry;
//     enssub -= 23.0f; enssub *= enssub;
//     prbqry -= 1.0f;
//     prbsub -= 1.0f;
//     return
//         (0.4f + 0.6f * __expf(-1.5e-5f *
//             (HDP1_ZZICOV * ensqry * ensqry + HDP1_ZZICOV * enssub * enssub +
//                 2.0f * HDP1_ZWICOV * ensqry * enssub)) )
//         *
//         (__expf( -7.5f *
//             (HDP1_XXICOV * prbqry * prbqry + HDP1_XXICOV * prbsub * prbsub +
//                 2.0f * HDP1_XYICOV * prbqry * prbsub)) )
//     ;
    //it seems that the compiler is able to optimize better
    return
        (0.4f + 0.6f * __expf(-1.5e-5f *(
            HDP1_ZZICOV * SLC_POW4(ensqry-23.0f) + HDP1_ZZICOV * SLC_POW4(enssub-23.0f) +
                2.0f * HDP1_ZWICOV * SLC_POW2(ensqry-23.0f) * SLC_POW2(enssub-23.0f))) )
        *
        (__expf( -7.5f * (
            HDP1_XXICOV * SLC_POW2(prbqry-1.0f) + HDP1_XXICOV * SLC_POW2(prbsub-1.0f) +
                2.0f * HDP1_XYICOV * (prbqry-1.0f)*(prbsub-1.0f))) )
    ;
}

// =========================================================================
// device functions for computing SSS scores
//
__global__ void CalcSM_HDP1S_SMEMUnroll2x(
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr attr,
    float hdp1swgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* outscores, 
    CUBSM_TYPE* outmodscores );

// -------------------------------------------------------------------------
// CalcHDP1ScoreSMEMUnroll2x: calculate two hdp1 scores according to the 
// SMEMUnroll2x approach
//
__device__ __forceinline__
void CalcHDP1ScoreSMEMUnroll2x( 
    const SerializedScoresAttr& attr,
    const float hdp1swgt,
    cudaTextureObject_t hdp1sTexo,
    const FPTYPE qrenoCache,
    const FPTYPE (* __restrict__ dbenoCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE* __restrict__ qrhdp1prbCache,
    const FPTYPE (* __restrict__ dbhdp1prbCache)[SMINIT_2DCACHE_DIM],
    const INTYPE* __restrict__ qrhdp1ndxCache,
    const INTYPE (* __restrict__ dbhdp1ndxCache)[SMINIT_2DCACHE_DIM],
    bool col2cond,
    float& score1, float& score2 )
{
    score1 = 0.0f, score2 = 0.0f;
    FPTYPE e1 = qrenoCache;
    INTYPE n1 = qrhdp1ndxCache[threadIdx.y];
    FPTYPE f1 = qrhdp1prbCache[threadIdx.y];
    //
    //using registers is a bit more effective although access to 
    //SMEM without bank conflicts (below) is almost as effective
    FPTYPE e21 = dbenoCache[0][threadIdx.x];
    FPTYPE e22 = dbenoCache[1][threadIdx.x];
    INTYPE n21 = dbhdp1ndxCache[0][threadIdx.x];
    INTYPE n22 = dbhdp1ndxCache[1][threadIdx.x];
    FPTYPE f21 = dbhdp1prbCache[0][threadIdx.x];
    FPTYPE f22 = dbhdp1prbCache[1][threadIdx.x];

#ifdef HDP1S_SCORESP1E2
    if( 0 <= n1 ) {
        if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
            score1 = SerializedScoresTM<CUBSM_TYPE>::GetScoreP1E2( 
                hdp1sTexo, attr.eth1_, attr.eth2_, attr.nelems_,
                n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
                e1, e21//dbenoCache[0][threadIdx.x]
            );
        if( col2cond && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
            score2 = SerializedScoresTM<CUBSM_TYPE>::GetScoreP1E2( 
                hdp1sTexo, attr.eth1_, attr.eth2_, attr.nelems_,
                n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
                e1, e22//dbenoCache[1][threadIdx.x]
            );
    }
#else
    if( 0 <= n1 ) {
        SerializedScoresTM<CUBSM_TYPE> hdp1s(
            hdp1sTexo, attr.szalloc_,
            attr.nplvs_, attr.nenos_, attr.card_,
            attr.ntbls_, attr.nelems_
        );
//         if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
//             score1 = hdp1s.GetScore( 
//                 n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
//                 f1, f21,//dbhdp1prbCache[0][threadIdx.x], 
//                 e1, e21//dbenoCache[0][threadIdx.x]
//             );
//         if( col2cond && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
//             score2 = hdp1s.GetScore( 
//                 n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
//                 f1, f22,//dbhdp1prbCache[1][threadIdx.x], 
//                 e1, e22//dbenoCache[1][threadIdx.x]
//             );
        if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
            score1 = hdp1s.GetScoreP1( 
                n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
                0.0f, 0.0f, 
                e1, e21//dbenoCache[0][threadIdx.x]
            );
        if( col2cond && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
            score2 = hdp1s.GetScoreP1( 
                n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
                0.0f, 0.0f, 
                e1, e22//dbenoCache[1][threadIdx.x]
            );
    }
#endif

    if( -99.0f < score1 ) {
        score1 *= GetFactorWrtENOAndHDPPrb(
                        e1, e21,//dbenoCache[0][threadIdx.x],
                        f1, f21);//dbhdp1prbCache[0][threadIdx.x]);
        score1 *= hdp1swgt;
    }
    else score1 = 0.0f;

    if( col2cond ) {
        if( -99.0f < score2 ) {
            score2 *= GetFactorWrtENOAndHDPPrb(
                            e1, e22,//dbenoCache[1][threadIdx.x],
                            f1, f22);//dbhdp1prbCache[1][threadIdx.x]);
            score2 *= hdp1swgt;
        }
        else score2 = 0.0f;
    }
}

#endif//__CuBatchScoreMatrix_hdp1s_h__
