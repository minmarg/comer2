/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_score_probs_h__
#define __CuBatchSS_score_probs_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_com.h"

//device functions for executing dynamic programming

__global__ void CalcScoreProbs(
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpss2datbuffers
);

// =========================================================================
// -------------------------------------------------------------------------
// VerifyScore: check whether a score is within the allowed interval of 
// scores
__device__ __forceinline__
CUBSM_TYPE VerifyScore( CUBSM_TYPE score )
{
    if( score < CUSS_SCORE_MIN_f )
        return CUSS_SCORE_MIN_f;
    else if( score > CUSS_SCORE_MAX_f )
        return CUSS_SCORE_MAX_f;
    return score;
}

#endif//__CuBatchSS_score_probs_h__
