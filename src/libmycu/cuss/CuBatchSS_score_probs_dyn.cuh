/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_score_probs_dyn_h__
#define __CuBatchSS_score_probs_dyn_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_com.h"

//device functions for calculating score probabilities using 
// dynamic parallelism

__global__ void CalcScoreProbsDynPlm(
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

#endif//__CuBatchSS_score_probs_dyn_h__
