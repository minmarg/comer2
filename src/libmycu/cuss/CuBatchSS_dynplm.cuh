/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_dynplm_h__
#define __CuBatchSS_dynplm_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_com.h"

//device functions for calculating alignment statistics using 
// dynamic parallelism

__global__ void CalcStatisticsDynPlm(
    bool Xuninf,
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    const int ssemodel,
    float qyeno,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpbuffer,
    float* __restrict__ tmpss2datbuffers,
    float* __restrict__ dp2alndatbuffers
);

// =========================================================================
// InitS2DatBufferForMAPDP: tmpss2datbuffers is not used from this point on, 
// but it is used in the MAPDP calculations;
// initialize tmpss2datbuffers for the total probability for each 
// pair of profiles calculated during the MAP dynamic programming
//
__device__ __forceinline__
void InitS2DatBufferForMAPDP(
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    if( threadIdx.x == 0 )
        tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES] = 0.0f;
}

// -------------------------------------------------------------------------

#endif//__CuBatchSS_dynplm_h__
