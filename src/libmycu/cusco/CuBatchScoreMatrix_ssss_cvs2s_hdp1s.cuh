/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ssss_cvs2s_hdp1s_h__
#define __CuBatchScoreMatrix_ssss_cvs2s_hdp1s_h__

#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "CuBatchScoreMatrix_com.h"

//device functions for computing SSS, CVS, and HDP scores at once

template <typename Func>
__global__ void CalcSM_SSSS_CVS2S_HDP1S_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr ssattr,
    float ssswgt,
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr cvattr,
    float cvswgt,
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr h1attr,
    float hdp1swgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    Func roundfunc,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores );

#include "CuBatchScoreMatrix_ssss_cvs2s_hdp1s.cu"

#endif//__CuBatchScoreMatrix_ssss_cvs2s_hdp1s_h__
