/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ssss_cvs2s_h__
#define __CuBatchScoreMatrix_ssss_cvs2s_h__

#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchScoreMatrix_com.h"

//device functions for computing SSS and CVS scores

//kernel declaration
__global__ void CalcSM_SSSS_CVS2S_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr ssattr,
    float ssswgt,
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr cvattr,
    float cvswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores, 
    CUBSM_TYPE* __restrict__ outmodscores );

#endif//__CuBatchScoreMatrix_ssss_cvs2s_h__
