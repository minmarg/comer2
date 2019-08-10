/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_Init_hdp1s_h__
#define __CuBatchScoreMatrix_Init_hdp1s_h__

#include "libmycu/cupro/SerializedScoresAttr.h"
#include "CuBatchScoreMatrix_com.h"

//device functions for computing initial and HDP scores at once

__global__ void CalcSM_Init_HDP1S_SMEMUnroll2x(
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr attr,
    float hdp1swgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores, 
    CUBSM_TYPE* __restrict__ outmodscores );

#endif//__CuBatchScoreMatrix_Init_hdp1s_h__
