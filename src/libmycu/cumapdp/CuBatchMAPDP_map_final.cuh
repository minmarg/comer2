/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_map_final_h__
#define __CuBatchMAPDP_map_final_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchMAPDP_com.h"

//device functions for finalizing MAP dynamic programming

__global__ void FinalizeMAPDP_MAP(
    bool printsss,
    float logevthld,
    uint ndb1pros,
    uint /*querprosOmtd*/, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    uint dblen2,
    uint dbxpad2,
    uint dbalnlen2,
    const CUBSM_TYPE* __restrict__ scores, 
    const float* __restrict__ tmpdpdiagbuffers,
    const uint* __restrict__ maxscoordsbuf,
    const char* __restrict__ btckdata,
    float* __restrict__ dp2alndatbuffers,
    char* __restrict__ outalns
);

// =========================================================================
// -------------------------------------------------------------------------

#endif//__CuBatchMAPDP_map_final_h__
