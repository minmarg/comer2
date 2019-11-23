/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_aln_final_h__
#define __CuBatchDP_aln_final_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"

//device functions for finalizing dynamic programming to get output alignments;
//analogous to FinalizeMAPDP_MAP

__global__ void FinalizeDP_ALN(
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
    const CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    const uint* __restrict__ maxscoordsbuf,
    const char* __restrict__ btckdata,
    float* __restrict__ dp2alndatbuffers,
    char* __restrict__ outalns
);

// =========================================================================
// -------------------------------------------------------------------------

#endif//__CuBatchDP_aln_final_h__
