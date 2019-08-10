/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_map_btck_h__
#define __CuBatchMAPDP_map_btck_h__

#include "extsp/psl.h"
#include "CuBatchMAPDP_com.h"

//device functions for executing dynamic programming to find MAP substitution
// probabilities

__global__ void ExecMAPDP_MAP_Btck_Unroll32x(
    const float alnextreg,
    float alnextgapreg,
    uint blkdiagnum,
    uint lastydiagnum,
    float logevthld,
    uint ndb1pros,
    uint /*querprosOmtd*/, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, //uint ndb1poss, uint ndbCposs, uint dbxpad,
    //uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    uint dblen2,
    uint dbxpad2,
    const float* __restrict__ tmpfwdprobs,
    const float* __restrict__ dp2alndatbuffers,
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpss2datbuffers,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata
);

#endif//__CuBatchMAPDP_map_btck_h__
