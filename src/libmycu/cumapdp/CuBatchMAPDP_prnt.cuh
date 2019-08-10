/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_prnt_h__
#define __CuBatchMAPDP_prnt_h__

//parent device functions for executing MAP dynamic programming 
//
__global__ void ExecMAPDP_Parent_DynPlm(
    const float alnextreg,
    const float alnextgapreg,
    //
    uint nblkdiags,
    dim3 nblcks, dim3 nthrds,
    //
    uint lastydiagnum,
    float logevthld,
    bool modscoresinuse,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    uint dblen2,
    uint dbxpad2,
    const CUBSM_TYPE* __restrict__ scores,
    const CUBSM_TYPE* __restrict__ modscores,
    const float* __restrict__ dp2alndatbuffers,
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpfwdprobs,
//float* __restrict__ tmpbwdprobs,
    float* __restrict__ tmpss2datbuffers,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata
);

#endif//__CuBatchMAPDP_prnt_h__
