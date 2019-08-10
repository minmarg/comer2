/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "CuBatchMAPDP_fwd.cuh"
#include "CuBatchMAPDP_bwd.cuh"
#include "CuBatchMAPDP_map_btck.cuh"

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
    char* __restrict__ btckdata )
{
    //launch blocks along block diagonals so that d=x+y-1
    for( uint d = 0; d < nblkdiags; d++) {
        ExecMAPDP_Fwd_Unroll32x<<<nblcks,nthrds>>>(
            d,
            lastydiagnum,
            logevthld,
            modscoresinuse,
            (uint)ndb1pros,
            (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            (uint)dblen2,
            dbxpad2,
            scores,//[in]
            modscores,//[in]
            dp2alndatbuffers,//[in]
            tmpdpdiagbuffers, 
            tmpdpbotbuffer,
            tmpfwdprobs
        );
        MYCUDACHECKLAST;
    }
    //launch blocks along block diagonals so that d=x+y-1
    for( int d = nblkdiags-1; 0 <= d; d--) {
        ExecMAPDP_Bwd_Unroll32x<<<nblcks,nthrds>>>( 
            d,
            lastydiagnum,
            logevthld,
            modscoresinuse,
            (uint)ndb1pros,
            (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            (uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            (uint)dblen2,
            dbxpad2,
            scores,//[in]
            modscores,//[in]
            dp2alndatbuffers,//[in]
            tmpdpdiagbuffers, 
            tmpdpbotbuffer,
            tmpfwdprobs,
//tmpbwdprobs,
            tmpss2datbuffers
        );
        MYCUDACHECKLAST;
    }
    //launch blocks along block diagonals so that d=x+y-1
    for( uint d = 0; d < nblkdiags; d++) {
        ExecMAPDP_MAP_Btck_Unroll32x<<<nblcks,nthrds>>>( 
            alnextreg,
            alnextgapreg,
            d,
            lastydiagnum,
            logevthld,
            (uint)ndb1pros,
            (uint)querprosOmtd, (uint)ndb1prosOmtd, (uint)ndbCprosOmtd,
            (uint)nqyposs, //(uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            //(uint)querposoffset, (uint)bdb1posoffset, (uint)bdbCposoffset,
            (uint)dblen2,
            dbxpad2,
            tmpfwdprobs,//[in]
            dp2alndatbuffers,//[in]
            tmpdpdiagbuffers, 
            tmpdpbotbuffer,
            tmpss2datbuffers,
            maxscoordsbuf,
            btckdata
        );
        MYCUDACHECKLAST;
    }
}
