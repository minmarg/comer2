/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_fwd_h__
#define __CuBatchMAPDP_fwd_h__

#include "extsp/psl.h"
#include "libmycu/cucom/myfastapproxmath.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init.cuh"
#include "CuBatchMAPDP_com.h"

//device functions for executing dynamic programming to find forward 
// probabilities

__global__ void ExecMAPDP_Fwd_Unroll32x(
    uint blkdiagnum,
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
    float* __restrict__ tmpfwdprobs
);

// =========================================================================

// =========================================================================
// -------------------------------------------------------------------------
__device__ __forceinline__
CUBDP_TYPE MAPDPLocGetProb0()
{
    return CUBDP_Q(0);
}

__device__ __forceinline__
CUBDP_TYPE MAPDPLocGetLogProb0()
{
    //return DPLocGetLog0();
    return CUBDP_Q(-99);
}

__device__ __forceinline__
CUBDP_TYPE MAPDPLocGetInitTrnProbValue()
{
    return MAPDPLocGetLogProb0();
}

// -------------------------------------------------------------------------
// MAPDPLocinitTrnProbs: initialize transition probabilities
//
__device__ __forceinline__
void MAPDPLocInitTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_D])
{
    smtrpCache[dptrMMp][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrMIc][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrMDp][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrIMp][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrIIc][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrDMp][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrDDp][threadIdx.x] = MAPDPLocGetInitTrnProbValue();
}

__device__ __forceinline__
void MAPDPLocInitTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_DpX],
    int ndxoff )
{
    smtrpCache[dptrMMp][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrMIc][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrMDp][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrIMp][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrIIc][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrDMp][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
    smtrpCache[dptrDDp][ndxoff+threadIdx.x] = MAPDPLocGetInitTrnProbValue();
}

// -------------------------------------------------------------------------
// MAPDPLocGetTrnProb: read transition probability
__device__ __forceinline__
FPTYPE MAPDPLocGetTrnProb( int fldsndx, int trn, int pos )
{
    return
        ((FPTYPE*)(dc_pm2dvfields_[fldsndx+ptr2DTrnPrbs/*ptr2DTrnPrbsExp*/+trn]))[pos];
}

// -------------------------------------------------------------------------
// MAPDPLocCacheTrnProbs: cache transition probabilities to SMEM at position 
// pos
//
__device__ __forceinline__
void MAPDPLocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_D],
    int fldsndx,
    int pos )
{
    smtrpCache[dptrMMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMMp, pos    );
    smtrpCache[dptrMIc][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMIc, pos+1  );
    smtrpCache[dptrMDp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMDp, pos    );
    smtrpCache[dptrIMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrIMp, pos    );
    smtrpCache[dptrIIc][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrIIc, pos+1  );
    smtrpCache[dptrDMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrDMp, pos    );
    smtrpCache[dptrDDp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrDDp, pos    );
}

// MAPDPLocCacheTrnProbs: overloaded for caching CUMAPDP_2DCACHE_DIMx2 
// positions
//
__device__ __forceinline__
void MAPDPLocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_DpX],
    int ndxoff,
    int fldsndx,
    int pos )
{
    int ndx = ndxoff+threadIdx.x;
    smtrpCache[dptrMMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMMp, pos    );
    smtrpCache[dptrMIc][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMIc, pos+1  );
    smtrpCache[dptrMDp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMDp, pos    );
    smtrpCache[dptrIMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrIMp, pos    );
    smtrpCache[dptrIIc][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrIIc, pos+1  );
    smtrpCache[dptrDMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrDMp, pos    );
    smtrpCache[dptrDDp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrDDp, pos    );
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MAPDPGetDynScoreParams: get dynamic score parameters based on the 
// expect of the alignment;
// epa, log of pair e-value;
// initialization is assumed to have been done;
//
__device__ __forceinline__
void MAPDPGetDynScoreParams( float epa, float* dyno, float* dmul )
{
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    *dyno = *dmul = 0.0f;
    epa += CUMAPDP_EPA_LOGFACTOR;
    if( SLC_LOG_SP_MIN < epa ) {
        *dyno = CUMAPDP_DYNS_OFFSET;
        if( epa < 0.0f ) {
            epa = __expf(epa);
            *dyno = CUMAPDP_DYNS_OFFSET * epa;
            *dmul = CUMAPDP_DYNS_MULTP * (1.0f - epa);
        }
    }
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MAPDPGetDynGapExtParams: get dynamic gap extension parameters 
// based on the expect of the alignment;
// epa, log of pair e-value;
// initialization is assumed to have been done;
//
__device__ __forceinline__
void MAPDPGetDynGapExtParams( float epa, float* dmul )
{
    *dmul = CUMAPDP_DYNGAPEXT_MULTP_LB;
    epa += CUMAPDP_EPAGAPEXT_LOGFACTOR;
    if( SLC_LOG_SP_MIN < epa ) {
        *dmul = CUMAPDP_DYNGAPEXT_MULTP_UB;
        if( epa < 0.0f ) {
            *dmul = CUMAPDP_DYNGAPEXT_MULTP_LB + 
                (CUMAPDP_DYNGAPEXT_MULTP_UB-CUMAPDP_DYNGAPEXT_MULTP_LB) *
                    (1.0f + 0.0113636f * epa);
        }
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
__device__ __forceinline__
float MAPDPfastapprox_expf( float arg )
{
    //NOTE: two comparisons are expensive when heavily accessed as here
    //return (-16.0f<arg)? (arg? myfastapexpf(arg): 1.0f): 0.0f;
//     return (-1.0f<arg)? 0.99f: 0.0f;//TEST DIFF in PERFORMANCE
    return (-16.0f<arg)? myfastapexpf(arg): 0.0f;
    //return (-64.0f<arg)? myfastapexpf(arg): 0.0f;
}
__device__ __forceinline__
float MAPDPregular_expf( float arg )
{
    return (-16.0f<arg)? __expf(arg): 0.0f;
    //return (-64.0f<arg)? (arg? __expf(arg): 1.0f): 0.0f;
}
// -------------------------------------------------------------------------
// MAPDPexpf: approximate exponential function for calculating accumulated 
// probabilities during dynamic programming
//
__device__ __forceinline__
float MAPDPexpf( float arg )
{
#ifdef CUMAPDP_USEFASTAPPROXMATH 
    return MAPDPfastapprox_expf(arg);
#else
    return MAPDPregular_expf(arg);
#endif
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
__device__ __forceinline__
float MAPDPfastapprox_logf( float arg )
{
//     return (0.3f<arg)? 0.0f: -99.0f;//TEST DIFF in PERFORMANCE
    return (3.8e-11f<arg)? myfastaplogf(arg): -99.0f;
    //return (1.6e-28f<arg)? myfastaplogf(arg): -99.0f;
}
__device__ __forceinline__
float MAPDPregular_logf( float arg )
{
    return (3.8e-11f<arg)? __logf(arg): -99.0f;
    //return (1.6e-28f<arg)? __logf(arg): -99.0f;
    //return (3.8e-11f<arg)? ((0.999999f<arg && arg<1.000001f)? 0.0f: __logf(arg)): -99.0f;
}
// -------------------------------------------------------------------------
// MAPDPlogf: approximate logarithmic function for calculating log of 
// accumulated probabilities during dynamic programming
//
__device__ __forceinline__
float MAPDPlogf( float arg )
{
#ifdef CUMAPDP_USEFASTAPPROXMATH 
    return MAPDPfastapprox_logf(arg);
#else
    return MAPDPregular_logf(arg);
#endif
}

#endif//__CuBatchMAPDP_fwd_h__
