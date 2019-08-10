/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_bwd_h__
#define __CuBatchMAPDP_bwd_h__

#include "extsp/psl.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init.cuh"
#include "CuBatchMAPDP_com.h"
#include "CuBatchMAPDP_fwd.cuh"

//device functions for executing dynamic programming to find backward 
// probabilities

__global__ void ExecMAPDP_Bwd_Unroll32x(
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
    float* __restrict__ tmpfwdprobs,
//float* __restrict__ tmpbwdprobs,
    float* __restrict__ tmpss2datbuffers
);

// =========================================================================
// -------------------------------------------------------------------------
// MAPDP_Bwd_LocCacheTrnProbs: cache transition probabilities to SMEM at 
// position pos+1 (current position)
//
__device__ __forceinline__
void MAPDP_Bwd_LocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_D],
    int fldsndx,
    int pos )
{
    smtrpCache[dptrMMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMMp, pos+1  );
    smtrpCache[dptrMIc][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMIc, pos+1  );
    smtrpCache[dptrMDp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrMDp, pos+1  );
    smtrpCache[dptrIMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrIMp, pos+1  );
    smtrpCache[dptrIIc][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrIIc, pos+1  );
    smtrpCache[dptrDMp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrDMp, pos+1  );
    smtrpCache[dptrDDp][threadIdx.x] = MAPDPLocGetTrnProb( fldsndx, dptrDDp, pos+1  );
}

// MAPDPLocCacheTrnProbs: overloaded for caching CUMAPDP_2DCACHE_DIMx2 
// positions
//
__device__ __forceinline__
void MAPDP_Bwd_LocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUMAPDP_2DCACHE_DIM_DpX],
    int ndxoff,
    int fldsndx,
    int pos )
{
    int ndx = ndxoff+threadIdx.x;
    smtrpCache[dptrMMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMMp, pos+1  );
    smtrpCache[dptrMIc][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMIc, pos+1  );
    smtrpCache[dptrMDp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrMDp, pos+1  );
    smtrpCache[dptrIMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrIMp, pos+1  );
    smtrpCache[dptrIIc][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrIIc, pos+1  );
    smtrpCache[dptrDMp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrDMp, pos+1  );
    smtrpCache[dptrDDp][ndx] = MAPDPLocGetTrnProb( fldsndx, dptrDDp, pos+1  );
}

// =========================================================================
// -------------------------------------------------------------------------
// MAPDP_ReduceMaxLogFwdProbability: reduce for the maximum of log forward 
// probability;
// pronr2, db profile serial number in phase 2;
// dbprolen, db profile length;
// dblen2, phase-2 db length in positions;
// dbbegpos2, beginning of the db profile in phase 2;
// tmpdpdiagbuffers, buffers containing maximum log probabilities to be 
// reduced;
// tmpss2datbuffers, temporary buffer to save the reduced value;
//
__device__ __forceinline__
float MAPDP_ReduceMaxLogFwdProbability(
    uint pronr2,
    int dbprolen,
    int dblen2,
    int dbbegpos2,
    const float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpss2datbuffers )
{
    float maxlogprob = 0.0f;
    int dbpos2off = 
        dpdssDiagM * nTDPDiagScoreSubsections * dblen2 + 
        dbbegpos2 + dbprolen-1-threadIdx.x;

    if( threadIdx.x < dbprolen )
        maxlogprob = tmpdpdiagbuffers[dbpos2off];
    if( threadIdx.x+CUMAPDP_2DCACHE_DIM_D < dbprolen )
        maxlogprob = 
                 myhdmax( maxlogprob, tmpdpdiagbuffers[dbpos2off-CUMAPDP_2DCACHE_DIM_D]);
    //warp reduce
    maxlogprob = myhdmax( maxlogprob, __shfl_down_sync(0xffffffff, maxlogprob, 16));
    maxlogprob = myhdmax( maxlogprob, __shfl_down_sync(0xffffffff, maxlogprob, 8));
    maxlogprob = myhdmax( maxlogprob, __shfl_down_sync(0xffffffff, maxlogprob, 4));
    maxlogprob = myhdmax( maxlogprob, __shfl_down_sync(0xffffffff, maxlogprob, 2));
    maxlogprob = myhdmax( maxlogprob, __shfl_down_sync(0xffffffff, maxlogprob, 1));

    //WRITE the maximum log probability in the second cell of tmpss2datbuffers for the 
    // given profile:
    if( threadIdx.x == 0 )
        tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+1] = maxlogprob;

    maxlogprob = __shfl_sync(0xffffffff, maxlogprob, 0);

    return maxlogprob;
}

// -------------------------------------------------------------------------
// MAPDP_ReadMaxLogFwdProbability: read the maximum of log forward 
// probability, which has been previously reduced;
// pronr2, db profile serial number in phase 2;
// tmpss2datbuffers, temporary buffer to read the reduced value from;
//
__device__ __forceinline__
float MAPDP_ReadMaxLogFwdProbability( 
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    float maxlogprob = 0.0f;

    //READ the maximum log probability from the second cell of tmpss2datbuffers of the 
    // given profile:
    if( threadIdx.x == 0 )
        maxlogprob = tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+1];

    maxlogprob = __shfl_sync(0xffffffff, maxlogprob, 0);

    return maxlogprob;
}

// -------------------------------------------------------------------------
// MAPDP_GetMaxLogFwdProbability: get the maximum of log forward 
// probability, which has been previously reduced;
// firstprocdiagonal, flag of the first block diagonal being processed;
// pronr2, db profile serial number in phase 2;
// tmpss2datbuffers, temporary buffer for the reduced value;
//
__device__ __forceinline__
float MAPDP_GetMaxLogFwdProbability(
    bool firstprocdiagonal,
    uint pronr2,
    int dbprolen,
    int dblen2,
    int dbbegpos2,
    const float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpss2datbuffers )
{
    return     
        firstprocdiagonal
        ?
        MAPDP_ReduceMaxLogFwdProbability(
            pronr2,
            dbprolen,
            dblen2,
            dbbegpos2,
            tmpdpdiagbuffers,
            tmpss2datbuffers )
        :
        MAPDP_ReadMaxLogFwdProbability( pronr2, tmpss2datbuffers );
}

// -------------------------------------------------------------------------
// MAPDP_UpdateForTotalProbability: update the sum of exponentials for 
// required computing log of the total probability
//
__device__ __forceinline__
void MAPDP_UpdateForTotalProbability(
    uint pronr2,
    int nqyposs,
    int dbprolen,
    int dblen2,
    int x, int y,
    int dbpos2,
    const float maxlogprob,
    const float* __restrict__ tmpfwdprobs,
    float* __restrict__ tmpss2datbuffers )
{
    float nprob;//normalized probability
    float sum = 0.0f;
    //
    //cache scores over the oblique diagonal block
    #pragma unroll 4
    for( int i = 0; i < CUMAPDP_2DCACHE_DIM_D; i++ ) {
        //going upwards
        nprob = 0.0f;
        //NOTE: cache a score block SHIFTED right and down
        if( 0 <= y-i && y-i < nqyposs && 0 <= x+i && x+i < dbprolen )
            //take account of the starting position of line i of the oblq. 
            // diagonal block: 
            nprob = MAPDPexpf(tmpfwdprobs[(y-i) * dblen2 + i  +  dbpos2] - maxlogprob);
        //reduce for the total sum over the block
        nprob += __shfl_down_sync(0xffffffff, nprob, 16);
        nprob += __shfl_down_sync(0xffffffff, nprob, 8);
        nprob += __shfl_down_sync(0xffffffff, nprob, 4);
        nprob += __shfl_down_sync(0xffffffff, nprob, 2);
        nprob += __shfl_down_sync(0xffffffff, nprob, 1);
        if( threadIdx.x == 0 )
            sum += nprob;
    }
    //atomic ADD of the sum of this block to the total sum to be converted to 
    // log total probability before using it
    if( threadIdx.x == 0 ) {
        //NOTE: tmpss2datbuffers has been initialized in CalcStatisticsDynPlm 
        // dynamic parallelism when calculating SS before beginning this stage of 
        // MAP dynamic programming
        atomicAdd( tmpss2datbuffers + pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES, sum );
    }
}



// =========================================================================
// -------------------------------------------------------------------------
// MAPDP_FinalizeLogTotalFwdProbability: finalize log of the total forward 
// probability by taking the logarithm of the calculated sum and adding the 
// max log term;
// save the result for future calls and return it;
// pronr2, db profile serial number in phase 2;
// tmpss2datbuffers, temporary buffer for the log total probabilities;
//
__device__ __forceinline__
float MAPDP_FinalizeLogTotalFwdProbability( 
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    float logtotprob = 0.0f;

    if( threadIdx.x == 0 ) {
        //READ sum of exponentials:
        logtotprob = tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES];
        //READ max log probability from the second cell of tmpss2datbuffers:
        float maxlogprob = tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+1];
        //finalize sum of exponentials to get log total probability
        logtotprob = MAPDPlogf(logtotprob) + maxlogprob;
        //WRITE back the result to the corresponding location
        tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES] = logtotprob;
    }

    logtotprob = __shfl_sync(0xffffffff, logtotprob, 0);

    return logtotprob;
}

// -------------------------------------------------------------------------
// MAPDP_ReadLogTotalFwdProbability: read log total forward probability;
// pronr2, db profile serial number in phase 2;
// tmpss2datbuffers, temporary buffer for the log total probabilities;
//
__device__ __forceinline__
float MAPDP_ReadLogTotalFwdProbability( 
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    float logtotprob = 0.0f;

    //READ the logarithm of the total forward probability for the given profile:
    if( threadIdx.x == 0 )
        logtotprob = tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES];

    logtotprob = __shfl_sync(0xffffffff, logtotprob, 0);

    return logtotprob;
}

// -------------------------------------------------------------------------
// MAPDP_GetTotalFwdProbability: get log of the total forward probability;
// if it is the first call (firstprocdiagonal), finalize and save the 
// probability for future calls;
// firstprocdiagonal, flag of the first block diagonal being processed;
// pronr2, db profile serial number in phase 2;
// tmpss2datbuffers, temporary buffer for the log total probabilities;
//
__device__ __forceinline__
float MAPDP_GetLogTotalFwdProbability( 
    bool firstprocdiagonal,
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    return     
        firstprocdiagonal
        ?
        MAPDP_FinalizeLogTotalFwdProbability( pronr2,tmpss2datbuffers )
        :
        MAPDP_ReadLogTotalFwdProbability( pronr2, tmpss2datbuffers );
}

// =========================================================================

#endif//__CuBatchMAPDP_bwd_h__
