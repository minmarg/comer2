/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_score_probs.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming;
// NOTE: Version for thread block of one warp!
// ndb1pros, number of profiles in the first profile data buffer db1;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// nqyposs, number of query positions to process;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// dbxpad, number of padded positions for memory alignment;
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
//

// -------------------------------------------------------------------------
// CalcScoreProbs: device code for calculating score probabilities using 
// shared memory for caching scores;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers that at this stage contain the 
// information of profiles passed;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// 
__global__ void CalcScoreProbs(
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint /*querposoffset*/, uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpss2datbuffers )
{
    __shared__ uint orgpronrCache;//original profile number/index
    __shared__ LNTYPE dbprodstCache;//distance in positions to the original db profile
    __shared__ INTYPE dbprolenCache;//length of the profile
    //__shared__ CUBSM_TYPE //cache of scores
    //        scoreCache[CUSS_2DCACHE_DIM][CUSS_2DCACHE_DIM];
    //__shared__ float probsCache[CUSS_N_DIFF_TOTAL_SCORES];//cache of score probabilities

    const int dblen = ndb1poss + ndbCposs + dbxpad;

    if( threadIdx.y == 0 && threadIdx.x == 0 ) {
        // blockIdx.z is the profile serial number in pass 2;
        orgpronrCache = tmpdpdiagbuffers[dblen+blockIdx.z];
        uint dbfldsndx;
        if( orgpronrCache < ndb1pros ) { orgpronrCache += ndb1prosOmtd;
                    dbfldsndx = pmv2DTotFlds;
        } else {    orgpronrCache += ndbCprosOmtd - ndb1pros;//jump to section ndbCposs
                    //bdb1posoffset = bdbCposoffset;
                    dbfldsndx = TIMES2(pmv2DTotFlds);
        }
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[orgpronrCache];
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[orgpronrCache];
    }

    //thread's linear number in the block
    //uint tidlin = blockIdx.y * blockDim.x + threadIdx.x;

    //initialize probabilities in SMEM
    //if( tidlin < CUSS_N_DIFF_TOTAL_SCORES )
    //    probsCache[tidlin] = 0.0f;

    __syncthreads();

    //if( dbprolenCache <= blockIdx.x * blockDim.x )
    //    //block does not participate: out of profile boundaries
    //    return;

    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    int dbpos;
    CUBSM_TYPE score;
    int intscore;
    int intscorescaled;
    //whether a thread is within boundaries
    bool inmtx = x < dbprolenCache && y < nqyposs;

    if( !inmtx )
        return;

    //if( inmtx )
    //{
        //db profile position in the dc_pm2dvfields_ buffers
        dbpos = x + dbprodstCache;
        //dbpos is now the x position in the score matrix plus the 
        //offset determined by thread id:
        dbpos = (orgpronrCache < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;

        //coalescent read
        score = scores[y*dblen + dbpos];
        score = VerifyScore( score );
        intscore = (int)rintf( score );
        intscorescaled = (int)rintf( score * CUSS_SCALE_FACTOR );

        //update counts in SMEM so that access to GMEM is minimized
//         probsCache[intscore-CUSS_SCORE_MIN] += 1.0f;
//         probsCache[CUSS_N_DIFF_SCORES + intscorescaled-CUSS_SCALED_SCORE_MIN] += 1.0f;
        //atomicAdd( probsCache + intscore-CUSS_SCORE_MIN, 1.0f );
        //atomicAdd( probsCache + CUSS_N_DIFF_SCORES + intscorescaled-CUSS_SCALED_SCORE_MIN, 1.0f );
    //}

    //__syncthreads();

    //if( !inmtx )
    //    return;

    atomicAdd( 
        (float*)(tmpss2datbuffers + (uint)(blockIdx.z*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES +
            (intscore - CUSS_SCORE_MIN))),
        1.0f );
    atomicAdd( 
        (float*)(tmpss2datbuffers + (uint)(blockIdx.z*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES +
            CUSS_N_DIFF_SCORES + (intscorescaled - CUSS_SCALED_SCORE_MIN))),
        1.0f );
    return;

    //update counts in GMEM in a coalescent way
    //if( tidlin < CUSS_N_DIFF_TOTAL_SCORES )
    //    atomicAdd( 
    //        tmpss2datbuffers + blockIdx.z*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+tidlin,
    //        probsCache[tidlin] );
}
