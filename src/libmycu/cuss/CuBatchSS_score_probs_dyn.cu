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
#include "CuBatchSS_minmax.cuh"
#include "CuBatchSS_expected.cuh"
#include "CuBatchSS_lambda.cuh"
#include "CuBatchSS_entropy.cuh"
#include "CuBatchSS_K.cuh"
#include "CuBatchSS_score_probs.cuh"
#include "CuBatchSS_score_probs_dyn.cuh"

#include "libmycu/cusf/CuSF_betainc.cuh"
#include "libmycu/cusf/CuSF_gammainc.cuh"

// #define CUSS_SCOREPROBS_DYNPLM_TESTPRINT 1286 //1768 //1304 //250 //488

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for calculating score probabilities;
// NOTE: Version for PARENT thread block of one thread!
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
// InitScoreProbsDynPlmProfile: device code for initializing score 
// probabilities for one given profile;
// pronr2, profile serial number in phase 2;
// NOTE: memory pointers should be aligned!
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// 
__global__ void InitScoreProbsDynPlmProfile(
    uint pronr2,
    float* __restrict__ tmpss2datbuffers )
{
    tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES + threadIdx.x] = 0.0f;
}

// -------------------------------------------------------------------------
// CalcScoreProbsDynPlmProfile: device code for calculating score 
// probabilities for one given profile;
// pronr2, profile serial number in phase 2;
// db1prolen, length of the (given) profile to be processed;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// NOTE: due to atomic operations with GMEM, this function takes the 
// most of the time to calculate statistical parameters;
// 
__global__ void CalcScoreProbsDynPlmProfile(
    uint pronr2,
    uint nqyposs, uint db1prolen, int dbpos, int dblen,
    CUBSM_TYPE* __restrict__ scores, 
    float* __restrict__ tmpss2datbuffers )
{
    //NOTE: using SMEM is inefficient here for multiple atomic 
    //operations and synchronization;

    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    //if outside the boundaries
    if( db1prolen <= x || nqyposs <= y )
        return;

    //coalescent read
    CUBSM_TYPE score = scores[y * dblen + dbpos + x];
    score = VerifyScore( score );

    //narrow range of scores induces delays for atomic operations;
    //calculate it later or along with statistical parameters;
//     int intscore = (int)rintf( score );
//     atomicAdd( 
//        tmpss2datbuffers + pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES +
//            intscore-CUSS_SCORE_MIN,
//        1.0f );

    int intscorescaled = (int)rintf( score * CUSS_SCALE_FACTOR );
    atomicAdd( 
        tmpss2datbuffers + pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES +
            CUSS_N_DIFF_SCORES + intscorescaled-CUSS_SCALED_SCORE_MIN,
        1.0f );
}

// -------------------------------------------------------------------------
// NormalizeScoreProbsDynPlmProfile: device code for calculating 
// probabilities of unscaled scores and normalizing probabilities; 
// pronr2, profile serial number in phase 2;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// NOTE: keep # registers <= 32
// 
__global__ void NormalizeScoreProbsDynPlmProfile(
    uint pronr2,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    float* __restrict__ tmpss2datbuffers )
{
    //cache for probabilities; add 1 to eliminate bank conflicts when accessed
    __shared__ float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1];
    __shared__ float unscldprobsCache[CUSS_N_DIFF_SCORES];
    __shared__ float tmpCache[TIMES2(CUSS_N_DIFF_SCORES)];
    __shared__ float totalsumshared;
    //threadIdx.y corresponds to an unormalized score;
    //threadIdx.x is a scaled score for the threadIdx.y score;
    //value of unscaled score:
    int unscaled = (int)threadIdx.y+CUSS_SCORE_MIN;
    //range of scaled scores for threadIdx.y:
    int leftndx = unscaled * CUSS_SCALE_FACTOR - (CUSS_SCALE_FACTOR>>1) - CUSS_SCALED_SCORE_MIN;
    int rghtndx = unscaled * CUSS_SCALE_FACTOR + (CUSS_SCALE_FACTOR>>1) - CUSS_SCALED_SCORE_MIN;

    //adjust the interval of the scaled score
    if( unscaled <= 0 ) {
        leftndx++;
        if( leftndx < 0 ) leftndx = 0;
    }
    if( 0 <= unscaled ) {
        rghtndx--;
        if( CUSS_N_DIFF_SCALED_SCORES-1 < rghtndx ) rghtndx = CUSS_N_DIFF_SCALED_SCORES-1;
    }

    //cache scaled probabilities (counts)
    if( threadIdx.x + leftndx <= rghtndx )
        scaledprobsCache[threadIdx.y][threadIdx.x] = 
            tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES + 
                CUSS_N_DIFF_SCORES + leftndx + threadIdx.x];
    else
        scaledprobsCache[threadIdx.y][threadIdx.x] = 0.0f;

    __syncthreads();

    float sumcnt = scaledprobsCache[threadIdx.y][threadIdx.x];//sum of counts for a given unscaled score
    float sumtot = 0.0f;//total sum of counts

    //warp reduce (sum) the counts of the scaled scores
    sumcnt += __shfl_xor_sync(0xffffffff, sumcnt, 16);
    sumcnt += __shfl_xor_sync(0xffffffff, sumcnt, 8);
    sumcnt += __shfl_xor_sync(0xffffffff, sumcnt, 4);
    sumcnt += __shfl_xor_sync(0xffffffff, sumcnt, 2);
    sumcnt += __shfl_xor_sync(0xffffffff, sumcnt, 1);

    //write the counts to SMEM
    if( threadIdx.x == 0 )
        unscldprobsCache[threadIdx.y] = sumcnt;
    __syncthreads();

    //reduce the partial counts of unscaled scores to get the total count
    if( threadIdx.y == 0 ) {
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
            sumtot = unscldprobsCache[threadIdx.x];
        sumtot += __shfl_xor_sync(0xffffffff, sumtot, 16);
        sumtot += __shfl_xor_sync(0xffffffff, sumtot, 8);
        sumtot += __shfl_xor_sync(0xffffffff, sumtot, 4);
        sumtot += __shfl_xor_sync(0xffffffff, sumtot, 2);
        sumtot += __shfl_xor_sync(0xffffffff, sumtot, 1);
        //write the (normalized) probabilities of unscaled scores
        if( threadIdx.x < CUSS_N_DIFF_SCORES )
        {
            unscldprobsCache[threadIdx.x] =
                __fdividef(unscldprobsCache[threadIdx.x], sumtot);
            //branching has a greater negative effect than additional writes
            //if( unscldprobsCache[threadIdx.x])
                tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES + threadIdx.x] = 
                    unscldprobsCache[threadIdx.x];
            if( threadIdx.x == 0 )
                totalsumshared = sumtot;
        }
    }

    __syncthreads();

// #ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
//     if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT) {
//         printf("Intervals: y=%u x=%u %d %d:  \n",threadIdx.y,threadIdx.x,leftndx,rghtndx);
//         if(threadIdx.y==0 && threadIdx.x==0) {
//             printf("\nCounts:\n");
//             for(int i=0;i<CUSS_N_DIFF_SCALED_SCORES;i++) printf(" %.1f,",
//                 tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+CUSS_N_DIFF_SCORES+i]);
//             printf("\n\n");
//             for(int i=0;i<CUSS_N_DIFF_SCORES;i++) {
//                 for(int j=0;j<CUSS_SCALE_FACTOR;j++) printf(" %.1f,",scaledprobsCache[i][j]);
//                 printf("\n");
//             }
//             printf("\n\n");
//             printf("Sums (total= %.1f %.1f):\n",sumtot,totalsumshared);
//             for(int i=0;i<CUSS_N_DIFF_SCORES;i++) printf(" %.1f,", unscldprobsCache[i]);
//             printf("\n\n");
//         }
//     }
//     __syncthreads();
// #endif

    //write the (normalized) probabilities of scaled scores;
    //NOTE: perform all calculations with scaled scores in cache;
    //if( threadIdx.x + leftndx <= rghtndx /*&&
    //    //branching has a greater negative effect than additional writes
    //    scaledprobsCache[threadIdx.y][threadIdx.x] */)
    //    tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES + 
    //        CUSS_N_DIFF_SCORES + leftndx + threadIdx.x] =
    //            __fdividef(scaledprobsCache[threadIdx.y][threadIdx.x], totalsumshared);

    //normalize probabilities only in SMEM
    scaledprobsCache[threadIdx.y][threadIdx.x] = 
        __fdividef(scaledprobsCache[threadIdx.y][threadIdx.x], totalsumshared);

    __syncthreads();

    //reuse register sumtot; it is now the expected score
    sumtot = CalcExpectedScoreProfile( tmpCache, scaledprobsCache, leftndx );
    GetScoreMinMaxProfile( tmpCache, scaledprobsCache, leftndx );

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT && threadIdx.y==0 && threadIdx.x==0) {
        printf("\n\n MIN= %.1f MAX= %.1f E= %.4f (%.4f)\n\n", 
            tmpCache[0], tmpCache[CUSS_N_DIFF_SCORES], 
            sumtot, sumtot/CUSS_SCALE_FACTOR );
    }
#endif

    float minscore = tmpCache[0];
    float maxscore = tmpCache[CUSS_N_DIFF_SCORES];

    float lambda = CalcLambdaProfile( 
        tmpCache,
        scaledprobsCache,
        leftndx,
        sumtot/*expected*/,
        /*tmpCache[0]*/minscore,
        /*tmpCache[CUSS_N_DIFF_SCORES]*/maxscore
    );

    float relent = CalcEntropyProfile(
        tmpCache,
        scaledprobsCache,
        leftndx,
        /*tmpCache[0]*/minscore,
        /*tmpCache[CUSS_N_DIFF_SCORES]*/maxscore,
        lambda
    );
//     float relent = CalcEntropyProfile(
//         tmpCache,
//         unscldprobsCache,
//         /*tmpCache[0]*/minscore,
//         /*tmpCache[CUSS_N_DIFF_SCORES]*/maxscore,
//         lambda*CUSS_SCALE_FACTOR
//     );

    lambda *= CUSS_SCALE_FACTOR;
    sumtot = __fdividef(sumtot, CUSS_SCALE_FACTOR);

    if( threadIdx.y == 0 )
        scaledprobsCache[0-CUSS_SCORE_MIN][threadIdx.x] = 1.0f;

    //NOTE: no need for __syncthreads(), which is the last operation in the next call
    GetScoreMinMaxProfile( tmpCache, unscldprobsCache );

    minscore = tmpCache[0];
    maxscore = tmpCache[CUSS_N_DIFF_SCORES];

    float K = CalcKProfile( 
        tmpCache,
        unscldprobsCache,
        scaledprobsCache,
        minscore,
        maxscore,
        maxscore-minscore+1.0f/*unscldscorerange*/,
        leftndx,
        rghtndx,
        sumtot/*expected*/,
        lambda,
        relent
    );

    if( threadIdx.y == 0 && threadIdx.x == 0)
       tmpss2datbuffers[pronr2 * CUSS_ALIGNED_N_DIFF_TOTAL_SCORES + 
           CUSS_N_DIFF_SCORES] = K;

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT && threadIdx.y==0 && threadIdx.x==0) {
        printf("\n MIN= %.1f MAX= %.1f Expected= %.4f LMB= %.4f H= %.4f K= %.4f\n\n", 
            minscore, maxscore, sumtot, lambda, relent, K );
    }
#endif

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT && threadIdx.y==0 && threadIdx.x==0) {
        for(int i=0;i<CUSS_N_DIFF_SCORES;i++) printf(" (%d)=%.4f,",
            i+CUSS_SCORE_MIN,
            tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+i]);
        printf("\n\n");
        for(int i=0;i<CUSS_N_DIFF_SCALED_SCORES;i++) printf(" (%d)=%.4f,",
            i+CUSS_SCALED_SCORE_MIN,
            tmpss2datbuffers[pronr2*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+CUSS_N_DIFF_SCORES+i]);
        printf("\n\n");
        //
//         for(int ppp=0;ppp<2000;ppp++){
//             float tmpres;
//             int locpsts=ppp;//(int)rintf( __fdividef((float)ppp*100000.0f, __powf((float)(23*52),1.5f)) );
//             tmpres = cusf_lbetaincf( locpsts/1000.0f, 0.1f, 0.757213f );
//             printf(" > psts= %d log(ib)= %g\n",locpsts,tmpres);
//         }
//         for(float ppp=0.0f;ppp<1000.0f;ppp+=1.f){
//             float tmpres;
//             tmpres = cusf_lgammainc_Qf( 1.346206f, ppp );
//             printf(" > p= %.1g log(ig)= %g\n",ppp,tmpres);
//         }
    }
#endif
}

// =========================================================================
// CalcScoreProbsDynPlm: device code for calculating score probabilities 
// using dynamic parallelism;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers that at this stage contain the 
// information of profiles passed;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// 
__global__ void CalcScoreProbsDynPlm(
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint /*querposoffset*/, uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpss2datbuffers )
{
    //NOTE: assuming only ONE thread in a block
    //
//     cudaStream_t streamss;
//     MYCUDACHECK( cudaStreamCreateWithFlags(&streamss, cudaStreamNonBlocking));

    //initialize
    InitScoreProbsDynPlmProfile<<<1,CUSS_ALIGNED_N_DIFF_TOTAL_SCORES/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
        tmpss2datbuffers
    );
    MYCUDACHECKLAST;


    uint serpronrSM;//profile serial number in the score matrix (current processing)
    uint orgpronrCache;//original profile number/index
    LNTYPE dbprodstCache;//distance in positions to the original db profile
    INTYPE dbprolenCache;//length of the profile

    const int dblen = ndb1poss + ndbCposs + dbxpad;

    // blockIdx.x is the profile serial number in phase 2;
    serpronrSM = tmpdpdiagbuffers[dblen+blockIdx.x];
    orgpronrCache = serpronrSM;
    uint dbfldsndx;
    if( orgpronrCache < ndb1pros ) { orgpronrCache += ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    orgpronrCache += ndbCprosOmtd - ndb1pros;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[orgpronrCache];
    dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[orgpronrCache];

    //db profile position in the dc_pm2dvfields_ buffers
    int dbpos = dbprodstCache;
    //dbpos is now the beginning x position in the score matrix:
    dbpos = (serpronrSM < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;


    //execution configuration for the profile this parent thread is responsible for:
    dim3 nthrds(CUSS_2DCACHE_DIM,CUSS_2DCACHE_DIM,1);
    dim3 nblcks((dbprolenCache+nthrds.x-1)/nthrds.x,
                (nqyposs+nthrds.y-1)/nthrds.y,
                1);

    //calculate counts
    CalcScoreProbsDynPlmProfile<<<nblcks,nthrds/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
        (uint)nqyposs, (uint)dbprolenCache, dbpos, dblen,
        scores,
        tmpss2datbuffers
    );
    MYCUDACHECKLAST;

    //normalize probabilities
    nthrds = dim3(CUSS_SCALE_FACTOR,CUSS_N_DIFF_SCORES,1);
    NormalizeScoreProbsDynPlmProfile<<<1,nthrds/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
        orgpronrCache,
#endif
        tmpss2datbuffers
    );
    MYCUDACHECKLAST;

//     MYCUDACHECK( cudaStreamDestroy( streamss ));
//     MYCUDACHECKLAST;
}
