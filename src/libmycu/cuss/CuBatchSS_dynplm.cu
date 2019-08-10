/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #define CUSS_SCOREPROBS_DYNPLM_TESTPRINT -1 //0 //1286 //1768 //1304 //250 //488

#include "liblib/mybase.h"

#include <math.h>

#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_score_probs.cuh"
#include "CuBatchSS_minmax.cuh"
#include "CuBatchSS_expected.cuh"
#include "CuBatchSS_lambda.cuh"
#include "CuBatchSS_entropy.cuh"
#include "CuBatchSS_K.cuh"
#include "CuBatchSS_NN.cuh"
#include "CuBatchSS_SS18.cuh"
#include "CuBatchSS_dynplm.cuh"

#include "libmycu/cusf/CuSF_betainc.cuh"
#include "libmycu/cusf/CuSF_gammainc.cuh"

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
            CUSS_N_DIFF_SCORES + (intscorescaled - CUSS_SCALED_SCORE_MIN),
        1.0f );
}

// -------------------------------------------------------------------------
// CalcStatParamsDynPlmProfile: device code for calculating 
// statistical parameters using the probabilities of scaled and unscaled 
// scores; the normalization of probabilities is performed inline; 
// pronr2, profile serial number in phase 2;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// NOTE: keep # registers <= 32
// 
__global__ void CalcStatParamsDynPlmProfile(
    uint pronr2,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    float* __restrict__ tmpss2datbuffers,
    float* __restrict__ dp2alndatbuffers )
{
    //cache for probabilities; add 1 to eliminate bank conflicts when accessed
    __shared__ float scaledprobsCache[CUSS_N_DIFF_SCORES][CUSS_SCALE_FACTOR+1];
    __shared__ float unscldprobsCache[CUSS_N_DIFF_SCORES];
    __shared__ float tmpCache[TIMES2(CUSS_N_DIFF_SCORES)];
    __shared__ float totalsumshared;
    __shared__ float statparsCache[nTDP2OutputPhase2_2-nTDP2OutputPhase2_1];
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

    //calculate expected score and rescale appropriately; synchronization inside;
    CalcExpectedScoreProfile( 
            &statparsCache[REFSS(dp2oadE)],
            tmpCache, scaledprobsCache, leftndx );
    //get the min and max scores; synchronization inside;
    GetScoreMinMaxProfile( 
            &statparsCache[REFSS(dp2oadMin)],
            &statparsCache[REFSS(dp2oadMax)],
            tmpCache, scaledprobsCache, leftndx );

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT && threadIdx.y==0 && threadIdx.x==0) {
        printf("\n\n MIN= %.1f MAX= %.1f E= %.4f (%.4f)\n\n", 
            statparsCache[REFSS(dp2oadMin)],
            statparsCache[REFSS(dp2oadMax)],
            statparsCache[REFSS(dp2oadE)],
            statparsCache[REFSS(dp2oadE)]/CUSS_SCALE_FACTOR );
    }
#endif

    //calculate lambda using scaled scores; synchronization inside;
    CalcLambdaProfile( 
        &statparsCache[REFSS(dp2oadLmbd)],
        tmpCache,
        scaledprobsCache,
        leftndx,
        statparsCache[REFSS(dp2oadE)],//scled expected score (used only for a check)
        statparsCache[REFSS(dp2oadMin)],//min
        statparsCache[REFSS(dp2oadMax)]//max
    );

    //calculate relent using scaled scores; synchronization inside;
    CalcEntropyProfile(
        &statparsCache[REFSS(dp2oadH)],
        tmpCache,
        scaledprobsCache,
        leftndx,
        statparsCache[REFSS(dp2oadMin)],//min
        statparsCache[REFSS(dp2oadMax)],//max
        statparsCache[REFSS(dp2oadLmbd)]//lambda
    );


    if( threadIdx.y == 0 ) {
        scaledprobsCache[0-CUSS_SCORE_MIN][threadIdx.x] = 1.0f;
        if( threadIdx.x == 0 ) {
            //rescale lambda and the expected score; synchronization below...
            statparsCache[REFSS(dp2oadLmbd)] *= CUSS_SCALE_FACTOR;
            statparsCache[REFSS(dp2oadE)] = 
                __fdividef(statparsCache[REFSS(dp2oadE)], CUSS_SCALE_FACTOR);
        }
    }

    //get the min and max scores for UNscaled scores; synchronization inside!
    GetScoreMinMaxProfile( 
            &statparsCache[REFSS(dp2oadMin)],
            &statparsCache[REFSS(dp2oadMax)],
            tmpCache, unscldprobsCache );

//     //calculate relent using unscaled scores; synchronization inside;
//     float relent = CalcEntropyProfile(
//         &statparsCache[REFSS(dp2oadH)],
//         tmpCache,
//         unscldprobsCache,
//         statparsCache[REFSS(dp2oadMin)],//min
//         statparsCache[REFSS(dp2oadMax)],//max
//         statparsCache[REFSS(dp2oadLmbd)],//rescaled lambda 
//     );

    CalcKProfile(
        &statparsCache[REFSS(dp2oadK)],
        tmpCache,
        unscldprobsCache,
        scaledprobsCache,
        statparsCache[REFSS(dp2oadMin)],
        statparsCache[REFSS(dp2oadMax)],
        //range of scores for unscaled scores unscldscorerange:
          statparsCache[REFSS(dp2oadMax)]-
          statparsCache[REFSS(dp2oadMin)]+1.0f,
        leftndx,
        rghtndx,
        statparsCache[REFSS(dp2oadE)],//expected
        statparsCache[REFSS(dp2oadLmbd)],
        statparsCache[REFSS(dp2oadH)]
    );

    // coalescent write of the results to GMEM
    if( threadIdx.y == 0 && threadIdx.x < REFSS(nTDP2OutputPhase2_2))
        dp2alndatbuffers[nTDP2OutputAlnData*pronr2 + nTDP2OutputPhase2_1 + threadIdx.x] = 
            statparsCache[threadIdx.x];


#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT && threadIdx.y==0 && threadIdx.x==0) {
        printf("\n MIN= %.1f MAX= %.1f Expected= %.4f LMB= %.4f H= %.4f K= %.4f\n\n", 
            statparsCache[REFSS(dp2oadMin)],
            statparsCache[REFSS(dp2oadMax)],
            statparsCache[REFSS(dp2oadE)],
            statparsCache[REFSS(dp2oadLmbd)],
            statparsCache[REFSS(dp2oadH)],
            statparsCache[REFSS(dp2oadK)] );
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

// -------------------------------------------------------------------------
// CalcSignificanceSS18DynPlmProfile: device code for calculating alignment 
// significance using the caclulated statistical parameters and based on the 
// profile attributes (model SS18); 
// pronr2, profile serial number in phase 2;
// nqyposs, number of query positions (the length of query);
// db1prolen, length of the (given db) target profile;
// qyeno, query ENO;
// db1proeno, ENO of the db target profile;
// searchspace, search space;
// reflambda, reference lambda (ungapped configuration);
// refK, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
// 
__global__ void CalcSignificanceSS18DynPlmProfile(
    uint pronr2,
    const int ssemodel,
    uint nqyposs, uint db1prolen,
    float qyeno, float db1proeno,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    CUBDP_TYPE* __restrict__ tmpdpbuffer,
    float* __restrict__ tmpss2datbuffers,
    float* __restrict__ dp2alndatbuffers )
{
    __shared__ float NNweightCache[sz_total_dc_cuss_NNwghts];//[read-only]
    __shared__ float lmbdKCache[2];//lambda and K [saved already]
    __shared__ float tmpdatCache[CUSS_2DCACHE_DIM];//[temporary]
    __shared__ float statsCache[nTDP2OutputPhase2_1];//[to be saved]
    float mu, scale;

    // coalescent read of statistics obtained during the DP finalization;
    // synchronization BELOW once all required data has been cached!
    if( threadIdx.x < nTDP2OutputPhase1 ) {
        statsCache[threadIdx.x] = 
            tmpdpbuffer[nTDP2OutputPhase1*pronr2 + threadIdx.x];
        // coalescent read of calculated lambda and K;
        if( threadIdx.x < 2 )
            lmbdKCache[threadIdx.x] = 
                dp2alndatbuffers[nTDP2OutputAlnData*pronr2 + dp2oadLmbd + threadIdx.x];
    }

    //cache NN weights
    NNweightCache[threadIdx.x] = dc_cuss_NN_weights_[threadIdx.x];
    if( blockDim.x + threadIdx.x < sz_total_dc_cuss_NNwghts )
        NNweightCache[blockDim.x + threadIdx.x] = dc_cuss_NN_weights_[blockDim.x + threadIdx.x];
    if( 2*blockDim.x + threadIdx.x < sz_total_dc_cuss_NNwghts )
        NNweightCache[2*blockDim.x + threadIdx.x] = dc_cuss_NN_weights_[2*blockDim.x + threadIdx.x];

    __syncthreads();

    //NOTE: should be called for MAPDP initialization 
    InitS2DatBufferForMAPDP( pronr2, tmpss2datbuffers );

    //estimate gapped lambda and K; NOTE: synchronization omitted
    EstimateGappedLambdaK(
        &statsCache[dp2oadLmbdEst],
        &statsCache[dp2oadKEst],
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK,
        expgappedlambda, expgappedK );

    //calculate e-value for the pair of profiles; synchronization inside
    GetPairEvalue(
        &statsCache[dp2oadEpA],
        nqyposs, db1prolen,
        statsCache[dp2oadScore], 
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK );


    //predict EVD mu and scale; 
    //NOTE: using registers for output values; only tid 0 will have the results;
    //synchronization before and after predictions inside!
    PredictEVDParametersSS18(
        &mu, &scale,
        ssemodel,
        nqyposs, db1prolen, qyeno, db1proeno, lmbdKCache[0],
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
        orgpronr,
#endif
        tmpdatCache, NNweightCache );

    //calculate e-value; synchronization inside
    CalculateEvalueSS18(
        &statsCache[dp2oadEvalue],
        ssemodel,
        nqyposs, db1prolen,
        searchspace,
        statsCache[dp2oadScore], 
        statsCache[dp2oadPstvs],
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK,
        expgappedlambda, expgappedK,
        mu, scale
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
        , orgpronr
#endif
        );


    // coalescent write of data and calculated measures
    if( threadIdx.x < nTDP2OutputPhase2_1 )
        dp2alndatbuffers[nTDP2OutputAlnData*pronr2 + threadIdx.x] =
            statsCache[threadIdx.x];

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT) {
        if(threadIdx.x==0) {
            printf("\n\n Q Len/Eno= %u/%.1f T Len/Eno= %u/%.1f Effss= %.1f "
                "lambda_ref= %.4f K_ref= %.4f lambda_exp= %.4f K_exp= %.4f\n",
                nqyposs,qyeno,db1prolen,db1proeno,searchspace,
                reflambda,refK,expgappedlambda,expgappedK);
            printf(" SS18 predicted mu= %.4f scale= %.4f\n", mu, scale);
            printf(" SS18 score= %.4f psts= %.1f lambda_est= %.4f K_est= %.4f "
                "LogEpA= %.4f LogE= %.4f\n\n", 
                statsCache[dp2oadScore], statsCache[dp2oadPstvs],
                statsCache[dp2oadLmbdEst], statsCache[dp2oadKEst], 
                statsCache[dp2oadEpA], statsCache[dp2oadEvalue]);
        }
//         for(int l1=123;l1<1234;l1+=220)
//             for(int l2=123;l2<1234;l2+=220)
//                 for(float e1=1.0f;e1<16.0f;e1+=3.0f)
//                     for(float e2=1.0f;e2<16.0f;e2+=3.0f) {
//                         FormatInputConfigurationE(l1,l2,e1,e2,inputCache);
//                         PredictParameterbyNNconfEProfile(
//                             &outparamsCache[cuss_NN_E_MU],tmpdatCache,inputCache,NNweightCache);
//                         if(threadIdx.x==0)
//                             printf(" l1= %d l2= %d e1= %.1f e2= %.1f predicted NN_E_MU= %.4f\n",
//                                 l1,l2,e1,e2,outparamsCache[cuss_NN_E_MU]);
//                     }
//         for(int l1=123;l1<1234;l1+=220)
//             for(int l2=123;l2<1234;l2+=220)
//                 for(float lm=0.1f;lm<6.0f;lm+=0.3f) {
//                     FormatInputConfigurationL(l1,l2,lm,inputCache);
//                     PredictParametersby2NNsConfLProfile( outparamsCache + cuss_NN_L_MU,
//                         tmpdatCache, inputCache, NNweightCache + n_dc_cuss_NNwghts_E_MU_ );
//                     if(threadIdx.x==0)
//                         printf(" l1= %d l2= %d lmb= %.1f predicted NN_L_MU= %.4f NN_L_SCALE= %.4f\n",
//                             l1,l2,lm,outparamsCache[cuss_NN_L_MU],outparamsCache[cuss_NN_L_SCALE]);
//                 }
    }
#endif
}

// -------------------------------------------------------------------------
// CalcSignificanceObsDynPlmProfile: device code for calculating the 
// significance of alignment score based on calculated lambda and K; 
// pronr2, profile serial number in phase 2;
// nqyposs, number of query positions (the length of query);
// db1prolen, length of the (given db) target profile;
// searchspace, search space;
// reflambda, reference lambda (ungapped configuration);
// refK, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
// 
__global__ void CalcSignificanceObsDynPlmProfile(
    uint pronr2,
    uint nqyposs, uint db1prolen,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    CUBDP_TYPE* __restrict__ tmpdpbuffer,
    float* __restrict__ tmpss2datbuffers,
    float* __restrict__ dp2alndatbuffers )
{
    __shared__ float lmbdKCache[2];//lambda and K [saved already]
    __shared__ float statsCache[nTDP2OutputPhase2_1];//[to be saved]

    // coalescent read of statistics obtained during the DP finalization;
    // synchronization BELOW once all required data has been cached!
    if( threadIdx.x < nTDP2OutputPhase1 ) {
        statsCache[threadIdx.x] = 
            tmpdpbuffer[nTDP2OutputPhase1*pronr2 + threadIdx.x];
        // coalescent read of calculated lambda and K;
        if( threadIdx.x < 2 )
            lmbdKCache[threadIdx.x] = 
                dp2alndatbuffers[nTDP2OutputAlnData*pronr2 + dp2oadLmbd + threadIdx.x];
    }

    __syncthreads();

    //NOTE: should be called for MAPDP initialization 
    InitS2DatBufferForMAPDP( pronr2, tmpss2datbuffers );

    //estimate gapped lambda and K; NOTE: synchronization omitted
    EstimateGappedLambdaK(
        &statsCache[dp2oadLmbdEst],
        &statsCache[dp2oadKEst],
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK,
        expgappedlambda, expgappedK );

    //calculate e-value for the pair of profiles; synchronization inside
    GetPairEvalue(
        &statsCache[dp2oadEpA],
        nqyposs, db1prolen,
        statsCache[dp2oadScore], 
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK );


    //calculate e-value; no synchronization!
    CalculateEvalueObs(
        &statsCache[dp2oadEvalue],
        searchspace,
        statsCache[dp2oadScore], 
        lmbdKCache[0]/*lambda*/, lmbdKCache[1]/*K*/,
        reflambda, refK,
        expgappedlambda, expgappedK );

    __syncthreads();


    // coalescent write of data and calculated measures
    if( threadIdx.x < nTDP2OutputPhase2_1 )
        dp2alndatbuffers[nTDP2OutputAlnData*pronr2 + threadIdx.x] =
            statsCache[threadIdx.x];

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT) {
        if(threadIdx.x==0) {
            printf("\n\n Q Len= %u T Len= %u Effss= %.1f "
                "lambda_ref= %.4f K_ref= %.4f lambda_exp= %.4f K_exp= %.4f\n",
                nqyposs,db1prolen,searchspace,
                reflambda,refK,expgappedlambda,expgappedK);
            printf(" SSObs score= %.4f lambda_est= %.4f K_est= %.4f LogEpA= %.4f LogE= %.4f\n\n", 
                statsCache[dp2oadScore],
                statsCache[dp2oadLmbdEst], statsCache[dp2oadKEst], 
                statsCache[dp2oadEpA], statsCache[dp2oadEvalue]);
        }
    }
#endif
}



// =========================================================================
// CalcStatisticsDynPlm: device code for calculating alignment statistics 
// using dynamic parallelism;
// NOTE: memory pointers should be aligned!
// ssemodel, number of the model for statistical significance estimation;
// qyeno, query ENO;
// searchspace, search space;
// reflambda, reference lambda (ungapped configuration);
// refK, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
// scores, calculated scores used as input;
// tmpdpbuffer, temporary buffer that at this stage contains the 
// information of profiles passed;
// tmpss2datbuffers, profile-specific temporary buffers for score 
// probabilities;
// dp2alndatbuffers, device memory for profile-profile alignment statistics;
// 
__global__ void CalcStatisticsDynPlm(
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint /*querposoffset*/, uint bdb1posoffset, uint /*bdbCposoffset*/,
    //
    const int ssemodel,
    float qyeno,
    float searchspace,
    float reflambda, float refK,
    float expgappedlambda, float expgappedK,
    //
    CUBSM_TYPE* __restrict__ scores,
    CUBDP_TYPE* __restrict__ tmpdpbuffer,
    float* __restrict__ tmpss2datbuffers,
    float* __restrict__ dp2alndatbuffers )
{
    //NOTE: assuming only ONE thread in a block
    //
//     cudaStream_t streamss;
//     MYCUDACHECK( cudaStreamCreateWithFlags(&streamss, cudaStreamNonBlocking));

    //initialize memory for calculating score probabilities
    InitScoreProbsDynPlmProfile<<<1,CUSS_ALIGNED_N_DIFF_TOTAL_SCORES/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
        tmpss2datbuffers
    );
    MYCUDACHECKLAST;


    uint serpronrSM;//profile serial number in the score matrix (current processing)
    uint orgpronrCache;//original profile number/index
    LNTYPE dbprodstCache;//distance in positions to the original db profile
    INTYPE dbprolenCache;//length of the profile
    FPTYPE dbproenoCache;//ENO of the profile

    const int dblen = ndb1poss + ndbCposs + dbxpad;

    // blockIdx.x is the profile serial number in phase 2;
    serpronrSM = *(uint*)(tmpdpbuffer + nTDP2OutputPhase1*blockIdx.x+dp2oadOrgProNo);
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
    dbproenoCache = ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO]))[orgpronrCache];

    //db profile position in the dc_pm2dvfields_ buffers
    int dbpos = dbprodstCache;
    //dbpos is now the beginning x position in the score matrix:
    dbpos = (serpronrSM < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;


    //execution configuration for the profile this parent thread is responsible for:
    dim3 nthrds(CUSS_2DCACHE_DIM,CUSS_2DCACHE_DIM,1);
    dim3 nblcks((dbprolenCache+nthrds.x-1)/nthrds.x,
                (nqyposs+nthrds.y-1)/nthrds.y,
                1);

    //calculate score probabilities (counts)
    CalcScoreProbsDynPlmProfile<<<nblcks,nthrds/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
        (uint)nqyposs, (uint)dbprolenCache, dbpos, dblen,
        scores,
        tmpss2datbuffers
    );
    MYCUDACHECKLAST;

    //calculate statistical parameters
    nthrds = dim3(CUSS_SCALE_FACTOR,CUSS_N_DIFF_SCORES,1);
    CalcStatParamsDynPlmProfile<<<1,nthrds/*,0,streamss*/>>>(
        blockIdx.x/*phase-2 profile number*/,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
        orgpronrCache,
#endif
        tmpss2datbuffers,
        dp2alndatbuffers
    );
    MYCUDACHECKLAST;


    //calculate alignment significance
    nthrds = dim3(CUSS_2DCACHE_DIM,1,1);
    if( ssemodel == 0 )
        CalcSignificanceObsDynPlmProfile<<<1,nthrds/*,0,streamss*/>>>(
            blockIdx.x,//phase-2 profile number
            (uint)nqyposs, (uint)dbprolenCache,
            searchspace,
            reflambda, refK,
            expgappedlambda, expgappedK,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
            orgpronrCache,
#endif
            tmpdpbuffer,
            tmpss2datbuffers,
            dp2alndatbuffers
        );
    else
        CalcSignificanceSS18DynPlmProfile<<<1,nthrds/*,0,streamss*/>>>(
            blockIdx.x,//phase-2 profile number
            ssemodel,
            (uint)nqyposs, (uint)dbprolenCache,
            qyeno, dbproenoCache,
            searchspace,
            reflambda, refK,
            expgappedlambda, expgappedK,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
            orgpronrCache,
#endif
            tmpdpbuffer,
            tmpss2datbuffers,
            dp2alndatbuffers
        );
    MYCUDACHECKLAST;

//     MYCUDACHECK( cudaStreamDestroy( streamss ));
//     MYCUDACHECKLAST;
}
