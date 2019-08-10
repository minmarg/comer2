/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_SS18_h__
#define __CuBatchSS_SS18_h__

#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusf/CuSF_betainc.cuh"
#include "libmycu/cusf/CuSF_gammainc.cuh"
#include "CuBatchSS_com.h"
#include "CuBatchSS_NN.cuh"

// minimum/maximum values of mu and scale
#define CUSS_SS18_MINMU 3.0f
#define CUSS_SS18_MINSCALE 0.1f
#define CUSS_SS18_MAXSCALE 5.0f

#define CUSS_SS18_LOG001 -4.605170186f

//inverse of the number of pairs used in experiments
#define CUSS_SS18_INVNPAIRS 0.0001f
#define CUSS_SS18_LOGINVNPAIRS -9.210340372f

//data of the derived statistic: const*psts/pow(l1*l2,1.5);
//factor of psts in the numerator:
#define CUSS_SS18_PSTS_FACTOR 100000.0f
//exponent in the denominator
#define CUSS_SS18_PSTS_EXP 1.5f
//shape parameter of an NBD of the derived statistic
#define CUSS_SS18_PSTS_NBD_SHAPE_EXP 0.1f
//probability parameter of an NBD of the derived statistic
#define CUSS_SS18_PSTS_NBD_PROB_EXP 0.757213f

//data for the empirical Brown's method;
//degrees of freedom * 0.5:
#define CUSS_SS18_EBM_DFo2 1.346206f
//scale factor for the dependent estimation
#define CUSS_SS18_EBM_C 1.485656f

// -------------------------------------------------------------------------

__device__ __forceinline__
void PredictEVDParametersSS18(
    float* __restrict__ out_mu,
    float* __restrict__ out_scale,
    int ssemodel,
    uint nqyposs, uint db1prolen,
    float qyeno, float db1proeno,
    float lambda,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    float* __restrict__ tmpdatbuf,
    const float* __restrict__ weights );

__device__ __forceinline__
void CalculateEvalueSS18(
    float* __restrict__ logevalue,
    int ssemodel,
    uint nqyposs, uint db1prolen,
    float searchspace,
    float alnscore, 
    float psts,
    float /*lambda*/, float K,
    float /*lambdaref*/, float Kref,
    float expgappedlambda, float expgappedK,
    float mu, float scale
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    ,uint orgpronr
#endif
    );

__device__ __forceinline__
void CalculateEvalueObs(
    float* __restrict__ logevalue,
    float searchspace,
    float alnscore, 
    float /*lambda*/, float /*K*/,
    float /*lambdaref*/, float /*Kref*/,
    float expgappedlambda, float expgappedK );


__device__ __forceinline__
void EstimateGappedLambdaK(
    float* __restrict__ drvdlambda,
    float* __restrict__ drvdK,
    float lambda, float K,
    float lambdaref, float Kref,
    float expgappedlambda, float expgappedK );

__device__ __forceinline__
void GetPairEvalue(
    float* __restrict__ logevalue,
    uint nqyposs, uint db1prolen,
    float alnscore, 
    float lambda, float K,
    float lambdaref, float Kref );

// =========================================================================
// -------------------------------------------------------------------------
// PredictEVDParametersSS18: predict the parameters of the EVD governing 
// alignment score distribution as described by Margelevicius (bioRxiv 2018 
// doi:10.1101/484485);
// out_mu, mu parameter to be predicted;
// out_scale, scale parameter to be predicted;
// ssemodel, variant index of model SS18;
// nqyposs, number of query positions (the length of query);
// db1prolen, length of the (given db) target profile;
// qyeno, query ENO;
// db1proeno, ENO of the db target profile;
// lambda, lambda obtained for a pair of profiles;
// tmpdatbuf, SMEM buffer for temporary data;
// weights, ANN weights;
//
__device__ __forceinline__
void PredictEVDParametersSS18(
    float* __restrict__ out_mu,
    float* __restrict__ out_scale,
    int ssemodel,
    uint nqyposs, uint db1prolen,
    float qyeno, float db1proeno,
    float lambda,
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    uint orgpronr,
#endif
    float* __restrict__ tmpdatbuf,
    const float* __restrict__ weights )
{
    __shared__ float inputCache[4];//[for local use]
    __shared__ float outparamsCache[nTcuss_NN_Ouputs];//[for local use]

    //format input for NN; synchronization inside!
    FormatInputConfigurationE(
        nqyposs, db1prolen, qyeno, db1proeno,
        inputCache );

    //predict mu based on profile attributes (configuration E); sync inside
    PredictParameterbyNNconfEProfile(
        &outparamsCache[cuss_NN_E_MU],
        tmpdatbuf, inputCache, weights );


    //format input for predicting parameters based on compositional similarity; 
    //synchronization inside!
    FormatInputConfigurationL(
        nqyposs, db1prolen, lambda,
        inputCache );

    //predict mu and scale SIMULTANEOUSLY based on compositional similarity (configuration L);
    //synchronization inside!
    PredictParametersby2NNsConfLProfile(
        outparamsCache + cuss_NN_L_MU,
        tmpdatbuf, inputCache, weights + n_dc_cuss_NNwghts_E_MU_ );


    // only tid 0 proceeds further!
    if( threadIdx.x == 0 ) 
    {
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
        if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT) {
            printf(" predicted NN_E_MU= %.4f NN_L_MU= %.4f NN_L_SCALE= %.4f\n\n",
                outparamsCache[cuss_NN_E_MU], outparamsCache[cuss_NN_L_MU], 
                outparamsCache[cuss_NN_L_SCALE]);
        }
#endif
        //adjustment parameters
        const float* adjparsSS18 = weights + n_dc_cuss_adjparamsSS18_start;

        switch( ssemodel ) {
            case SS18variant1: adjparsSS18 = weights + n_dc_cuss_adjparamsSS18_start; break;
            case SS18variant2: adjparsSS18 = weights + n_dc_cuss_adjparamsSS18_var1; break;
        }

        if( outparamsCache[cuss_NN_E_MU] < CUSS_SS18_MINMU ) 
            outparamsCache[cuss_NN_E_MU] = CUSS_SS18_MINMU;
        if( outparamsCache[cuss_NN_L_MU] < CUSS_SS18_MINMU ) 
            outparamsCache[cuss_NN_L_MU] = CUSS_SS18_MINMU;


        if( outparamsCache[cuss_NN_L_SCALE] < CUSS_SS18_MINSCALE ) 
            outparamsCache[cuss_NN_L_SCALE] = CUSS_SS18_MINSCALE;

        outparamsCache[cuss_NN_L_SCALE] = adjparsSS18[ss18scaleR_F] * 
                (__expf(outparamsCache[cuss_NN_L_SCALE]) - 1.0f);

        if( outparamsCache[cuss_NN_L_SCALE] < CUSS_SS18_MINSCALE ) 
            outparamsCache[cuss_NN_L_SCALE] = CUSS_SS18_MINSCALE;
        if( CUSS_SS18_MAXSCALE < outparamsCache[cuss_NN_L_SCALE] ) 
            outparamsCache[cuss_NN_L_SCALE] = CUSS_SS18_MAXSCALE;


        float scaleE = adjparsSS18[ss18scaleE_S] * outparamsCache[cuss_NN_E_MU] + 
                       adjparsSS18[ss18scaleE_I];

        if( scaleE < CUSS_SS18_MINSCALE ) scaleE = CUSS_SS18_MINSCALE;
        if( CUSS_SS18_MAXSCALE < scaleE ) scaleE = CUSS_SS18_MAXSCALE;


        outparamsCache[cuss_NN_E_MU] = 
            adjparsSS18[ss18muE_S] * outparamsCache[cuss_NN_E_MU] + adjparsSS18[ss18muE_I];
        outparamsCache[cuss_NN_L_MU] = 
            adjparsSS18[ss18muR_S] * outparamsCache[cuss_NN_L_MU] + adjparsSS18[ss18muR_I];

//         if( outparamsCache[cuss_NN_E_MU] < CUSS_SS18_MINMU ) 
//             outparamsCache[cuss_NN_E_MU] = CUSS_SS18_MINMU;
//         if( outparamsCache[cuss_NN_L_MU] < CUSS_SS18_MINMU ) 
//             outparamsCache[cuss_NN_L_MU] = CUSS_SS18_MINMU;


        float alpha = adjparsSS18[ss18alpha];
        float ibeta = adjparsSS18[ss18ibeta];
        if( 0.0f < lambda ) {
            float rcplmbd = __fdividef(1.0f, lambda);
            //alpha += (1.0f - alpha) * 0.1f * red;
            alpha -= alpha * 0.1f * SLC_MIN(10.0f,rcplmbd);
            //
            //ibeta += (1.0f - ibeta) * 0.1f * red;
            ibeta += (1.0f - ibeta) * 0.1f * SLC_MIN(10.0f,rcplmbd);
        }

        //WRITE conditional mean estimate of mu
        *out_mu = alpha * outparamsCache[cuss_NN_E_MU] + 
            (1.0f - alpha) * outparamsCache[cuss_NN_L_MU];
        //WRITE conditional mean estimate of scale
        *out_scale = (1.0f - ibeta) * scaleE + 
            ibeta * outparamsCache[cuss_NN_L_SCALE];
    }

    __syncthreads();
}

// -------------------------------------------------------------------------
// CalculateEvalueSS18: calculate E-value based on statistical model SS18 as 
// described by Margelevicius (bioRxiv 2018 doi:10.1101/484485);
// logevalue, log of resulting E-value;
// ssemodel, variant index of model SS18;
// nqyposs, number of query positions (the length of query);
// db1prolen, length of the (given db) target profile;
// searchspace, search space;
// alnscore, profile-profile alignment score;
// psts, # positive substitution scores in the profile-profile alignment;
// lambda, lambda calculated for a pair of profiles;
// K, K calculated for a pair of profiles;
// lambdaref, reference lambda (ungapped configuration);
// Kref, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
// mu, predicted mu parameter of an EVD;
// scale, predicted scale parameter of an EVD;
//
__device__ __forceinline__
void CalculateEvalueSS18(
    float* __restrict__ logevalue,
    int ssemodel,
    uint nqyposs, uint db1prolen,
    float searchspace,
    float alnscore, 
    float psts,
    float /*lambda*/, float K,
    float /*lambdaref*/, float Kref,
    float expgappedlambda, float expgappedK,
    float mu, float scale
#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
    ,uint orgpronr
#endif
    )
{
    // only tid 0 proceeds!
    if( threadIdx.x == 0 ) 
    {
        if( mu <= 0.0f || scale <= 0.0f ) {
            //fallback estimation...
            //NOTE: no synchronization inside the following function
            CalculateEvalueObs(
                logevalue,
                searchspace, alnscore, 
                1.0f/*lambda*/, K,
                1.0f/*lambdaref*/, Kref,
                expgappedlambda, expgappedK );
        }
        else {
            float l1l2 = nqyposs * db1prolen;
            float logKmn_der = __fdividef(mu, scale);
            if( K > Kref && 0.0f < K && K < 1.0f && 0.0f < Kref && Kref < 1.0f )
                logKmn_der += __logf(__fdividef(K, Kref));//K_d=K*K_e/K_u

            //significance of alignment score
            float logexpect = logKmn_der - __fdividef(alnscore, scale);

            if( ssemodel == SS18variant2 )
                logexpect += CUSS_SS18_LOGINVNPAIRS;

            float logpvalue = logexpect;

            if( CUSS_SS18_LOG001 <= logexpect ) {
                if( SLC_LOG_SP_MAX - 3.0f < logexpect )
                    logpvalue = 0.0f;
                else {
                    logpvalue = __expf(logexpect);
                    if( SLC_LOG_SP_MAX - 3.0f < logpvalue )
                        logpvalue = 0.0f;
                    else
                        logpvalue = __logf(1.0f - __expf(-logpvalue));
                }
            }

            //significance of # positive scores;
            //PSTS normalized, general case
            psts = rintf( 
                __fdividef(CUSS_SS18_PSTS_FACTOR*psts, 
                        __powf( myhdmax(l1l2,CUSS_PSTS_MIN_LENPRODUCT),
                                CUSS_SS18_PSTS_EXP )) );

            float logpvaluepsts = cusf_lbetaincf( psts+1.0f, 
                        CUSS_SS18_PSTS_NBD_SHAPE_EXP, CUSS_SS18_PSTS_NBD_PROB_EXP );

            //combining dependent significance by the empirical Brown's method;
            //combined log of p-value based on the score and psts measures:
            logpvalue = cusf_lgammainc_Qf( CUSS_SS18_EBM_DFo2, 
                        __fdividef(-logpvalue-logpvaluepsts, CUSS_SS18_EBM_C));

#ifdef CUSS_SCOREPROBS_DYNPLM_TESTPRINT
            if(orgpronr==CUSS_SCOREPROBS_DYNPLM_TESTPRINT) {
                printf(" SS18 psts_norm= %.4f logpvaluepsts= %.4f logpvalue= %.4f\n\n",
                    psts, logpvaluepsts, logpvalue );
            }
#endif
            //transform to the log of e-value
            logexpect = logpvalue;

            if( CUSS_SS18_LOG001 <= logpvalue ) {
                logpvalue = __expf(logpvalue);
                logexpect = 77.0f;
                if( logpvalue < 1.0f ) {
                    logexpect = -__logf(1.0f - logpvalue);
                    logexpect = __logf(logexpect);
                }
            }

            //take account of the db size
            logexpect += __logf(__fdividef(searchspace, l1l2));

            *logevalue = logexpect;
        }
    }

    __syncthreads();
}


// -------------------------------------------------------------------------
// CalculateEvalueObs: calculate E-value based on solved statistical 
// parameter;
// logevalue, log of resulting E-value;
// searchspace, search space;
// alnscore, profile-profile alignment score;
// lambda, lambda calculated for a pair of profiles;
// K, K calculated for a pair of profiles;
// lambdaref, reference lambda (ungapped configuration);
// Kref, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
//
__device__ __forceinline__
void CalculateEvalueObs(
    float* __restrict__ logevalue,
    float searchspace,
    float alnscore, 
    float /*lambda*/, float /*K*/,
    float /*lambdaref*/, float /*Kref*/,
    float expgappedlambda, float expgappedK )
{
    // tid 0 calculates
    if( threadIdx.x == 0 ) 
    {
        *logevalue = __logf(expgappedK * searchspace) - expgappedlambda * alnscore;
    }

    //NOTE: do not synchronize here
    //__syncthreads();
}


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// EstimateGappedLambdaK: calculate the gapped parameters lambda and K;
// drvdlambda, address to hold calculated value of lambda;
// drvdK, address to hold calculated value of K;
// lambda, lambda calculated for a pair of profiles;
// K, K calculated for a pair of profiles;
// lambdaref, reference lambda (ungapped configuration);
// Kref, reference K (ungapped configuration);
// expgappedlambda, empirically determined lambda for gapped alignments;
// expgappedK, empirically determined K for gapped alignments;
//
__device__ __forceinline__
void EstimateGappedLambdaK(
    float* __restrict__ drvdlambda,
    float* __restrict__ drvdK,
    float lambda, float K,
    float lambdaref, float Kref,
    float expgappedlambda, float expgappedK )
{
    // tid 0 calculates
    if( threadIdx.x == 0 ) 
    {
        //(computed ungapped) * (empirically determined gapped) / (reference computed ungapped)
        if( 0.0f < lambda )
            *drvdlambda = lambda * __fdividef(expgappedlambda, lambdaref);
        if( 0.0f < K )
            *drvdK = K * __fdividef(expgappedK, Kref);
    }

    //NOTE: sync can be omitted as long as another synchronizing call follows
    //__syncthreads();
}

// -------------------------------------------------------------------------
// GetPairEvalue: calculate E-value for the pair of profiles given their 
// alignment score;
// logevalue, log of resulting E-value;
// nqyposs, number of query positions (the length of query);
// db1prolen, length of the (given db) target profile;
// alnscore, profile-profile alignment score;
// lambda, lambda calculated for a pair of profiles;
// K, K calculated for a pair of profiles;
// lambdaref, reference lambda (ungapped configuration);
// Kref, reference K (ungapped configuration);
//
__device__ __forceinline__
void GetPairEvalue(
    float* __restrict__ logevalue,
    uint nqyposs, uint db1prolen,
    float alnscore, 
    float lambda, float K,
    float lambdaref, float Kref )
{
    // tid 0 calculates
    if( threadIdx.x == 0 ) 
    {
        float l1l2 = nqyposs * db1prolen;

        if( lambda <= 0.0f || K <= 0.0f )
            *logevalue = __logf(Kref * l1l2) - lambdaref * alnscore;
        else
            *logevalue = __logf(K * l1l2) - lambda * alnscore;
    }

    __syncthreads();
}

#endif//__CuBatchSS_SS18_h__
