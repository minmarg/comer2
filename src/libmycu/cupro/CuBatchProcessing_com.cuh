/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchProcessing_com_h__
#define __CuBatchProcessing_com_h__

//calculate DP correlation scores inline
// #define CUDP_CALC_CORRSCORES_INLINE

//expected length of the proteins in the database;
//NOTE: increasing it reduces memory requirements, but 
// misprediction will cost additional computation time;
// CHANGED with an option
// #define CUDP_EXPECTED_DBPROLEN 50

// -------------------------------------------------------------------------
// general constant memory details
// 
//total number of fields obtained after packing:
const size_t fsg_sztotflds2 = TIMES2(pmv2DTotFlds);
const size_t fsg_sztotflds3 = TIMES3(pmv2DTotFlds);

//enum for indexing constant memory arrays:
enum {
    ndx_dc_pm2dvfields_,
    ndx_dc_newpm2dvfields_
};
//device addresses of vectors representing query AND cached database profile model data plus
//device addresses of vectors representing new database profile model data:
extern __constant__ char* dc_pm2dvfields_[fsg_sztotflds3];


// -------------------------------------------------------------------------
// constant memory details for significance estimation
// 

//number of inputs for ANNs based on profile attributes
#define CUSS_E_NINPUT 4
//number of inputs for ANNs based on profile compositional similarity;
//actually, it is 3, but add 1 for convenient parallel calculation
#define CUSS_L_NINPUT 4 //3+1
//number of hidden units in each of ANNs
#define CUSS_NHIDDEN 4
//number of layers in each of ANNs
#define CUSS_NLAYERS 1
//number of outputs in each of ANNs
#define CUSS_NOUTPUT 1

// log base 2 of CUSS_E_NINPUT, CUSS_L_NINPUT, and CUSS_NHIDDEN
#define CUSS_LOG_E_NINPUT 2
#define CUSS_LOG_L_NINPUT 2
#define CUSS_LOG_NHIDDEN 2

//adjustment parameters for model SS18
enum TSS18hypers {
    ss18muE_S,//slope for muE
    ss18muE_I,//intercept for muE
    ss18muR_S,//slope for muR
    ss18muR_I,//intercept for muR
    ss18scaleE_S,//slope for scaleE
    ss18scaleE_I,//intercept for scaleE
    ss18scaleR_F,//factor for scaleR
    ss18alpha,//contribution weight of muE to mu
    ss18ibeta,//contribution weight of scaleR to scale
    nTSS18hypers
};
enum TSS18models {
    SS18variant1 = 1,
    SS18variant2 = 2,
};

// cumulative number of weights for each ANN and parameters:
enum {
    //+1 for the bias unit in each layer (including the output layer)
    n_dc_cuss_NNwghts_E_MU_ = 
            (CUSS_E_NINPUT+1)*CUSS_NHIDDEN + (CUSS_NHIDDEN+1)*CUSS_NOUTPUT,
    n_dc_cuss_NNwghts_L_MU_ = n_dc_cuss_NNwghts_E_MU_ + 
            (CUSS_L_NINPUT+1)*CUSS_NHIDDEN + (CUSS_NHIDDEN+1)*CUSS_NOUTPUT,
    n_dc_cuss_NNwghts_L_SCALE_ = n_dc_cuss_NNwghts_L_MU_ + 
            (CUSS_L_NINPUT+1)*CUSS_NHIDDEN + (CUSS_NHIDDEN+1)*CUSS_NOUTPUT,
    //
    n_dc_cuss_adjparamsSS18_start = n_dc_cuss_NNwghts_L_SCALE_,
    n_dc_cuss_adjparamsSS18_var1 = n_dc_cuss_adjparamsSS18_start + nTSS18hypers,
    n_dc_cuss_adjparamsSS18_var2 = n_dc_cuss_adjparamsSS18_var1 + nTSS18hypers,
    sz_total_dc_cuss_NNwghts = n_dc_cuss_adjparamsSS18_var2
};
//weights of the ANNs for etimating statistical significance reside in constant memeory
extern __constant__ float dc_cuss_NN_weights_[sz_total_dc_cuss_NNwghts];

// -------------------------------------------------------------------------

#endif//__CuBatchProcessing_com_h__
