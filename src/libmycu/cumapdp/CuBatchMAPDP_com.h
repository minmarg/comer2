/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_com_h__
#define __CuBatchMAPDP_com_h__

// Dimension of 2D cache of shared memory for performing dynamic programming;
// it also represents 2D dimensions for the CUDA thread block

//NOTE: CUMAPDP_2DCACHE_DIM_D and CUMAPDP_2DCACHE_DIM_X should be EQUAL!
//cuurently these values are the only option:
//diagonal size (warp size)
#define CUMAPDP_2DCACHE_DIM_D 32
//size along x dimension
#define CUMAPDP_2DCACHE_DIM_X 32
//CUMAPDP_2DCACHE_DIM_D + CUMAPDP_2DCACHE_DIM_X
#define CUMAPDP_2DCACHE_DIM_DpX 64


//log factor for pair e-value
#define CUMAPDP_EPA_LOGFACTOR 22.0f //26.0f
//dynamic score offset
#define CUMAPDP_DYNS_OFFSET 0.05f //0.07f
//dynamic score multiplier
#define CUMAPDP_DYNS_MULTP 0.9f //0.7f

//log factor for pair e-value used in gap extension
#define CUMAPDP_EPAGAPEXT_LOGFACTOR 30.0f
//upper bound of the dynamic multiplier for the MAP gap extension threshold
#define CUMAPDP_DYNGAPEXT_MULTP_UB 1.0f
//lower bound dynamic multiplier for the MAP gap extension threshold
#define CUMAPDP_DYNGAPEXT_MULTP_LB 0.3f


//use extremely fast but approximate exponential and logarithmic functions:
//NOTE: making this an option induces huge occupation of registers;
// it is better to control this MACRO through compiler options (at compile time)
#define CUMAPDP_USEFASTAPPROXMATH 


//MAPDP probability scales at each query position
enum TMAPDPProbScales {
    mapdppsPosScale,//scale at a particular position
    mapdppsTotScale,//total log scale obtained up to a particular position
    nTMAPDPProbScales
};

#endif//__CuBatchMAPDP_com_h__
