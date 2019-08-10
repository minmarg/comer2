/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_com_h__
#define __CuBatchDP_com_h__

#include "libpro/srcpro/TRANSPROBS.h"

#ifndef CUBATCHDP_FLOAT
#define CUBATCHDP_FLOAT
#define CUBDP_TYPE  float
#define CUBDP_Q(CNST) (CNST ## . ## f)
#else
#define CUBATCHDP_INT
#define CUBDP_TYPE  int
#define CUBDP_Q(CNST) (CNST)
#endif//CUBATCHDP_FLOAT


// Dimension of 2D cache of shared memory for performing dynamic programming;
// it also represents 2D dimensions for the CUDA thread block

#define CUDP_2DCACHE_DIM_DequalsX
//diagonal size (warp size)
#define CUDP_2DCACHE_DIM_D 32
//size along x dimension (cannot be greater than CUDP_2DCACHE_DIM_D)
#define CUDP_2DCACHE_DIM_X 32 //32 //24 //16 //8
//CUDP_2DCACHE_DIM_D + CUDP_2DCACHE_DIM_X
#define CUDP_2DCACHE_DIM_DpX 64 //64 //56 //48 //40


//make sure it equals to MAIN_TRNPRB_LOGFCT
#define CUDP_TRNPRB_LOGFCT MAIN_TRNPRB_LOGFCT

//number of scores used to calculated correlation score
#define CUDP_CORR_NSCORES 4

//delimiter for floating point variables
#define CUDP_FLOAT_DELM 64.0f


//maximum edge length of a bottom-left and top-right triangular of the DP matrix
//that can be ignored
#define CUDP_MAXCORNERLEN 60
//percentage of profile length used to calculate the maximum edge length of corner 
//triangulars
#define CUDP_CORNERLENPRC 0.34f


enum TDPDiagScoreSections {
    //DP scores
    dpdssDiag1,//buffer of the first diagonal
    dpdssDiag2,//buffer of the second diagonal
    dpdssDiagM,//buffer of the maximum scores found up to the second diagonal
    //
    //match and correlation scores 
    dpdssCorrDiagM,//buffer of correlation scores corresponding to max DP aln scores (dpdssDiagM)
    dpdssCorrDiag1,//buffer of the scores of the first diagonal for calculating correlations 
    dpdssCorrDiag2,//buffer of the scores of the second diagonal for calculating correlations 
    dpdssCorrDiagS1,//buffer of correlation scores of the first diagonal
    dpdssCorrDiagS2,//buffer of correlation scores of the second diagonal
    nTDPDiagScoreSections
};

enum TDPBottomScoreSections {
    //DP scores
    dpbssBottm,//buffer of the bottom DP scores
    //
    //match and correlation scores 
    dpbssCorrBottm,//buffer of the bottom scores for calculating correlations
    dpbssCorrBottmS,//buffer of correlation scores of the blocks' bottom
    nTDPBottomScoreSections
};

enum TDPDiagScoreSubsections {
    dpdsssStateMM,//match-to-match state
    dpdsssStateMI,//match-to-insert
    dpdsssStateIM,//insert-to-match
    dpdsssStateDN,//delete-to-no-state(gap)
    dpdsssStateND,//no-state-to-delete
    nTDPDiagScoreSubsections
};

//used transitions in the profile transition model while executing DP;
//NOTE: the order must match that of the TPTRANS enum
enum TDPUsedProTrn {
    dptrMMp,//match-match, previous position
    dptrMIc,//match-insert, current position
    dptrMDp,//match-delete, previous position
    dptrIMp,//insert-match, previous position
    dptrIIc,//insert-insert, current position
    dptrDMp,//delete-match, previous position
    dptrDDp,//delete-delete, previous position
    nTDPUsedProTrn
};

//backtracking direction constants
enum TDPBtckDirection {
    dpbtckSTOP = 0,
    dpbtckLEFT = 1,//left, gap in query
    dpbtckUP = 2,//up, gap in db profile
    dpbtckDIAG = 3,//diagonal, match between query and db target
};

#endif//__CuBatchDP_com_h__
