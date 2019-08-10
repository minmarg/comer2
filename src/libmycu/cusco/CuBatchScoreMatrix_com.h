/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_com_h__
#define __CuBatchScoreMatrix_com_h__

#ifndef CUBATCHSCOREMATRIX_FLOAT
#define CUBATCHSCOREMATRIX_FLOAT
#define CUBSM_TYPE  float
#define CUBSM_Q(CNST) (CNST ## . ## f)
#else
#define CUBATCHSCOREMATRIX_INT
#define CUBSM_TYPE  int
#define CUBSM_Q(CNST) (CNST)
#endif//CUBATCHSCOREMATRIX_FLOAT

// Dimension of 2D cache of shared memory for calculating score matrix
#define SMINIT_2DCACHE_DIM 32

#endif//__CuBatchScoreMatrix_com_h__
