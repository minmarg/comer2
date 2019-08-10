/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mylimits_h__
#define __mylimits_h__

#define LOG_PROB_MIN ( -32767 )
#define SCORE_MIN ( -32767 )

#ifdef UINT_MAX
#   define MYSIZE_MAX UINT_MAX
#else
#   define MYSIZE_MAX 4294967295
#endif

//1024,1024*1024,1024*1024*1024
#define ONEK       1024
#define ONEM    1048576
#define ONEG 1073741824UL

#ifndef INT_MAX
#   define INT_MAX 2147483647
#endif

#ifndef UINT_MAX
#   define UINT_MAX 4294967295U
#endif

// default filename of scores to be read 
#define DEFAULT_PSCORES_FILENAME    "PSCORES"

// default filename of mixture parameters to be read
#define DEFAULT_MIXTURES_FILENAME   "MIXTURES"

// initial hit reserve
#define ALLOCHITS   1024

// default size for a short string
#define BUF_MAX     256
#define KBYTE       1024





// maximum number of columns a sequence or profile can have;
// if exceeded the program terminates
#define MAXCOLUMNS  50000

// alignment output width per line
#define OUTPUTWIDTH 60

// output indent
#define OUTPUTINDENT 13

// length of annotations
#define ANNOTATIONLEN 75

// maximum description length in the program's output
#define MAX_DESCRIPTION_LENGTH 4096

// maximum number of iterations before terminating scaling of scores
// to obtain value for the statistical parameter lambda
#define MAX_SCALE_ITERATIONS 10

// whether to use binary search in scaling of scores to
// obtain lambda
#define SCALE_BINARY_SEARCH

// maximum number of iterations allowed to perform before
// the root of conservation function is found
#define MAXIT 100

// accuracy to which lambda parameter should be calculated
#define LAMBDA_ACCURACY ( 1.e-5f )

// small value above zero
#define EPSILON ( 1.e-6f )



// tolerance for optimization of target frequencies by the Newton's method 
#define TARGET_FREQ_OPT_TOLERANCE   ( 1.e-6f )

// maximum number of Newton's method iterations in optimizing the target frequencies
#define TARGET_FREQ_OPT_MAXITERATIONS   ( 1000 )


// number of COLUMN residues to take into window 
#define MAC_SEGSEQ_WIN_LENGTH       ( 30 )
// log2( 30/2 ) *.85, *.90
#define MAC_SEGSEQ_LOW_ENTROPY      ( 3.315f )
#define MAC_SEGSEQ_HIGH_ENTROPY     ( 3.510f )

// minimum required window size in positions of extent
#define MAC_EXTENT_MINWIN           ( 7 )
// minimum required sequence percentage an extent must cover
#define MAC_EXTENT_SEQ_COVER        ( 0.05f )

// relative weight for residue pseudocount frequencies, beta: 7, 9, 13
#define MAC_WEIGHT_PSEUDO_COUNTS    ( 7 )


// value to which frequencies always sum
#define FREQUENCY_SUM 100

#ifndef SCALE_PROFILES
#define SCALE_PROFILES
#endif


// scaling constant used to scale profile-profile scores before
// using them to align two profiles
#define PP_SCALE_CONSTANT 32

// maximum range of scores allowed
#define MAX_RANGE_OF_SCORES 2048

// maximum number of iterations allowed to perform while
// computing infinite sum of the parameter K
#define PARAMETER_K_MAXIT 100

// accuracy to compute alignment probability expressed in
// infinite sum to; the probability is used to compute K
#define PARAMETER_K_ACCURACY ( 1.e-4f )

// maximum number of times to iterate searching for fixed
// value of length adjustment expression
#define LENGTH_ADJUSTMENT_MAXIT 20



#endif//__mylimits_h__

