/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_com_h__
#define __CuBatchSS_com_h__

// 2D shared memory cache for calculating score probabilities;
// this also corresponds to the block dimension of the execution 
// configuration
#define CUSS_2DCACHE_DIM 32

// minimum and maximum scores assumed when calculating 
// score probabilities
#define CUSS_SCORE_MIN   (-10)
#define CUSS_SCORE_MIN_f (-10.0f)
#define CUSS_SCORE_MAX    10
#define CUSS_SCORE_MAX_f  10.0f
// maximum number of different scores (max-min+1)
#define CUSS_N_DIFF_SCORES 21

// scale factor of scores for calculating statistical parameters;
// it is important that it equals the warp size
#define CUSS_SCALE_FACTOR 32
// minimum and maximum of scaled scores
#define CUSS_SCALED_SCORE_MIN (-320)
#define CUSS_SCALED_SCORE_MAX  320
// maximum number of different scaled scores (scale*(max-min)+1)
#define CUSS_N_DIFF_SCALED_SCORES 641
// maximum number of total (scaled and not) scores:
#define CUSS_N_DIFF_TOTAL_SCORES 662
// aligned (CUL2CLINESIZE=32 and warpSize=32) maximum number of 
// total (scaled and not) scores:
#define CUSS_ALIGNED_N_DIFF_TOTAL_SCORES 672

// minimum length product for the significance obtained from 
// #positive substitution scores to take effect
#define CUSS_PSTS_MIN_LENPRODUCT 625.0f

// referencing statistical parameters in the TDP2OutputAlnData array
#define REFSS(ndx) (ndx-nTDP2OutputPhase2_1)

// alignment data calculated as output for host
enum TDP2OutputAlnData {
    dp2oadScore,//alignment score
    dp2oadOrgProNo,//original profile serial number
    dp2oadProNewDst,//profile distance in the new buffer composed during phase 2
    dp2oadPstvs,//number of positive position-specific matches in the alignment (phase 1-2)
    dp2oadBegCoords,//beginning coordinates of the phase-1-2 alignment
    dp2oadEndCoords,//end coordinates of the phase-1-2 alignment
    nTDP2OutputPhase1,//NOTE:it is NOT allowed to exceed dpdssDiagM * nTDPDiagScoreSubsections!
    // --- phase 2 ---
    dp2oadLmbdEst=nTDP2OutputPhase1,//estimated gapped lambda based on estimates of ungapped configuration
    dp2oadKEst,//estimated K based on estimates of ungapped configuration
    dp2oadEvalue,//(log of) e-value
    dp2oadEpA,//e-value for the pair of profiles/alignment length on exit of phase 2
    nTDP2OutputPhase2_1,
    //
    dp2oadMin=nTDP2OutputPhase2_1,//minimum substitution score
    dp2oadMax,//maximum substitution score
    dp2oadE,//expected value of scores
    dp2oadH,//relative entropy of scores
    dp2oadLmbd,//calculated stat. param. lambda
    dp2oadK,//calculated K
    nTDP2OutputPhase2_2,
    //
    dp2oadIdnts=nTDP2OutputPhase2_2,//number of identical amino acids in the alignment (phase 2)
    dp2oadNGaps,//number of gaps in the alignment (phase 2)
    nTDP2OutputAlnData
};

// sections for the alignment itself
enum TDP2OutputAlignment {
    dp2oaQuery,//query alignment sequence
    dp2oaMiddle,//middle (information) sequence of the alignment 
    dp2oaTarget,//target alignment sequence
    nTDP2OutputAlignment,
    dp2oaQuerySSS=nTDP2OutputAlignment,//query SS state sequence
    dp2oaTargetSSS,//target SS state sequence
    nTDP2OutputAlignmentSSS
};

#endif//__CuBatchSS_com_h__
