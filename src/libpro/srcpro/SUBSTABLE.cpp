/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
// #include <math.h>
#include <cmath>
#include <ctype.h>

#include "extsp/psl.h"
#include "liblib/msg.h"
#include "SUBSTABLE.h"

_TCOMPUTED_BLOSUM80 COMPUTED_BLOSUM80;
_TCOMPUTED_BLOSUM62 COMPUTED_BLOSUM62;
_TCOMPUTED_BLOSUM45 COMPUTED_BLOSUM45;
_TCOMPUTED_PSCORES_ COMPUTED_PSCORES_;
_TCOMPUTED_GONNET_  COMPUTED_GONNET;

SUBSTABLE   STABLE;

// =========================================================================
//
void SetSTABLE( const mystring& value, const char* )
{
    mystring msg = "Reading ";
    mystring lcval = value;
    lcval.lower();
    if( lcval == "blosum80")
        STABLE.SetClass( SUBSTABLE::Class80 );
    else if( lcval == "blosum62")
        STABLE.SetClass( SUBSTABLE::Class62 );
    else if( lcval == "blosum45")
        STABLE.SetClass( SUBSTABLE::Class45 );
    else if( lcval == "gonnet")
        STABLE.SetClass( SUBSTABLE::ClassGn );
    else if( lcval == "pscores") {
        msg = "PSCORES not included in the implementation. Using default matrix.";
        message( msg.c_str());
    }
    else {
        warning("Unknown substitution matrix. Using default.");
    }
}

// =========================================================================
// SUBSTABLE: constructor
//
SUBSTABLE::SUBSTABLE()
:   auxprobabs_( NULL )
{
    SetClass( ClassGn );
    ComputeLogProbabilities();
    InitLogProbabs( logprobs_1_ );
    InitLogProbabs( logprobs_2_ );
    InitLogProbabs( auxlogps_ );
}

// StoreProbabilities_1: store probabilities of the first profile
void SUBSTABLE::StoreProbabilities_1( const float* probs )
{
    probs_1_ = probs;
    ComputeLogProbabilities( probs_1_, logprobs_1_ );
}
// RestoreProbabilities: restore probabilities
void SUBSTABLE::RestoreProbabilities_1()
{
    probs_1_ = NULL;
    InitLogProbabs( logprobs_1_ );
}
// StoreProbabilities_2: store probabilities of the second profile
void SUBSTABLE::StoreProbabilities_2( const float* probs )
{
    probs_2_ = probs;
    ComputeLogProbabilities( probs_2_, logprobs_2_ );
}
// RestoreProbabilities: restore probabilities
void SUBSTABLE::RestoreProbabilities_2()
{
    probs_2_ = NULL;
    InitLogProbabs( logprobs_2_ );
}
// StoreProbabilities: temporalily store probabilities
void SUBSTABLE::StoreProbabilities( const float* probs )
{
    auxprobabs_ = probs;
    ComputeLogProbabilities( auxprobabs_, auxlogps_ );
}
// RestoreProbabilities: restore probabilities
void SUBSTABLE::RestoreProbabilities()
{
    auxprobabs_ = NULL;
    InitLogProbabs( auxlogps_ );
}

// InitLogProbabs: reset logs of tempralily stored probabilities
//
void SUBSTABLE::InitLogProbabs( float* logprobs )
{
    if( !logprobs )
        return;
    for( int a = 0; a < NUMAA; a++ )
        logprobs[a] = ( float )LOG_PROB_MIN;
}

// PROBABility: background probability
//
float SUBSTABLE::PROBABility( int a )
{
#if 1//def __DEBUG__
    if( a < 0 || NUMAA <= a )
        return 0.0;
        //throw MYRUNTIME_ERROR("SUBSTABLE::PROBABility: Memory access error.");
#endif
    if( GetProbabilities())
        return GetProbabilities()[a];

    switch( GetClass()) {
        case Class80:
        case Class62:
        case Class45: return Robinson_PROBS[a];
        case ClassGn: return Robinson_PROBS[a];//Gonnet_PROBS[a];
        case ClassPS: ;
        default:
            throw MYRUNTIME_ERROR("SUBSTABLE::PROBABility: Invalid substitution matrix.");
    };
    return 0.0f;
}

// LogPROBABility: log of background probability
//
float SUBSTABLE::LogPROBABility( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a )
        throw MYRUNTIME_ERROR("SUBSTABLE::LogPROBABility: Memory access error.");
#endif
    if( GetProbabilities())
        return auxlogps_[a];
    return logprobs_[a];
}

// PROBABILITY_1: background probability of the first profile
float SUBSTABLE::PROBABILITY_1( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a )
        throw MYRUNTIME_ERROR("SUBSTABLE::PROBABILITY_1: Memory access error.");
#endif
    if( GetProbabilities_1())
        return probs_1_[a];
    return PROBABility( a );
}
// LogPROBABILITY_1: log of background probability of the first profile
float SUBSTABLE::LogPROBABILITY_1( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a )
        throw MYRUNTIME_ERROR("SUBSTABLE::LogPROBABILITY_1: Memory access error.");
#endif
    if( GetProbabilities_1())
        return logprobs_1_[a];
    return LogPROBABility( a );
}

// PROBABILITY_2: background probability of the second profile
float SUBSTABLE::PROBABILITY_2( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a )
        throw MYRUNTIME_ERROR("SUBSTABLE::PROBABILITY_2: Memory access error.");
#endif
    if( GetProbabilities_2())
        return probs_2_[a];
    return PROBABility( a );
}
// LogPROBABILITY_2: log of background probability of the second profile
float SUBSTABLE::LogPROBABILITY_2( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a )
        throw MYRUNTIME_ERROR("SUBSTABLE::LogPROBABILITY_2: Memory access error.");
#endif
    if( GetProbabilities_2())
        return logprobs_2_[a];
    return LogPROBABility( a );
}

// ComputeLogProbabilities: recompute the logs of probabilities
//
void SUBSTABLE::ComputeLogProbabilities()
{
    for( int a = 0; a < NUMAA; a++ ) {
        if( 0.0 < PROBABility( a ))
            logprobs_[a] = logf( PROBABility(a));
        else
            logprobs_[a] = ( float )LOG_PROB_MIN;
    }
}

// ComputeLogProbabilities: compute the logs of tempralily stored probabilities
//
void SUBSTABLE::ComputeLogProbabilities( const float* probs, float* logprobs )
{
    if( !probs || !logprobs )
        throw MYRUNTIME_ERROR("SUBSTABLE::ComputeLogProbabilities: Memory access error.");

    for( int a = 0; a < NUMAA; a++ ) {
        if( 0.0 < probs[a] )
            logprobs[a] = logf( probs[a] );
        else
            logprobs[a] = (float)LOG_PROB_MIN;
    }
}

// FreqRatio: return appropriate frequency ratio entry related to the 
//  subs. matrix used
//
float SUBSTABLE::FreqRatio( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a ||
        b < 0 || NUMAA <= b )
        throw MYRUNTIME_ERROR("SUBSTABLE::FreqRatio: Memory access error.");
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "SUBSTABLE: Restore probabilities first before calling FreqRatio." ));

    switch( GetClass()) {
        case Class80: return BLOSUM80_FREQRATIOS[a][b];
        case Class62: return BLOSUM62_FREQRATIOS[a][b];
        case Class45: return BLOSUM45_FREQRATIOS[a][b];
        case ClassGn: return GONNET_FREQRATIOS[a][b];
        case ClassPS: ;
        default:
            throw MYRUNTIME_ERROR("SUBSTABLE::FreqRatio: Invalid substitution matrix.");
    };
    return -1.0;
}

// PrecomputedEntry: return precomputed entry of the corresponding table
//
float SUBSTABLE::PrecomputedEntry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a ||
        b < 0 || NUMAA <= b )
        throw MYRUNTIME_ERROR("SUBSTABLE::PrecomputedEntry: Memory access error.");
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "SUBSTABLE: Restore probabilities first before calling PrecomputedEntry." ));

    switch( GetClass()) {
        case Class80: return COMPUTED_BLOSUM80( a, b );
        case Class62: return COMPUTED_BLOSUM62( a, b );
        case Class45: return COMPUTED_BLOSUM45( a, b );
        case ClassGn: return COMPUTED_GONNET( a, b );
        case ClassPS: ;
        default:
            throw MYRUNTIME_ERROR("SUBSTABLE::PrecomputedEntry: Invalid substitution matrix.");
    };
    return (float)SCORE_MIN;
}

// Entry: return entry rounded to the nearest integer
//
int SUBSTABLE::Entry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMAA <= a ||
        b < 0 || NUMAA <= b )
        throw MYRUNTIME_ERROR("SUBSTABLE::Entry: Memory access error.");
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "SUBSTABLE: Restore probabilities first before calling Entry." ));

    switch( GetClass()) {
        case Class80: return BLOSUM80[a][b];
        case Class62: return BLOSUM62[a][b];
        case Class45: return BLOSUM45[a][b];
        case ClassGn: throw MYRUNTIME_ERROR("SUBSTABLE::Entry: Precomputed values should be used for Gonnet matrix.");
        case ClassPS: ;
        default:
            throw MYRUNTIME_ERROR("SUBSTABLE::Entry: Invalid substitution matrix.");
    };
    return SCORE_MIN;
}

// StatisParam: return a statistical parameter of the substitution matrix, 
//  indexed by gap scheme and field 
//
float SUBSTABLE::StatisParam( int scheme, int field )
{
#ifdef __DEBUG__
    if( field < 0 || NumFields <= field )
        throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Memory access error.");
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "SUBSTABLE: Restore probabilities first before calling StatisParam." ));

    switch( GetClass()) {
        case Class80:
            if( scheme < 0 || Num80Entries <= scheme )
                throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Memory access error.");
            return BLOSUM80_VALUES[scheme][field];
        case Class62:
            if( scheme < 0 || NumEntries <= scheme )
                throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Memory access error.");
            return BLOSUM62_VALUES[scheme][field];
        case Class45:
            if( scheme < 0 || Num45Entries <= scheme )
                throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Memory access error.");
            return BLOSUM45_VALUES[scheme][field];
        case ClassGn:
            if( scheme < 0 || NumGnEntries <= scheme )
                throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Memory access error.");
            return GONNET_VALUES[scheme][field];
        case ClassPS: ;
        default:
            throw MYRUNTIME_ERROR("SUBSTABLE::StatisParam: Invalid substitution matrix.");
    };
    return -1.0;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM80: precompute BLOSUM80 values given frequency ratios
//
_TCOMPUTED_BLOSUM80::_TCOMPUTED_BLOSUM80()
{
    for( int a = 0; a < NUMAA; a++ )
        for( int b = 0; b < NUMAA; b++ )
            if( BLOSUM80_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = logf( BLOSUM80_FREQRATIOS[a][b] ) / SLC_LN2 * Blosum80ScalingConstant;
            else
                data[a][b] = (float)SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM62: precompute BLOSUM62 values given frequency ratios
//
_TCOMPUTED_BLOSUM62::_TCOMPUTED_BLOSUM62()
{
    for( int a = 0; a < NUMAA; a++ )
        for( int b = 0; b < NUMAA; b++ )
            if( BLOSUM62_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = logf( BLOSUM62_FREQRATIOS[a][b] ) / SLC_LN2 * Blosum62ScalingConstant;
            else
                data[a][b] = (float)SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM45: precompute BLOSUM45 values given frequency ratios
//
_TCOMPUTED_BLOSUM45::_TCOMPUTED_BLOSUM45()
{
    for( int a = 0; a < NUMAA; a++ )
        for( int b = 0; b < NUMAA; b++ )
            if( BLOSUM45_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = logf( BLOSUM45_FREQRATIOS[a][b] ) / SLC_LN2 * Blosum45ScalingConstant;
            else
                data[a][b] = (float)SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_GONNET_: precompute scaled Gonnet matrix values given
//  frequency ratios
//
_TCOMPUTED_GONNET_::_TCOMPUTED_GONNET_()
{
    int a, b;
    for( a = 0; a < NUMAA; a++ )
        for( b = 0; b < NUMAA; b++ )
            if( 0.0 < GONNET_FREQRATIOS[a][b])
                data[a][b] = logf( GONNET_FREQRATIOS[a][b] ) / SLC_LN2 * Gonnet_ScalingConstant;
            else
                data[a][b] = (float)SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_PSCORES_: PSCORES_ constructor: initialize scores
//
_TCOMPUTED_PSCORES_::_TCOMPUTED_PSCORES_()
{
    for( int a = 0; a < NUMAA; a++ )
        for( int b = 0; b < NUMAA; b++ )
            data[a][b] = (float)SCORE_MIN;
}

// Compute: recomputes scores by using COUNTS information
//
void _TCOMPUTED_PSCORES_::Compute()
{
    throw MYRUNTIME_ERROR("_TCOMPUTED_PSCORES_::Compute: Invalid substitution matrix.");
}

// =========================================================================
// Constants ---------------------------------------------------------------
//
// AB Robinson and LR Robinson background probabilities
// (PNAS USA 88, 1991, 8880-4)
const float Robinson_PROBS[NUMAA] = {
//   A          R          N          D          C          Q          E          G          H          I 
//   L          K          M          F          P          S          T          W          Y          V 
  0.078050f, 0.051290f, 0.044870f, 0.053640f, 0.019250f, 0.042640f, 0.062950f, 0.073770f, 0.021990f, 0.051420f,
  0.090190f, 0.057440f, 0.022430f, 0.038560f, 0.052030f, 0.071200f, 0.058410f, 0.013300f, 0.032160f, 0.064410f
};

// Background probabilities by pscores
const float Pscores_PROBS[NUMAA] = {
//    A          R          N          D          C          Q          E          G          H          I 
//    L          K          M          F          P          S          T          W          Y          V 
  0.094028f, 0.053829f, 0.039187f, 0.060046f, 0.013052f, 0.032910f, 0.060580f, 0.079290f, 0.024236f, 0.062608f,
  0.105052f, 0.045876f, 0.019318f, 0.037945f, 0.043358f, 0.060625f, 0.051353f, 0.012555f, 0.029567f, 0.074585f
};

// Frequency ratios for BLOSUM80
const float BLOSUM80_FREQRATIOS[NUMAA][NUMAA] = {
// A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V
{4.773f,  0.555f,  0.510f,  0.451f,  0.732f,  0.696f,  0.703f,  0.957f,  0.514f,  0.543f,  0.505f,  0.723f,  0.625f,  0.397f,  0.771f,  1.535f,  0.980f,  0.309f,  0.436f,  0.866f}, // A
{0.555f,  8.245f,  0.773f,  0.477f,  0.233f,  1.394f,  0.832f,  0.377f,  0.925f,  0.299f,  0.363f,  2.192f,  0.506f,  0.287f,  0.446f,  0.695f,  0.598f,  0.294f,  0.418f,  0.354f}, // R
{0.510f,  0.773f,  8.963f,  1.584f,  0.300f,  0.958f,  0.811f,  0.761f,  1.124f,  0.258f,  0.250f,  0.938f,  0.382f,  0.273f,  0.398f,  1.165f,  0.908f,  0.221f,  0.385f,  0.297f}, // N
{0.451f,  0.477f,  1.584f,  9.106f,  0.214f,  0.763f,  1.635f,  0.541f,  0.594f,  0.214f,  0.197f,  0.677f,  0.252f,  0.234f,  0.452f,  0.774f,  0.611f,  0.145f,  0.245f,  0.245f}, // D
{0.732f,  0.233f,  0.300f,  0.214f, 20.702f,  0.295f,  0.180f,  0.272f,  0.221f,  0.581f,  0.493f,  0.241f,  0.499f,  0.395f,  0.269f,  0.576f,  0.602f,  0.302f,  0.308f,  0.634f}, // C
{0.696f,  1.394f,  0.958f,  0.763f,  0.295f,  8.340f,  1.906f,  0.425f,  1.316f,  0.309f,  0.407f,  1.524f,  0.887f,  0.285f,  0.538f,  0.859f,  0.724f,  0.408f,  0.462f,  0.411f}, // Q
{0.703f,  0.832f,  0.811f,  1.635f,  0.180f,  1.906f,  6.995f,  0.399f,  0.901f,  0.264f,  0.276f,  1.195f,  0.429f,  0.249f,  0.581f,  0.845f,  0.685f,  0.241f,  0.333f,  0.369f}, // E
{0.957f,  0.377f,  0.761f,  0.541f,  0.272f,  0.425f,  0.399f,  7.882f,  0.387f,  0.184f,  0.210f,  0.483f,  0.286f,  0.249f,  0.347f,  0.784f,  0.492f,  0.264f,  0.230f,  0.251f}, // G
{0.514f,  0.925f,  1.124f,  0.594f,  0.221f,  1.316f,  0.901f,  0.387f, 16.070f,  0.258f,  0.314f,  0.740f,  0.432f,  0.572f,  0.420f,  0.661f,  0.540f,  0.390f,  1.819f,  0.289f}, // H
{0.543f,  0.299f,  0.258f,  0.214f,  0.581f,  0.309f,  0.264f,  0.184f,  0.258f,  4.868f,  1.665f,  0.313f,  1.512f,  0.841f,  0.286f,  0.379f,  0.701f,  0.343f,  0.539f,  2.496f}, // I
{0.505f,  0.363f,  0.250f,  0.197f,  0.493f,  0.407f,  0.276f,  0.210f,  0.314f,  1.665f,  4.463f,  0.357f,  2.123f,  1.114f,  0.303f,  0.368f,  0.561f,  0.439f,  0.581f,  1.220f}, // L
{0.723f,  2.192f,  0.938f,  0.677f,  0.241f,  1.524f,  1.195f,  0.483f,  0.740f,  0.313f,  0.357f,  6.326f,  0.534f,  0.283f,  0.597f,  0.820f,  0.736f,  0.241f,  0.408f,  0.370f}, // K
{0.625f,  0.506f,  0.382f,  0.252f,  0.499f,  0.887f,  0.429f,  0.286f,  0.432f,  1.512f,  2.123f,  0.534f,  8.883f,  0.893f,  0.362f,  0.498f,  0.758f,  0.561f,  0.550f,  1.224f}, // M
{0.397f,  0.287f,  0.273f,  0.234f,  0.395f,  0.285f,  0.249f,  0.249f,  0.572f,  0.841f,  1.114f,  0.283f,  0.893f,  9.486f,  0.237f,  0.369f,  0.445f,  1.089f,  2.780f,  0.649f}, // F
{0.771f,  0.446f,  0.398f,  0.452f,  0.269f,  0.538f,  0.581f,  0.347f,  0.420f,  0.286f,  0.303f,  0.597f,  0.362f,  0.237f, 15.155f,  0.652f,  0.560f,  0.178f,  0.258f,  0.370f}, // P
{1.535f,  0.695f,  1.165f,  0.774f,  0.576f,  0.859f,  0.845f,  0.784f,  0.661f,  0.379f,  0.368f,  0.820f,  0.498f,  0.369f,  0.652f,  5.106f,  1.663f,  0.271f,  0.462f,  0.494f}, // S
{0.980f,  0.598f,  0.908f,  0.611f,  0.602f,  0.724f,  0.685f,  0.492f,  0.540f,  0.701f,  0.561f,  0.736f,  0.758f,  0.445f,  0.560f,  1.663f,  6.205f,  0.285f,  0.474f,  0.891f}, // T
{0.309f,  0.294f,  0.221f,  0.145f,  0.302f,  0.408f,  0.241f,  0.264f,  0.390f,  0.343f,  0.439f,  0.241f,  0.561f,  1.089f,  0.178f,  0.271f,  0.285f, 41.552f,  2.036f,  0.342f}, // W
{0.436f,  0.418f,  0.385f,  0.245f,  0.308f,  0.462f,  0.333f,  0.230f,  1.819f,  0.539f,  0.581f,  0.408f,  0.550f,  2.780f,  0.258f,  0.462f,  0.474f,  2.036f, 12.194f,  0.489f}, // Y
{0.866f,  0.354f,  0.297f,  0.245f,  0.634f,  0.411f,  0.369f,  0.251f,  0.289f,  2.496f,  1.220f,  0.370f,  1.224f,  0.649f,  0.370f,  0.494f,  0.891f,  0.342f,  0.489f,  4.584f}, // V
};

// Frequency ratios for BLOSUM62 as determined by Stephen Altschul;
//  Stephen and Jorja Henikoff used different number for B, Z, X.
//  Each entry in the table equals to substitution frequency qij
//  devided by the product of background probabilities pi * pj
// const double BLOSUM62_FREQRATIOS[ NUMALPH ][ NUMALPH ] = {
// A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V       B       Z       X       *      -
// {3.903,  0.613,  0.588,  0.545,  0.868,  0.757,  0.741,  1.057,  0.569,  0.632,  0.602,  0.775,  0.723,  0.465,  0.754,  1.472,  0.984,  0.416,  0.543,  0.936,  0.565,  0.747,  0.750,  0.250,  0.000}, // A
// {0.613,  6.666,  0.859,  0.573,  0.309,  1.406,  0.961,  0.450,  0.917,  0.355,  0.474,  2.077,  0.623,  0.381,  0.481,  0.767,  0.678,  0.395,  0.556,  0.420,  0.703,  1.133,  0.750,  0.250,  0.000}, // R
// {0.588,  0.859,  7.094,  1.554,  0.398,  1.001,  0.911,  0.864,  1.222,  0.328,  0.310,  0.940,  0.474,  0.354,  0.500,  1.232,  0.984,  0.278,  0.486,  0.369,  4.071,  0.946,  0.750,  0.250,  0.000}, // N
// {0.545,  0.573,  1.554,  7.398,  0.301,  0.897,  1.688,  0.634,  0.679,  0.339,  0.287,  0.784,  0.346,  0.299,  0.599,  0.913,  0.695,  0.232,  0.346,  0.337,  4.743,  1.382,  0.750,  0.250,  0.000}, // D
// {0.868,  0.309,  0.398,  0.301, 19.577,  0.366,  0.286,  0.420,  0.355,  0.653,  0.642,  0.349,  0.611,  0.439,  0.380,  0.738,  0.741,  0.450,  0.434,  0.756,  0.345,  0.317,  0.750,  0.250,  0.000}, // C
// {0.757,  1.406,  1.001,  0.897,  0.366,  6.244,  1.902,  0.539,  1.168,  0.383,  0.477,  1.554,  0.864,  0.334,  0.641,  0.966,  0.791,  0.509,  0.611,  0.467,  0.944,  3.582,  0.750,  0.250,  0.000}, // Q
// {0.741,  0.961,  0.911,  1.688,  0.286,  1.902,  5.470,  0.481,  0.960,  0.331,  0.373,  1.308,  0.500,  0.331,  0.679,  0.950,  0.741,  0.374,  0.496,  0.429,  1.335,  4.090,  0.750,  0.250,  0.000}, // E
// {1.057,  0.450,  0.864,  0.634,  0.420,  0.539,  0.481,  6.876,  0.493,  0.275,  0.284,  0.589,  0.396,  0.341,  0.477,  0.904,  0.579,  0.422,  0.349,  0.337,  0.739,  0.503,  0.750,  0.250,  0.000}, // G
// {0.569,  0.917,  1.222,  0.679,  0.355,  1.168,  0.960,  0.493, 13.506,  0.326,  0.381,  0.779,  0.584,  0.652,  0.473,  0.737,  0.557,  0.444,  1.798,  0.339,  0.925,  1.040,  0.750,  0.250,  0.000}, // H
// {0.632,  0.355,  0.328,  0.339,  0.653,  0.383,  0.331,  0.275,  0.326,  3.998,  1.694,  0.396,  1.478,  0.946,  0.385,  0.443,  0.780,  0.409,  0.630,  2.417,  0.334,  0.351,  0.750,  0.250,  0.000}, // I
// {0.602,  0.474,  0.310,  0.287,  0.642,  0.477,  0.373,  0.284,  0.381,  1.694,  3.797,  0.428,  1.994,  1.155,  0.371,  0.429,  0.660,  0.568,  0.692,  1.314,  0.297,  0.413,  0.750,  0.250,  0.000}, // L
// {0.775,  2.077,  0.940,  0.784,  0.349,  1.554,  1.308,  0.589,  0.779,  0.396,  0.428,  4.764,  0.625,  0.344,  0.704,  0.932,  0.793,  0.359,  0.532,  0.457,  0.855,  1.403,  0.750,  0.250,  0.000}, // K
// {0.723,  0.623,  0.474,  0.346,  0.611,  0.864,  0.500,  0.396,  0.584,  1.478,  1.994,  0.625,  6.481,  1.004,  0.424,  0.599,  0.794,  0.610,  0.708,  1.269,  0.405,  0.641,  0.750,  0.250,  0.000}, // M
// {0.465,  0.381,  0.354,  0.299,  0.439,  0.334,  0.331,  0.341,  0.652,  0.946,  1.155,  0.344,  1.004,  8.129,  0.287,  0.440,  0.482,  1.374,  2.769,  0.745,  0.324,  0.332,  0.750,  0.250,  0.000}, // F
// {0.754,  0.481,  0.500,  0.599,  0.380,  0.641,  0.679,  0.477,  0.473,  0.385,  0.371,  0.704,  0.424,  0.287, 12.838,  0.755,  0.689,  0.282,  0.363,  0.443,  0.554,  0.664,  0.750,  0.250,  0.000}, // P
// {1.472,  0.767,  1.232,  0.913,  0.738,  0.966,  0.950,  0.904,  0.737,  0.443,  0.429,  0.932,  0.599,  0.440,  0.755,  3.843,  1.614,  0.385,  0.557,  0.565,  1.058,  0.956,  0.750,  0.250,  0.000}, // S
// {0.984,  0.678,  0.984,  0.695,  0.741,  0.791,  0.741,  0.579,  0.557,  0.780,  0.660,  0.793,  0.794,  0.482,  0.689,  1.614,  4.832,  0.431,  0.573,  0.981,  0.826,  0.761,  0.750,  0.250,  0.000}, // T
// {0.416,  0.395,  0.278,  0.232,  0.450,  0.509,  0.374,  0.422,  0.444,  0.409,  0.568,  0.359,  0.610,  1.374,  0.282,  0.385,  0.431, 38.108,  2.110,  0.374,  0.253,  0.426,  0.750,  0.250,  0.000}, // W
// {0.543,  0.556,  0.486,  0.346,  0.434,  0.611,  0.496,  0.349,  1.798,  0.630,  0.692,  0.532,  0.708,  2.769,  0.363,  0.557,  0.573,  2.110,  9.832,  0.658,  0.409,  0.541,  0.750,  0.250,  0.000}, // Y
// {0.936,  0.420,  0.369,  0.337,  0.756,  0.467,  0.429,  0.337,  0.339,  2.417,  1.314,  0.457,  1.269,  0.745,  0.443,  0.565,  0.981,  0.374,  0.658,  3.692,  0.351,  0.444,  0.750,  0.250,  0.000}, // V
// {0.565,  0.703,  4.071,  4.743,  0.345,  0.944,  1.335,  0.739,  0.925,  0.334,  0.297,  0.855,  0.405,  0.324,  0.554,  1.058,  0.826,  0.253,  0.409,  0.351,  4.438,  1.184,  0.750,  0.250,  0.000}, // B
// {0.747,  1.133,  0.946,  1.382,  0.317,  3.582,  4.090,  0.503,  1.040,  0.351,  0.413,  1.403,  0.641,  0.332,  0.664,  0.956,  0.761,  0.426,  0.541,  0.444,  1.184,  3.893,  0.750,  0.250,  0.000}, // Z
// {0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.750,  0.250,  0.000}, // X
// {0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  1.333,  0.250}, // *
// {0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.250,  0.000}, // -
// };

// Frequency ratios for BLOSUM62 computed purely from BLOSUM62 scores;
//  Each entry in the table equals to target probability qij
//  devided by the product of marginal probabilities pi * pj
//  (computed originaly)
const float BLOSUM62_FREQRATIOS[NUMAA][NUMAA] = {
// A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V
{3.9030f, 0.6127f, 0.5883f, 0.5446f, 0.8680f, 0.7568f, 0.7413f, 1.0569f, 0.5694f, 0.6325f, 0.6019f, 0.7754f, 0.7232f, 0.4649f, 0.7541f, 1.4721f, 0.9844f, 0.4166f, 0.5426f, 0.9365f}, // A
{0.6127f, 6.6654f, 0.8586f, 0.5732f, 0.3089f, 1.4058f, 0.9608f, 0.4500f, 0.9171f, 0.3548f, 0.4739f, 2.0769f, 0.6226f, 0.3807f, 0.4815f, 0.7672f, 0.6777f, 0.3951f, 0.5560f, 0.4201f}, // R
{0.5883f, 0.8586f, 7.0940f, 1.5538f, 0.3978f, 1.0006f, 0.9113f, 0.8637f, 1.2220f, 0.3279f, 0.3100f, 0.9398f, 0.4745f, 0.3543f, 0.4999f, 1.2316f, 0.9842f, 0.2778f, 0.4860f, 0.3690f}, // N
{0.5446f, 0.5732f, 1.5538f, 7.3978f, 0.3015f, 0.8971f, 1.6879f, 0.6343f, 0.6785f, 0.3390f, 0.2866f, 0.7841f, 0.3464f, 0.2990f, 0.5987f, 0.9135f, 0.6948f, 0.2321f, 0.3457f, 0.3365f}, // D
{0.8680f, 0.3089f, 0.3978f, 0.3015f,19.5772f, 0.3658f, 0.2859f, 0.4204f, 0.3551f, 0.6535f, 0.6423f, 0.3491f, 0.6114f, 0.4390f, 0.3796f, 0.7384f, 0.7406f, 0.4500f, 0.4342f, 0.7559f}, // C
{0.7568f, 1.4058f, 1.0006f, 0.8971f, 0.3658f, 6.2446f, 1.9017f, 0.5386f, 1.1680f, 0.3829f, 0.4773f, 1.5543f, 0.8642f, 0.3340f, 0.6413f, 0.9655f, 0.7913f, 0.5094f, 0.6111f, 0.4668f}, // Q
{0.7413f, 0.9608f, 0.9113f, 1.6879f, 0.2859f, 1.9017f, 5.4695f, 0.4813f, 0.9601f, 0.3305f, 0.3729f, 1.3083f, 0.5003f, 0.3307f, 0.6792f, 0.9503f, 0.7414f, 0.3743f, 0.4965f, 0.4290f}, // E
{1.0569f, 0.4500f, 0.8637f, 0.6343f, 0.4204f, 0.5386f, 0.4813f, 6.8761f, 0.4930f, 0.2750f, 0.2845f, 0.5889f, 0.3955f, 0.3406f, 0.4774f, 0.9036f, 0.5793f, 0.4217f, 0.3487f, 0.3369f}, // G
{0.5694f, 0.9171f, 1.2220f, 0.6785f, 0.3551f, 1.1680f, 0.9601f, 0.4930f,13.5057f, 0.3263f, 0.3807f, 0.7789f, 0.5841f, 0.6520f, 0.4729f, 0.7367f, 0.5575f, 0.4441f, 1.7979f, 0.3395f}, // H
{0.6325f, 0.3548f, 0.3279f, 0.3390f, 0.6535f, 0.3829f, 0.3305f, 0.2750f, 0.3263f, 3.9981f, 1.6944f, 0.3964f, 1.4777f, 0.9458f, 0.3847f, 0.4432f, 0.7798f, 0.4089f, 0.6304f, 2.4175f}, // I
{0.6019f, 0.4739f, 0.3100f, 0.2866f, 0.6423f, 0.4773f, 0.3729f, 0.2845f, 0.3807f, 1.6944f, 3.7966f, 0.4283f, 1.9943f, 1.1546f, 0.3711f, 0.4289f, 0.6603f, 0.5681f, 0.6921f, 1.3142f}, // L
{0.7754f, 2.0769f, 0.9398f, 0.7841f, 0.3491f, 1.5543f, 1.3083f, 0.5889f, 0.7789f, 0.3964f, 0.4283f, 4.7644f, 0.6253f, 0.3441f, 0.7038f, 0.9319f, 0.7929f, 0.3589f, 0.5322f, 0.4565f}, // K
{0.7232f, 0.6226f, 0.4745f, 0.3464f, 0.6114f, 0.8642f, 0.5003f, 0.3955f, 0.5841f, 1.4777f, 1.9943f, 0.6253f, 6.4814f, 1.0044f, 0.4239f, 0.5986f, 0.7938f, 0.6103f, 0.7084f, 1.2689f}, // M
{0.4649f, 0.3807f, 0.3543f, 0.2990f, 0.4390f, 0.3340f, 0.3307f, 0.3406f, 0.6520f, 0.9458f, 1.1546f, 0.3441f, 1.0044f, 8.1286f, 0.2875f, 0.4400f, 0.4817f, 1.3744f, 2.7695f, 0.7451f}, // F
{0.7541f, 0.4815f, 0.4999f, 0.5987f, 0.3796f, 0.6413f, 0.6792f, 0.4774f, 0.4729f, 0.3847f, 0.3711f, 0.7038f, 0.4239f, 0.2875f,12.8376f, 0.7555f, 0.6889f, 0.2818f, 0.3635f, 0.4431f}, // P
{1.4721f, 0.7672f, 1.2316f, 0.9135f, 0.7384f, 0.9655f, 0.9503f, 0.9036f, 0.7367f, 0.4432f, 0.4289f, 0.9319f, 0.5986f, 0.4400f, 0.7555f, 3.8429f, 1.6140f, 0.3853f, 0.5575f, 0.5652f}, // S
{0.9844f, 0.6777f, 0.9842f, 0.6948f, 0.7406f, 0.7913f, 0.7414f, 0.5793f, 0.5575f, 0.7798f, 0.6603f, 0.7929f, 0.7938f, 0.4817f, 0.6889f, 1.6140f, 4.8323f, 0.4309f, 0.5732f, 0.9809f}, // T
{0.4166f, 0.3951f, 0.2778f, 0.2321f, 0.4500f, 0.5094f, 0.3743f, 0.4217f, 0.4441f, 0.4089f, 0.5681f, 0.3589f, 0.6103f, 1.3744f, 0.2818f, 0.3853f, 0.4309f,38.1074f, 2.1098f, 0.3745f}, // W
{0.5426f, 0.5560f, 0.4860f, 0.3457f, 0.4342f, 0.6111f, 0.4965f, 0.3487f, 1.7979f, 0.6304f, 0.6921f, 0.5322f, 0.7084f, 2.7695f, 0.3635f, 0.5575f, 0.5732f, 2.1098f, 9.8321f, 0.6580f}, // Y
{0.9365f, 0.4201f, 0.3690f, 0.3365f, 0.7559f, 0.4668f, 0.4290f, 0.3369f, 0.3395f, 2.4175f, 1.3142f, 0.4565f, 1.2689f, 0.7451f, 0.4431f, 0.5652f, 0.9809f, 0.3745f, 0.6580f, 3.6922f}, // V
};

// Frequency ratios for BLOSUM45
const float BLOSUM45_FREQRATIOS[NUMAA][NUMAA] = {
// A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V
{2.950f,  0.700f,  0.789f,  0.689f,  0.800f,  0.867f,  0.825f,  1.080f,  0.654f,  0.747f,  0.712f,  0.786f,  0.821f,  0.587f,  0.709f,  1.300f,  1.001f,  0.565f,  0.639f,  1.010f}, // A
{0.700f,  4.747f,  0.898f,  0.771f,  0.472f,  1.329f,  1.011f,  0.570f,  0.973f,  0.488f,  0.601f,  1.943f,  0.776f,  0.590f,  0.582f,  0.799f,  0.715f,  0.580f,  0.807f,  0.578f}, // R
{0.789f,  0.898f,  4.478f,  1.502f,  0.680f,  1.032f,  0.893f,  1.059f,  1.220f,  0.564f,  0.484f,  1.119f,  0.682f,  0.572f,  0.640f,  1.200f,  1.115f,  0.390f,  0.606f,  0.552f}, // N
{0.689f,  0.771f,  1.502f,  5.356f,  0.533f,  0.958f,  1.643f,  0.740f,  0.976f,  0.440f,  0.463f,  0.942f,  0.494f,  0.431f,  0.724f,  0.929f,  0.876f,  0.373f,  0.645f,  0.494f}, // D
{0.800f,  0.472f,  0.680f,  0.533f, 17.090f,  0.486f,  0.545f,  0.557f,  0.491f,  0.543f,  0.673f,  0.547f,  0.604f,  0.602f,  0.411f,  0.797f,  0.822f,  0.334f,  0.489f,  0.715f}, // C
{0.867f,  1.329f,  1.032f,  0.958f,  0.486f,  4.407f,  1.531f,  0.687f,  1.151f,  0.578f,  0.642f,  1.330f,  0.941f,  0.444f,  0.716f,  1.092f,  0.781f,  0.645f,  0.829f,  0.547f}, // Q
{0.825f,  1.011f,  0.893f,  1.643f,  0.545f,  1.531f,  3.873f,  0.576f,  0.962f,  0.485f,  0.571f,  1.277f,  0.615f,  0.498f,  0.911f,  0.912f,  0.833f,  0.519f,  0.617f,  0.555f}, // E
{1.080f,  0.570f,  1.059f,  0.740f,  0.557f,  0.687f,  0.576f,  5.071f,  0.662f,  0.416f,  0.450f,  0.678f,  0.585f,  0.480f,  0.702f,  1.058f,  0.693f,  0.591f,  0.549f,  0.479f}, // G
{0.654f,  0.973f,  1.220f,  0.976f,  0.491f,  1.151f,  0.962f,  0.662f,  9.512f,  0.453f,  0.670f,  0.890f,  0.918f,  0.679f,  0.661f,  0.854f,  0.706f,  0.452f,  1.472f,  0.457f}, // H
{0.747f,  0.488f,  0.564f,  0.440f,  0.543f,  0.578f,  0.485f,  0.416f,  0.453f,  3.233f,  1.596f,  0.532f,  1.455f,  1.064f,  0.610f,  0.618f,  0.848f,  0.565f,  0.906f,  2.176f}, // I
{0.712f,  0.601f,  0.484f,  0.463f,  0.673f,  0.642f,  0.571f,  0.450f,  0.670f,  1.596f,  2.997f,  0.554f,  1.731f,  1.303f,  0.478f,  0.556f,  0.781f,  0.671f,  0.965f,  1.334f}, // L
{0.786f,  1.943f,  1.119f,  0.942f,  0.547f,  1.330f,  1.277f,  0.678f,  0.890f,  0.532f,  0.554f,  3.327f,  0.738f,  0.529f,  0.781f,  0.890f,  0.885f,  0.562f,  0.737f,  0.592f}, // K
{0.821f,  0.776f,  0.682f,  0.494f,  0.604f,  0.941f,  0.615f,  0.585f,  0.918f,  1.455f,  1.731f,  0.738f,  4.114f,  1.063f,  0.644f,  0.660f,  0.860f,  0.634f,  1.023f,  1.236f}, // M
{0.587f,  0.590f,  0.572f,  0.431f,  0.602f,  0.444f,  0.498f,  0.480f,  0.679f,  1.064f,  1.303f,  0.529f,  1.063f,  5.748f,  0.451f,  0.610f,  0.716f,  1.355f,  2.185f,  0.953f}, // F
{0.709f,  0.582f,  0.640f,  0.724f,  0.411f,  0.716f,  0.911f,  0.702f,  0.661f,  0.610f,  0.478f,  0.781f,  0.644f,  0.451f,  8.819f,  0.750f,  0.856f,  0.525f,  0.479f,  0.540f}, // P
{1.300f,  0.799f,  1.200f,  0.929f,  0.797f,  1.092f,  0.912f,  1.058f,  0.854f,  0.618f,  0.556f,  0.890f,  0.660f,  0.610f,  0.750f,  2.782f,  1.472f,  0.428f,  0.706f,  0.728f}, // S
{1.001f,  0.715f,  1.115f,  0.876f,  0.822f,  0.781f,  0.833f,  0.693f,  0.706f,  0.848f,  0.781f,  0.885f,  0.860f,  0.716f,  0.856f,  1.472f,  3.139f,  0.454f,  0.744f,  1.040f}, // T
{0.565f,  0.580f,  0.390f,  0.373f,  0.334f,  0.645f,  0.519f,  0.591f,  0.452f,  0.565f,  0.671f,  0.562f,  0.634f,  1.355f,  0.525f,  0.428f,  0.454f, 29.702f,  1.801f,  0.473f}, // W
{0.639f,  0.807f,  0.606f,  0.645f,  0.489f,  0.829f,  0.617f,  0.549f,  1.472f,  0.906f,  0.965f,  0.737f,  1.023f,  2.185f,  0.479f,  0.706f,  0.744f,  1.801f,  5.753f,  0.809f}, // Y
{1.010f,  0.578f,  0.552f,  0.494f,  0.715f,  0.547f,  0.555f,  0.479f,  0.457f,  2.176f,  1.334f,  0.592f,  1.236f,  0.953f,  0.540f,  0.728f,  1.040f,  0.473f,  0.809f,  2.871f}, // V
};


// Frequency ratios obtained by program pscores from this (or previous version of) software package
const float PSCORES_FREQRATIOS[NUMAA][NUMAA] = {
// A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V
{2.6399f, 0.7772f, 0.6046f, 0.6396f, 1.3343f, 0.8524f, 0.8788f, 0.9549f, 0.6629f, 0.7821f, 0.7160f, 0.7947f, 0.8863f, 0.6782f, 0.8406f, 1.0633f, 0.9068f, 0.5385f, 0.7350f, 1.0114f}, // A
{0.7772f, 3.2222f, 1.1060f, 0.9739f, 0.5769f, 1.7280f, 1.4248f, 0.7038f, 1.1446f, 0.5082f, 0.5734f, 2.0486f, 0.7491f, 0.5213f, 0.8564f, 0.9549f, 0.9844f, 0.6488f, 0.8224f, 0.5188f}, // R
{0.6046f, 1.1060f, 4.6090f, 1.7032f, 1.0571f, 1.1988f, 1.0683f, 1.0617f, 1.2828f, 0.3409f, 0.3881f, 1.1640f, 0.6028f, 0.5028f, 0.7141f, 1.2564f, 1.0116f, 0.4786f, 0.7474f, 0.3642f}, // N
{0.6396f, 0.9739f, 1.7032f, 4.5842f, 0.5278f, 1.2140f, 1.6612f, 0.8334f, 0.9902f, 0.2692f, 0.2782f, 1.1105f, 0.4417f, 0.3635f, 0.8345f, 1.1700f, 0.8626f, 0.3788f, 0.5435f, 0.3105f}, // D
{1.3343f, 0.5769f, 1.0571f, 0.5278f, 6.5544f, 0.6014f, 0.4600f, 0.8342f, 0.7558f, 1.0468f, 0.9810f, 0.5346f, 1.1102f, 0.9978f, 0.5809f, 1.1723f, 1.0666f, 0.7596f, 0.9339f, 1.2918f}, // C
{0.8524f, 1.7280f, 1.1988f, 1.2140f, 0.6014f, 2.5555f, 1.8075f, 0.7634f, 1.1916f, 0.4924f, 0.5738f, 1.7922f, 0.8599f, 0.5273f, 0.8579f, 1.0435f, 1.0210f, 0.5745f, 0.7825f, 0.5117f}, // Q
{0.8788f, 1.4248f, 1.0683f, 1.6612f, 0.4600f, 1.8075f, 3.1430f, 0.7491f, 0.9492f, 0.4201f, 0.4277f, 1.6446f, 0.5885f, 0.4155f, 1.0072f, 0.9867f, 0.9604f, 0.4420f, 0.5855f, 0.4509f}, // E
{0.9549f, 0.7038f, 1.0617f, 0.8334f, 0.8342f, 0.7634f, 0.7491f, 4.9690f, 0.7570f, 0.2915f, 0.2879f, 0.8177f, 0.4809f, 0.4218f, 0.7177f, 1.0752f, 0.7110f, 0.3942f, 0.4658f, 0.3672f}, // G
{0.6629f, 1.1446f, 1.2828f, 0.9902f, 0.7558f, 1.1916f, 0.9492f, 0.7570f, 6.8668f, 0.5083f, 0.5619f, 1.0191f, 1.0884f, 0.9062f, 0.8053f, 1.0680f, 0.8575f, 1.0543f, 1.5969f, 0.5232f}, // H
{0.7821f, 0.5082f, 0.3409f, 0.2692f, 1.0468f, 0.4924f, 0.4201f, 0.2915f, 0.5083f, 2.8700f, 1.8051f, 0.4448f, 1.4890f, 1.4123f, 0.5529f, 0.4636f, 0.7253f, 0.9336f, 0.8864f, 2.2487f}, // I
{0.7160f, 0.5734f, 0.3881f, 0.2782f, 0.9810f, 0.5738f, 0.4277f, 0.2879f, 0.5619f, 1.8051f, 2.8619f, 0.4879f, 1.6867f, 1.5207f, 0.4791f, 0.4363f, 0.6036f, 0.9125f, 1.1057f, 1.3686f}, // L
{0.7947f, 2.0486f, 1.1640f, 1.1105f, 0.5346f, 1.7922f, 1.6446f, 0.8177f, 1.0191f, 0.4448f, 0.4879f, 3.1664f, 0.6585f, 0.4404f, 0.9710f, 0.9900f, 1.0097f, 0.4682f, 0.6541f, 0.4636f}, // K
{0.8863f, 0.7491f, 0.6028f, 0.4417f, 1.1102f, 0.8599f, 0.5885f, 0.4809f, 1.0884f, 1.4890f, 1.6867f, 0.6585f, 3.7274f, 1.4213f, 0.5855f, 0.6614f, 0.8331f, 1.1025f, 1.1362f, 1.2605f}, // M
{0.6782f, 0.5213f, 0.5028f, 0.3635f, 0.9978f, 0.5273f, 0.4155f, 0.4218f, 0.9062f, 1.4123f, 1.5207f, 0.4404f, 1.4213f, 4.3794f, 0.5739f, 0.5469f, 0.6617f, 2.7991f, 2.4194f, 1.1852f}, // F
{0.8406f, 0.8564f, 0.7141f, 0.8345f, 0.5809f, 0.8579f, 1.0072f, 0.7177f, 0.8053f, 0.5529f, 0.4791f, 0.9710f, 0.5855f, 0.5739f, 7.2858f, 0.9525f, 0.8545f, 0.5244f, 0.5648f, 0.6290f}, // P
{1.0633f, 0.9549f, 1.2564f, 1.1700f, 1.1723f, 1.0435f, 0.9867f, 1.0752f, 1.0680f, 0.4636f, 0.4363f, 0.9900f, 0.6614f, 0.5469f, 0.9525f, 2.6125f, 1.6190f, 0.5226f, 0.6688f, 0.5819f}, // S
{0.9068f, 0.9844f, 1.0116f, 0.8626f, 1.0666f, 1.0210f, 0.9604f, 0.7110f, 0.8575f, 0.7253f, 0.6036f, 1.0097f, 0.8331f, 0.6617f, 0.8545f, 1.6190f, 3.1768f, 0.6107f, 0.7246f, 0.9277f}, // T
{0.5385f, 0.6488f, 0.4786f, 0.3788f, 0.7596f, 0.5745f, 0.4420f, 0.3942f, 1.0543f, 0.9336f, 0.9125f, 0.4682f, 1.1025f, 2.7991f, 0.5244f, 0.5226f, 0.6107f,16.6611f, 2.9300f, 0.8216f}, // W
{0.7350f, 0.8224f, 0.7474f, 0.5435f, 0.9339f, 0.7825f, 0.5855f, 0.4658f, 1.5969f, 0.8864f, 1.1057f, 0.6541f, 1.1362f, 2.4194f, 0.5648f, 0.6688f, 0.7246f, 2.9300f, 5.3274f, 0.8269f}, // Y
{1.0114f, 0.5188f, 0.3642f, 0.3105f, 1.2918f, 0.5117f, 0.4509f, 0.3672f, 0.5232f, 2.2487f, 1.3686f, 0.4636f, 1.2605f, 1.1852f, 0.6290f, 0.5819f, 0.9277f, 0.8216f, 0.8269f, 2.8866f}, // V
};

// -------------------------------------------------------------------------
// (Gonnet marginal probabilities solved from O p = 1, where 
//  O is odds matrix, and p is probability vector to solve for)
const float Gonnet_PROBS[NUMAA] = {
// A           R           N           D           C           Q           E           G
// H           I           L           K           M           F           P           S
// T           W           Y           V
0.07691477f, 0.05208473f, 0.04450198f, 0.05403615f, 0.01915176f, 0.03930480f, 0.05939026f, 0.07545474f,
0.02418457f, 0.05054209f, 0.09926734f, 0.06011962f, 0.02195156f, 0.04124713f, 0.04543059f, 0.04717934f,
0.07169904f, 0.01292610f, 0.03165732f, 0.07295611f
};

// Gonnet frequency ratios;
// characteristics of the bit scaled matrix
//             Scale: 1/4 (ln2/4); 1/3 (ln2/3)
//           Entropy: 0.2465141
//     ExpectedScore: -0.2065916
//      HighestScore: 19; 14 (after rescaling)
//       LowestScore: -7; -5 (after rescaling)
const float GONNET_FREQRATIOS[NUMAA][NUMAA] = {
//  A         R         N         D         C         Q         E         G         H         I         L         K         M         F         P         S         T         W         Y         V
{1.737801f, 0.870964f, 0.933254f, 0.933254f, 1.122018f, 0.954993f, 1.000000f, 1.122018f, 0.831764f, 0.831764f, 0.758578f, 0.912011f, 0.851138f, 0.588844f, 1.071519f, 1.288250f, 1.148154f, 0.436516f, 0.602560f, 1.023293f}, // A 
{0.870964f, 2.951209f, 1.071519f, 0.933254f, 0.602560f, 1.412538f, 1.096478f, 0.794328f, 1.148154f, 0.575440f, 0.602560f, 1.862087f, 0.676083f, 0.478630f, 0.812831f, 0.954993f, 0.954993f, 0.691831f, 0.660693f, 0.630957f}, // R 
{0.933254f, 1.071519f, 2.398833f, 1.659587f, 0.660693f, 1.174898f, 1.230269f, 1.096478f, 1.318257f, 0.524807f, 0.501187f, 1.202264f, 0.602560f, 0.489779f, 0.812831f, 1.230269f, 1.122018f, 0.436516f, 0.724436f, 0.602560f}, // N 
{0.933254f, 0.933254f, 1.659587f, 2.951209f, 0.478630f, 1.230269f, 1.862087f, 1.023293f, 1.096478f, 0.416869f, 0.398107f, 1.122018f, 0.501187f, 0.354813f, 0.851138f, 1.122018f, 1.000000f, 0.301995f, 0.524807f, 0.512861f}, // D 
{1.122018f, 0.602560f, 0.660693f, 0.478630f,14.125375f, 0.575440f, 0.501187f, 0.630957f, 0.741310f, 0.776247f, 0.707946f, 0.524807f, 0.812831f, 0.831764f, 0.489779f, 1.023293f, 0.891251f, 0.794328f, 0.891251f, 1.000000f}, // C 
{0.954993f, 1.412538f, 1.174898f, 1.230269f, 0.575440f, 1.862087f, 1.479108f, 0.794328f, 1.318257f, 0.645654f, 0.691831f, 1.412538f, 0.794328f, 0.549541f, 0.954993f, 1.047129f, 1.000000f, 0.537032f, 0.676083f, 0.707946f}, // Q 
{1.000000f, 1.096478f, 1.230269f, 1.862087f, 0.501187f, 1.479108f, 2.290868f, 0.831764f, 1.096478f, 0.537032f, 0.524807f, 1.318257f, 0.630957f, 0.407380f, 0.891251f, 1.047129f, 0.977237f, 0.371535f, 0.537032f, 0.645654f}, // E 
{1.122018f, 0.794328f, 1.096478f, 1.023293f, 0.630957f, 0.794328f, 0.831764f, 4.570882f, 0.724436f, 0.354813f, 0.363078f, 0.776247f, 0.446684f, 0.301995f, 0.691831f, 1.096478f, 0.776247f, 0.398107f, 0.398107f, 0.467735f}, // G 
{0.831764f, 1.148154f, 1.318257f, 1.096478f, 0.741310f, 1.318257f, 1.096478f, 0.724436f, 3.981072f, 0.602560f, 0.645654f, 1.148154f, 0.741310f, 0.977237f, 0.776247f, 0.954993f, 0.933254f, 0.831764f, 1.659587f, 0.630957f}, // H 
{0.831764f, 0.575440f, 0.524807f, 0.416869f, 0.776247f, 0.645654f, 0.537032f, 0.354813f, 0.602560f, 2.511886f, 1.905461f, 0.616595f, 1.778279f, 1.258925f, 0.549541f, 0.660693f, 0.870964f, 0.660693f, 0.851138f, 2.041738f}, // I 
{0.758578f, 0.602560f, 0.501187f, 0.398107f, 0.707946f, 0.691831f, 0.524807f, 0.363078f, 0.645654f, 1.905461f, 2.511886f, 0.616595f, 1.905461f, 1.584893f, 0.588844f, 0.616595f, 0.741310f, 0.851138f, 1.000000f, 1.513561f}, // L 
{0.912011f, 1.862087f, 1.202264f, 1.122018f, 0.524807f, 1.412538f, 1.318257f, 0.776247f, 1.148154f, 0.616595f, 0.616595f, 2.089296f, 0.724436f, 0.467735f, 0.870964f, 1.023293f, 1.023293f, 0.446684f, 0.616595f, 0.676083f}, // K 
{0.851138f, 0.676083f, 0.602560f, 0.501187f, 0.812831f, 0.794328f, 0.630957f, 0.446684f, 0.741310f, 1.778279f, 1.905461f, 0.724436f, 2.691535f, 1.445440f, 0.575440f, 0.724436f, 0.870964f, 0.794328f, 0.954993f, 1.445440f}, // M 
{0.588844f, 0.478630f, 0.489779f, 0.354813f, 0.831764f, 0.549541f, 0.407380f, 0.301995f, 0.977237f, 1.258925f, 1.584893f, 0.467735f, 1.445440f, 5.011872f, 0.416869f, 0.524807f, 0.602560f, 2.290868f, 3.235937f, 1.023293f}, // F 
{1.071519f, 0.812831f, 0.812831f, 0.851138f, 0.489779f, 0.954993f, 0.891251f, 0.691831f, 0.776247f, 0.549541f, 0.588844f, 0.870964f, 0.575440f, 0.416869f, 5.754399f, 1.096478f, 1.023293f, 0.316228f, 0.489779f, 0.660693f}, // P 
{1.288250f, 0.954993f, 1.230269f, 1.122018f, 1.023293f, 1.047129f, 1.047129f, 1.096478f, 0.954993f, 0.660693f, 0.616595f, 1.023293f, 0.724436f, 0.524807f, 1.096478f, 1.659587f, 1.412538f, 0.467735f, 0.645654f, 0.794328f}, // S 
{1.148154f, 0.954993f, 1.122018f, 1.000000f, 0.891251f, 1.000000f, 0.977237f, 0.776247f, 0.933254f, 0.870964f, 0.741310f, 1.023293f, 0.870964f, 0.602560f, 1.023293f, 1.412538f, 1.778279f, 0.446684f, 0.645654f, 1.000000f}, // T 
{0.436516f, 0.691831f, 0.436516f, 0.301995f, 0.794328f, 0.537032f, 0.371535f, 0.398107f, 0.831764f, 0.660693f, 0.851138f, 0.446684f, 0.794328f, 2.290868f, 0.316228f, 0.467735f, 0.446684f,26.302680f, 2.570396f, 0.549541f}, // W 
{0.602560f, 0.660693f, 0.724436f, 0.524807f, 0.891251f, 0.676083f, 0.537032f, 0.398107f, 1.659587f, 0.851138f, 1.000000f, 0.616595f, 0.954993f, 3.235937f, 0.489779f, 0.645654f, 0.645654f, 2.570396f, 6.025596f, 0.776247f}, // Y 
{1.023293f, 0.630957f, 0.602560f, 0.512861f, 1.000000f, 0.707946f, 0.645654f, 0.467735f, 0.630957f, 2.041738f, 1.513561f, 0.676083f, 1.445440f, 1.023293f, 0.660693f, 0.794328f, 1.000000f, 0.549541f, 0.776247f, 2.187762f}, // V 
};


// -------------------------------------------------------------------------
//Substitution table: 'BLOSUM80'
//             Scale: 1/2
//           Entropy: 0.9868
//     ExpectedScore: -0.7442
//      HighestScore: 11
//       LowestScore: -6
const int BLOSUM80[NUMAA][NUMAA] =
{// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
{   5,  -2,  -2,  -2,  -1,  -1,  -1,   0,  -2,  -2,  -2,  -1,  -1,  -3,  -1,   1,   0,  -3,  -2,   0}, // A 
{  -2,   6,  -1,  -2,  -4,   1,  -1,  -3,   0,  -3,  -3,   2,  -2,  -4,  -2,  -1,  -1,  -4,  -3,  -3}, // R 
{  -2,  -1,   6,   1,  -3,   0,  -1,  -1,   0,  -4,  -4,   0,  -3,  -4,  -3,   0,   0,  -4,  -3,  -4}, // N
{  -2,  -2,   1,   6,  -4,  -1,   1,  -2,  -2,  -4,  -5,  -1,  -4,  -4,  -2,  -1,  -1,  -6,  -4,  -4}, // D
{  -1,  -4,  -3,  -4,   9,  -4,  -5,  -4,  -4,  -2,  -2,  -4,  -2,  -3,  -4,  -2,  -1,  -3,  -3,  -1}, // C
{  -1,   1,   0,  -1,  -4,   6,   2,  -2,   1,  -3,  -3,   1,   0,  -4,  -2,   0,  -1,  -3,  -2,  -3}, // Q
{  -1,  -1,  -1,   1,  -5,   2,   6,  -3,   0,  -4,  -4,   1,  -2,  -4,  -2,   0,  -1,  -4,  -3,  -3}, // E
{   0,  -3,  -1,  -2,  -4,  -2,  -3,   6,  -3,  -5,  -4,  -2,  -4,  -4,  -3,  -1,  -2,  -4,  -4,  -4}, // G
{  -2,   0,   0,  -2,  -4,   1,   0,  -3,   8,  -4,  -3,  -1,  -2,  -2,  -3,  -1,  -2,  -3,   2,  -4}, // H
{  -2,  -3,  -4,  -4,  -2,  -3,  -4,  -5,  -4,   5,   1,  -3,   1,  -1,  -4,  -3,  -1,  -3,  -2,   3}, // I
{  -2,  -3,  -4,  -5,  -2,  -3,  -4,  -4,  -3,   1,   4,  -3,   2,   0,  -3,  -3,  -2,  -2,  -2,   1}, // L
{  -1,   2,   0,  -1,  -4,   1,   1,  -2,  -1,  -3,  -3,   5,  -2,  -4,  -1,  -1,  -1,  -4,  -3,  -3}, // K
{  -1,  -2,  -3,  -4,  -2,   0,  -2,  -4,  -2,   1,   2,  -2,   6,   0,  -3,  -2,  -1,  -2,  -2,   1}, // M
{  -3,  -4,  -4,  -4,  -3,  -4,  -4,  -4,  -2,  -1,   0,  -4,   0,   6,  -4,  -3,  -2,   0,   3,  -1}, // F
{  -1,  -2,  -3,  -2,  -4,  -2,  -2,  -3,  -3,  -4,  -3,  -1,  -3,  -4,   8,  -1,  -2,  -5,  -4,  -3}, // P
{   1,  -1,   0,  -1,  -2,   0,   0,  -1,  -1,  -3,  -3,  -1,  -2,  -3,  -1,   5,   1,  -4,  -2,  -2}, // S
{   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -2,  -1,  -1,  -2,  -2,   1,   5,  -4,  -2,   0}, // T
{  -3,  -4,  -4,  -6,  -3,  -3,  -4,  -4,  -3,  -3,  -2,  -4,  -2,   0,  -5,  -4,  -4,  11,   2,  -3}, // W
{  -2,  -3,  -3,  -4,  -3,  -2,  -3,  -4,   2,  -2,  -2,  -3,  -2,   3,  -4,  -2,  -2,   2,   7,  -2}, // Y
{   0,  -3,  -4,  -4,  -1,  -3,  -3,  -4,  -4,   3,   1,  -3,   1,  -1,  -3,  -2,   0,  -3,  -2,   4}, // V
};

//Substitution table: 'BLOSUM62'
//             Scale: 1/2
//           Entropy: 0.6979
//     ExpectedScore: -0.5209
//      HighestScore: 11
//       LowestScore: -4
const int BLOSUM62[NUMAA][NUMAA] =
{// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
{   4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0}, // A 
{  -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3}, // R 
{  -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3}, // N 
{  -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3}, // D 
{   0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1}, // C 
{  -1,   1,   0,   0,  -3,   5,   2,  -2,   0 , -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2}, // Q 
{  -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2}, // E 
{   0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3}, // G 
{  -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3}, // H 
{  -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,   0,  -3,  -2,  -1,  -3,  -1,   3}, // I 
{  -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4 , -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1}, // L 
{  -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2}, // K
{  -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1}, // M 
{  -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1}, // F 
{  -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2}, // P 
{   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2}, // S 
{   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0}, // T
{  -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3}, // W 
{  -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1}, // Y 
{   0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4}, // V 
};

//Substitution table: 'BLOSUM45'
//             Scale: 1/3
//           Entropy: 0.3795
//     ExpectedScore: -0.2789
//      HighestScore: 15
//       LowestScore: -5
const int BLOSUM45[NUMAA][NUMAA] =
{// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
{   5,  -2,  -1,  -2,  -1,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -2,  -2,   0}, // A
{  -2,   7,   0,  -1,  -3,   1,   0,  -2,   0,  -3,  -2,   3,  -1,  -2,  -2,  -1,  -1,  -2,  -1,  -2}, // R
{  -1,   0,   6,   2,  -2,   0,   0,   0,   1,  -2,  -3,   0,  -2,  -2,  -2,   1,   0,  -4,  -2,  -3}, // N
{  -2,  -1,   2,   7,  -3,   0,   2,  -1,   0,  -4,  -3,   0,  -3,  -4,  -1,   0,  -1,  -4,  -2,  -3}, // D
{  -1,  -3,  -2,  -3,  12,  -3,  -3,  -3,  -3,  -3,  -2,  -3,  -2,  -2,  -4,  -1,  -1,  -5,  -3,  -1}, // C
{  -1,   1,   0,   0,  -3,   6,   2,  -2,   1,  -2,  -2,   1,   0,  -4,  -1,   0,  -1,  -2,  -1,  -3}, // Q
{  -1,   0,   0,   2,  -3,   2,   6,  -2,   0,  -3,  -2,   1,  -2,  -3,   0,   0,  -1,  -3,  -2,  -3}, // E
{   0,  -2,   0,  -1,  -3,  -2,  -2,   7,  -2,  -4,  -3,  -2,  -2,  -3,  -2,   0,  -2,  -2,  -3,  -3}, // G
{  -2,   0,   1,   0,  -3,   1,   0,  -2,  10,  -3,  -2,  -1,   0,  -2,  -2,  -1,  -2,  -3,   2,  -3}, // H
{  -1,  -3,  -2,  -4,  -3,  -2,  -3,  -4,  -3,   5,   2,  -3,   2,   0,  -2,  -2,  -1,  -2,   0,   3}, // I
{  -1,  -2,  -3,  -3,  -2,  -2,  -2,  -3,  -2,   2,   5,  -3,   2,   1,  -3,  -3,  -1,  -2,   0,   1}, // L
{  -1,   3,   0,   0,  -3,   1,   1,  -2,  -1,  -3,  -3,   5,  -1,  -3,  -1,  -1,  -1,  -2,  -1,  -2}, // K
{  -1,  -1,  -2,  -3,  -2,   0,  -2,  -2,   0,   2,   2,  -1,   6,   0,  -2,  -2,  -1,  -2,   0,   1}, // M
{  -2,  -2,  -2,  -4,  -2,  -4,  -3,  -3,  -2,   0,   1,  -3,   0,   8,  -3,  -2,  -1,   1,   3,   0}, // F
{  -1,  -2,  -2,  -1,  -4,  -1,   0,  -2,  -2,  -2,  -3,  -1,  -2,  -3,   9,  -1,  -1,  -3,  -3,  -3}, // P
{   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -3,  -1,  -2,  -2,  -1,   4,   2,  -4,  -2,  -1}, // S
{   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,   2,   5,  -3,  -1,   0}, // T
{  -2,  -2,  -4,  -4,  -5,  -2,  -3,  -2,  -3,  -2,  -2,  -2,  -2,   1,  -3,  -4,  -3,  15,   3,  -3}, // W
{  -2,  -1,  -2,  -2,  -3,  -1,  -2,  -3,   2,   0,   0,  -1,   0,   3,  -3,  -2,  -1,   3,   8,  -1}, // Y
{   0,  -2,  -3,  -3,  -1,  -3,  -3,  -3,  -3,   3,   1,  -2,   1,   0,  -3,  -1,   0,  -3,  -1,   5}, // V
};

//Substitution table produced by pscores
//             Scale: 1/3
//           Entropy: 0.3360
//     ExpectedScore: -0.2945
//      HighestScore: 12
//       LowestScore: -6
const int PSCORES_[NUMAA][NUMAA] =
{// A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
{   4,  -1,  -2,  -2,   1,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   0,   0,  -3,  -1,   0}, // A
{  -1,   5,   0,   0,  -2,   2,   2,  -2,   1,  -3,  -2,   3,  -1,  -3,  -1,   0,   0,  -2,  -1,  -3}, // R
{  -2,   0,   7,   2,   0,   1,   0,   0,   1,  -5,  -4,   1,  -2,  -3,  -1,   1,   0,  -3,  -1,  -4}, // N
{  -2,   0,   2,   7,  -3,   1,   2,  -1,   0,  -6,  -6,   0,  -4,  -4,  -1,   1,  -1,  -4,  -3,  -5}, // D
{   1,  -2,   0,  -3,   8,  -2,  -3,  -1,  -1,   0,   0,  -3,   0,   0,  -2,   1,   0,  -1,   0,   1}, // C
{  -1,   2,   1,   1,  -2,   4,   3,  -1,   1,  -3,  -2,   3,  -1,  -3,  -1,   0,   0,  -2,  -1,  -3}, // Q
{  -1,   2,   0,   2,  -3,   3,   5,  -1,   0,  -4,  -4,   2,  -2,  -4,   0,   0,   0,  -4,  -2,  -3}, // E
{   0,  -2,   0,  -1,  -1,  -1,  -1,   7,  -1,  -5,  -5,  -1,  -3,  -4,  -1,   0,  -1,  -4,  -3,  -4}, // G
{  -2,   1,   1,   0,  -1,   1,   0,  -1,   8,  -3,  -2,   0,   0,   0,  -1,   0,  -1,   0,   2,  -3}, // H
{  -1,  -3,  -5,  -6,   0,  -3,  -4,  -5,  -3,   5,   3,  -4,   2,   1,  -3,  -3,  -1,   0,  -1,   4}, // I
{  -1,  -2,  -4,  -6,   0,  -2,  -4,  -5,  -2,   3,   5,  -3,   2,   2,  -3,  -4,  -2,   0,   0,   1}, // L
{  -1,   3,   1,   0,  -3,   3,   2,  -1,   0,  -4,  -3,   5,  -2,  -4,   0,   0,   0,  -3,  -2,  -3}, // K
{  -1,  -1,  -2,  -4,   0,  -1,  -2,  -3,   0,   2,   2,  -2,   6,   2,  -2,  -2,  -1,   0,   1,   1}, // M
{  -2,  -3,  -3,  -4,   0,  -3,  -4,  -4,   0,   1,   2,  -4,   2,   6,  -2,  -3,  -2,   4,   4,   1}, // F
{  -1,  -1,  -1,  -1,  -2,  -1,   0,  -1,  -1,  -3,  -3,   0,  -2,  -2,   9,   0,  -1,  -3,  -2,  -2}, // P
{   0,   0,   1,   1,   1,   0,   0,   0,   0,  -3,  -4,   0,  -2,  -3,   0,   4,   2,  -3,  -2,  -2}, // S
{   0,   0,   0,  -1,   0,   0,   0,  -1,  -1,  -1,  -2,   0,  -1,  -2,  -1,   2,   5,  -2,  -1,   0}, // T
{  -3,  -2,  -3,  -4,  -1,  -2,  -4,  -4,   0,   0,   0,  -3,   0,   4,  -3,  -3,  -2,  12,   5,  -1}, // W
{  -1,  -1,  -1,  -3,   0,  -1,  -2,  -3,   2,  -1,   0,  -2,   1,   4,  -2,  -2,  -1,   5,   7,  -1}, // Y
{   0,  -3,  -4,  -5,   1,  -3,  -3,  -4,  -3,   4,   1,  -3,   1,   1,  -2,  -2,   0,  -1,  -1,   5}, // V
};

// -------------------------------------------------------------------------
//
const float Blosum80ScalingConstant = 2.0f;
const float Blosum62ScalingConstant = 2.0f;
const float Blosum45ScalingConstant = 3.0f;
const float PScores_ScalingConstant = 3.0f;
// const float Gonnet_ScalingConstant = 4.0f;
const float Gonnet_ScalingConstant = 3.0f;

// -------------------------------------------------------------------------
//

#define FINF ((float)(32767))

// (Data from the NCBI BLAST toolkit)
const float BLOSUM80_VALUES[Num80Entries][NumFields] =
{//Open Extend  Lambda  K       H       alpha    beta
{ FINF, FINF,   0.3430f, 0.177f,  0.6568f, 0.5222f, -1.6f   },
{ 25.f,   2.f,      0.342f,  0.17f,   0.66f,   0.52f,   -1.6f   },
{ 13.f,   2.f,      0.336f,  0.15f,   0.57f,   0.59f,   -3.0f   },
{ 9.f,    2.f,      0.319f,  0.11f,   0.42f,   0.76f,   -6.0f   },
{ 8.f,    2.f,      0.308f,  0.090f,  0.35f,   0.89f,   -9.0f   },
{ 7.f,    2.f,      0.293f,  0.070f,  0.27f,   1.10f,   -14.0f  },
{ 6.f,    2.f,      0.268f,  0.045f,  0.19f,   1.40f,   -19.0f  },
{ 11.f,   1.f,      0.314f,  0.095f,  0.35f,   0.90f,   -9.0f   },
{ 10.f,   1.f,      0.299f,  0.071f,  0.27f,   1.10f,   -14.0f  },
{ 9.f,    1.f,      0.279f,  0.048f,  0.20f,   1.40f,   -19.0f  }
};

// (Data from the NCBI BLAST toolkit)
const float BLOSUM62_VALUES[NumEntries][NumFields] =
{//Open Extend  Lambda  K       H       alpha    beta
{ FINF, FINF,   0.3176f, 0.134f,  0.4012f, 0.7916f, -3.2f   },
{ 11.f,   2.f,      0.297f,  0.082f,  0.27f,   1.1f,    -10.0f  },
{ 10.f,   2.f,      0.291f,  0.075f,  0.23f,   1.3f,    -15.0f  },
{ 9.f,    2.f,      0.279f,  0.058f,  0.19f,   1.5f,    -19.0f  },
{ 8.f,    2.f,      0.264f,  0.045f,  0.15f,   1.8f,    -26.0f  },
{ 7.f,    2.f,      0.239f,  0.027f,  0.10f,   2.5f,    -46.0f  },
{ 6.f,    2.f,      0.201f,  0.012f,  0.061f,  3.3f,    -58.0f  },
{ 13.f,   1.f,      0.292f,  0.071f,  0.23f,   1.2f,    -11.0f  },
{ 12.f,   1.f,      0.283f,  0.059f,  0.19f,   1.5f,    -19.0f  },
{ 11.f,   1.f,      0.267f,  0.041f,  0.14f,   1.9f,    -30.0f  },
{ 10.f,   1.f,      0.243f,  0.024f,  0.10f,   2.5f,    -44.0f  },
{ 9.f,    1.f,      0.206f,  0.010f,  0.052f,  4.0f,    -87.0f  }
};

// (Data from the NCBI BLAST toolkit)
const float BLOSUM45_VALUES[Num45Entries][NumFields] =
{//Open Extend  Lambda  K       H       alpha    beta
{ FINF, FINF,   0.2291f, 0.0924f, 0.2514f, 0.9113f, -5.7f    },
{ 13.f,   3.f,      0.207f,  0.049f,  0.14f,   1.5f,    -22.0f   },
{ 12.f,   3.f,      0.199f,  0.039f,  0.11f,   1.8f,    -34.0f   },
{ 11.f,   3.f,      0.190f,  0.031f,  0.095f,  2.0f,    -38.0f   },
{ 10.f,   3.f,      0.179f,  0.023f,  0.075f,  2.4f,    -51.0f   },
{ 16.f,   2.f,      0.210f,  0.051f,  0.14f,   1.5f,    -24.0f   },
{ 15.f,   2.f,      0.203f,  0.041f,  0.12f,   1.7f,    -31.0f   },
{ 14.f,   2.f,      0.195f,  0.032f,  0.10f,   1.9f,    -36.0f   },
{ 13.f,   2.f,      0.185f,  0.024f,  0.084f,  2.2f,    -45.0f   },
{ 12.f,   2.f,      0.171f,  0.016f,  0.061f,  2.8f,    -65.0f   },
{ 19.f,   1.f,      0.205f,  0.040f,  0.11f,   1.9f,    -43.0f   },
{ 18.f,   1.f,      0.198f,  0.032f,  0.10f,   2.0f,    -43.0f   },
{ 17.f,   1.f,      0.189f,  0.024f,  0.079f,  2.4f,    -57.0f   },
{ 16.f,   1.f,      0.176f,  0.016f,  0.063f,  2.8f,    -67.0f   }
};

// -------------------------------------------------------------------------
//

const float PSCORES_VALUES[NumPSEntries][NumFields] =
{//Open Extend  Lambda  K       H       alpha    beta
{ FINF, FINF,   0.2218f, 0.0927f, 0.2199f, 0.9500f, -6.0f    },
};

const float GONNET_VALUES[NumGnEntries][NumFields] =
{//Open Extend  Lambda  K       H       alpha    beta
// { FINF, FINF,   0.1730, 0.0744, 0.2002, 0.5000, -1.0    },//scale==4
{ FINF, FINF,   0.2305f, 0.0781f, 0.2143f, 0.5000f, -1.0f    },//scale==3
};
