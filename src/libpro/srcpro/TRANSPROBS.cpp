/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <string.h>

#include "TRANSPROBS.h"

// =========================================================================
// Globals -----------------------------------------------------------------
//
_TRANSPROBS TRANSPROBS;

const char* gTPTRANS_NAMES[P_NSTATES] = {
    "MM", "MI", "MD",
    "IM", "II", "ID",
    "DM", "DI", "DD"
};

const char* gTPSTATES_NAMES[PS_NSTATES] = {
    "M", "I", "D"
};

//number of transitions per state
const int gTPTRANS_NTPS = 3;

//Dirichlet distribution parameters, alphas, for transition priors
// * G.Mitchison, S.R.Eddy data
// float gTPTRANS_ALPHAS[P_NSTATES] = {
//     0.7939, 0.0278, 0.0135,   //MM, MI, MD,
//     0.1551, 0.1331, 0.0000,   //IM, II, ID,
//     0.9002, 0.0000, 0.5630    //DM, DI, DD,
// };
float gTPTRANS_ALPHAS[P_NSTATES] = {
    0.7939f, 0.0234f, 0.0234f,   //MM, MI, MD,
    0.1551f, 0.4653f, 0.0000f,   //IM, II, ID,
    0.3001f, 0.0000f, 0.9002f    //DM, DI, DD,
};

// =========================================================================
// _TRANSPROBS: constructor
//
_TRANSPROBS::_TRANSPROBS()
{
    effnos_[PS_M] = 1.0f;
    effnos_[PS_I] = 0.0f;
    effnos_[PS_D] = 0.0f;
    InitOnce();
    Initialize();
    SetClass( Dirichlet );
}

// InitOnce: initialize data once
//
void _TRANSPROBS::InitOnce()
{
    int st, states, n;
    float sum;

    memset( priors_, 0, sizeof(float) * P_NSTATES );

    for( n = 0, states = 0; states < P_NSTATES; states += gTPTRANS_NTPS, n++ ) {
        sum = 0.0;
        for( st = states; st < states + gTPTRANS_NTPS && st < P_NSTATES; st++ )
            sum += gTPTRANS_ALPHAS[st];
        if( n < PS_NSTATES )
            sumalphas_[n] = sum;
    }
}

// Initialize: initialize or reset private data
//
void _TRANSPROBS::Initialize()
{
    memset( pmestimators_, 0, sizeof(float) * P_NSTATES );
    memset( logpmestimators_, 0, sizeof(float) * P_NSTATES );
}

// SetClass: set a class for priors
//
void _TRANSPROBS::SetClass( TPriorClass value )
{
    float dummy[P_NSTATES];
    class_ = value;

    memset( dummy, 0, sizeof(float) * P_NSTATES );
    PME( &dummy );
    if( GetPMEstimators())
        memcpy( priors_, *GetPMEstimators(), sizeof(float) * P_NSTATES );
}

// GetEffNoSequences: get effective no. observations for state `mid'
//
float _TRANSPROBS::GetEffNoSequences( int mid ) const
{
    if( mid < 0 || PS_NSTATES <= mid )
        throw MYRUNTIME_ERROR("_TRANSPROBS::GetEffNoSequences: Invalid state.");
    return effnos_[mid];
}

// SetEffNoSequences: set effective no. observations for all states
//
void _TRANSPROBS::SetEffNoSequences( float nm, float ni, float nd )
{
    effnos_[PS_M] = nm;
    effnos_[PS_I] = ni;
    effnos_[PS_D] = nd;
}

// PME: calculate posterior mean estimates of probabilities
//
void _TRANSPROBS::PME( const float (*obsfreqs)[P_NSTATES])
{
    Initialize();
    switch( GetClass()) {
        case Dirichlet: DirPME( obsfreqs ); break;
        default:
            throw MYRUNTIME_ERROR("_TRANSPROBS::PME: Invalid class of priors.");
    };
    return;
}

// DirPME: apply Dirichlet priors to posterior mean estimate
// NOTE: have to have counts rather than frequencies to compute properly;
//  if not, effective number of sequences from which the observations were 
//  derived should be provided
//
void _TRANSPROBS::DirPME( const float (*obsfreqs)[P_NSTATES])
{
    int n, st, states;
    const int maxtr = gTPTRANS_NTPS;
    float est, logest;
    float sum, tsum;
    float mid[PS_NSTATES];
    float /*midobs, */nobs;
    const float accuracy = 1.0e-6f;

    if( obsfreqs == NULL )
        throw MYRUNTIME_ERROR("_TRANSPROBS::DirPME: Memory access error.");

    mid[PS_M] = GetEffNoSequences( PS_M );
    mid[PS_I] = GetEffNoSequences( PS_I );
    mid[PS_D] = GetEffNoSequences( PS_D );
//     midobs = mid[PS_M] + mid[PS_I] + mid[PS_D];

    for( n = 0, states = 0; states < P_NSTATES; states += maxtr, n++ ) {
        sum = 0.0f;
        for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
//             if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
            sum += (*obsfreqs)[st];
        }
        if( PS_NSTATES <= n ) {
            throw MYRUNTIME_ERROR( "_TRANSPROBS::DirPME: Invalid number of states." );
            break;
        }
        nobs = mid[n];
        tsum = nobs * sum + sumalphas_[n];
        if( 0.0f <= sum && 0.0f < tsum ) {
            //to remove g++ warning about possible overflow; this is not the case here:
            //for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
            for( st = states; st - maxtr < states && st < P_NSTATES; st++ ) {
                est = 0.0f;
                logest = (float)LOG_PROB_MIN;
//                 if( st != P_ID || st != P_DI ) {//IGNORE states P_ID and P_DI
                    est = ( nobs * (*obsfreqs)[st] + gTPTRANS_ALPHAS[st]) / tsum;
                    if( 0.0f < est ) 
                        logest = logf( est );
//                 }
                pmestimators_[st] = est;
                logpmestimators_[st] = logest;
            }
        }
    }

    //verify probability conservation
    for( n = 0, states = 0; states < P_NSTATES; states += maxtr, n++ ) {
        sum = 0.0f;
        for( st = states; st < states + maxtr && st < P_NSTATES; st++ )
//             if( st != P_ID && st != P_DI )//omit to-be-ignored states
                sum += pmestimators_[st];
        if( sum < 1.0f - accuracy || sum > 1.0f + accuracy )
            throw MYRUNTIME_ERROR("_TRANSPROBS::DirPME: "
                    "Transition probability estimates are invalid.");
    }
}
