/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/pslerror.h"
#include "extsp/pslvector.h"
#include "liblib/alpha.h"
#include "libpro/srcpro/SUBSTABLE.h"
#include "libpro/srcpro/Configuration.h"
// #include "libpro/tfopt/TargetFreqOptimizerH.h"
#include "AbstractScoreMatrix.h"
#include "IntegerScoreMatrix.h"
#include "ScoresAttr.h"
#include "ProfileMatrix.h"

template<typename TScore>
float ProfileMatrix<TScore>::position_[NUMALPH];

//score for masked positions
const float g_scoreX = 0.0f;

// -------------------------------------------------------------------------
// Constructor
//
template<typename TScore>
ProfileMatrix<TScore>::ProfileMatrix(
    const float (*pssmscores)[NUMALPH],
    const char*     ress,
    int             length,
    Configuration   config[NoCTypes],
    int             precscale )
:
    IntegerScoreMatrix<TScore>( 
        AbstractScoreMatrix<TScore>::PositionSpecific, 
        config, 
        AbstractScoreMatrix<TScore>::ComputeStatistics,
        precscale ),
    residues_( ress ),
    scores_( length * NUMALPH ),
    rprobs_( length ),
    cprobs_( NUMALPH )
{
    IntegerScoreMatrix<TScore>::Init( length, NUMALPH );
    FillMatrix( pssmscores );
    AbstractScoreMatrix<TScore>::SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// Default construction is invalid
//
template<typename TScore>
ProfileMatrix<TScore>::ProfileMatrix()
:
    IntegerScoreMatrix<TScore>(),
    residues_( NULL ),
    scores_( 1 ),
    rprobs_( 1 ),
    cprobs_( 1 )
{
}

// -------------------------------------------------------------------------
// Destructor
//
template<typename TScore>
ProfileMatrix<TScore>::~ProfileMatrix()
{
}

// -------------------------------------------------------------------------
// FillMatrix: calculate profile scores
//
template<typename TScore>
void ProfileMatrix<TScore>::FillMatrix( const float (*pssmscores)[NUMALPH] )
{
    if( GetQuerySize() < 1 || GetSubjectSize() != NUMALPH )
        throw MYRUNTIME_ERROR( "ProfileMatrix::FillMatrix: Invalid matrix dimensions." );

    const int fct = IntegerScoreMatrix<TScore>::GetScoresFactor();
    int     negatives = 0;
    bool    allnegats = true;
    float   value;
    TScore  score;

    //fill matrix with values
    for( int m = 0; m < GetSubjectSize(); m++ ) {
        negatives = 0;

        for( int n = 0; n < GetQuerySize(); n++ )
        {
            value = pssmscores[n][m];
            score = (TScore)rintf( value * (float)fct );

            IntegerScoreMatrix<TScore>::SetScore( m, n, score );

            //update counter on each negative score found
            if( allnegats && score < 0 ) negatives++;
        }

        if( allnegats && negatives != GetQuerySize())
            allnegats = false;
    }

    IntegerScoreMatrix<TScore>::SetAllNegatives( allnegats );
}

// -------------------------------------------------------------------------
// GetVectorAt: get the vecor of values at the given position of the matrix
//
template<typename TScore>
const float ( *ProfileMatrix<TScore>::GetVectorAt( int n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( n < 0 || GetQuerySize() <= n || GetSubjectSize() != NUMALPH )
        throw MYRUNTIME_ERROR( "ProfileMatrix::GetVectorAt: Memory access error." );
#endif

    for( int r = 0; r < NUMALPH; r++ )
        position_[r] = 
            IntegerScoreMatrix<TScore>::GetFinalScore( IntegerScoreMatrix<TScore>::GetScore(r,n));

    return &position_;//( CONST_SCORE_ARRAY )( GetScores() + n );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: callback method: calculate probabilities of 
// scores observed at positions (i,j). 
// Probability is computed as
//
//      1     _
//    ----   \   Pa
//   length  /_ 
//           i,j:
//         sij=sk
//
//  where Pa are background probabilities; sij are score at (i,j), sk is a
//  discrete score value; length is for query;
//
//  This method allocates space for probabilities!
//
template<typename TScore>
void ProfileMatrix<TScore>::ComputeScoreProbabilities( ScoresAttr<TScore>* PATTR )
{
    const mystring preamb = "ProfileMatrix::ComputeScoreProbabilities: ";

    if( !GetResidues() || GetSubjectSize() != NUMALPH )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile matrix." );

    if( PATTR == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null attributes." );

    const float accuracy = 1.e-4f;
    char    strbuf[BUF_MAX];
    TScore  score = 0;
    float   avgscore = 0.0f;
    float   consv = 0.0f;
    float   prob = 0.0f;
    size_t  length = 0;
    char    res;
    int     n, m, nn, mm;

    static int  symB = HashAlphSymbol('B');//N,D
    static int  symZ = HashAlphSymbol('Z');//Q,E
    static int  symJ = HashAlphSymbol('J');//I,L
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');
    static int  resI = HashAlphSymbol('I');
    static int  resL = HashAlphSymbol('L');

    scores_.Clear();
    rprobs_.Clear();
    cprobs_.Clear();


    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
        res = GetResidueAt( n );
//         if( res == X )
//             continue;
        prob = STABLE.PROBABility( res );

        if( res == X )
            prob = 0.05f;
        else if( res == symB )
            prob = ( STABLE.PROBABility( resN ) + STABLE.PROBABility( resD )) * 0.5f;
        else if( res == symZ )
            prob = ( STABLE.PROBABility( resQ ) + STABLE.PROBABility( resE )) * 0.5f;
        else if( res == symJ )
            prob = ( STABLE.PROBABility( resI ) + STABLE.PROBABility( resL )) * 0.5f;

        if( prob <= 0.0f )
            continue;
        rprobs_.AddValueAt( nn++, prob );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
        prob = STABLE.PROBABility( m );
        if( prob <= 0.0f )
            continue;
        cprobs_.AddValueAt( mm++, prob );
    }


    for( n = 0; n < GetQuerySize(); n++ ) {
//         if( GetResidueAt( n ) == X )
//             continue;

//         if( STABLE.PROBABility( n ) <= 0.0f )
//             continue;

        length++;
        for( m = 0; m < GetSubjectSize(); m++ )
        {
            prob = STABLE.PROBABility( m );
            score = IntegerScoreMatrix<TScore>::GetScore( m, n );

            if( prob <= 0.0f )
                continue;

            if( score <= SCORE_MIN )
                    scores_.Push( g_scoreX );
            else    scores_.Push((float)score );

            if( score <= SCORE_MIN )
                continue;

            PATTR->IncProbabilityOf( score, prob );
        }
    }

    //normalize probabilities
    //(min/max scores have been found in ScoresAttr<TScore>::ComputeScoreProbabilities,
    // where this member function is called back from)
    for( TScore sc = PATTR->GetMinScore(); sc <= PATTR->GetMaxScore(); sc++ ) {
        PATTR->DivideProbabilityOf( sc, (float)length );
        avgscore += sc * PATTR->GetProbabilityOf( sc );
    }


    consv = rprobs_.Sum();
    if( 0.0 < consv ) {
        rprobs_.MultiplyBy( 1.0f / consv );
        consv = rprobs_.Sum();
        if( consv < 1.0f - accuracy || consv > 1.0f + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw MYRUNTIME_ERROR( preamb + strbuf );
        }
    }
    consv = cprobs_.Sum();
    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy ) {
        sprintf( strbuf, "Background probabilities not conserved: %f", consv );
        throw MYRUNTIME_ERROR( preamb + strbuf );
    }


    PATTR->SetExpectedScore( avgscore );
}


// =========================================================================
// OPTIMIZATION OF TARGET FREQUENCIES
//
// OptimizeTargetFrequencies: Optimize target frequencies
//
template<typename TScore>
void ProfileMatrix<TScore>::OptimizeTargetFrequencies()
{
//     TODO:
//     float lambda;
// 
//     if( IntegerScoreMatrix<TScore>::GetAllNegatives()) {
//         warning( "Unable to optimize target frequencies." );
//         return;
//     }
// 
//     IntegerScoreMatrix<TScore>::ComputeStatisticalParameters();
// 
//     if( scores_.GetSize() < 1 || rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 )
//         return;
// 
//     if( 0.0f <= IntegerScoreMatrix<TScore>::GetExpectedScore() || 
//         IntegerScoreMatrix<TScore>::GetLambda() <= 0.0f ) {
//         warning( "Unable to optimize target frequencies." );
//         return;
//     }
// 
//     const int fct = IntegerScoreMatrix<TScore>::GetScoresFactor();
//     lambda = IntegerScoreMatrix<TScore>::GetLambda();
// 
//     scores_.MultiplyBy( lambda / (float)fct );
// 
// 
//     TargetFreqOptimizerH    optimizer( scores_, rprobs_, cprobs_ );
//     int                     status;
//     int                     n, m, ind = 0;
// 
// //{{TESTING..
// // PrintScoreMatrix( stderr );
// // fprintf( stderr, "\n\n" );
// // for( n = 0; n < rprobs_.GetSize(); n++ ) {
// //     for( m = 0; m < cprobs_.GetSize(); m++ )
// //         fprintf( stderr, "%12f ", scores_.GetValueAt( ind++ ));
// //     fprintf( stderr, "\n");
// // }
// //}}
// 
//     //no need to set lambda if scores were previously multiplied by
// //     optimizer.SetLambda( lambda );
//     optimizer.SetConstrainedH( IntegerScoreMatrix<TScore>::GetRefH());
//     optimizer.SetNMResidTol( TARGET_FREQ_OPT_TOLERANCE );
//     optimizer.SetMaxNoNMIterations( TARGET_FREQ_OPT_MAXITERATIONS );
// 
//     status = optimizer.Optimize();
// 
//     if( status != PSL_SUCCESS /*&& status != PSL_MAXITERATS */) {
// #ifdef TFOPTESTPRINT
//         warning( extspsl::TranslatePSLError( status ));
// #else
//         warning( "Target frequencies not optimized." );
// #endif
//         return;
//     }
// 
//     if( status != PSL_SUCCESS )
//         optimizer.Normalize();
// 
//     if( optimizer.Negatives()) {
//         warning( "Invalid target frequencies." );
//         return;
//     }
// 
//     IncorporateTargetFrequencies( optimizer.GetOptTargetFreqs());
}

// -------------------------------------------------------------------------
// IncorporateTargetFrequencies: calculate new scores based on optimized
// target frequencies
//
template<typename TScore>
void ProfileMatrix<TScore>::IncorporateTargetFrequencies( const extspsl::Pslvector& optfreqs )
{
    if( rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 ||
        optfreqs.GetSize() != rprobs_.GetSize() * cprobs_.GetSize())
        return;

    const float accuracy = 1.0e-5f;
    const int fct = IntegerScoreMatrix<TScore>::GetScoresFactor();
    int     n, m, nn, mm, ind;
    float   tfval = 0.0f;//target frequency value
    float   proquery = 0.0f;//background probability at query position
    float   prosbjct = 0.0f;//background probability at subject position
    float   loc_score = 0.0f;
    float   consv;
    float   value;
    TScore  score;
    //char    res;

    consv = optfreqs.Sum();
    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy ) {
        warning( "Target frequencies not conserved." );
        return;
    }

    ind = 0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
//         res = GetResidueAt( n );
//         if( res == X )
//             continue;
//         if( ASTERISK <= res )
//             continue;

        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ )
        {
            if( STABLE.PROBABility( m ) <= 0.0f )
                continue;

            prosbjct = cprobs_.GetValueAt( mm++ );

            tfval = optfreqs.GetValueAt( ind++ );

            if( proquery <= 0.0f || prosbjct <= 0.0f )
                    throw MYRUNTIME_ERROR(
                        "ProfileMatrix::IncorporateTargetFrequencies: Invalid probabilities.");

            value = logf( tfval / ( proquery * prosbjct ));

            loc_score = (float)IntegerScoreMatrix<TScore>::GetScore( m, n );

            if( loc_score <= SCORE_MIN )
                continue;

            score = (TScore)rintf( value * (float)fct );
            IntegerScoreMatrix<TScore>::SetScore( m, n, score );
        }
    }
}



// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: print table of computed statistical parameter
//     values to a stream
//
template<typename TScore>
void ProfileMatrix<TScore>::PrintParameterTable( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    float expscore = IntegerScoreMatrix<TScore>::GetExpectedScore();
    float parK = IntegerScoreMatrix<TScore>::GetK();
    float lmbd = IntegerScoreMatrix<TScore>::GetLambda();
    float refK = IntegerScoreMatrix<TScore>::GetRefK();
    float reflmbd = IntegerScoreMatrix<TScore>::GetRefLambda();
    float parH = IntegerScoreMatrix<TScore>::GetEntropy();
    float parE = IntegerScoreMatrix<TScore>::GetExpectedScore();

    if( 0.0f <= expscore ) {
        print_func( vpn, "Expected score per position is non-negative, %.4f!%s%s", 
                    expscore, NL, NL );
        return;
    }

    char kinf[BUF_MAX], lmbdinf[BUF_MAX];

    if( 0.0f < parK )   sprintf( kinf, "%6.4f", parK );
        else            sprintf( kinf, "n/a" );

    if( 0.0 < lmbd )    sprintf( lmbdinf, "%6.4f", lmbd );
        else            sprintf( lmbdinf, "n/a" );

    print_func( vpn, "%-25s  %-6s   %-6s%s", " ", "K", "Lambda", NL );
    print_func( vpn, "%-25s  %6s   %6s%s", "Computed  ungapped,", kinf, lmbdinf, NL );
    print_func( vpn, "%-25s  %6.4f   %6.4f%s", "Reference ungapped,", refK, reflmbd, NL );

    print_func( vpn, "Entropy, %6.4f; Expected, %6.4f%s%s", parH, parE, NL, NL );
}

// // PrintFinal: nothing to print
// //
// 
// void ProfileMatrix::PrintFinal( TPrintFunction, void* ) const
// {
// }

// -------------------------------------------------------------------------
// PrintScoreMatrix: print profile matrix
//
template<typename TScore>
void ProfileMatrix<TScore>::PrintScoreMatrix( FILE* fp )
{
    int l = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific score matrix%s", 32, NL );

    fprintf( fp, "%9c", 32 );

    for( int m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, "%4c", DehashCode( m ));

    for( int n = 0; n < GetQuerySize(); n++ ) {
        fprintf( fp, "%s%5d %c   ", NL, ++l, DehashCode( GetResidueAt(n)));

        for( int m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%3d ", (int)IntegerScoreMatrix<TScore>::GetScore(m,n));
    }
    fprintf( fp, "%s", NL );
}

// ==== Functions ==========================================================
// ComputedSubMatrixWithParams: output computed substitution matrix with
// statistical parameters calculated for it
//
void ComputedSubMatrixWithParams()
{
    //filename is not set, so that an error occurs when attempting to read or write 
    Configuration   config[NoCTypes];
    SetUngappedParams( config[CTUngapped] );//ungapped configuration
    SetUngappedParams( config[CTGapped] );//not using gapped configuration, make it identical to ungapped configuration

    const char resids[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

    float locgonnet[NUMAA][NUMALPH];
    float (*p)[NUMALPH] = locgonnet;
    //make a gonnet matrix suitable for this class
    for( int a = 0; a < NUMAA; a++, p++ ) {
        memset( p, 0, NUMALPH * sizeof(float));
        memcpy( p, COMPUTED_GONNET.data + a, NUMAA * sizeof(float));
    }

    ProfileMatrix<int> submat( locgonnet, resids, NUMAA, config, 1/*precscale*/ );
    submat.ComputeStatisticalParameters();
//     submat.ScaleScoreMatrix(); //if scaling needed
    submat.PrintScoreMatrix( stderr );
    dynamic_cast<AbstractScoreMatrix<int>&>(submat).PrintParameterTable( stderr );
}

// =========================================================================

template class ProfileMatrix<int>;
