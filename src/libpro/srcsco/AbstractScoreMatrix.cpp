/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/psl.h"
#include "libpro/srcpro/Configuration.h"
#include "AbstractScoreMatrix.h"

// template<typename TScore>
// size_t AbstractScoreMatrix<TScore>::deltaLength_ = 0;

// template<typename TScore>
// uint64_mt AbstractScoreMatrix<TScore>::searchSpace_ = 0ULL;

// -------------------------------------------------------------------------
// constructor
//
template<typename TScore>
AbstractScoreMatrix<TScore>::AbstractScoreMatrix(
        TType           type,
        Configuration   config[],
        TStats          stats,
        int             precscale )
:
    type_( type ),
    stats_( stats ),
    name_( NULL ),

    usemodimage_( false ),

    scaled_by_factor_( precscale ),
    score_multiplier_( 1.0f ),

    scadjrelent_( false ),
    scadjcorrsc_( false ),

    min_score_( 0 ),
    max_score_( 0 ),

    expscore_( 0.0f ),

    referenceLambda_( -1.0f ),
    referenceH_( 0.0f ),
    referenceK_( -1.0f ),
    experimentalGappedLambda_( -1.0f ),
    experimentalGappedK_( -1.0f ),
    derivedGappedLambda_( -1.0f ),
    derivedGappedK_( -1.0f ),

    lambda_( -1.0f ),
    entropy_( 0.0f ),
    parameterK_( -1.0f ),

    queryLen_( 0 ),
    subjectLen_( 0 ),

    supportoptfreq_( false ),

    configuration_( config )
{
    InitializeStatParams();
}

// -------------------------------------------------------------------------
// Init: data initialization
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::Init( int querylen, int sbjctlen )
{
    SetQuerySize( querylen );
    SetSubjectSize( sbjctlen );
}

// -------------------------------------------------------------------------
// default constructor
//
template<typename TScore>
AbstractScoreMatrix<TScore>::AbstractScoreMatrix()
:
    type_( NoType ),
    stats_( StatisticsGiven ),
    name_( NULL ),

    usemodimage_( false ),

    scaled_by_factor_( 1 ),
    score_multiplier_( 1.0f ),

    scadjrelent_( false ),
    scadjcorrsc_( false ),

    min_score_( 0 ),
    max_score_( 0 ),

    expscore_( 0.0f ),

    referenceLambda_( -1.0f ),
    referenceH_( 0.0f ),
    referenceK_( -1.0f ),
    experimentalGappedLambda_( -1.0f ),
    experimentalGappedK_( -1.0f ),
    derivedGappedLambda_( -1.0f ),
    derivedGappedK_( -1.0f ),

    lambda_( -1.0f ),
    entropy_( 0.0f ),
    parameterK_( -1.0f ),

    queryLen_( 0 ),
    subjectLen_( 0 ),

    supportoptfreq_( false ),

    configuration_( NULL )
{
    throw MYRUNTIME_ERROR("AbstractScoreMatrix::AbstractScoreMatrix: "
                "Default initialization is prohibited.");
}

// -------------------------------------------------------------------------
// destructor
//
template<typename TScore>
AbstractScoreMatrix<TScore>::~AbstractScoreMatrix()
{
}

// -------------------------------------------------------------------------
// InitializeStatParams: initialize statistical parameters
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::InitializeStatParams()
{
    const Configuration& ungapped_config = GetConfiguration( CTUngapped );
    const Configuration& gapped_config = GetConfiguration( CTGapped );

    SetRefLambda(   ungapped_config.GetLambda());
    SetRefH(        ungapped_config.GetH());
    SetRefK(        ungapped_config.GetK());
    SetExpGappedLambda( gapped_config.GetLambda());
    SetExpGappedK(      gapped_config.GetK());

    if( GetStats() == StatisticsGiven ) {
        //all necessary parameters are assumed to have been
        //precomputed earlier and saved to configuration
        SetLambda(  ungapped_config.GetLambda());
        SetEntropy( ungapped_config.GetH());
        SetK(       ungapped_config.GetK());
        //set expected score to validate statistics
        SetExpectedScore( -1.0f );

        SetDerivedGappedLambda( gapped_config.GetLambda());
        SetDerivedGappedK(      gapped_config.GetK());

        SetMultiplier( ungapped_config.GetScaleFactor());
    }
}

// -------------------------------------------------------------------------
// DeriveGappedParameters: estimate gapped statistical parameters given
//     experimentally defined values and computed ungapped ones
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::DeriveGappedParameters()
{
    //(computed ungapped for a score system) * 
    //(gapped defined experimentally) / (ungapped computed analitically)
    if( 0.0f < GetLambda())
        //SetDerivedGappedLambda( GetExpGappedLambda());
        SetDerivedGappedLambda( GetLambda() * GetExpGappedLambda() / GetRefLambda());
    if( 0.0f < GetK())
        SetDerivedGappedK( GetK() * GetExpGappedK() / GetRefK());
}

// -------------------------------------------------------------------------
// ComputeExpectation: calculate e-value given statistical parameters
//
template<typename TScore>
float AbstractScoreMatrix<TScore>::ComputeExpectation(
        float score,
        float* ref_expect,
        float* pure_expect,
        float* pair_expect,
        float* bitscore ) const
{
    float expect = -1.0f;
    float pairsspace = 0.0f;

    if( GetLambda() <= 0.0f || GetK() < 0.0f ) {
        if( ref_expect )
            //calculate expectation for ungapped alignment:
            // when expected score per position is positive, 
            // predict approx. e-value of alignment score
            *ref_expect = GetExpGappedK() * GetSearchSpace() * expf( -GetExpGappedLambda() * score );

        if( pair_expect ) {
            if( (int)GetDeltaLength() < GetQuerySize() && (int)GetDeltaLength() < GetSubjectSize())
                pairsspace = (float)(( GetQuerySize() - GetDeltaLength()) * ( GetSubjectSize() - GetDeltaLength()));
            else
                pairsspace = (float)(GetQuerySize() * GetSubjectSize());
            *pair_expect = GetRefK() * pairsspace * expf( -GetRefLambda() * score );
        }

        if( pure_expect )
            *pure_expect = GetRefK() * GetSearchSpace() * expf( -GetRefLambda() * score );

        if( bitscore && 0.0f < GetExpGappedK())
            *bitscore = ( GetExpGappedLambda() * score - logf( GetExpGappedK())) / SLC_LN2;

        return expect;
    }

    if( GetStats() == StatisticsGiven )
        expect = GetExpGappedK() * GetSearchSpace() * expf( -GetExpGappedLambda() * score );
    else
        expect = GetDerivedGappedK() * GetSearchSpace() * expf( -GetDerivedGappedLambda() * score );

    if( bitscore && 0.0f < GetExpGappedK())
        *bitscore = ( GetDerivedGappedLambda() * score - logf( GetDerivedGappedK())) / SLC_LN2;

    if( pair_expect ) {
        if( (int)GetDeltaLength() < GetQuerySize() && (int)GetDeltaLength() < GetSubjectSize())
            pairsspace = (float)(( GetQuerySize() - GetDeltaLength()) * ( GetSubjectSize() - GetDeltaLength()));
        else
            pairsspace = (float)(GetQuerySize() * GetSubjectSize());
        *pair_expect = GetK() * pairsspace * expf( -GetLambda() * score );
    }
    if( pure_expect )
        *pure_expect = GetK() * GetSearchSpace() * expf( -GetLambda() * score );

    return expect;
}

// -------------------------------------------------------------------------
// ComputeLengthAdjustment: helper method to calculate correction for length
//
template<typename TScore>
bool AbstractScoreMatrix<TScore>::ComputeLengthAdjustment(
    size_t query_len, uint64_mt db_len, size_t no_sequences )
{
    const Configuration& config = GetConfiguration(CTGapped);
    return ComputeLengthAdjustment(
        config.GetLambda(), config.GetK(), config.GetAlpha(), config.GetBeta(),
        query_len, db_len, no_sequences );
}

// -------------------------------------------------------------------------
// ComputeLengthAdjustment: compute length adjustment for edge-effect
//     correction;
// Altschul et al. in their paper in Nucleic Acids Res. 29 (2001) showed
// that expected alignment length can be linearly expressed in scores s: 
// l(s) = alpha s + beta. Substituting this expression into the
// e-value expression K (m-l) (n-Nl) exp( -lambda s ) and solving it for l
// given e-value equals 1 gives
//
//  l = alpha / lambda ln (K (m-l)(n-Nl)) + beta,
//
// where m is the query length, n is the DB length and N is the number of
// sequences in the DB. The equation is constrained by 
// K (m-l)(n-Nl) >= max{m,n}; also, the solution has to be an integer 
// number.
//
// Choosing E-value = 1: The greater this number is, the smaller 
// value of l is obtained and consequently, greater e-value in result is 
// obtained. Setting e-value to one corresponds to a p-value of 0.6.
// The algorithm is as given by E. Gertz (2005)
// -------------------------------------------------------------------------
//
template<typename TScore>
bool AbstractScoreMatrix<TScore>::ComputeLengthAdjustment(
    float lambda, float K, float alpha, float beta,
    size_t m/*query_len*/, uint64_mt n/*db_len*/, size_t N/*no_sequences*/ )
{
    MYMSG( "AbstractScoreMatrix::ComputeLengthAdjustment", 5 );

    SetSearchSpace( m * n );
    SetDeltaLength( 0 );

    if( lambda <= 0.0f || K <= 0.0f || N < 1 )
        return false;

    const int       maxit = LENGTH_ADJUSTMENT_MAXIT;
    const float     logK = logf( K );
    const float     alratio = alpha / lambda;
    const float     max_mn = (float)(( m > n )? m: n);

    float       length = 0.0f;      //expected alignment length
    float       min_length = 0.0f;  //lower bound of interval of alignment length
    float       max_length = 0.0f;  //upper bound of interval of alignment length
    float       nxt_length = 0.0f;  //alignment length at the next iteration
    float       sspace = 0.0f;      //effective search space
    float       space2 = 0.0f;      //effective search space
    bool        converged = false;  //iteration converged

    //upper bound of interval: length satisfies the constraint;
    //we have quadratic equation for l: mn - max{m,n}/K - (n + mN)l + Nl*l
    float       a = (float)N;		//coefficient a of the quadratic equation
    float       b = (float)(n + m * N);//coefficient -b
    float       c = (float)(m * n) - max_mn / K;

    if( c < 0.0f ) {
        MYMSGBEGl(3)
            char msgbuf[ONEK];
            sprintf( msgbuf, "AbstractScoreMatrix::ComputeLengthAdjustment: "
                    "search_space= %llu", GetSearchSpace());
            MYMSG( msgbuf, 3 );
        MYMSGENDl
        return false;
    }

    //take the smaller root of the equation
    //max_length = ( b - sqrtf( SQUARE(b) - 4.0f*a*c )) / ( 2.0f*a );
    //since *b>>4ac is possible, use alternative form for a solution
    max_length = 2.0f*c / ( b + sqrtf( SQUARE(b) - 4.0f*a*c ));

    for( int j = 0; j < maxit; j++ ) {
        length = nxt_length;
        sspace = ( m - length ) * ( n - N * length );
        nxt_length = alratio * ( logK + logf( sspace )) + beta;

        if( length <= nxt_length ) {
            min_length = length;
            if( nxt_length - length <= 1.0f ) {
                converged = true;
                break;
            }
            if( min_length == max_length )
                break;
        } else
            max_length = length;

        if( nxt_length < min_length || max_length < nxt_length )
            //outside the range
            nxt_length = ( !j )? max_length : ( min_length + max_length ) * 0.5f;
    }

    if( converged ) {
        //make sure that floor(min_length) + 1 != floor(min_length)
        length = ceilf(min_length);
        if( length <= max_length ) {
            space2 = ( m - length ) * ( n - N * length );
            nxt_length = alratio * ( logK + logf(space2)) + beta;

            if( length <= nxt_length ) {
                min_length = length;
                sspace = space2;
            }
        }
    }
    //if not converged, save the closest value to the solution
    SetSearchSpace((uint64_mt)sspace );
    SetDeltaLength((size_t)min_length );

    MYMSGBEGl(3)
        char msgbuf[ONEK];
        sprintf( msgbuf, "AbstractScoreMatrix::ComputeLengthAdjustment: "
                "search_space= %llu delta_length= %zu", GetSearchSpace(), GetDeltaLength());
        MYMSG( msgbuf, 3 );
    MYMSGENDl

    return converged;
}



// =========================================================================
// PRINT ROUTINES
//
// PrintReferenceParameterTable: print refrence parameter table to string
// stream; space for stream must be PRE-ALLOCATED!
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::PrintReferenceParameterTable( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//to ensure appending to the stream end
    PrintReferenceParameterTable( &string_print, sp );
}

// PrintReferenceParameterTable: print reference parameter table to file
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::PrintReferenceParameterTable( FILE* fp ) const
{
    PrintReferenceParameterTable( &file_print, fp );
}

// PrintReferenceParameterTable: print reference values of statistical
// pareameters
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::PrintReferenceParameterTable( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;
    print_func( vpn, "%-20s  %-6s   %-6s%s", "Reference values of", "K", "Lambda", NL );
    print_func( vpn, "%-20s  %6.4f   %6.4f%s", "Ungapped", GetRefK(), GetRefLambda(), NL );
    print_func( vpn, "%-20s  %6.4f   %6.4f%s", "Gapped",   GetExpGappedK(), GetExpGappedLambda(), NL );
}

// -------------------------------------------------------------------------
// PrintParameterTable: print parameter table to string stream; 
// space for stream must be PRE-ALLOCATED!
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::PrintParameterTable( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//printing at the stream end
    PrintParameterTable( &string_print, sp );
}

// PrintParameterTable: print parameter table to file
//
template<typename TScore>
void AbstractScoreMatrix<TScore>::PrintParameterTable( FILE* fp ) const
{
    PrintParameterTable( &file_print, fp );
}

// // -------------------------------------------------------------------------
// // PrintProbabilities: print score probabilities if available
// //
// // PrintFinal: final print to string stream; space for stream must be
// //     PRE-ALLOCATED before!
// //
// 
// void AbstractScoreMatrix::PrintFinal( char* sp ) const
// {
//     if( !sp )
//         return;
//     *sp = 0;  //to ensure printing to the end of the stream
//     PrintFinal( &string_print, sp );
// }
// 
// // PrintFinal: final print to file
// //
// 
// void AbstractScoreMatrix::PrintFinal( FILE* fp ) const
// {
//     PrintFinal( &file_print, fp );
// }
// 

// =========================================================================

template class AbstractScoreMatrix<int>;
template class AbstractScoreMatrix<float>;
