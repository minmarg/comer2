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
#include "liblib/msg.h"
#include "libpro/srcpro/Configuration.h"
#include "ScoresAttr.h"
#include "AbstractScoreMatrix.h"
#include "IntegerScoreMatrix.h"

// -------------------------------------------------------------------------
// constructor
//
template<typename TScore>
IntegerScoreMatrix<TScore>::IntegerScoreMatrix(
        typename AbstractScoreMatrix<TScore>::TType type,
        Configuration config[],
        typename AbstractScoreMatrix<TScore>::TStats stats,
        int precscale )
:
    AbstractScoreMatrix<TScore>( type, config, stats, precscale ),
    scores_( NULL ),
    modscores_( NULL ),
    allnegatives_( false ),
    attr_( NULL )
{
}

// -------------------------------------------------------------------------
// Init: data initialization
//
template<typename TScore>
void IntegerScoreMatrix<TScore>::Init( int querylen, int sbjctlen )
{
    AbstractScoreMatrix<TScore>::Init( querylen, sbjctlen );

    scores_ = ( TScore** )malloc( sizeof(TScore*) * GetSubjectSize());
    if( !scores_ )
        throw MYRUNTIME_ERROR("IntegerScoreMatrix::Init: Not enough memory.");

    for( int m = 0; m < GetSubjectSize(); m++ ) {
        scores_[m] = ( TScore* )malloc( sizeof(TScore) * GetQuerySize());
        if( !scores_[m])
            throw MYRUNTIME_ERROR("IntegerScoreMatrix::Init: Not enough memory.");
        memset( scores_[m], 0, sizeof(TScore) * GetQuerySize());
    }

    if( AbstractScoreMatrix<TScore>::GetUseModScores())
    {
        modscores_ = ( TScore** )malloc( sizeof(TScore*) * GetSubjectSize());
        if( !modscores_ )
            throw MYRUNTIME_ERROR("IntegerScoreMatrix::Init: Not enough memory.");

        for( int m = 0; m < GetSubjectSize(); m++ ) {
            modscores_[m] = ( TScore* )malloc( sizeof(TScore) * GetQuerySize());
            if( !modscores_[m])
                throw MYRUNTIME_ERROR("IntegerScoreMatrix::Init: Not enough memory.");
            memset( modscores_[m], 0, sizeof(TScore) * GetQuerySize());
        }
    }

    NewAttr();
}

// -------------------------------------------------------------------------
// default constructor
//
template<typename TScore>
IntegerScoreMatrix<TScore>::IntegerScoreMatrix()
:
    AbstractScoreMatrix<TScore>(),
    scores_( NULL ),
    modscores_( NULL ),
    allnegatives_( false ),
    attr_( NULL )
{
    throw MYRUNTIME_ERROR("IntegerScoreMatrix::IntegerScoreMatrix: "
                "Default initialization is prohibited.");
}

// -------------------------------------------------------------------------
// destructor
//
template<typename TScore>
IntegerScoreMatrix<TScore>::~IntegerScoreMatrix()
{
    DeleteAttr();
    if( scores_ ) {
        for( int m = 0; m < GetSubjectSize(); m++ )
            if( scores_[m] )
                free( scores_[m] );
        free( scores_ );
        scores_ = NULL;
    }
    if( modscores_ ) {
        for( int m = 0; m < GetSubjectSize(); m++ )
            if( modscores_[m] )
                free( modscores_[m] );
        free( modscores_ );
        modscores_ = NULL;
    }
}

// -------------------------------------------------------------------------
// ScanForHSPs: perform search of multiple high-scoring pairs (HSPs)
//     in the same diagonal of the score system;
// minhspscore, minimum HSP score;
// hsplen, HSP length;
// nohsps, minimum number of HSPs required;
// maxdist, max distance between HSPs;
// possbjct, posquery, subject and query positions of HSPs found in the
//     same diagonal (to be returned);
// return true if such HSPs have been found
//
template<typename TScore>
bool IntegerScoreMatrix<TScore>::ScanForHSPs(
    float minhspscore, int hsplen, int nohsps, int maxdist,
    int* possbjct, int* posquery )
{
    bool bret = false;

    if( !GetAttr())
        throw MYRUNTIME_ERROR( "IntegerScoreMatrix::ScanForHSPs: Null scores attributes." );

    bret = GetAttr()->SearchForHSPs(
        (TScore)rint( minhspscore * GetScoresFactor()),
        hsplen, nohsps, maxdist,
        possbjct, posquery );

    return bret;
}

// -------------------------------------------------------------------------
// ComputeStatisticalParameters: calculate statistical parameters
//
template<typename TScore>
void IntegerScoreMatrix<TScore>::ComputeStatisticalParameters( bool computelambda )
{
    int fct = GetScoresFactor();

    if( !GetAttr())
        throw MYRUNTIME_ERROR(
            "IntegerScoreMatrix::ComputeStatisticalParameters: Null scores attributes." );

    GetAttr()->ComputeStatisticalParameters( computelambda, true );

    if( 0.0f < GetAttr()->GetLambda())
        SetLambda(      GetAttr()->GetLambda() * (float)fct );
    else
        SetLambda(      GetAttr()->GetLambda());
    SetEntropy(         GetAttr()->GetH());
    SetK(               GetAttr()->GetK());
    SetExpectedScore(   GetAttr()->GetExpectedScore() / (float)fct );

    SetMinScore(( TScore )rintf((float)GetAttr()->GetMinScore() / (float)fct ));
    SetMaxScore(( TScore )rintf((float)GetAttr()->GetMaxScore() / (float)fct ));

    AbstractScoreMatrix<TScore>::DeriveGappedParameters();
}

// -------------------------------------------------------------------------
// ScaleScoreMatrix: scale score matrix to have reference composition 
// between profiles
//
template<typename TScore>
void IntegerScoreMatrix<TScore>::ScaleScoreMatrix()
{
    if( AbstractScoreMatrix<TScore>::GetStats() == 
        AbstractScoreMatrix<TScore>::StatisticsGiven ) {
        ComputeStatisticalParameters();//NOTE!!
        //Nothing doing since all statistical parameters are assumed to be 
        // known from configuration
        return;
    }

    //NOTE: score_multiplier_ is assumed to be set appropriately
    ScaleScoreMatrixHelper();
    ComputeStatisticalParameters( false/*computelambda*/);
// PrintScoringMatrix( stderr );//***TEST***
}

// -------------------------------------------------------------------------
// ScaleScoreMatrixHelper: helper method to scale matrix; 
// perform iterative scaling until the lambda of the score matrix is 
// equal to the reference value;
// binary search for lambda is much more precise but much slower; 
// otherwise, operating with integer scores induces errors, and the 
// multiplication of scores by the ratio of lambdas rarely leads to an exact 
// solution.
//
template<typename TScore>
void IntegerScoreMatrix<TScore>::ScaleScoreMatrixHelper()
{
    const mystring preamb = "IntegerScoreMatrix::ScaleScoreMatrixHelper: ";
    const int fct = GetScoresFactor();

    if( GetQuerySize() < 1 || GetSubjectSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Zero matrix size." );

    if( GetAllNegatives())
        throw MYRUNTIME_ERROR2( preamb + "Scores are all negative.", SCALING );

    if( !GetAttr())
        throw MYRUNTIME_ERROR( preamb + "Null scores attributes." );

    if( fct < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid scores factor.");

#ifdef SCALE_BINARY_SEARCH
    const bool lbUSE_BOUND_PRECOMPUTATION = false;//true;

    float   prev_lambda;
    float   diff_target = 1.0f;//difference of lambdas
    float   min_diff_target = 10.0f;//minimum value of the difference
    float   multiplier = 1.0f;
    float   prev_multiplier = 1.0f;//previous value of multiplier

    GetAttr()->ComputeProbabilitiesAndLambda();

    float       right = 1.0f;//right interval bound
    float       left  = 1.0f;//left interval bound
    float       step = 0.05f;//initial step for narrowing the interval
    float       max_multiplier = 1.0f;
    const float accuracy = 1.e-4f / (float)fct;
    int i = 0;

    bool greater = GetRefLambda() < GetLambda();

    if( !lbUSE_BOUND_PRECOMPUTATION ) {

//        if( greater )   right = 10.0;
//        else            left = 0.0;
        //float lmb_l = GetRefLambda()/1.1f;
        float lmb_u = GetRefLambda()/1.0f;//0.98;
        if( greater ) {
            if( lmb_u < GetLambda())
                return;
            right = 1.02f;
        } else {
            //if( GetLambda() < lmb_l )
            //    return;
            left = 0.9f;//0.2f;
        }
//        if( greater )   return;
//        else if(1.f < multiplier) right = multiplier;
//             else left = multiplier;

    }
    else {//lbUSE_BOUND_PRECOMPUTATION
// fprintf(stderr,"init: fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());
        for( i = 0; 0.0f < GetLambda() && i < MAX_SCALE_ITERATIONS; i++, step += step ) {
            if( greater ) {
                    if( GetLambda() <= GetRefLambda())   break;
            } else  if( GetLambda() >= GetRefLambda())   break;

            if( greater ) {
                left   = right;
                right += step;
                prev_multiplier = multiplier;
                multiplier = right;
            } else {
                right = left;
                left -= step;
                prev_multiplier = multiplier;
                multiplier = left;
                if( left < 0.0f ) {
                    left = 0.0f;
                    break;
                }
            }

            prev_lambda = GetLambda();
            MultiplyScoresBy( multiplier / prev_multiplier );//NOTE
            GetAttr()->ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"bound: fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());

            if( GetLambda() < 0.0f )
                break;

            diff_target = fabsf( GetRefLambda() - GetLambda());

            if( diff_target < min_diff_target ) {
                min_diff_target = diff_target;
                max_multiplier = multiplier;
            }
            if( diff_target < accuracy )
                break;
            if( fabsf( GetLambda() - prev_lambda ) < accuracy )
                break;
        }//for
        if( 0.0f < GetLambda()) {
            if( greater ) {
                    if( GetLambda() > GetRefLambda())   right += step + step;
            } else  if( GetLambda() < GetRefLambda())   left = 0.0f;
        }
    }//if USE_BOUND_PRECOMPUTATION

    if( 0.0 < GetLambda())
        for( ; i < MAX_SCALE_ITERATIONS && accuracy < diff_target; i++ ) {
            prev_multiplier = multiplier;
            multiplier = ( left + right ) * 0.5f;

            prev_lambda = GetLambda();
            MultiplyScoresBy( multiplier / prev_multiplier );//NOTE
            GetAttr()->ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());

            if( GetLambda() < 0.0f )
                break;

            diff_target = fabsf( GetRefLambda() - GetLambda());

            if( diff_target < min_diff_target ) {
                min_diff_target = diff_target;
                max_multiplier = multiplier;
            }
            if( fabsf( GetLambda() - prev_lambda ) < accuracy )
                break;

            if( GetRefLambda() < GetLambda())
                left  = multiplier;
            else
                right = multiplier;
        }//for

    if( GetLambda() < 0.0f || min_diff_target < diff_target ) {
        MultiplyScoresBy( multiplier = max_multiplier / multiplier );//NOTE
        GetAttr()->ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"adjust: fct=%f, l=%g ref=%g\n", multiplier,GetLambda(),GetRefLambda());
    }

#else //if not defined SCALE_BINARY_SEARCH

    //sometimes, scaling floating point numbers and rounding them to nearest
    //integer can produce positive expected score per position;
    //this may happen when expected score is nearly zero. 
    //In this case, approximate integer scaling is used
    multiplier = IntegerScaleMatrix();

#endif

    SetMultiplier( multiplier );
}

// -------------------------------------------------------------------------
// IntegerScaleMatrix: perform iterative scaling until lambda equals the 
// reference value
//
template<typename TScore>
float IntegerScoreMatrix<TScore>::IntegerScaleMatrix()
{
    float   prev_lambda;
    float   diff_lambda;
    float   multiplier = 1.0f;
    int     iter = 0;
    const float least_difference = GetRefLambda() * 0.001f;

    if( !GetAttr())
        throw MYRUNTIME_ERROR("IntegerScoreMatrix::IntegerScaleMatrix: Null scores attributes.");

    GetAttr()->ComputeProbabilitiesAndLambda();

    while( 0.0f < GetLambda() && iter++ < MAX_SCALE_ITERATIONS )
    {
        prev_lambda = GetLambda();

        MultiplyScoresBy( multiplier = GetLambda() / GetRefLambda() );

        GetAttr()->ComputeProbabilitiesAndLambda( false/*wrn*/);
        diff_lambda = GetLambda() - prev_lambda;

        if( -least_difference < diff_lambda && diff_lambda < least_difference )
            break;
    }

    return multiplier;
}

// -------------------------------------------------------------------------
// MultiplyScoresBy: multiply scores in the matrix by the given factor
//
template<typename TScore>
void IntegerScoreMatrix<TScore>::MultiplyScoresBy( float multiplier )
{
    //const int fct = GetScoresFactor();
    TScore val, newval;
    int m, n;

    if( multiplier == 1.0f )
        return;

    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            val = GetScore( m, n );
            if( SCORE_MIN < val ) {
                newval = (TScore)rintf( multiplier * (float)val );
                SetScore( m, n, newval );
            }
        }
    }
}

// =========================================================================

template class IntegerScoreMatrix<int>;
