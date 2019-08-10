/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __ScoresAttr_h__
#define __ScoresAttr_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "libpro/srcpro/PMTransModel.h"

template<typename TScore>
class AbstractScoreMatrix;

////////////////////////////////////////////////////////////////////////////
// CLASS ScoresAttr
// Implementation of score attributes associated with a score table
//
template<typename TScore>
class ScoresAttr
{
public:
    typedef void    ( AbstractScoreMatrix<TScore>::*TProbabilityFunction )( ScoresAttr<TScore>* );
    typedef void    ( ScoresAttr<TScore>::*TConservationFunction )( float, float*, float*, int, void* );

    ScoresAttr(
            AbstractScoreMatrix<TScore>* prnt,
            TProbabilityFunction probfun,
            float reflambda,
            float refH,
            float refK,
            const TScore** scores,
            int scaled_by_factor = 1 );

    virtual ~ScoresAttr();

//     bool                IsValid() const;//whether score matrix is valid

    TScore              GetParentScore/*GetImageScore*/( int m, int n ) const;//score at specified positions

    //non-scaled scores for the H and K calculation
    bool                GetUsingNSscoresforHK() const { return calculate_HK_ns_; }
    void                SetUsingNSscoresforHK( bool value ) { calculate_HK_ns_ = value; }

    int                 GetQuerySize() const                        { return queryLen_;   }
    int                 GetSubjectSize() const                      { return subjectLen_; }

    float               GetRefLambda() const                        { return referenceLambda_; }
    float               GetRefH() const                             { return referenceH_; }
    float               GetRefK() const                             { return referenceK_; }

    float               GetLambda() const                           { return lambda_; }//computed ungapped Lambda
    float               GetH() const                                { return entropy_; }//computed relative entropy
    float               GetK() const                                { return parameterK_; }//computed ungapped K
    float               GetExpectedScore() const                    { return expscore_; }//expected score
    void                SetExpectedScore( float  E )                { expscore_ = E; }

    bool                GetAllNegatives() const                     { return allnegatives_; }
    void                SetAllNegatives( bool value )               { allnegatives_ = value; }

    virtual void        Init( int querylen, int sbjctlen );

    bool SearchForHSPs( TScore minhspscore, int hsplen, int nohsps, int maxdist, int* = NULL, int* = NULL );

    void                ComputeStatisticalParameters( bool = true, bool wrn = true);
    void                ComputeProbabilitiesAndLambda( bool wrn = true);

    int                 GetScoresFactor/*GetAutoScalingFactor*/() const { return scaled_by_factor_; }

//     double              GetPrivateMultiplier() const            { return private_multiplier; }
//     void                SetPrivateMultiplier( double value )    { private_multiplier = value; }

    TScore              GetMinScore() const { return min_score_; }
    TScore              GetMaxScore() const { return max_score_; }

    void                SetMinMaxScores( TScore min, TScore max );
    float               GetProbabilityOf( TScore ) const;
    void                SetProbabilityOf( TScore, float value );
    void                IncProbabilityOf( TScore, float value );
    void                DivideProbabilityOf( TScore, float value );

protected:
    explicit ScoresAttr();

    void                SetQuerySize( int value ) { queryLen_ = value;   }
    void                SetSubjectSize( int value ) { subjectLen_ = value; }

    void                Init(){}//IMPORTANT: query and subject sizes are supposed to be initialized

    void                SetRefLambda( float value )         { referenceLambda_ = value; }
    void                SetRefH( float value )              { referenceH_ = value; }
    void                SetRefK( float value )              { referenceK_ = value; }

    void                SetLambda( float newlambda )        { lambda_ = newlambda; }
    void                SetH( float H )                     { entropy_ = H; }
    void                SetK( float K )                     { parameterK_ = K; }

    //finding root and conservation equation
    virtual float       findLambdaRoot( int, float, float, float, int );
    static void         lmbd_FdFfunc( float, float*, float*, void* = NULL );
    //

    void                ComputeGCD();//greatest common divisor

    virtual void        ComputeScoreProbabilities();
    virtual void        ComputeScoreProbabilitiesHelper();
    void                ComputeScoreProbabilitiesNS();
    void                ComputeLambda( bool wrn = true );//compute scaling parameter lambda
    void                ComputeKHGivenLambda();
    virtual void        ComputeEntropyGivenLambda();//compute relative entropy given lambda
    virtual void        ComputeKarlinsK();//compute parameter K

    void                NewProbabilities( size_t );//allocate memory for probabilities
    void                DestroyProbabilities();//deallocate memory for probabilities

    TScore              GetGCD() const { return score_gcd_; }
    void                SetGCD( TScore value ) { score_gcd_ = value; }

    void                SetMinScore( TScore value ) { min_score_ = value; }
    void                SetMaxScore( TScore value ) { max_score_ = value; }

    void                PrintProbabilities( FILE* );

    //{{Member functions for non-scaled scores...
    void                NewProbabilitiesNS( size_t );
    void                DestroyProbabilitiesNS();
    float               GetProbabilityOfNS( TScore ) const;
    void                IncProbabilityOfNS( TScore, float value );

    void                ComputeGCDNS();
    TScore              GetGCDNS() const { return score_gcd_ns_; }
    void                SetGCDNS( TScore value ) { score_gcd_ns_ = value; }

    void                SetMinMaxScoresNS( TScore min, TScore max );
    TScore              GetMinScoreNS() const { return min_score_ns_; }
    void                SetMinScoreNS( TScore value ) { min_score_ns_ = value; }

    TScore              GetMaxScoreNS() const { return max_score_ns_; }
    void                SetMaxScoreNS( TScore value ) { max_score_ns_ = value; }
    //}}

    //Private routines for the manipulation of probability vector used to compute K
    void                NewPrivateProbVector( size_t maxit, size_t range );
    void                DestroyPrivateProbVector();
    void                ResetPrivateProbVector( size_t = 0 );
    float*              GetPrivateProbVector() const        { return priv_prob_vector_; }
    size_t              GetPrivateProbVectorSize() const    { return prob_vector_size_; }
    //...non-scaled scores...
    void                NewPrivateProbVectorNS( size_t maxit, size_t range );
    void                DestroyPrivateProbVectorNS();
    void                ResetPrivateProbVectorNS( size_t = 0 );
    float*              GetPrivateProbVectorNS() const { return priv_prob_vector_ns_; }
    size_t              GetPrivateProbVectorSizeNS() const { return prob_vector_size_ns_; }
    //End of private routines.

protected:
    AbstractScoreMatrix<TScore>*        GetParent()             { return parent_; }
    const AbstractScoreMatrix<TScore>*  GetParent() const       { return parent_; }
    TProbabilityFunction    GetParentProbFunction()     { return parent_probfunction_; }

    const TScore**      GetParentScores() const         { return parent_scores_; }

private:
    AbstractScoreMatrix<TScore>* parent_;//Parent score matrix
    TProbabilityFunction    parent_probfunction_;//Parent probability function
    const TScore**          parent_scores_;     //parent scores
//     float           private_multiplier;  
    const int               scaled_by_factor_;
    bool                    calculate_HK_ns_;   //calculate parameters H and K using non-scaled scores

    float*                  probabilities_;     //score probabilities
    TScore                  min_score_;         //minimum score
    TScore                  max_score_;         //maximum score
    TScore                  score_gcd_;         //the greatest common divisor of scores

    float*                  priv_prob_vector_;  //private probability vector
    size_t                  prob_vector_size_;  //size of private probability vector

    //{{attributes for non-scaled scores
    float*                  probabilities_ns_;  //probabilities of non-scaled scores 
    TScore                  min_score_ns_;      //minimum of non-scaled scores
    TScore                  max_score_ns_;      //maximum of non-scaled scores
    TScore                  score_gcd_ns_;      //the greatest common divisor of non-scaled scores

    float*                  priv_prob_vector_ns_;//private probability vector for non-scaled scores
    size_t                  prob_vector_size_ns_;//size of private probability vector for non-scaled scores
    //}}

    float                   expscore_;          //expected score per column pair
    bool                    allnegatives_;      //all scores negative

    float                   referenceLambda_;   //reference lambda parameter
    float                   referenceH_;        //reference parameter H
    float                   referenceK_;        //reference parameter K

    float                   lambda_;            //scaling parameter lambda
    float                   entropy_;           //entropy describing information per aligned pair of positions (parameter H)
    float                   parameterK_;        //parameter K

    int                     queryLen_;          //length of query profile
    int                     subjectLen_;        //length of subject profile
};

// -------------------------------------------------------------------------
// INLINES ...

// inline bool ScoresAttr::IsValid() const
// {
//     if( GetKeepInMemory())
//         return GetScores() != NULL && 0 < GetQuerySize() && 0 < GetSubjectSize();
//     return true;
// }

// -------------------------------------------------------------------------
// GetParentScore: get parent score at profile positions
//
template<typename TScore>
inline
TScore ScoresAttr<TScore>::GetParentScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !parent_scores_ || subjectLen_ <= m || m < 0 || queryLen_ <= n || n < 0 )
        throw MYRUNTIME_ERROR("ScoresAttr::GetParentScore: Memory access error.");
#endif

    return parent_scores_[m][n];
}

// -------------------------------------------------------------------------
// SetMinMaxScores: set min max scores and allocates space for probability
//     vector if needed
//
template<typename TScore>
inline
void ScoresAttr<TScore>::SetMinMaxScores( TScore min, TScore max )
{
    if( max < min )
        throw MYRUNTIME_ERROR( "ScoresAttr::SetMinMaxScores: max score < min score." );

    NewProbabilities( max - min + 1 );

    SetMinScore( min );
    SetMaxScore( max );
}


// -------------------------------------------------------------------------
// NewProbabilities: allocate memory for probabilities
//
template<typename TScore>
inline
void ScoresAttr<TScore>::NewProbabilities( size_t size )
{
    DestroyProbabilities();

    const mystring preamb = "ScoresAttr::NewProbabilities: ";

    if( prob_vector_size_ < size ) {
        char strbuf[BUF_MAX];
        sprintf( strbuf, "%zu > %zu", size, prob_vector_size_ );
        throw MYRUNTIME_ERROR( preamb + "Size exceeds the maximum allowed: " + strbuf );
    }

    probabilities_ = ( float* )malloc( sizeof(float) * size );

    if( !probabilities_ ) 
        throw MYRUNTIME_ERROR( preamb + "Not enough memory." );

    memset( probabilities_, 0, sizeof(float) * size );
}

// -------------------------------------------------------------------------
// DestroyProbabilities: deallocate memory allocated for probabilities
//
template<typename TScore>
inline
void ScoresAttr<TScore>::DestroyProbabilities()
{
    if( probabilities_ ) {
        free( probabilities_ );
        probabilities_ = NULL;
    }
}

// -------------------------------------------------------------------------
// GetProbabilityOf: get score probability
//
template<typename TScore>
inline
float ScoresAttr<TScore>::GetProbabilityOf( TScore score ) const
{
#ifdef __DEBUG__
    if( !probabilities_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::GetProbabilityOf: Memory access error." );

    if( score < min_score_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::GetProbabilityOf: Memory access error." );
#endif
    return probabilities_[score-min_score_];
}

// -------------------------------------------------------------------------
// SetProbabilityOf: set probability of score
//
template<typename TScore>
inline
void ScoresAttr<TScore>::SetProbabilityOf( TScore score, float value )
{
#ifdef __DEBUG__
    if( !probabilities_ )
        throw myruntime_error( "ScoresAttr::SetProbabilityOf: Memory access error." );

    if( score < min_score_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::SetProbabilityOf: Memory access error." );
#endif
    probabilities_[score-min_score_] = value;
}

// -------------------------------------------------------------------------
// IncProbabilityAt: increment score probability by value
//
template<typename TScore>
inline
void ScoresAttr<TScore>::IncProbabilityOf( TScore score, float value )
{
#ifdef __DEBUG__
    if( !probabilities_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::IncProbabilityOf: Memory access error." );

    if( score < min_score_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::IncProbabilityOf: Memory access error." );
#endif
    probabilities_[score-min_score_] += value;
}

// -------------------------------------------------------------------------
// DivideProbabilityAt: divide score probability by value
//
template<typename TScore>
inline
void ScoresAttr<TScore>::DivideProbabilityOf( TScore score, float value )
{
#ifdef __DEBUG__
    if( !probabilities_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::DivideProbabilityOf: Memory access error." );

    if( score < min_score_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::DivideProbabilityOf: Memory access error." );
#endif
    if( value == 0.0f )
        throw MYRUNTIME_ERROR( "ScoresAttr::DivideProbabilityOf: Illegal operation." );

    if( probabilities_[score-min_score_])
        probabilities_[score-min_score_] /= value;
}


// -------------------------------------------------------------------------
// SetMinMaxScoresNS: set min max non-scaled scores and allocates space for 
// probability vector
//
template<typename TScore>
inline
void ScoresAttr<TScore>::SetMinMaxScoresNS( TScore min, TScore max )
{
    if( max < min )
        throw MYRUNTIME_ERROR( "ScoresAttr::SetMinMaxScoresNS: max score < min score." );

    NewProbabilitiesNS( max - min + 1 );

    SetMinScoreNS( min );
    SetMaxScoreNS( max );
}

// -------------------------------------------------------------------------
// NewProbabilitiesNS: allocate memory for probabilities for non-scaled 
// scores
//
template<typename TScore>
inline
void ScoresAttr<TScore>::NewProbabilitiesNS( size_t size )
{
    DestroyProbabilitiesNS();

    const mystring preamb = "ScoresAttr::NewProbabilitiesNS: ";

    if( prob_vector_size_ns_ < size ) {
        char strbuf[BUF_MAX];
        sprintf( strbuf, "%zu > %zu", size, prob_vector_size_ns_ );
        throw MYRUNTIME_ERROR( preamb + "Size exceeds the maximum allowed: " + strbuf );
    }

    probabilities_ns_ = ( float* )malloc( sizeof(float) * size );

    if( !probabilities_ns_ ) 
        throw MYRUNTIME_ERROR( preamb + "Not enough memory." );

    memset( probabilities_ns_, 0, sizeof(float) * size );
}

// -------------------------------------------------------------------------
// DestroyProbabilitiesNS: deallocate memory allocated for probabilities 
// used for non-scaled scores
//
template<typename TScore>
inline
void ScoresAttr<TScore>::DestroyProbabilitiesNS()
{
    if( probabilities_ns_ ) {
        free( probabilities_ns_ );
        probabilities_ns_ = NULL;
    }
}

// -------------------------------------------------------------------------
// GetProbabilityOfNS: get probability of a non-scaled score
//
template<typename TScore>
inline
float ScoresAttr<TScore>::GetProbabilityOfNS( TScore score ) const
{
#ifdef __DEBUG__
    if( !probabilities_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::GetProbabilityOfNS: Memory access error." );

    if( score < min_score_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::GetProbabilityOfNS: Memory access error." );
#endif
    return probabilities_ns_[score-min_score_ns_];
}

// -------------------------------------------------------------------------
// IncProbabilityOfNS: increment the probability of a non-scaled score by 
// value
//
template<typename TScore>
inline
void ScoresAttr<TScore>::IncProbabilityOfNS( TScore score, float value )
{
#ifdef __DEBUG__
    if( !probabilities_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::IncProbabilityOfNS: Memory access error." );

    if( score < min_score_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::IncProbabilityOfNS: Memory access error." );
#endif
    probabilities_ns_[score-min_score_ns_] += value;
}


// -------------------------------------------------------------------------
// InitializePrivateProbVector: allocate memory for private probability 
// vector
//
template<typename TScore>
inline
void ScoresAttr<TScore>::NewPrivateProbVector( size_t maxit, size_t range )
{
    DestroyPrivateProbVector();

    size_t size = maxit * range + 1;//plus to include the zero value

    priv_prob_vector_ = ( float* )malloc( sizeof(float) * size );

    if( !priv_prob_vector_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::NewPrivateProbVector: Not enough memory." );

    prob_vector_size_ = size;
    ResetPrivateProbVector();
}

// -------------------------------------------------------------------------
// DestroyPrivateProbVector: deallocate memory allocated for private 
// probability vector
//
template<typename TScore>
inline
void ScoresAttr<TScore>::DestroyPrivateProbVector()
{
    if( priv_prob_vector_ ) {
        free( priv_prob_vector_ );
        priv_prob_vector_ = NULL;
    }
    prob_vector_size_ = 0;
}

// -------------------------------------------------------------------------
// ResetPrivateProbVector: reset values of private probability vector
//
template<typename TScore>
inline
void ScoresAttr<TScore>::ResetPrivateProbVector( size_t size )
{
#ifdef __DEBUG__
    if( !priv_prob_vector_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::ResetPrivateProbVector: Memory access error." );

    if( prob_vector_size_ < size ) {
        warning( "ScoresAttr::ResetPrivateProbVector: Wrong size." );
        size = prob_vector_size_;
    }
#endif
    if( size == 0 )
        size = prob_vector_size_;
    memset( priv_prob_vector_, 0, sizeof(float) * size );
}


// -------------------------------------------------------------------------
// NewPrivateProbVectorNS: allocate memory for private probability 
// vector for non-scaled scores
//
template<typename TScore>
inline
void ScoresAttr<TScore>::NewPrivateProbVectorNS( size_t maxit, size_t range )
{
    DestroyPrivateProbVectorNS();

    size_t size = maxit * range + 1;//plus to include the zero value

    priv_prob_vector_ns_ = ( float* )malloc( sizeof(float) * size );

    if( !priv_prob_vector_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::NewPrivateProbVectorNS: Not enough memory." );

    prob_vector_size_ns_ = size;
    ResetPrivateProbVectorNS();
}

// -------------------------------------------------------------------------
// DestroyPrivateProbVectorNS: deallocate memory allocated for private 
// probability vector for non-scaled scores
//
template<typename TScore>
inline
void ScoresAttr<TScore>::DestroyPrivateProbVectorNS()
{
    if( priv_prob_vector_ns_ ) {
        free( priv_prob_vector_ns_ );
        priv_prob_vector_ns_ = NULL;
    }
    prob_vector_size_ns_ = 0;
}

// -------------------------------------------------------------------------
// ResetPrivateProbVectorNS: reset values of private probability vector 
// used for non-scaled scores
//
template<typename TScore>
inline
void ScoresAttr<TScore>::ResetPrivateProbVectorNS( size_t size )
{
#ifdef __DEBUG__
    if( !priv_prob_vector_ns_ )
        throw MYRUNTIME_ERROR( "ScoresAttr::ResetPrivateProbVectorNS: Memory access error." );

    if( prob_vector_size_ns_ < size ) {
        warning( "ScoresAttr::ResetPrivateProbVectorNS: Wrong size." );
        size = prob_vector_size_ns_;
    }
#endif
    if( size == 0 )
        size = prob_vector_size_ns_;
    memset( priv_prob_vector_ns_, 0, sizeof(float) * size );
}

#endif//__ScoresAttr_h__
