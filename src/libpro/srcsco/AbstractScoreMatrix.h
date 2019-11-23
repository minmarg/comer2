/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __AbstractScoreMatrix_h__
#define __AbstractScoreMatrix_h__

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <map>
#include <utility>

#include "libpro/srcpro/Configuration.h"
// #include "libpro/srcpro/GapScheme.h"
// #include "libpro/srcpro/DistributionMatrix.h"

template<typename TScore>
class ScoresAttr;

////////////////////////////////////////////////////////////////////////////
// CLASS AbstractScoreMatrix
// Base class for implementations of Score Matrices
//
template<typename TScore>
class AbstractScoreMatrix
{
public:
    enum TType {
            PositionSpecific,
            ProfileSpecific,
            ProfSpecLSO,
            HDPProfileSpecific,
            HDP0CtxProfileSpecific,
            CuBatchLSO,
            NoType
    };

    enum TStats {
            StatisticsGiven,
            ComputeStatistics
    };

//     enum TScaling {
//             NoScaling = 0,
//             AutoScalling = 1,
//             FPScaling = 2
//     };

public:
    AbstractScoreMatrix( TType, Configuration config[NoCTypes], TStats, int precscale = 1 );
    virtual ~AbstractScoreMatrix();

    TType               GetType() const { return type_; }
    TStats              GetStats/*GetBehaviour*/() const { return stats_; }
//     TMask               GetMaskingApproach() const      { return commonmasking; }

//     virtual const char* GetMethodName() const = 0;

    const char*         GetName() const { return name_; }
    void                SetName( const char* value ) { name_ = value; }

//     double              GetInformationThreshold() const { return information_content; }
//     void                SetInfoThresholdByEval( double eval );

//     bool                GetAutocorrectionPositional() const             { return  positcorrects; }
//     void                SetAutocorrectionPositional( bool value )       { positcorrects = value; }

//     double              GetInfoCorrectionUpperBound2nd() const          { return ainfoupper2nd; }
//     void                SetInfoCorrectionUpperBound2nd( double value )  { ainfoupper2nd = value; }

//     double              GetInfoCorrectionNumerator2nd() const           { return ainfonumerator2nd; }
//     void                SetInfoCorrectionNumerator2nd( double value )   { ainfonumerator2nd = value; }

//     double              GetInfoCorrectionScale2nd() const               { return ainfoscale2nd; }
//     void                SetInfoCorrectionScale2nd( double value )       { ainfoscale2nd = value; }

//     double              GetInfoCorrectionNumeratorAlt() const           { return ainfonumeratoralt; }
//     void                SetInfoCorrectionNumeratorAlt( double value )   { ainfonumeratoralt = value; }

//     double              GetInfoCorrectionScaleAlt() const               { return ainfoscalealt; }
//     void                SetInfoCorrectionScaleAlt( double value )       { ainfoscalealt = value; }

    bool                GetSupportOptimFreq() const { return supportoptfreq_; }
    virtual void        OptimizeTargetFrequencies() { return; }

    bool                GetScoreAdjRelEnt() const { return scadjrelent_; }
    void                SetScoreAdjRelEnt( bool value ) { scadjrelent_ = value; }

    bool                GetScoreAdjCorrel() const { return scadjcorrsc_; }
    void                SetScoreAdjCorrel( bool value ) { scadjcorrsc_ = value; }


//     bool                GetFPScaling() const            { return auto_scaling == FPScaling; }
//     bool                GetAutoScaling() const          { return auto_scaling == AutoScalling; }
    virtual int         GetScoresFactor/*GetAutoScalingFactor*/() const { return scaled_by_factor_; }

                                                                    //returns score at specified profile positions
    virtual TScore      GetScore( int m, int n ) const = 0;// { return GetScoreBase( m, n ); }
//     double              GetImageScore( int m, int n ) const;        //returns image score at specified profile positions
//     double              GetModImageScore( int m, int n ) const;
    virtual TScore      GetModScore( int m, int n ) const = 0;

    bool                GetUseModScores/*GetUseModImage*/() const { return usemodimage_; }
    void                SetUseModScores/*SetUseModImage*/( bool value ) { usemodimage_ = value; }

//     TMask               GetMasked( int m, int n  ) const;           //get mask value at

    virtual int         GetQuerySize() const { return queryLen_; }
    virtual int         GetSubjectSize() const { return subjectLen_; }

    virtual bool        ScanForHSPs( float /*minhspscore*/, int /*hsplen*/, int /*nohsps*/, int /*maxdist*/, 
                                     int* = NULL, int* = NULL ) { return false; }

    float               GetRefLambda() const            { return referenceLambda_; }
    float               GetRefH() const                 { return referenceH_; }
    float               GetRefK() const                 { return referenceK_; }
    float               GetExpGappedLambda() const      { return experimentalGappedLambda_; }
    float               GetExpGappedK() const           { return experimentalGappedK_; }
    float               GetDerivedGappedLambda() const  { return derivedGappedLambda_; }
    float               GetDerivedGappedK() const       { return derivedGappedK_; }

    float               GetLambda() const               { return lambda_; }     //computed ungapped Lambda
    float               GetEntropy() const              { return entropy_; }    //computed relative entropy
    float               GetK() const                    { return parameterK_; } //computed ungapped parameter K
    float               GetExpectedScore() const        { return expscore_; }   //expected score

    TScore              GetMinScore() const             { return min_score_; }
    TScore              GetMaxScore() const             { return max_score_; }

    float               GetMultiplier() const           { return score_multiplier_; }

    virtual void        ComputeScoreMatrix/*ComputeProfileScoringMatrix*/( bool final = false ) = 0;

//     void                PostScalingProc(
//                                 const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
//                                 GapScheme& qg, GapScheme& sg,
//                                 bool autoc, int acwindow, bool firstpass );

    virtual float       GetFinalScore( TScore value ) const = 0;

    virtual void        ScaleScoreMatrix/*ScaleScoringMatrix*/() = 0;
    virtual void        ComputeStatisticalParameters( bool = true ) = 0;
    virtual float       ComputeExpectation( float, float* = NULL, float* = NULL, float* = NULL, float* = NULL ) const;

    virtual void        ComputeScoreProbabilities( ScoresAttr<TScore>* ) = 0;

    void                PrintParameterTable( char* ) const;//print calculated statistical parameters
    void                PrintParameterTable( FILE* ) const;//print calculated statistical parameters
    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const = 0;
//     void                PrintFinal( char* ) const;
//     void                PrintFinal( FILE* ) const;
//     virtual void        PrintFinal( TPrintFunction, void* vpn ) const = 0;

    virtual void        PrintScoreMatrix/*PrintScoringMatrix*/( FILE* ) {;};

    //public methods...
    bool        ComputeLengthAdjustment(//compute length adjustment to correct for edge-effects
                size_t query_len, uint64_mt db_len, size_t no_sequences );

    size_t      GetDeltaLength() const { return deltaLength_; }//correction for length
    uint64_mt   GetSearchSpace() const { return searchSpace_; }

    void                PrintReferenceParameterTable( char* ) const;
    void                PrintReferenceParameterTable( FILE* ) const;
    void                PrintReferenceParameterTable( TPrintFunction, void* vpn ) const;

protected:
    explicit AbstractScoreMatrix();

    virtual void        Init( int querylen, int sbjctlen );
    void                InitializeStatParams/*InitializeSSParameters*/();

//     virtual double      GetScoreBase( int m, int n ) const;         //returns score at specified profile positions

    virtual void        SetScore( int m, int n, TScore ) = 0;
    virtual void        SetModScore( int m, int n, TScore ) = 0;

//     void                SetInformationThreshold( double value )     { information_content = value; }

//     bool                GetMaskedUnmasked( int m, int n  ) const    { return GetMasked( m, n  ) == Unmasked; }
//     bool                GetMaskedToIgnore( int m, int n  ) const    { return GetMasked( m, n  ) == MaskToIgnore; }
//     bool                GetMaskedToConsider( int m, int n  ) const  { return GetMasked( m, n  ) == MaskToConsider; }

//     void                SetMasked( TMask, int m, int n  );          //set mask value at

    const Configuration& GetConfiguration( TConfigType ps ) const;

    virtual bool        GetAllNegatives() const = 0;
    virtual void        SetAllNegatives( bool value ) = 0;

    void                SetMinScore( TScore value ) { min_score_ = value; }
    void                SetMaxScore( TScore value ) { max_score_ = value; }

    void                SetMultiplier( float value ) { score_multiplier_ = value; }

    void                SetSupportOptimFreq( bool value ) { supportoptfreq_ = value; }

private:
    void                SetRefLambda( float value )             { referenceLambda_ = value; }
    void                SetRefH( float value )                  { referenceH_ = value; }
    void                SetRefK( float value )                  { referenceK_ = value; }

public:
    void                SetExpGappedLambda( float value )       { experimentalGappedLambda_ = value; }
    void                SetExpGappedK( float value )            { experimentalGappedK_ = value; }

private:
    void                SetDerivedGappedLambda( float value )   { derivedGappedLambda_ = value; }
    void                SetDerivedGappedK( float value )        { derivedGappedK_ = value; }

protected:
    void                SetLambda( float newlambda )            { lambda_ = newlambda; }
    void                SetEntropy( float H )                   { entropy_ = H; }
    void                SetK( float K )                         { parameterK_ = K; }
    void                SetExpectedScore( float E )             { expscore_ = E; }

    void                DeriveGappedParameters();


private:
    //protected methods...
    bool        ComputeLengthAdjustment(//compute length adjustment to correct for edge-effects
                float lambda, float K, float alpha, float beta,
                size_t query_len, uint64_mt db_len, size_t no_sequences );

    void        SetDeltaLength( size_t value ) { deltaLength_ = value; }
    void        SetSearchSpace( uint64_mt value ) { searchSpace_ = value; }

// protected:
//     const AttributableScores*   GetCorresScores() const             { return corres_scores; }
//     AttributableScores*         GetCorresScores()                   { return corres_scores; }
//     const AttributableScores*   GetScaledScores() const             { return scaled_scores; }
//     AttributableScores*         GetScaledScores()                   { return scaled_scores; }

private:
    virtual const TScore** GetScores/*GetImage*/() const = 0;
//     const TMask**       GetMask() const                     { return ( const TMask** )mask; }

    void                SetQuerySize( int len ) { queryLen_ = len; }
    void                SetSubjectSize( int len ) { subjectLen_ = len; }

private:
    const TType         type_;//scoring method and policy
    const TStats        stats_;//approach by which statistical parameters are obtained
//     const TMask         commonmasking;              //common masking approach
    const char*         name_;

    bool                usemodimage_;//modular scores are in use

//     TScaling            auto_scaling;               //whether to use auto scaling to increase precision of scores
protected:
    const int           scaled_by_factor_;//factor by which scores are scaled to increase precision

private:
    float               score_multiplier_;//multiplier of scores associated with used in scaling of scores

    //     double              information_content;        //information content threshold
//     double              ainfoupper2nd;              //upper bound of information content threshold used in 2nd-pass computations
//     double              ainfonumerator2nd;          //numerator of expression to compute 2nd-pass information content threshold
//     double              ainfoscale2nd;              //logarithmic scale to compute 2nd-pass inf. content threshold
//     double              ainfonumeratoralt;          //numerator of alternative expression to compute inf. content threshold
//     double              ainfoscalealt;              //logarithmic scale to alternatively compute inf. content threshold
//     bool                positcorrects;              //whether corrections computed positionally by entropies are in use

    bool                scadjrelent_;//relative entropy score adjustment
    bool                scadjcorrsc_;//adjustment by correlation of scores 

    TScore              min_score_;//minimum score value; just for information
    TScore              max_score_;//maximum score value; just for information

    float               expscore_;//expected score per column pair

    float               referenceLambda_;//reference lambda parameter
    float               referenceH_;//reference parameter H
    float               referenceK_;//reference parameter K
    float               experimentalGappedLambda_;//experimental gapped lambda
    float               experimentalGappedK_;//experimental gapped K
    float               derivedGappedLambda_;//gapped lambda derived for the score matrix under consideration
    float               derivedGappedK_;//gapped K derived for the score matrix under consideration

    float               lambda_;//computed Lambda
    float               entropy_;//computed relative entropy
    float               parameterK_;//computed parameter K

    int                 queryLen_;//length of query profile
    int                 subjectLen_;//length of subject profile

    bool                supportoptfreq_;//optimization of target frequencies

    const Configuration* configuration_;//parameter configuration

    size_t          deltaLength_;//correction for length
    uint64_mt       searchSpace_;//search space
    //map of length correction values already calculated:
    std::map<size_t,std::pair<size_t,uint64_mt>> map_lencrtd_;
};

// INLINES ...
//
// // -------------------------------------------------------------------------
// // PostScalingProc: performs post scaling procedure
// // -------------------------------------------------------------------------
// inline
// void AbstractScoreMatrix::PostScalingProc(
//     const LogOddsMatrix&    querylogo,
//     const LogOddsMatrix&    sbjctlogo,
//     GapScheme&              querygaps,
//     GapScheme&              sbjctgaps,
//     bool    autoc,
//     int     acwindow,
//     bool    firstpass )
// {
// #ifdef __DEBUG__
//     if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
//         throw myruntime_error(
//             mystring( "AbstractScoreMatrix: Memory access error." ));
// #endif
// 
//     if( GetAutoScaling())
//         GetScaledScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
//     else
//         if( autoc )
//             GetCorresScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
// }

// -------------------------------------------------------------------------
// GetConfiguration: get configuration of statistical parameters
//
template<typename TScore>
inline
const Configuration& AbstractScoreMatrix<TScore>::GetConfiguration( TConfigType ps ) const
{
    if( !configuration_ )
        throw MYRUNTIME_ERROR( "AbstractScoreMatrix::GetConfiguration: Null configuration." );

    if( NoCTypes <= ps )
        throw MYRUNTIME_ERROR( "AbstractScoreMatrix::GetConfiguration: Invalid configuration index." );

    return configuration_[ps];
}


#endif//__AbstractScoreMatrix_h__
