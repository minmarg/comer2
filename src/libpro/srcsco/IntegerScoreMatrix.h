/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __IntegerScoreMatrix_h__
#define __IntegerScoreMatrix_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>

#include "liblib/msg.h"
#include "libpro/srcpro/Configuration.h"
#include "AbstractScoreMatrix.h"
#include "ScoresAttr.h"

////////////////////////////////////////////////////////////////////////////
// CLASS IntegerScoreMatrix
// Integer implementations of Score Matrices
//
template<typename TScore>
class IntegerScoreMatrix: public AbstractScoreMatrix<TScore>
{
public:
    IntegerScoreMatrix( typename AbstractScoreMatrix<TScore>::TType, 
                        Configuration config[NoCTypes], 
                        typename AbstractScoreMatrix<TScore>::TStats,
                        int precscale = 1 );
    virtual ~IntegerScoreMatrix();

    virtual TScore      GetScore( int m, int n ) const;
    virtual TScore      GetModScore( int m, int n ) const;

    virtual int         GetScoresFactor() const { return AbstractScoreMatrix<TScore>::GetScoresFactor(); }
    virtual int         GetQuerySize() const { return AbstractScoreMatrix<TScore>::GetQuerySize(); }
    virtual int         GetSubjectSize() const { return AbstractScoreMatrix<TScore>::GetSubjectSize(); }

    virtual bool        ScanForHSPs( float minhspscore, int hsplen, int nohsps, int maxdist, 
                                     int* = NULL, int* = NULL );

    float               GetRefLambda() const            { return AbstractScoreMatrix<TScore>::GetRefLambda(); }
    float               GetRefH() const                 { return AbstractScoreMatrix<TScore>::GetRefH(); }
    float               GetRefK() const                 { return AbstractScoreMatrix<TScore>::GetRefK(); }

    float               GetLambda() const               { return AbstractScoreMatrix<TScore>::GetLambda(); }
    float               GetEntropy() const              { return AbstractScoreMatrix<TScore>::GetEntropy(); }
    float               GetK() const                    { return AbstractScoreMatrix<TScore>::GetK(); }
    float               GetExpectedScore() const        { return AbstractScoreMatrix<TScore>::GetExpectedScore(); }

    void                SetLambda( float newlambda )    { AbstractScoreMatrix<TScore>::SetLambda( newlambda ); }
    void                SetEntropy( float H )           { AbstractScoreMatrix<TScore>::SetEntropy( H ); }
    void                SetK( float K )                 { AbstractScoreMatrix<TScore>::SetK( K ); }
    void                SetExpectedScore( float E )     { AbstractScoreMatrix<TScore>::SetExpectedScore( E ); }

    void                SetMinScore( TScore value )     { AbstractScoreMatrix<TScore>::SetMinScore( value ); }
    void                SetMaxScore( TScore value )     { AbstractScoreMatrix<TScore>::SetMaxScore( value ); }
    void                SetMultiplier( float value )    { AbstractScoreMatrix<TScore>::SetMultiplier( value ); }

    
    virtual void        ComputeScoreMatrix( bool final = false ) = 0;

//     void                PostScalingProc(
//                                 const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
//                                 GapScheme& qg, GapScheme& sg,
//                                 bool autoc, int acwindow, bool firstpass );

    float               GetFinalScore( TScore value ) const;

    virtual void        ScaleScoreMatrix();
    void                ScaleScoreMatrixHelper();
    float               IntegerScaleMatrix();
    void                MultiplyScoresBy( float fact );
    virtual void        ComputeStatisticalParameters( bool computelambda = true );

//     virtual void        ComputeScoreProbabilities( ScoresAttr<TScore>* ) = 0;

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const = 0;
    virtual void        PrintScoreMatrix( FILE* ) {;};

protected:
    explicit IntegerScoreMatrix();

    virtual void        Init( int querylen, int sbjctlen );

    void                SetScore( int m, int n, TScore );
    void                SetModScore( int m, int n, TScore );

    bool                GetAllNegatives() const { return allnegatives_; }
    void                SetAllNegatives( bool value );

protected:
    const ScoresAttr<TScore>*   GetAttr() const { return attr_; }
    ScoresAttr<TScore>*         GetAttr() { return attr_; }

    void    DeleteAttr();
    void    NewAttr();

private:
    virtual const TScore** GetScores() const { return (const TScore**)scores_; }

private:
    TScore**            scores_;//scores
    TScore**            modscores_;//modular scores
    bool                allnegatives_;//whether scores are all negative
    ScoresAttr<TScore>* attr_;//attributes of scores
};

// INLINES ...
//
// // -------------------------------------------------------------------------
// // PostScalingProc: performs post scaling procedure
// // -------------------------------------------------------------------------
// inline
// void IntegerScoreMatrix::PostScalingProc(
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
//             mystring( "IntegerScoreMatrix: Memory access error." ));
// #endif
// 
//     if( GetAutoScaling())
//         GetScaledScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
//     else
//         if( autoc )
//             GetCorresScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
// }

// -------------------------------------------------------------------------
// GetScore: get the score at profile positions
//
template<typename TScore>
inline
TScore IntegerScoreMatrix<TScore>::GetScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !scores_ || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw MYRUNTIME_ERROR( "IntegerScoreMatrix::GetScore: Memory access error." );
#endif
    return scores_[m][n];
}

// -------------------------------------------------------------------------
// GetModScore: get the modular score at profile positions
//
template<typename TScore>
inline
TScore IntegerScoreMatrix<TScore>::GetModScore(int m, int n) const
{
#ifdef __DEBUG__
    if( !modscores_ || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw MYRUNTIME_ERROR("IntegerScoreMatrix::GetModScore: Memory access error.");
#endif
    return modscores_[m][n];
}

// -------------------------------------------------------------------------
// GetFinalScore: get final score which is the score divided by the factor 
// used to scale scores to increase precision
//
template<typename TScore>
inline
float IntegerScoreMatrix<TScore>::GetFinalScore( TScore value ) const
{
    float val = (float)value;
    int fct = GetScoresFactor();

    if( val && ( fct <= -1 || fct > 1 ))
        val /= (float)fct;

    return val;
}
// -------------------------------------------------------------------------
// SetScore: set a score at profile positions
//
template<typename TScore>
inline
void IntegerScoreMatrix<TScore>::SetScore( int m, int n, TScore value )
{
#ifdef __DEBUG__
    if( !scores_ || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw MYRUNTIME_ERROR( "IntegerScoreMatrix::SetScore: Memory access error." );
#endif
    scores_[m][n] = value;
}

// -------------------------------------------------------------------------
// SetModScore: set a modular score at specified positions
//
template<typename TScore>
inline
void IntegerScoreMatrix<TScore>::SetModScore( int m, int n, TScore value )
{
#ifdef __DEBUG__
    if( !modscores_ || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw MYRUNTIME_ERROR( "IntegerScoreMatrix::SetModScore: Memory access error." );
#endif
    modscores_[m][n] = value;
}

// -------------------------------------------------------------------------
// SetAllNegatives: set a flag indicating the presence or absence of all 
// negative scores
//
template<typename TScore>
inline
void IntegerScoreMatrix<TScore>::SetAllNegatives( bool value )
{
    if( GetAttr()) GetAttr()->SetAllNegatives( value );
    allnegatives_ = value;
}

// -------------------------------------------------------------------------
// DeleteAttr: destroy the scores attribute object
//
template<typename TScore>
inline
void IntegerScoreMatrix<TScore>::DeleteAttr()
{
    if( attr_ ) {
        delete attr_;
        attr_ = NULL;
    }
}

// -------------------------------------------------------------------------
// NewAttr: create a new scores attribute object
//
template<typename TScore>
inline
void IntegerScoreMatrix<TScore>::NewAttr()
{
    DeleteAttr();

    attr_ = new ScoresAttr<TScore>(
        this,
        &AbstractScoreMatrix<TScore>::ComputeScoreProbabilities,
        GetRefLambda() / GetScoresFactor(),//divide lambda by a factor and...
        GetRefH(),
        GetRefK(),
        GetScores(),
        GetScoresFactor()//provide the factor used to scale scores
    );

    if( !attr_ )
        throw MYRUNTIME_ERROR( "IntegerScoreMatrix::NewAttr: Not enough memory." );

    attr_->Init( GetQuerySize(), GetSubjectSize());
}

#endif//__IntegerScoreMatrix_h__
