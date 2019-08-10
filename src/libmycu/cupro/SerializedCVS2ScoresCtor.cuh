/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedCVS2ScoresCtor_h__
#define __SerializedCVS2ScoresCtor_h__

#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cucom/myassert.h"
#include "SerializedCVS2Scores.cuh"

// -------------------------------------------------------------------------
// class SerializedCVS2ScoresCtor
//
// Implementation of construction of serialized CVS2Scores for parallel 
// processing;
// serialization is transforming a number of linear score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedCVS2ScoresCtor: public SerializedCVS2Scores<TScore>
{
    using SerializedCVS2Scores<TScore>::STEP;
    using SerializedCVS2Scores<TScore>::h_scores_;
    using SerializedCVS2Scores<TScore>::szalloc_;
    using SerializedCVS2Scores<TScore>::nenos_;
    using SerializedCVS2Scores<TScore>::ntbls_;

public:
    SerializedCVS2ScoresCtor( const CVS2Scores& cvs2s );
    virtual __host__ __device__ ~SerializedCVS2ScoresCtor();

protected:
    void NewScores( int sztotal );
    __host__ __device__ void DestroyScores();
};

// -------------------------------------------------------------------------
// INLINES
//
// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
__host__ __device__
inline
SerializedCVS2ScoresCtor<TScore>::~SerializedCVS2ScoresCtor()
{
    CUMSG( "SerializedCVS2ScoresCtor::~SerializedCVS2ScoresCtor", 3 );
    DestroyScores();
}

// -------------------------------------------------------------------------
// constructor: serialize a number of linear score tables
//
template <typename TScore>
inline
SerializedCVS2ScoresCtor<TScore>::SerializedCVS2ScoresCtor( 
    const CVS2Scores& cvs2s )
:   SerializedCVS2Scores<TScore>()
{
    CUMSG( "SerializedCVS2ScoresCtor::SerializedCVS2ScoresCtor", 3 );
    const mystring preamb = "SerializedCVS2ScoresCtor::SerializedCVS2ScoresCtor: ";
#ifdef __DEBUG__
    if( cvs2s.GetNoTables() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid argument.");
#endif
    nenos_ = cvs2s.GetLevels().GetSize();
    ntbls_ = ( nenos_ * ( nenos_ + 1 )) >> 1;

#ifdef __DEBUG__
    if( ntbls_ != cvs2s.GetNoTables() || STEP != cvs2s.GetStep())
        throw MYRUNTIME_ERROR( preamb + "Inconsistent argument.");
#endif

    //indices, shift, and total number of entries over all tables:
    int n, i, j, totentrs;
#ifdef __DEBUG__
    const extspsl::Pslvector* keys = NULL;
#endif
    const extspsl::Pslvector* scos = NULL;

    for( i = 0, totentrs = 0; i < ntbls_; i++ ) {
        scos = cvs2s.GetScores(i);
#ifdef __DEBUG__
        keys = cvs2s.GetKeys(i);
        if( !keys || !scos || keys->GetSize() != scos->GetSize())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent scores.");
        if( scos->GetSize() < 1 )
            throw MYRUNTIME_ERROR( preamb + "No scores.");
        if( cvs2s.GetShift(i) < 0 )
            throw MYRUNTIME_ERROR( preamb + "Negative shift.");
        if( cvs2s.GetScale(i) <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid scale.");
#endif
        totentrs += scos->GetSize();
    }

    //size calculated for scores: levels, cardinality, shift, and scale for each 
    // table, and tables themselves;
    //each table represents a pair of key and value linear arrays
    int szscores = ( nenos_ + ntbls_*3 + totentrs*2 ) * sizeof(TScore);

    NewScores( szscores );

    n = 0;

    //write enos levels
    for( i = 0; i < nenos_; i++ )
        *(h_scores_ + n++) = (TScore)cvs2s.GetLevels().GetValueAt(i);

    //write cardinality for each score table
    for( i = 0; i < ntbls_; i++ )
        *(h_scores_ + n++) = (TScore)cvs2s.GetScores(i)->GetSize();

    //write shift for each score table
    for( i = 0; i < ntbls_; i++ )
        *(h_scores_ + n++) = (TScore)cvs2s.GetShift(i);

    //write scale for each score table
    for( i = 0; i < ntbls_; i++ )
        *(h_scores_ + n++) = (TScore)cvs2s.GetScale(i);

    //write all score tables
    for( i = 0; i < ntbls_; i++ ) {
        for( j = 0; j < cvs2s.GetKeys(i)->GetSize(); j++ )
            *(h_scores_ + n++) = (TScore)cvs2s.GetKeys(i)->GetValueAt(j);
        for( j = 0; j < cvs2s.GetScores(i)->GetSize(); j++ )
            *(h_scores_ + n++) = (TScore)cvs2s.GetScores(i)->GetValueAt(j);
    }
}

// -------------------------------------------------------------------------
// DestroyScores: destroy scores
// 
template <typename TScore>
__host__ __device__
inline
void SerializedCVS2ScoresCtor<TScore>::DestroyScores()
{
    CUMSG( "SerializedCVS2ScoresCtor::DestroyScores", 3 );
    if( h_scores_ ) {
        free( h_scores_ );
        h_scores_ = NULL;
    }
    szalloc_ = 0;
}

// -------------------------------------------------------------------------
// NewScores: allocate memory for scores
// 
template <typename TScore>
inline
void SerializedCVS2ScoresCtor<TScore>::NewScores( int sztotal )
{
    CUMSG( "SerializedCVS2ScoresCtor::NewScores", 3 );
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ )
        throw MYRUNTIME_ERROR("SerializedCVS2ScoresCtor::NewScores: Not enough memory.");
    szalloc_ = sztotal;
}

#endif//__SerializedCVS2ScoresCtor_h__
