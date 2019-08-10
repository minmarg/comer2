/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedCVS2Scores_h__
#define __SerializedCVS2Scores_h__

#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cucom/cucommon.h"

// -------------------------------------------------------------------------
// class SerializedCVS2Scores
//
// Implementation of serialized CVS2Scores for parallel processing;
// serialization is transforming a number of linear score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedCVS2Scores
{
    enum {
        //number of intermediate values +1 between two adjacent integers
        STEP = 2
    };
    enum TSCSInfOrder{
        scsioCARDINALITY, 
        scsioSHIFT, 
        scsioSCALE
    };

public:
    SerializedCVS2Scores( const CVS2Scores& cvs2s );
    SerializedCVS2Scores();
    ~SerializedCVS2Scores();

    void GetScore( TScore*, TScore cvscore, float fstens, float secens ) const;
    void GetScoreE1S1Step2( TScore* score, TScore cvscore, 
        int card, int shft, 
        float, float ) const;

protected:
    void NewScores( int sztotal );
    void DestroyScores();

protected:
    //{{host data
    TScore* h_scores_;//multi-dimensional scores
    //}}
    int szalloc_;//size (bytes) allocated for scores
    char nenos_;//number of levels for eff. number of observations
    int ntbls_;//number of tables per a pair of probability level values
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore>
inline
SerializedCVS2Scores<TScore>::SerializedCVS2Scores()
:   h_scores_( NULL ),
    szalloc_( 0 ),
    nenos_( 0 ),
    ntbls_( 0 )
{
    CUMSG( "SerializedCVS2Scores::SerializedCVS2Scores", 3 );
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
inline
SerializedCVS2Scores<TScore>::~SerializedCVS2Scores()
{
    CUMSG( "SerializedCVS2Scores::~SerializedCVS2Scores", 3 );
    DestroyScores();
}

// -------------------------------------------------------------------------
// constructor: serialize a number of linear score tables
//
template <typename TScore>
inline
SerializedCVS2Scores<TScore>::SerializedCVS2Scores( 
    const CVS2Scores& cvs2s )
:   h_scores_( NULL ),
    szalloc_( 0 ),
    nenos_( 0 ),
    ntbls_( 0 )
{
    CUMSG( "SerializedCVS2Scores::SerializedCVS2Scores", 3 );
    const mystring preamb = "SerializedCVS2Scores::SerializedCVS2Scores: ";
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
    const extspsl::Pslvector* keys = NULL;
    const extspsl::Pslvector* scos = NULL;

    for( i = 0, totentrs = 0; i < ntbls_; i++ ) {
        keys = cvs2s.GetKeys(i);
        scos = cvs2s.GetScores(i);
#ifdef __DEBUG__
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
inline
void SerializedCVS2Scores<TScore>::DestroyScores()
{
    CUMSG( "SerializedCVS2Scores::DestroyScores", 3 );
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
void SerializedCVS2Scores<TScore>::NewScores( int sztotal )
{
    CUMSG( "SerializedCVS2Scores::NewScores", 3 );
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ )
        throw MYRUNTIME_ERROR("SerializedCVS2Scores::NewScores: Not enough memory.");
    szalloc_ = sztotal;
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified/mapped by cvscore; 
// table is selected according to eff. no. observations given by 
// `fstens' and `secens' 
// 
template <typename TScore>
inline
void SerializedCVS2Scores<TScore>::GetScore( 
    TScore* score, TScore cvscore, 
    float fstens, float secens ) const
{
    CUMSG( "SerializedCVS2Scores::GetScore", 5 );
    int e, ee;//level indices

    MYASSERT( h_scores_ && score, "SerializedCVS2Scores::GetScore: Memory access error.");

    if( fstens < (float)h_scores_[0] || secens < (float)h_scores_[0]) {
        *score = TScore(0);
        return;
    }

    CNDSWAP( float, secens<fstens, secens, fstens );

    //get the corresponding score table index;
    //the index is determined by the pair of enos levels
    for( e = 1; e < nenos_; e++ ) {
        if( fstens < h_scores_[e] )
            break;
    }
    for( ee = e; ee < nenos_; ee++ ) {
        if( secens < h_scores_[ee] )
            break;
    }
    e--;
    ee--;

//     MYASSERT( e < nenos_ && ee < nenos_ && e <= ee, 
//         "SerializedCVS2Scores::GetScore: Wrong levels for effective number of observations.");

    //pair index calculated as sum of arithm. series
    ee = (( e * ( TIMES2(nenos_)-e+1 ))>>1 ) + ee-e;

    int card = (int)h_scores_[nenos_+ee];
    int shft = (int)h_scores_[nenos_+ntbls_+ee];
    float scale = (float)h_scores_[nenos_+ntbls_+ntbls_+ee];

    if( scale != 1.0f )
        cvscore *= scale;

    int beg = nenos_+ntbls_*3;//beginning index for scores
    for( int i = 0; i < ee; i++ )
        //there are two series of scores: keys and values
        beg += (int)h_scores_[nenos_+i]*2;

    //if cvscore <= the fist key value
    if( cvscore <= (float)h_scores_[beg]) {
        //return the first value from the scores
        *score = h_scores_[beg+card];
        return;
    }
    //if the last key value <= cvscore
    if( (float)h_scores_[beg+card-1] <= cvscore ) {
        //return the last value from the scores
        *score = h_scores_[beg+card+card-1];
        return;
    }

    //calculate index within a linear score table;
    //the index is written to beg
    if( STEP == 2 )
        beg += shft + (int)rintf(TIMES2(cvscore));
    else if( STEP == 1 )
        beg += shft + (int)rintf(cvscore);
    else
        beg += shft + (int)rintf((float)STEP*cvscore);

    MYASSERT((beg+card)*(int)sizeof(TScore) < szalloc_,
        "SerializedCVS2Scores::GetScore: Index out of range.");

    //return the score found at the index from the scores (+card)
    *score = h_scores_[beg+card];
}

// -------------------------------------------------------------------------
// GetScoreE1S1Step2: get score at table position identified/mapped by 
// cvscore; 
// NOTE: a single level for ENOS is assumed to be present, covering the 
// full range of values; 
// hence, only one (the first) score table is referenced;
// the cardinality and score shift are given by parameters card and shift
// 
template <typename TScore>
inline
void SerializedCVS2Scores<TScore>::GetScoreE1S1Step2( 
    TScore* score, TScore cvscore, 
    int card, int shft, 
    float, float ) const
{
    CUMSG( "SerializedCVS2Scores::GetScoreE1S1Step2", 5 );
    MYASSERT( h_scores_ && score, 
        "SerializedCVS2Scores::GetScoreE1S1Step2: Memory access error.");

    //beginning index for scores: 1 for #levels, 3 for 
    // cardinality, shift, and scale
    //int beg = 1+3;

    //if cvscore <= the fist key value
    if( cvscore <= (float)h_scores_[4/*beg*/]) {
        //return the first value from the scores
        *score = h_scores_[4/*beg*/+card];
        return;
    }
    //if the last key value <= cvscore
    if( (float)h_scores_[4/*beg*/+card-1] <= cvscore ) {
        //return the last value from the scores
        *score = h_scores_[4/*beg*/+card+card-1];
        return;
    }

    //calculate index within a linear score table;
    //the index is overwritten to shft;
    //( STEP == 2 )
    shft += 4/*beg*/ + (int)rintf(TIMES2(cvscore));

    MYASSERT((shft+card)*(int)sizeof(TScore) < szalloc_,
        "SerializedCVS2Scores::GetScoreE1S1Step2: Index out of range.");

    //return the score found at the index from the scores (+card)
    *score = h_scores_[shft+card];
}

#endif//__SerializedCVS2Scores_h__
