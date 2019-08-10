/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScores_h__
#define __SerializedScores_h__

#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"

// -------------------------------------------------------------------------
// class SerializedScores
//
// Implementation of serialized VirtScores for parallel processing;
// serialization is transforming N-dimensional score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedScores
{
    enum {
        NAVAL = -9999,
    };
public:
    SerializedScores( 
        const extspsl::Pslvector*** scores, 
        int npairplvs, int ntbls, int card, 
        const extspsl::Pslvector& prblvs, const extspsl::Pslvector& levels
    );
    SerializedScores();
    ~SerializedScores();

    int GetSizeAlloc() const { return szalloc_; }
    const TScore* GetScores() const { return h_scores_; }

    void GetScore( TScore*, int row, int col, 
        float fstprb, float secprb, float fstens, float secens ) const;
    void GetScoreP1( TScore*, int row, int col, 
        float, float, float fstens, float secens ) const;
    void GetScoreP1E1( TScore*, int row, int col, 
        float, float, float, float ) const;

    bool ScoreNA( TScore value ) const { return value <= TScore(NAVAL); }

protected:
    void NewScores( int sztotal );
    void DestroyScores();

protected:
    //{{host data
    TScore* h_scores_;//multi-dimensional scores
    //}}
    int szalloc_;//size (bytes) allocated for scores
    char nplvs_;//number of probability levels
    char nenos_;//number of levels for eff. number of observations
    int card_;//cardinality determinator, # rows of a square score table
    int ntbls_;//number of tables per a pair of probability level values
    int nelems_;//number of entries in one table
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore>
inline
SerializedScores<TScore>::SerializedScores()
:   h_scores_( NULL ),
    szalloc_( 0 ),
    nplvs_( 0 ),
    nenos_( 0 ),
    card_( 0 ),
    ntbls_( 0 ),
    nelems_( 0 )
{
    CUMSG( "SerializedScores::SerializedScores", 3 );
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
inline
SerializedScores<TScore>::~SerializedScores()
{
    CUMSG( "SerializedScores::~SerializedScores", 3 );
    DestroyScores();
}

// -------------------------------------------------------------------------
// constructor: serialize multi-dimensional scores;
// scores, scores to be serialized;
// npairplvs, number of different pair values for probability levels;
// ntbls, number of tables at each pair of values of probability levels;
// card, number of rows of a square score table for given probability and 
//  enos levels;
// prblvs, probability level values;
// levels, enos level values;
//
template <typename TScore>
inline
SerializedScores<TScore>::SerializedScores( 
        const extspsl::Pslvector*** scores, 
        int nopairplvs, int notbls, int card, 
        const extspsl::Pslvector& prblvs, const extspsl::Pslvector& levels )
:   h_scores_( NULL ),
    szalloc_( 0 ),
    nplvs_( 0 ),
    nenos_( 0 ),
    card_( 0 ),
    ntbls_( 0 ),
    nelems_( 0 )
{
    CUMSG( "SerializedScores::SerializedScores", 3 );
#ifdef __DEBUG__
    if( !scores || nopairplvs < 0 || notbls < 1 || card < 1 )
        throw MYRUNTIME_ERROR("SerializedScores::SerializedScores: Invalid arguments.");
#endif
    nplvs_ = SLC_MAX( 1, prblvs.GetSize());
    nenos_ = levels.GetSize();
    card_ = card;

    int npairplvscalc = ( nplvs_ * ( nplvs_ + 1 )) >> 1;
    ntbls_ = ( nenos_ * ( nenos_ + 1 )) >> 1;
    nelems_ = ( card_ * ( card_+1 )) >> 1;//number of entries in one table

#ifdef __DEBUG__
    if( npairplvscalc != nopairplvs || ntbls_ != notbls )
        throw MYRUNTIME_ERROR("SerializedScores::SerializedScores: Inconsistent argument values.");
#endif

    //size calculated for scores:
    int szscores = ( nplvs_ + nenos_ + npairplvscalc*ntbls_*nelems_ ) * sizeof(TScore);
    int n, i, j, k;

    NewScores( szscores );

    n = 0;

    //write probability levels
    if( 0 < prblvs.GetSize())
        for( i = 0; i < nplvs_; i++ )
            *(h_scores_ + n++) = (TScore)prblvs.GetValueAt(i);
    else
        h_scores_[n++] = TScore(0);

    //write enos levels
    for( i = 0; i < nenos_; i++ )
        *(h_scores_ + n++) = (TScore)levels.GetValueAt(i);

    //write all score tables
    for( i = 0; i < npairplvscalc; i++ )
        for( j = 0; j < ntbls_; j++ )
            for( k = 0; k < nelems_; k++ )
                *(h_scores_ + n++) = (TScore)scores[i][j]->GetValueAt(k);
}

// -------------------------------------------------------------------------
// DestroyScores: destroy scores
// 
template <typename TScore>
inline
void SerializedScores<TScore>::DestroyScores()
{
    CUMSG( "SerializedScores::DestroyScores", 3 );
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
void SerializedScores<TScore>::NewScores( int sztotal )
{
    CUMSG( "SerializedScores::NewScores", 3 );
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ )
        throw MYRUNTIME_ERROR("SerializedScores::NewScores: Not enough memory.");
    szalloc_ = sztotal;
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by row, col; table is 
//  selected according to probability levels `fstprb' and `secprb' and eff. 
//  no. observations given by `fstens' and `secens' 
// 
template <typename TScore>
inline
void SerializedScores<TScore>::GetScore( TScore* score, int row, int col, 
    float fstprb, float secprb, float fstens, float secens ) const
{
    CUMSG( "SerializedScores::GetScore", 5 );
    int e, ee, pi, ppi;//level indices

    MYASSERT( h_scores_ && score, "SerializedScores::GetScore: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col && row < card_ && col < card_, 
        "SerializedScores::GetScore: Invalid indices.");

    if( fstprb < (float)h_scores_[0] || secprb < (float)h_scores_[0]) {
        *score = TScore(0);
        return;
    }

    CNDSWAP( float, secprb<fstprb, secprb, fstprb );

    if( fstens < (float)h_scores_[nplvs_+0] || secens < (float)h_scores_[nplvs_+0]) {
        *score = TScore(0);
        return;
    }

    CNDSWAP( float, secens<fstens, secens, fstens );

    //get the corresponding score table index;
    //the first index is determined by the pair of probability levels
    for( pi = 1; pi < nplvs_; pi++ ) {
        if( fstprb < h_scores_[pi] )
            break;
    }
    for( ppi = pi; ppi < nplvs_; ppi++ ) {
        if( secprb < h_scores_[ppi] )
            break;
    }
    pi--;
    ppi--;

    MYASSERT( pi < nplvs_ && ppi < nplvs_ && pi <= ppi, 
        "SerializedScores::GetScore: Wrong probability level values.");

    //pair index calculated as sum of arithm. series
    ppi = (( pi * ( TIMES2(nplvs_)-pi+1 ))>>1 ) + ppi-pi;

    //the second index is determined by the pair of enos levels
    for( e = nplvs_+1; e < nplvs_+nenos_; e++ ) {
        if( fstens < h_scores_[e] )
            break;
    }
    for( ee = e; ee < nplvs_+nenos_; ee++ ) {
        if( secens < h_scores_[ee] )
            break;
    }
    e -= nplvs_ + 1;
    ee -= nplvs_ + 1;

//     MYASSERT( e < nplvs_+nenos_ && ee < nplvs_+nenos_ && e <= ee, 
//         "SerializedScores::GetScore: Wrong levels for effective number of observations.");

    //pair index calculated as sum of arithm. series
    ee = (( e * ( TIMES2(nenos_)-e+1 ))>>1 ) + ee-e;

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    col += nplvs_+nenos_ + ( ppi*ntbls_ + ee ) * nelems_;//global index in h_scores_ now

    MYASSERT( col*(int)sizeof(TScore) < szalloc_,
        "SerializedScores::GetScore: Index out of range.");

    *score = h_scores_[col];
}

// -------------------------------------------------------------------------
// GetScoreP1: get score at table position identified by row, col; table is 
//  selected according to eff. no. observations given by `fstens' and 
// `secens';
// NOTE: a single probability level is assumed to be present, covering the 
//  full range of values 
// 
template <typename TScore>
inline
void SerializedScores<TScore>::GetScoreP1( TScore* score, int row, int col, 
    float, float, float fstens, float secens ) const
{
    CUMSG( "SerializedScores::GetScoreP1", 5 );
    int e, ee;//level indices

    MYASSERT( h_scores_ && score, "SerializedScores::GetScoreP1: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col && row < card_ && col < card_, 
        "SerializedScores::GetScoreP1: Invalid indices.");

    //the first address is for a probability level by definition
    if( fstens < (float)h_scores_[1] || secens < (float)h_scores_[1]) {
        *score = TScore(0);
        return;
    }

    CNDSWAP( float, secens<fstens, secens, fstens );

    //the index is determined by the pair of enos levels
    //note: 1 cell for a prob level plus 1 gives 2
    //TODO: USE SHARED MEM for getting e and ee!
    for( e = 2; e < 1+nenos_; e++ ) {
        if( fstens < h_scores_[e] )
            break;
    }
    for( ee = e; ee < 1+nenos_; ee++ ) {
        if( secens < h_scores_[ee] )
            break;
    }
    e -= 2;
    ee -= 2;

//     MYASSERT( e < nplvs_+nenos_ && ee < nplvs_+nenos_ && e <= ee, 
//         "SerializedScores::GetScoreP1: Wrong levels for effective number of observations.");

    //pair index calculated as sum of arithm. series
    ee = (( e * ( TIMES2(nenos_)-e+1 ))>>1 ) + ee-e;

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    col += 1+nenos_ + ee * nelems_;//global index in h_scores_ now

    MYASSERT( col*(int)sizeof(TScore) < szalloc_,
        "SerializedScores::GetScoreP1: Index out of range.");

    *score = h_scores_[col];
}

// -------------------------------------------------------------------------
// GetScoreP1E1: get score at table position identified by row, col; 
// NOTE: a single probability level and a single level for ENOS is 
// assumed to be present, covering the full range of values; 
// hence, only one (the first) score table is referenced;
// 
template <typename TScore>
inline
void SerializedScores<TScore>::GetScoreP1E1( 
    TScore* score, int row, int col, 
    float, float, float, float ) const
{
    CUMSG( "SerializedScores::GetScoreP1E1", 5 );
    MYASSERT( h_scores_ && score, "SerializedScores::GetScoreP1E1: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col && row < card_ && col < card_, 
        "SerializedScores::GetScoreP1E1: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    //note: 1 cell for a prob level plus 1 cell for a level of ENOs gives 2
    col += 2;//global index in h_scores_ now

    MYASSERT( col*(int)sizeof(TScore) < szalloc_,
        "SerializedScores::GetScoreP1E1: Index out of range.");

    *score = h_scores_[col];
}

#endif//__SerializedScores_h__
