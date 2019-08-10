/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScores_h__
#define __SerializedScores_h__

#include <stdlib.h>

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
        //NAVAL = -9999,
        NAVAL = -999,//reduce original constant to avoid fp representation issues
    };
public:
    __host__ __device__ SerializedScores(
        int szalloc,
        char nplvs, char nenos, int card,
        int ntbls, int nelems )
    :   szalloc_(szalloc), 
        nplvs_(nplvs), nenos_(nenos), card_(card),
        ntbls_(ntbls), nelems_(nelems) 
    {};

    __host__ __device__ SerializedScores();
    virtual __host__ __device__ ~SerializedScores();

    virtual __host__ __device__ const TScore* GetScores() const = 0;
    __host__ __device__ int GetSizeAlloc() const { return szalloc_; }
    __host__ __device__ char GetNProbLevs() const { return nplvs_; }
    __host__ __device__ char GetNENOLevs() const { return nenos_; }
    __host__ __device__ int GetCardinality() const { return card_; }
    __host__ __device__ int GetNTables() const { return ntbls_; }
    __host__ __device__ int GetNEntries() const { return nelems_; }

    __host__ __device__ void GetEth12( float* eth1, float* eth2 ) const;
    __host__ __device__ TScore GetScore( int row, int col, 
        float fstprb, float secprb, float fstens, float secens ) const;
    __host__ __device__ TScore GetScoreP1( int row, int col, 
        float, float, float fstens, float secens ) const;
    __host__ __device__ TScore GetScoreP1E1( int row, int col, 
        float, float, float, float ) const;

    static __host__ __device__ bool ScoreNA( TScore value ) { return value <= TScore(NAVAL); }

protected:
    virtual __host__ __device__ TScore AccessScore( int ndx ) const = 0;

protected:
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
__host__ __device__
inline
SerializedScores<TScore>::SerializedScores()
:   szalloc_( 0 ),
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
__host__ __device__
inline
SerializedScores<TScore>::~SerializedScores()
{
    CUMSG( "SerializedScores::~SerializedScores", 3 );
}

// -------------------------------------------------------------------------
// GetEth12: read the first two ENO thresholds in the buffer of scores
// 
template <typename TScore>
__host__ __device__ 
inline
void SerializedScores<TScore>::GetEth12( float* eth1, float* eth2 ) const
{
    CUMSG( "SerializedScores::GetEth12", 5 );
    MYASSERT( GetScores(), "SerializedScores::GetEth12: Memory access error.");
    if( eth1 )
        *eth1 = (float)AccessScore(nplvs_);
    if( eth2 ) {
        *eth2 = 0.0f;
        if( 1 < nenos_ )
            *eth2 = (float)AccessScore(nplvs_+1);
    }
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by row, col; table is 
//  selected according to probability levels `fstprb' and `secprb' and eff. 
//  no. observations given by `fstens' and `secens' 
// 
template <typename TScore>
__host__ __device__
inline
TScore SerializedScores<TScore>::GetScore( int row, int col, 
    float fstprb, float secprb, float fstens, float secens ) const
{
    CUMSG( "SerializedScores::GetScore", 5 );
    int e, ee, pi, ppi;//level indices

    //make compiler remove variable card_; it's used only in assert statements
    MYASSERT( GetScores(), "SerializedScores::GetScore: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedScores::GetScore: Invalid indices.");

    if( fstprb < (float)AccessScore(0) || secprb < (float)AccessScore(0))
        return TScore(0);

    CNDSWAP( float, secprb<fstprb, secprb, fstprb );

    if( fstens < (float)AccessScore(nplvs_+0) || secens < (float)AccessScore(nplvs_+0))
        return TScore(0);

    CNDSWAP( float, secens<fstens, secens, fstens );

    //get the corresponding score table index;
    //the first index is determined by the pair of probability levels
    for( pi = 1; pi < nplvs_; pi++ ) {
        if( fstprb < AccessScore(pi) )
            break;
    }
    for( ppi = pi; ppi < nplvs_; ppi++ ) {
        if( secprb < AccessScore(ppi) )
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
        if( fstens < AccessScore(e) )
            break;
    }
    for( ee = e; ee < nplvs_+nenos_; ee++ ) {
        if( secens < AccessScore(ee) )
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

    return AccessScore(col);
}

// -------------------------------------------------------------------------
// GetScoreP1: get score at table position identified by row, col; table is 
//  selected according to eff. no. observations given by `fstens' and 
// `secens';
// NOTE: a single probability level is assumed to be present, covering the 
//  full range of values 
// 
template <typename TScore>
__host__ __device__
inline
TScore SerializedScores<TScore>::GetScoreP1( int row, int col, 
    float, float, float fstens, float secens ) const
{
    CUMSG( "SerializedScores::GetScoreP1", 5 );
    int e, ee;//level indices

    MYASSERT( GetScores(), "SerializedScores::GetScoreP1: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedScores::GetScoreP1: Invalid indices.");

    //the first address is for a probability level by definition
    if( fstens < (float)AccessScore(1) || secens < (float)AccessScore(1))
        return TScore(0);

    CNDSWAP( float, secens<fstens, secens, fstens );

    //the index is determined by the pair of enos levels
    //note: 1 cell for a prob level plus 1 gives 2
    //TODO: USE SHARED MEM for getting e and ee!
    for( e = 2; e < 1+nenos_; e++ ) {
        if( fstens < AccessScore(e) )
            break;
    }
    for( ee = e; ee < 1+nenos_; ee++ ) {
        if( secens < AccessScore(ee) )
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

    return AccessScore(col);
}

// -------------------------------------------------------------------------
// GetScoreP1E1: get score at table position identified by row, col; 
// NOTE: a single probability level and a single level for ENOS is 
// assumed to be present, covering the full range of values; 
// hence, only one (the first) score table is referenced;
// 
template <typename TScore>
__host__ __device__
inline
TScore SerializedScores<TScore>::GetScoreP1E1( 
    int row, int col, 
    float, float, float, float ) const
{
    CUMSG( "SerializedScores::GetScoreP1E1", 5 );
    MYASSERT( GetScores(), "SerializedScores::GetScoreP1E1: Memory access error.");
    MYASSERT( 0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedScores::GetScoreP1E1: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    //note: 1 cell for a prob level plus 1 cell for a level of ENOs gives 2
    col += 2;//global index in h_scores_ now

    MYASSERT( col*(int)sizeof(TScore) < szalloc_,
        "SerializedScores::GetScoreP1E1: Index out of range.");

    return AccessScore(col);
}

#endif//__SerializedScores_h__
