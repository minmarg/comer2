/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScoresSM_h__
#define __SerializedScoresSM_h__

#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"
#include "SerializedScores.cuh"

// -------------------------------------------------------------------------
// class SerializedScoresSM
//
// Implementation of serialized VirtScores for parallel processing 
// (shared/general memory version);
// serialization is transforming N-dimensional score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedScoresSM: public SerializedScores<TScore>
{
public:
    __host__ __device__ SerializedScoresSM(
        TScore* scores, int szalloc,
        char nplvs, char nenos, int card,
        int ntbls, int nelems )
    :   SerializedScores<TScore>(szalloc,nplvs,nenos,card,ntbls,nelems),
        h_scores_(scores)
    {};

    __host__ __device__ SerializedScoresSM();
    virtual __host__ __device__ ~SerializedScoresSM();

    virtual __host__ __device__ const TScore* GetScores() const { return h_scores_; }

    __host__ __device__ TScore GetScore( int row, int col, 
        float fstprb, float secprb, float fstens, float secens ) const
    {
        return SerializedScores<TScore>::GetScore( row, col, fstprb, secprb, fstens, secens );
    }
    __host__ __device__ TScore GetScoreP1( int row, int col, 
        float dum1, float dum2, float fstens, float secens ) const
    {
        return SerializedScores<TScore>::GetScoreP1( row, col, dum1, dum2, fstens, secens );
    }
    __host__ __device__ TScore GetScoreP1E1( int row, int col, 
        float dum1, float dum2, float dum3, float dum4 ) const
    {
        return SerializedScores<TScore>::GetScoreP1E1( row, col, dum1, dum2, dum3, dum4 );
    }

    static __host__ __device__ TScore GetScoreP1E1( const TScore* scores, int row, int col );

protected:
    virtual __host__ __device__ TScore AccessScore( int ndx ) const { return h_scores_[ndx]; }

protected:
    //{{host/device data
    TScore* h_scores_;//multi-dimensional scores
    //}}
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore>
__host__ __device__
inline
SerializedScoresSM<TScore>::SerializedScoresSM()
:   SerializedScores<TScore>(),
    h_scores_( NULL )
{
    CUMSG( "SerializedScoresSM::SerializedScoresSM", 3 );
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
__host__ __device__
inline
SerializedScoresSM<TScore>::~SerializedScoresSM()
{
    CUMSG( "SerializedScoresSM::~SerializedScoresSM", 3 );
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
TScore SerializedScoresSM<TScore>::GetScoreP1E1( 
    const TScore* scores, int row, int col )
{
    CUMSG( "SerializedScoresSM::GetScoreP1E1 [static]", 5 );
    //MYASSERT( scores, "SerializedScoresSM::GetScoreP1E1: Memory access error.");
    //MYASSERT( 0 <= row && 0 <= col, 
    //    "SerializedScoresSM::GetScoreP1E1: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    //note: 1 cell for a prob level plus 1 cell for a level of ENOs gives 2
    col += 2;//global index in scores now

    return scores[col];
}

#endif//__SerializedScoresSM_h__
