/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScoresTM_h__
#define __SerializedScoresTM_h__

#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"
#include "SerializedScores.cuh"

// -------------------------------------------------------------------------
// class SerializedScoresTM
//
// Implementation of serialized VirtScores for parallel processing (CUDA 
// Texture memory version);
// serialization is transforming N-dimensional score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedScoresTM: public SerializedScores<TScore>
{
public:
    __host__ __device__ SerializedScoresTM(
        cudaTextureObject_t texObj, int szalloc,
        char nplvs, char nenos, int card,
        int ntbls, int nelems )
    :   SerializedScores<TScore>(szalloc,nplvs,nenos,card,ntbls,nelems),
        texObj_(texObj)
    {};

    __host__ __device__ SerializedScoresTM();
    virtual __host__ __device__ ~SerializedScoresTM();

    __host__ __device__ const TScore* GetScores() const { return (TScore*)texObj_; }

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

    static __device__ TScore GetScoreP1E2( 
        cudaTextureObject_t texo, float eth1, float eth2, int nelems,
        int row, int col, float fstens, float secens );

protected:
    virtual __host__ __device__ TScore AccessScore( int ndx ) const { return tex1Dfetch<TScore>(texObj_,ndx); }
    static __host__ __device__ TScore AccessScore( cudaTextureObject_t texo, int ndx ) { 
        return tex1Dfetch<TScore>(texo,ndx); 
    }

protected:
    //{{device data
    cudaTextureObject_t texObj_;//identifier of texture object of multi-dimensional scores
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
SerializedScoresTM<TScore>::SerializedScoresTM()
:   SerializedScores<TScore>(),
    texObj_( 0 )
{
    CUMSG( "SerializedScoresTM::SerializedScoresTM", 3 );
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
__host__ __device__
inline
SerializedScoresTM<TScore>::~SerializedScoresTM()
{
    CUMSG( "SerializedScoresTM::~SerializedScoresTM", 3 );
}

// -------------------------------------------------------------------------
// GetScoreP1E2: get score at table position identified by row, col; 
// NOTE: a SINGLE probability level and TWO levels for ENOS, 
// specified by the thresholds eth1 and eth2, are assumed, covering the full 
// range of values; 
// eff. no. observations `fstens' and `secens' pick a table;
// nelems, number of entries in one table;
// 
template <typename TScore>
__device__
__forceinline__
TScore SerializedScoresTM<TScore>::GetScoreP1E2( 
    cudaTextureObject_t texo, float eth1, float eth2, int nelems,
    int row, int col, float fstens, float secens)
{
    CUMSG( "SerializedScoresTM::GetScoreP1E2", 5 );
    //MYASSERT( texo, "SerializedScoresTM::GetScoreP1E2: Memory access error.");

    //if( row < 0 || col < 0 )
    //   return TScore(0);

    //the first address is for a probability level by definition
    //if( fstens < (float)AccessScore(1) || secens < (float)AccessScore(1))
    if( fstens < eth1 || secens < eth1 )
        return TScore(0);

    CNDSWAP( float, secens<fstens, secens, fstens );

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + (( row * ( row-1 )) >> 1 );//index of the element now

    //calculate table index
    //reuse registers: row is now pair index (calculated at once)
    if( eth2 <= fstens )
        row = 2;//[last table]
    else if( eth2 <= secens )
        row = 1;//[intermediate]
    else
        row = 0;//[first table]

    //col += 1+nenos_ + ee * nelems_;//global index in h_scores_ now
    col += 3 + row * nelems;//global index in h_scores_ now

    //MYASSERT( col*(int)sizeof(TScore) < szalloc_,
    //    "SerializedScoresTM::GetScoreP1E2: Index out of range.");

    return AccessScore(texo, col);
//     TScore sc = AccessScore(texo, col);
//     if( SerializedScores<TScore>::ScoreNA(sc))
//         return TScore(0);
//     return sc;
}

#endif//__SerializedScoresTM_h__
