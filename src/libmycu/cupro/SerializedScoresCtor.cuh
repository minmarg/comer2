/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedScoresCtor_h__
#define __SerializedScoresCtor_h__

#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"
#include "SerializedScoresSM.cuh"

// -------------------------------------------------------------------------
// class SerializedScoresCtor
//
// Implementation of construction of serialized VirtScores for parallel 
// processing;
// serialization is transforming N-dimensional score tables to a 
// 1-dimensional array;
//
template <typename TScore>
class SerializedScoresCtor: public SerializedScoresSM<TScore>
{
    using SerializedScoresSM<TScore>::h_scores_;
    using SerializedScoresSM<TScore>::szalloc_;
    using SerializedScoresSM<TScore>::nplvs_;
    using SerializedScoresSM<TScore>::nenos_;
    using SerializedScoresSM<TScore>::card_;
    using SerializedScoresSM<TScore>::ntbls_;
    using SerializedScoresSM<TScore>::nelems_;

public:
    SerializedScoresCtor( 
        const extspsl::Pslvector*** scores, 
        int npairplvs, int ntbls, int card, 
        const extspsl::Pslvector& prblvs, const extspsl::Pslvector& levels
    );

    __host__ __device__ void GetEth12( float* eth1, float* eth2 ) const
    {
        SerializedScoresSM<TScore>::GetEth12( eth1, eth2 );
    }

    virtual __host__ __device__ ~SerializedScoresCtor();

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
SerializedScoresCtor<TScore>::~SerializedScoresCtor()
{
    CUMSG( "SerializedScoresCtor::~SerializedScoresCtor", 3 );
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
SerializedScoresCtor<TScore>::SerializedScoresCtor( 
        const extspsl::Pslvector*** scores,
        int
#ifdef __DEBUG__
        nopairplvs
#endif
        , int
#ifdef __DEBUG__
        notbls
#endif
        , int card,
        const extspsl::Pslvector& prblvs, const extspsl::Pslvector& levels )
:   SerializedScoresSM<TScore>()
{
    CUMSG( "SerializedScoresCtor::SerializedScoresCtor", 3 );
#ifdef __DEBUG__
    if( !scores || nopairplvs < 0 || notbls < 1 || card < 1 )
        throw MYRUNTIME_ERROR("SerializedScoresCtor::SerializedScoresCtor: Invalid arguments.");
#endif
    nplvs_ = SLC_MAX( 1, prblvs.GetSize());
    nenos_ = levels.GetSize();
    card_ = card;

    int npairplvscalc = ( nplvs_ * ( nplvs_ + 1 )) >> 1;
    ntbls_ = ( nenos_ * ( nenos_ + 1 )) >> 1;
    nelems_ = ( card_ * ( card_+1 )) >> 1;//number of entries in one table

#ifdef __DEBUG__
    if( npairplvscalc != nopairplvs || ntbls_ != notbls )
        throw MYRUNTIME_ERROR("SerializedScoresCtor::SerializedScoresCtor: Inconsistent argument values.");
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
__host__ __device__
inline
void SerializedScoresCtor<TScore>::DestroyScores()
{
    CUMSG( "SerializedScoresCtor::DestroyScores", 3 );
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
void SerializedScoresCtor<TScore>::NewScores( int sztotal )
{
    CUMSG( "SerializedScoresCtor::NewScores", 3 );
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ )
        throw MYRUNTIME_ERROR("SerializedScoresCtor::NewScores: Not enough memory.");
    szalloc_ = sztotal;
}

#endif//__SerializedScoresCtor_h__
