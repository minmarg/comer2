/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedCVS2Scores_h__
#define __SerializedCVS2Scores_h__

#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"

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
protected:
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
    __host__ __device__ SerializedCVS2Scores(
        TScore* scores, int szalloc, char nenos, int ntbls )
    :   h_scores_(scores), szalloc_(szalloc), 
        nenos_(nenos), ntbls_(ntbls)
    {};

    __host__ __device__ SerializedCVS2Scores();
    virtual __host__ __device__ ~SerializedCVS2Scores();

    __host__ __device__ const TScore* GetScores() const { return h_scores_; }
    __host__ __device__ int GetSizeAlloc() const { return szalloc_; }
    __host__ __device__ char GetNENOLevs() const { return nenos_; }
    __host__ __device__ int GetNTables() const { return ntbls_; }

    __host__ __device__ void GetCardShiftScale( 
        float fstens, float secens, int* card, int* shft, float* scale,
        TScore* keyfirst, TScore* valuefirst, TScore* keylast, TScore* valuelast) const;
    __host__ __device__ void Get1stCardShiftScale( 
        int* card, int* shft, float* scale,
        TScore* key1first = NULL, TScore* value1first = NULL, 
        TScore* key1last = NULL, TScore* value1last = NULL) const;

    __host__ __device__ TScore GetScore( TScore cvscore, float fstens, float secens ) const;
    __host__ __device__ TScore GetScoreE1S1Step2( 
        TScore cvscore, int card, int shft, float, float ) const;
    static __host__ __device__ TScore GetScoreE1S1Step2( 
        const TScore* scores, TScore cvscore, int card, int shft );
    static __host__ __device__ TScore GetScoreE1S1Step2Boundary( 
        const TScore* scores, TScore cvscore, int card, int shft,
        TScore key1first, TScore value1first, 
        TScore key1last, TScore value1last);

protected:
    //{{host/device data
    TScore* h_scores_;//multi-dimensional scores
    //}}
    int szalloc_;//size (bytes) allocated for scores
    char nenos_;//number of levels for eff. number of observations
    int ntbls_;//number of tables per a pair of eno level values
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore>
__host__ __device__
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
__host__ __device__
inline
SerializedCVS2Scores<TScore>::~SerializedCVS2Scores()
{
    CUMSG( "SerializedCVS2Scores::~SerializedCVS2Scores", 3 );
}

// -------------------------------------------------------------------------
// GetCardShiftScale: read cardinality, shift, and scale written in the 
// buffer of scores; 
// 
template <typename TScore>
__host__ __device__ 
inline
void SerializedCVS2Scores<TScore>::GetCardShiftScale( 
    float fstens, float secens, int* card, int* shft, float* scale,
    TScore* keyfirst, TScore* valuefirst, TScore* keylast, TScore* valuelast) const
{
    CUMSG( "SerializedCVS2Scores::GetCardShiftScale", 5 );
    int e, ee;//level indices

    MYASSERT( h_scores_, "SerializedCVS2Scores::GetCardShiftScale: Memory access error.");

    if( card ) *card = 0;
    if( shft ) *shft = 0;
    if( scale ) *scale = 0.0f;
    if( keyfirst ) *keyfirst = TScore(0);
    if( valuefirst ) *valuefirst = TScore(0);
    if( keylast ) *keylast = TScore(0);
    if( valuelast ) *valuelast = TScore(0);

    if( fstens < (float)h_scores_[0] || secens < (float)h_scores_[0])
        return;

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

    if( card ) *card = (int)h_scores_[nenos_+ee];
    if( shft ) *shft = (int)h_scores_[nenos_+ntbls_+ee];
    if( scale ) *scale = (float)h_scores_[nenos_+ntbls_+ntbls_+ee];

    int crd = (int)h_scores_[nenos_+ee];
    int beg = nenos_+ntbls_*3;//beginning index for scores
    for( int i = 0; i < ee; i++ )
        //there are two series of scores: keys and values
        beg += (int)h_scores_[nenos_+i]*2;

    //the fist key value
    if( keyfirst ) *keyfirst = h_scores_[beg];
    //the first value from the score map of interest
    if( valuefirst ) *valuefirst = h_scores_[beg+crd];

    //the last key value 
    if( keylast ) *keylast = h_scores_[beg+crd-1];
    //the last value from the map of interest
    if( valuelast ) *valuelast = h_scores_[beg+crd+crd-1];
}

// -------------------------------------------------------------------------
// GetCardShiftScale: read cardinality, shift, and scale of the first table 
// found in the buffer of scores; 
// 
template <typename TScore>
__host__ __device__ 
inline
void SerializedCVS2Scores<TScore>::Get1stCardShiftScale( 
    int* card, int* shft, float* scale,
    TScore* key1first, TScore* value1first, TScore* key1last, TScore* value1last) const
{
    CUMSG( "SerializedCVS2Scores::Get1stCardShiftScale", 5 );

    MYASSERT( h_scores_, "SerializedCVS2Scores::Get1stCardShiftScale: Memory access error.");

    float fstens = (float)h_scores_[0]; 
    float secens = (float)h_scores_[0];
    GetCardShiftScale( 
        fstens, secens, card, shft, scale, 
        key1first, value1first, key1last, value1last );
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified/mapped by cvscore; 
// table is selected according to eff. no. observations given by 
// `fstens' and `secens' 
// 
template <typename TScore>
__host__ __device__
inline
TScore SerializedCVS2Scores<TScore>::GetScore( 
    TScore cvscore, float fstens, float secens ) const
{
    CUMSG( "SerializedCVS2Scores::GetScore", 5 );
    int e, ee;//level indices

    MYASSERT( h_scores_, "SerializedCVS2Scores::GetScore: Memory access error.");

    if( fstens < (float)h_scores_[0] || secens < (float)h_scores_[0])
        return TScore(0);

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
        return h_scores_[beg+card];
    }
    //if the last key value <= cvscore
    if( (float)h_scores_[beg+card-1] <= cvscore ) {
        //return the last value from the scores
        return h_scores_[beg+card+card-1];
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
    return h_scores_[beg+card];
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
__host__ __device__
inline
TScore SerializedCVS2Scores<TScore>::GetScoreE1S1Step2( 
    TScore cvscore, int card, int shft, float, float ) const
{
    CUMSG( "SerializedCVS2Scores::GetScoreE1S1Step2", 5 );
    MYASSERT( h_scores_,
        "SerializedCVS2Scores::GetScoreE1S1Step2: Memory access error.");

    //beginning index for scores: 1 for #levels, 3 for 
    // cardinality, shift, and scale
    //int beg = 1+3;

    //if cvscore <= the fist key value
    if( cvscore <= (float)h_scores_[4/*beg*/]) {
        //return the first value from the scores
        return h_scores_[4/*beg*/+card];
    }
    //if the last key value <= cvscore
    if( (float)h_scores_[4/*beg*/+card-1] <= cvscore ) {
        //return the last value from the scores
        return h_scores_[4/*beg*/+card+card-1];
    }

    //calculate index within a linear score table;
    //the index is overwritten to shft;
    //( STEP == 2 )
    shft += 4/*beg*/ + (int)rintf(TIMES2(cvscore));

    MYASSERT((shft+card)*(int)sizeof(TScore) < szalloc_,
        "SerializedCVS2Scores::GetScoreE1S1Step2: Index out of range.");

    //return the score found at the index from the scores (+card)
    return h_scores_[shft+card];
}

// -------------------------------------------------------------------------
// GetScoreE1S1Step2 (static version): get score at table position 
// identified/mapped by cvscore; 
// NOTE: a single level for ENOS is assumed to be present, covering the 
// full range of values; 
// hence, only one (the first) score table is referenced;
// the cardinality and score shift are given by parameters card and shift;
// 
template <typename TScore>
__host__ __device__
__forceinline__
TScore SerializedCVS2Scores<TScore>::GetScoreE1S1Step2( 
    const TScore* scores, TScore cvscore, int card, int shft )
{
    CUMSG( "SerializedCVS2Scores::GetScoreE1S1Step2 [static]", 5 );
    //MYASSERT(scores,"SerializedCVS2Scores::GetScoreE1S1Step2: Memory access error.");

    //beginning index for scores: 1 for #levels, 3 for 
    // cardinality, shift, and scale
    //int beg = 1+3;

    //if cvscore <= the fist key value
    if( cvscore <= (float)scores[4/*beg*/]) {
        //return the first value from the scores
        return scores[4/*beg*/+card];
    }
    //if the last key value <= cvscore
    if( (float)scores[4/*beg*/+card-1] <= cvscore ) {
        //return the last value from the scores
        return scores[4/*beg*/+card+card-1];
    }

    //calculate index within a linear score table;
    //the index is overwritten to shft;
    //( STEP == 2 )
    shft += 4/*beg*/ + (int)rintf(TIMES2(cvscore));

    //MYASSERT((shft+card)*(int)sizeof(TScore) < szalloc_,
    //    "SerializedCVS2Scores::GetScoreE1S1Step2: Index out of range.");

    //return the score found at the index from the scores (+card)
    return scores[shft+card];
}

// -------------------------------------------------------------------------
// GetScoreE1S1Step2Boundary: get score at table position 
// identified/mapped by cvscore given boundary keys and values of the first 
// map; 
// NOTE: a single level for ENOS is assumed to be present, covering the 
// full range of values; 
// hence, only one (the first) score table is referenced;
// the cardinality and score shift are given by parameters card and shift;
// 
template <typename TScore>
__host__ __device__
__forceinline__
TScore SerializedCVS2Scores<TScore>::GetScoreE1S1Step2Boundary( 
    const TScore* scores, TScore cvscore, int card, int shft,
    TScore key1first, TScore value1first, 
    TScore key1last, TScore value1last)
{
    CUMSG( "SerializedCVS2Scores::GetScoreE1S1Step2Boundary [static]", 5 );
    //MYASSERT(scores,"SerializedCVS2Scores::GetScoreE1S1Step2: Memory access error.");

    //beginning index for scores: 1 for #levels, 3 for 
    // cardinality, shift, and scale
    //int beg = 1+3;

    //if cvscore <= the fist key value
    if( cvscore <= key1first ) {
        //return the first value from the scores
        return value1first;
    }
    //if the last key value <= cvscore
    if( key1last <= cvscore ) {
        //return the last value from the scores
        return value1last;
    }

    //calculate index within a linear score table;
    //the index is overwritten to shft;
    //( STEP == 2 )
    shft += 4/*beg*/ + (int)rintf(TIMES2(cvscore));

    //MYASSERT((shft+card)*(int)sizeof(TScore) < szalloc_,
    //    "SerializedCVS2Scores::GetScoreE1S1Step2: Index out of range.");

    //return the score found at the index from the scores (+card)
    return scores[shft+card];
}

#endif//__SerializedCVS2Scores_h__
