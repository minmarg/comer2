/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_h__
#define __CuBatchScoreMatrix_h__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"
#include "libHDP/HDPscores.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "CuBatchScoreMatrix_com.h"

template<typename TScore>
class ScoresAttr;

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchScoreMatrix
// Batch computation of profile-profile substitution scores using the CUDA 
// architecture or CPU threads
//
class CuBatchScoreMatrix: AbstractScoreMatrix<CUBSM_TYPE>
{
public:
    CuBatchScoreMatrix( Configuration config[NoCTypes], int precscale = 1 );
    virtual ~CuBatchScoreMatrix();

    virtual CUBSM_TYPE  GetScore( int, int ) const { return CUBSM_Q(0); }
    virtual CUBSM_TYPE  GetModScore( int, int ) const { return CUBSM_Q(0); }


    virtual void ComputeScoreMatrix( bool final = false );

    virtual float       GetFinalScore( CUBSM_TYPE value ) const { return value; }

    virtual void        ScaleScoreMatrix/*ScaleScoringMatrix*/() {;};
    virtual void        ComputeStatisticalParameters( bool = true ) {;};

    virtual void        ComputeScoreProbabilities( ScoresAttr<CUBSM_TYPE>* ) {;};

    virtual void        PrintParameterTable( TPrintFunction, void* /*vpn*/ ) const {;};


    float GetExpGappedLambda() const { return AbstractScoreMatrix<CUBSM_TYPE>::GetExpGappedLambda(); }
    float GetExpGappedK() const { return AbstractScoreMatrix<CUBSM_TYPE>::GetExpGappedK(); }

    float GetRefLambda() const { return AbstractScoreMatrix<CUBSM_TYPE>::GetRefLambda(); }
    float GetRefK() const { return AbstractScoreMatrix<CUBSM_TYPE>::GetRefK(); }

    size_t GetDeltaLength() const {return AbstractScoreMatrix<CUBSM_TYPE>::GetDeltaLength();}
    uint64_mt GetSearchSpace() const {return AbstractScoreMatrix<CUBSM_TYPE>::GetSearchSpace();}
    bool ComputeLengthAdjustment( size_t query_len, uint64_mt db_len, size_t no_sequences ) {
        return AbstractScoreMatrix<CUBSM_TYPE>::ComputeLengthAdjustment( 
            query_len, db_len, no_sequences 
        );
    }

    void ComputeScoreMatrixDevice(
        cudaStream_t streamproc,
        CUBSM_TYPE* sssscores,
        SerializedScoresAttr ssssattr,
        CUBSM_TYPE* cvs2scores,
        SerializedCVS2ScoresAttr cvs2sattr,
        cudaTextureObject_t hdp1sTexo,
        SerializedScoresAttr hdpsattr,
        bool hdp1scoresinuse,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        CUBSM_TYPE* outscores,
        CUBSM_TYPE* outmodscores
    );

    void TESTPrintProProScores1(
        char** cachedbdb1pmbeg, 
        char** cachedbdb1pmend, 
        char** querypmbeg, char** querypmend,
        char** bdb1pmbeg, char** bdb1pmend, char** bdbCpmbeg, char** bdbCpmend,
        unsigned int dbxpad,
        size_t szscores,
        CUBSM_TYPE* outscores
    );

protected:
    explicit CuBatchScoreMatrix();

    virtual void Init( int querylen, int sbjctlen );
    virtual void        SetScore( int, int, CUBSM_TYPE ) {;};
    virtual void        SetModScore( int, int, CUBSM_TYPE ) {;};

    virtual bool        GetAllNegatives() const { return false;}
    virtual void        SetAllNegatives( bool ) {;};

private:
    virtual const CUBSM_TYPE** GetScores/*GetImage*/() const { return NULL; }

private:
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__CuBatchScoreMatrix_h__
