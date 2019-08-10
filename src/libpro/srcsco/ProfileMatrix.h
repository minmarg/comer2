/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __ProfileMatrix_h__
#define __ProfileMatrix_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include "extsp/pslvector.h"
#include "liblib/alpha.h"
#include "libpro/srcpro/Configuration.h"
#include "IntegerScoreMatrix.h"
#include "ScoresAttr.h"

////////////////////////////////////////////////////////////////////////////
// CLASS ProfileMatrix
// Profile matrix template for scaling scores and calculating statistical 
// parameters
//
template<typename TScore>
class ProfileMatrix: public IntegerScoreMatrix<TScore>
{
public:
    ProfileMatrix(
            const float   (*pssmscores)[NUMALPH],
            const char*     ress,
            int             length,
            Configuration   config[NoCTypes],
            int precscale = PP_SCALE_CONSTANT
    );
    virtual ~ProfileMatrix();

    virtual const char* GetMethodName() const { return "Position-specific"; }

    virtual int         GetQuerySize() const { return IntegerScoreMatrix<TScore>::GetQuerySize(); }
    virtual int         GetSubjectSize() const { return IntegerScoreMatrix<TScore>::GetSubjectSize(); }

    const float      ( *GetVectorAt( int n ) const )[NUMALPH];

    virtual void        ComputeScoreMatrix( bool = false );

    virtual void        ComputeScoreProbabilities( ScoresAttr<TScore>* );

    virtual void        OptimizeTargetFrequencies();

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintScoreMatrix( FILE* );
//     virtual void        PrintFinal( TPrintFunction, void* vpn ) const;

protected:
    explicit ProfileMatrix();

    const char*         GetResidues() const { return residues_; }
    char                GetResidueAt( int n ) const;

    void                FillMatrix( const float (*pssmscores)[NUMALPH] );

    //helper routines for optimization of target frequencies
    void                IncorporateTargetFrequencies( const extspsl::Pslvector& );

private:
    const char*         residues_;//vector of residues
    static float        position_[NUMALPH];

    extspsl::Pslvector  scores_;//scores for optimization of target frequencies
    extspsl::Pslvector  rprobs_;//row background probabilities of the score system
    extspsl::Pslvector  cprobs_;//column background probabilities of the score system
};

void ComputedSubMatrixWithParams();

// INLINES ...
// -------------------------------------------------------------------------
// GetResidueAt: get residue at position n
//
template<typename TScore>
inline
char ProfileMatrix<TScore>::GetResidueAt( int n ) const
{
#ifdef __DEBUG__
    if( n < 0 || GetQuerySize() <= n )
        throw MYRUNTIME_ERROR( "ProfileMatrix::GetResidueAt: Memory access error." );
#endif
    return residues_[n];
}

// -------------------------------------------------------------------------
// ComputeScoreMatrix is not to be called from a class object
//
template<typename TScore>
inline
void ProfileMatrix<TScore>::ComputeScoreMatrix( bool )
{
    throw MYRUNTIME_ERROR( "ProfileMatrix::ComputeScoreMatrix: Should not be called." );
}

#endif//__ProfileMatrix_h__
