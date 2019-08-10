/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __DbProProbs_h__
#define __DbProProbs_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "extsp/ivector.h"


// -------------------------------------------------------------------------
// class DbProProbs
// 
// Implementation of probabilities for database profiles 
//
class DbProProbs
{
public:
    DbProProbs();
    ~DbProProbs();

    float           GetProbAt( int n ) const;
    float           GetProbAt( int row, int col ) const;
    float           GetProbAt( int e1, int l1, int e2, int l2 ) const;
    float           GetProb( float E1, int L1, float E2, int L2 ) const;
    bool            ProbsReady() const { return probs_ && probs_->GetSize(); }

    const extspsl::Pslvector& GetEffLvs() const { return efflvs_; }
    const extspsl::Ivector& GetLenLvs() const { return lenlvs_; }

    int             GetIndex( float E1, int L1 ) const;
    void            GetLevIndices( float E1, int L1, int* e1, int* l1 ) const;

    int             GetCardinality() const { return card_; }
    bool            Prob0( float value ) const { return value <= 0.0f; }

    void            ReadProbs( const char* filename );
    void            WriteProbs( const char* filename );

    void            ClearMargs();
    void            UpdateMargs( float E1, int L1 );
    void            CalcProbs();

    const extspsl::Pslvector* GetProbs() const;

protected:
    void            ReadProbsHelper( FILE*, extspsl::Pslvector* probs );
    void            WriteProbsHelper( FILE*, extspsl::Pslvector* probs );
    void            SetCardinality( int value ) { card_ = value; }

    extspsl::Pslvector* GetProbs();
    void            NewProbs();
    void            DestroyProbs();
    void            VerifyProbs( extspsl::Pslvector* probs, float acc = 1.e-5f );

protected:
    extspsl::Pslvector* probs_;//probabilities
    extspsl::Ivector    margs_;//marginal pair counts
    extspsl::Pslvector  efflvs_;//levels of eff. no. sequences
    extspsl::Ivector    lenlvs_;//distinct values of profile lengths
    extspsl::Pslvector  midefflvs_;//intermediate levels of eff. no. sequences
    extspsl::Ivector    midlenlvs_;//intermediate values of profile lengths
    int                 card_;//cardinality
};

// -------------------------------------------------------------------------
// INLINES
//
// ClearMargs: clear marginal counts
inline
void DbProProbs::ClearMargs()
{
    margs_.Clear();
    margs_.Reserve( GetCardinality());
}
// DestroyProbs: destroy probabilities
inline
void DbProProbs::DestroyProbs()
{
    if( probs_ ) {
        delete probs_;
        probs_ = NULL;
    }
}
// NewProbs: allocate new probabilities table
inline
void DbProProbs::NewProbs()
{
    DestroyProbs();
    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        return;
    probs_ = new extspsl::Pslvector();
    if( probs_ == NULL )
        throw MYRUNTIME_ERROR("DbProProbs::NewProbs: Not enough memory.");
}

// -------------------------------------------------------------------------
// GetProbs: get probabilities table
inline
extspsl::Pslvector* DbProProbs::GetProbs()
{
    return probs_;
}

// -------------------------------------------------------------------------
// GetProbs: get probabilities table
inline
const extspsl::Pslvector* DbProProbs::GetProbs() const
{
    return probs_;
}

// -------------------------------------------------------------------------
// GetIndex: get nearest level indices for the given E1 and L1
inline
void DbProProbs::GetLevIndices( float E1, int L1, int* e1, int* l1 ) const
{
    if( e1 == NULL || l1 == NULL )
        throw MYRUNTIME_ERROR("DbProProbs::GetLevIndices: Null parameters.");

    for( *e1 = 0; *e1 < midefflvs_.GetSize(); (*e1)++ )//FIXED:(*e1)++
        if( E1 < midefflvs_.GetValueAt(*e1))
            break;
    for( *l1 = 0; *l1 < midlenlvs_.GetSize(); (*l1)++ )//FIXED:(*l1)++
        if( L1 < midlenlvs_.GetValueAt(*l1))
            break;
}

// -------------------------------------------------------------------------
// GetIndex: get index (matrix row/col) given E1 and L1
inline
int DbProProbs::GetIndex( float E1, int L1 ) const
{
    int     e1, l1;//level indices
    int     row;

    for( e1 = 0; e1 < midefflvs_.GetSize(); e1++ )
        if( E1 < midefflvs_.GetValueAt(e1))
            break;
    for( l1 = 0; l1 < midlenlvs_.GetSize(); l1++ )
        if( L1 < midlenlvs_.GetValueAt(l1))
            break;

    row = e1 * lenlvs_.GetSize() + l1;
    return row;
}

// -------------------------------------------------------------------------
// GetProbAt: get probability at the specified location 
inline
float DbProProbs::GetProbAt( int n ) const
{
    if( probs_ == NULL || n < 0 || probs_->GetSize() <= n )
        throw MYRUNTIME_ERROR("DbProProbs::GetProbAt: Memory access error.");
    return probs_->GetValueAt(n);
}

// -------------------------------------------------------------------------
// GetProb: get probability at a table (matrix) position given by 
// row and col
inline
float DbProProbs::GetProbAt( int row, int col ) const
{
    int nn;

    if( probs_ == NULL || probs_->GetSize() < 1 )
        return 0.0f;

    if( row < col )
        throw MYRUNTIME_ERROR("DbProProbs::GetProbAt: Invalid indices.");
    if( row < 0 )
        throw MYRUNTIME_ERROR("DbProProbs::GetProbAt: Negative indices.");

    nn = col;
    if( row )
        nn += row + (( row * ( row-1 )) >> 1 );
    if( probs_->GetSize() <= nn )
        throw MYRUNTIME_ERROR("DbProProbs::GetProbAt: Memory access error.");
    return probs_->GetValueAt(nn);
}

// -------------------------------------------------------------------------
// GetProbAt: get probability at a table (matrix) position identified by 
// level indices e1,l1,e2,l2
inline
float DbProProbs::GetProbAt( int e1, int l1, int e2, int l2 ) const
{
    int row, col;
    row = e1 * lenlvs_.GetSize() + l1;
    col = e2 * lenlvs_.GetSize() + l2;
    return GetProbAt( row, col );
}

// -------------------------------------------------------------------------
// GetProb: get probability at a table (matrix) position identified by 
// E1,L1,E2,L2
inline
float DbProProbs::GetProb( float E1, int L1, float E2, int L2 ) const
{
    int     row, col;
    int     tmp;
//     static char locbuf[KBYTE];

    if( probs_ == NULL || probs_->GetSize() < 1 )
        return 0.0f;

    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        throw MYRUNTIME_ERROR("DbProProbs::GetProb: No probability data.");

    row = GetIndex( E1, L1 );
    col = GetIndex( E2, L2 );
    if( row < col ) {
        tmp = col; col = row; row = tmp;
    }

    return GetProbAt( row, col );
}

#endif//__DbProProbs_h__
