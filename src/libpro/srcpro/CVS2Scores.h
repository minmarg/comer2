/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CVS2Scores_h__
#define __CVS2Scores_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <cmath>

#include "extsp/psl.h"
#include "extsp/pslvector.h"

//global object declaration
class CVS2Scores;
extern CVS2Scores CVS2SCORES;

// -------------------------------------------------------------------------
// class CVS2Scores: implementation of the translation of context normal 
// vector scores to alignment scores
//
class CVS2Scores
{
public:
    CVS2Scores();
    virtual ~CVS2Scores();

    float           GetScore( float cvscore, float fstens, float secens ) const;
    void            ReadScores( const char* filename );

    float           GetCVSWeight() const { return cvsweight_; }
    void            SetCVSWeight( float value ) { cvsweight_ = value; }

    int             GetNoTables() const { return notbls_; }

    bool            ScoreNA( float value ) const { return value <= naval_; }
    bool            ScorePlus( float value ) const { return plusval_ <= value; }

    int             GetStep() const { return step_; }
    float           GetNAvalue() const { return naval_; }
    float           GetValuePlus() const { return plusval_; }

    const extspsl::Pslvector& GetLevels() const { return levels_; }

    const extspsl::Pslvector* GetKeys( int lvl ) const;
    const extspsl::Pslvector* GetScores( int lvl ) const;

    int             GetShift( int lvl ) const;
    float           GetScale( int lvl ) const;

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( FILE* fp );
    void            ReadScoresHelper( FILE* fp, extspsl::Pslvector* keys, extspsl::Pslvector* scos, 
                                      int* shift, int* card, float* scale );
    void            SmoothScores( extspsl::Pslvector* keys, extspsl::Pslvector* scos, 
                                  int* shift, int* card );

    extspsl::Pslvector* GetKeys( int lvl );
    extspsl::Pslvector* GetScores( int lvl );

    void            NewScores( int nolevels );
    void            DestroyScores();

private:
    float       cvsweight_;//weight of translated vector scores
    extspsl::Pslvector** keys_;//scores of normal vectors
    extspsl::Pslvector** scores_;//translated scores
    int*        shifts_;//number of negative keys for each table
    extspsl::Pslvector levels_;//levels of eff. no. sequences
    int         notbls_;//number of score tables 
    int*        cards_;//cardinalities for each table
    float*      scale_;//scale factors for each table
    const int   step_;//number of intermediate values +1 between two adjacent integers
    const float naval_;//NA value
    const float plusval_;//value for '+'
};

// -------------------------------------------------------------------------
// INLINES
//
// DestroyScores: destroy all score tables
inline
void CVS2Scores::DestroyScores()
{
    int m;
    if( shifts_ ) {
        free( shifts_ );
        shifts_ = NULL;
    }
    if( cards_ ) {
        free( cards_ );
        cards_ = NULL;
    }
    if( scale_ ) {
        free( scale_ );
        scale_ = NULL;
    }
    if( keys_ ) {
        for( m = 0; m < notbls_; m++ ) {
            if( keys_[m] == NULL )
                continue;
            delete keys_[m];
            keys_[m] = NULL;
        }
        free( keys_ );
        keys_ = NULL;
    }
    if( scores_ ) {
        for( m = 0; m < notbls_; m++ ) {
            if( scores_[m] == NULL )
                continue;
            delete scores_[m];
            scores_[m] = NULL;
        }
        free( scores_ );
        scores_ = NULL;
    }
    notbls_ = 0;
}

// NewScores: allocate new score table for each level of eff. no sequences;
//  nolevels: number of level values for eff. no. sequences
inline
void CVS2Scores::NewScores( int nolevels )
{
    int m;
    DestroyScores();
    if( nolevels < 1 )
        return;
    keys_ = ( extspsl::Pslvector** )malloc( nolevels * sizeof( extspsl::Pslvector* ));
    scores_ = ( extspsl::Pslvector** )malloc( nolevels * sizeof( extspsl::Pslvector* ));
    shifts_ = ( int* )malloc( nolevels * sizeof( int ));
    cards_ = ( int* )malloc( nolevels * sizeof( int ));
    scale_ = ( float* )malloc( nolevels * sizeof( float ));
    if( keys_ == NULL || scores_ == NULL || shifts_ == NULL || cards_ == NULL || scale_ == NULL )
        throw MYRUNTIME_ERROR("CVS2Scores::NewScores: Not enough memory.");
    memset( shifts_, 0, nolevels * sizeof( int ));
    memset( cards_, 0, nolevels * sizeof( int ));
    notbls_ = nolevels;

    for( m = 0; m < notbls_; m++ )
        scale_[m] = 1.0f;

    for( m = 0; m < notbls_; m++ ) {
        keys_[m] = new extspsl::Pslvector();
        scores_[m] = new extspsl::Pslvector();
        if( keys_[m] == NULL || scores_[m] == NULL )
            throw MYRUNTIME_ERROR("CVS2Scores::NewScores: Not enough memory.");
    }
}

// -------------------------------------------------------------------------
// GetKeys: get keys corresponding to the level of eff. no sequences
inline
extspsl::Pslvector* CVS2Scores::GetKeys( int lvl )
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetKeys: Memory access error.");
#endif
    return keys_[lvl];
}

// GetKeys: get keys corresponding to the level of eff. no sequences
inline
const extspsl::Pslvector* CVS2Scores::GetKeys( int lvl ) const
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetKeys: Memory access error.");
#endif
    return keys_[lvl];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the level of eff. no 
// sequences
inline
extspsl::Pslvector* CVS2Scores::GetScores( int lvl )
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetScores: Memory access error.");
#endif
    return scores_[lvl];
}

// GetScores: get score table corresponding to the level of eff. no 
// sequences
inline
const extspsl::Pslvector* CVS2Scores::GetScores( int lvl ) const
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetScores: Memory access error.");
#endif
    return scores_[lvl];
}

// -------------------------------------------------------------------------
// GetShift: get shift value for the table corresponding to the level of 
//  eff. no sequences
inline
int CVS2Scores::GetShift( int lvl ) const
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetShift: Memory access error.");
#endif
    return shifts_[lvl];
}

// -------------------------------------------------------------------------
// GetScale: get scale factor for keys corresponding to the level of 
//  eff. no sequences
inline
float CVS2Scores::GetScale( int lvl ) const
{
#ifdef __DEBUG__
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("CVS2Scores::GetScale: Memory access error.");
#endif
    return scale_[lvl];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by log-odds score 
//  `cvscore' of normal vectors
inline
float CVS2Scores::GetScore( float cvscore, float fstens, float secens ) const
{
    const mystring preamb = "CVS2Scores::GetScore: ";
    int     e, ee, ln;//level indices
    int     n, shft, size;
    int     step = GetStep();
    int     nolvs = levels_.GetSize();
    float   dtmp, scl;
    float   e1, e2, ee1, ee2;
    const extspsl::Pslvector* keys = NULL;
    const extspsl::Pslvector* scos = NULL;
    char    locbuf[BUF_MAX];

    if( step < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid step.");
    if( nolvs < 1 )
        throw MYRUNTIME_ERROR( preamb + "Levels of eff. number of sequences are not given.");
    if( fstens < levels_.GetValueAt(0) || secens < levels_.GetValueAt(0))
        return 0.0f;
    if( secens < fstens ) {
        dtmp = secens; secens = fstens; fstens = dtmp; 
    }

    //get corresponding score table index
    // corresponding to level of eff. no. sequences
    for( e = 0; e < nolvs; e++ ) {
        e1 = levels_.GetValueAt(e);
        e2 = 20.01f;
        if( e + 1 < nolvs )
            e2 = levels_.GetValueAt(e+1);
        if( e1 <= fstens && fstens < e2 )
            break;
    }
    for( ee = e; ee < nolvs; ee++ ) {
        ee1 = levels_.GetValueAt(ee);
        ee2 = 20.01f;
        if( ee + 1 < nolvs )
            ee2 = levels_.GetValueAt(ee+1);
        if( ee1 <= secens && secens < ee2 )
            break;
    }
    if( nolvs <= e || nolvs <= ee || ee < e ) {
        sprintf( locbuf, 
          "Segment covering eff. number of sequences not found: %g vs %g.", fstens, secens );
        throw MYRUNTIME_ERROR( preamb + locbuf );
    }
    //index calculated as sum of arithm. series
    ln = (( e*( TIMES2(nolvs)-e+1 ))>>1 ) + ee-e;

    keys = GetKeys(ln);
    scos = GetScores(ln);
    if( scos == NULL || keys == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null score table.");
    shft = GetShift(ln);
    if( shft < 0 )
        throw MYRUNTIME_ERROR( preamb + "Negative shift.");
    scl = GetScale(ln);
    if( scl <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale.");

    size = scos->GetSize();

    if( size < 1 )
        return 0.0f;

    if( scl != 1.0f )
        cvscore *= scl;

    if( cvscore <= keys->GetValueAt(0))
        return scos->GetValueAt(0);
    if( keys->GetValueAt(size-1) <= cvscore )
        return scos->GetValueAt(size-1);

    //calculate index within score table
    if( step == 2 )
        n = shft + (int)rintf(TIMES2(cvscore));
    else if( step == 1 )
        n = shft + (int)rintf(cvscore);
    else
        n = shft + (int)rintf((float)step*cvscore);

    if( size <= n || n < 0 )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    return scos->GetValueAt(n);
}

#endif//__CVS2Scores_h__
