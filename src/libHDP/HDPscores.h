/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __HDPscores__
#define __HDPscores__

#include <stdio.h>
#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/mybase.h"

//global object for application
class HDPscores;
extern HDPscores HDPSCORES;
// extern HDPscores HDPctSCORES;

// -------------------------------------------------------------------------
// class HDPscores: Hierarchical Dirichlet Process model-based (derived) 
// scores
//
class HDPscores
{
public:
    HDPscores();
    ~HDPscores();

    int             GetCardinality() const { return card_; }
    int             GetNoPrbLvs() const { return noplvs_; }
    int             GetNoTables() const { return notbls_; }

    const extspsl::Pslvector*** GetScores() const { return const_cast<const extspsl::Pslvector***>(scores_); }
    const extspsl::Pslvector& GetPrbLevels() const { return const_cast<const extspsl::Pslvector&>(prblvs_); }
    const extspsl::Pslvector& GetLevels() const { return const_cast<const extspsl::Pslvector&>(levels_); }

    float           GetScore( int row, int col, 
                      float fstprb, float secprb, float fstens, float secens ) const;
    void            ReadScores( const char* filename );

    bool            ScoreNA( float value ) const { return value <= naval_; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( const char* filename );
    void            ReadScoresHelper( FILE*, extspsl::Pslvector* scos );

    void            SetCardinality( int value ) { card_ = value; }
    float           GetNAvalue() const { return naval_; }

    extspsl::Pslvector* GetScores( int plvl, int lvl );
    const extspsl::Pslvector* GetScores( int plvl, int lvl ) const;
    void            NewScores( int noplvs, int nolevels );
    void            DestroyScores();

private:
    extspsl::Pslvector*** scores_;//scores
    int         noplvs_;//number of probability levels
    int         notbls_;//number of score tables for each prob. level
    extspsl::Pslvector prblvs_;//probability (of clusters) levels 
    extspsl::Pslvector levels_;//levels of eff. no. sequences
    const float naval_;//NA value
    int         card_;//cardinality
};

// -------------------------------------------------------------------------
// INLINES
//
// DestroyScores: destroy all score tables
inline
void HDPscores::DestroyScores()
{
    int n, m;
    if( scores_ ) {
        for( m = 0; m < noplvs_; m++ ) {
            if( scores_[m] == NULL )
                continue;
            for( n = 0; n < notbls_; n++ ) {
                if( scores_[m][n])
                    delete scores_[m][n];
                scores_[m][n] = NULL;
            }
            free( scores_[m] );
            scores_[m] = NULL;
        }
        free( scores_ );
        scores_ = NULL;
        notbls_ = 0;
    }
    noplvs_ = 0;
}
// NewScores: allocate new score table for each level of eff. no sequences
//  noplvs: no. probability level values
//  nolevels: no. level values of eff. no. sequences
inline
void HDPscores::NewScores( int noplvs, int nolevels )
{
    int n, m;
    DestroyScores();
    if( noplvs < 1 || nolevels < 1 )
        return;
    scores_ = (extspsl::Pslvector***)malloc( noplvs * sizeof(extspsl::Pslvector**));
    if( scores_ == NULL )
        throw MYRUNTIME_ERROR("HDPscores::NewScores: Not enough memory.");
    noplvs_ = noplvs;
    notbls_ = nolevels;

    for( m = 0; m < noplvs_; m++ ) {
        scores_[m] = (extspsl::Pslvector**)malloc( nolevels * sizeof(extspsl::Pslvector*));
        if( scores_[m] == NULL )
            throw MYRUNTIME_ERROR("HDPscores::NewScores: Not enough memory.");
        for( n = 0; n < notbls_; n++ ) {
            scores_[m][n] = new extspsl::Pslvector();
            if( scores_[m][n] == NULL )
                throw MYRUNTIME_ERROR("HDPscores::NewScores: Not enough memory.");
        }
    }
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of probability and 
//  eff. no seqns
inline
extspsl::Pslvector* HDPscores::GetScores( int plvl, int lvl )
{
    if( plvl < 0 || noplvs_ <= plvl )
        throw MYRUNTIME_ERROR("HDPscores::GetScores: Memory access error.");
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("HDPscores::GetScores: Memory access error.");
    return scores_[plvl][lvl];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of probability and 
//  eff. no seqns
inline
const extspsl::Pslvector* HDPscores::GetScores( int plvl, int lvl ) const
{
    if( plvl < 0 || noplvs_ <= plvl )
        throw MYRUNTIME_ERROR("HDPscores::GetScores: Memory access error.");
    if( lvl < 0 || notbls_ <= lvl )
        throw MYRUNTIME_ERROR("HDPscores::GetScores: Memory access error.");
    return scores_[plvl][lvl];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by row, col; table is 
//  selected according to probability levels `fstprb' and `secprb' and eff. 
//  no. sequences given by `fstens' and `secens' 
inline
float HDPscores::GetScore( int row, int col, 
    float fstprb, float secprb, float fstens, float secens ) const
{
    int     e, ee, ln, pi, ppi, plvl;//level indices
    int     n, tmp;
    int     noplvs = SLC_MAX( 1, prblvs_.GetSize());
    int     nolvs = levels_.GetSize();
    float   dtmp;
    float   pi1, pi2, ppi1, ppi2;
    float   e1, e2, ee1, ee2;
    const extspsl::Pslvector* scos = NULL;
    static char locbuf[KBYTE];

    if( row < 0 || col < 0 || card_ <= row || card_ <= col )
        throw MYRUNTIME_ERROR("HDPscores::GetScore: Memory access error.");

    if( noplvs < 1 )
        throw MYRUNTIME_ERROR("HDPscores::GetScore: No probability levels.");
    if( prblvs_.GetSize())
        if( fstprb < prblvs_.GetValueAt(0) || secprb < prblvs_.GetValueAt(0))
            return 0.0f;
    if( secprb < fstprb ) {
        dtmp = secprb; secprb = fstprb; fstprb = dtmp; 
    }

    if( nolvs < 1 )
        throw MYRUNTIME_ERROR("HDPscores::GetScore: No levels of eff. no. sqns.");
    if( fstens < levels_.GetValueAt(0) || secens < levels_.GetValueAt(0))
        return 0.0f;
    if( secens < fstens ) {
        dtmp = secens; secens = fstens; fstens = dtmp; 
    }

    //get corresponding score table
    //index corresponding to probability level
    plvl = 0;
    if( prblvs_.GetSize()) {
        for( pi = 0; pi < noplvs; pi++ ) {
            pi1 = prblvs_.GetValueAt(pi);
            pi2 = 1.01f;
            if( pi + 1 < noplvs )
                pi2 = prblvs_.GetValueAt(pi+1);
            if( pi1 <= fstprb && fstprb < pi2 )
                break;
        }
        for( ppi = pi; ppi < noplvs; ppi++ ) {
            ppi1 = prblvs_.GetValueAt(ppi);
            ppi2 = 1.01f;
            if( ppi + 1 < noplvs )
                ppi2 = prblvs_.GetValueAt(ppi+1);
            if( ppi1 <= secprb && secprb < ppi2 )
                break;
        }
        if( noplvs <= pi || noplvs <= ppi || ppi < pi ) {
            sprintf( locbuf, "Segment covering probability level value not found: %g vs %g.", 
                     fstprb, secprb );
            throw MYRUNTIME_ERROR("HDPscores::GetScore: " + mystring( locbuf ));
        }
        //index calculated as sum of arithm. series
        plvl = (( pi*( TIMES2(noplvs)-pi+1 ))>>1 ) + ppi-pi;
    }
    //index corresponding to level of eff. no. sequences
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
        sprintf( locbuf, "Segment covering eff. no. sequences not found: %g vs %g.", fstens, secens );
        throw MYRUNTIME_ERROR("HDPscores::GetScore: " + mystring( locbuf ));
    }
    //index calculated as sum of arithm. series
    ln = (( e*( TIMES2(nolvs)-e+1 ))>>1 ) + ee-e;

    scos = GetScores( plvl, ln );
    if( scos == NULL )
        throw MYRUNTIME_ERROR("HDPscores::GetScore: Null score table.");

    //calculate index within score table
    if( row < col ) {
        tmp = row; row = col; col = tmp; 
    }
    n = col;
    if( row )
        n += row + (( row * ( row-1 )) >> 1 );
    if( scos->GetSize() <= n )
        throw MYRUNTIME_ERROR("HDPscores::GetScore: Memory access error.");
    return scos->GetValueAt(n);
}

#endif//__HDPscores__
