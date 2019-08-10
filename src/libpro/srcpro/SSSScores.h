/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SSSScores_h__
#define __SSSScores_h__

#include "liblib/mybase.h"
#include "VirtScores.h"

//global object declaration
class SSSScores;
extern SSSScores SSSSCORES;

// -------------------------------------------------------------------------
// class SSSScores: Secondary structure state scores
//
class SSSScores: public VirtScores
{
public:
    SSSScores();
    virtual ~SSSScores();

    float           GetSSScore( char fstsss, char secsss, 
                            float fstprb, float secprb, float fstens, float secens ) const;
    void            ReadScores( const char* filename );

    float           GetSSSWeight() const { return sssweight_; }
    void            SetSSSWeight( float value ) { sssweight_ = value; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( FILE* fp );

private:
    float   sssweight_;//weight of SS state scores
};

// -------------------------------------------------------------------------
// INLINES
//
// -------------------------------------------------------------------------
// GetScore: get score at table position identified by SS states 
//  `fstsss' and `secsss' and probabilities `fstprb' and `secprb'
inline
float SSSScores::GetSSScore( char fstsss, char secsss, 
    float fstprb, float secprb, float fstens, float secens ) const
{
    MYMSG( "SSSScores::GetScore", 5 );
    int row = ( int )fstsss * 10 + int( fstprb * 10.0f );
    int col = ( int )secsss * 10 + int( secprb * 10.0f );
    return VirtScores::GetScore( row, col, fstprb, secprb, fstens, secens );
}

#endif//__SSSScores_h__
