/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <math.h>
#include <cmath>

#include "SSSScores.h"

// =========================================================================
//global object definition
SSSScores SSSSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
SSSScores::SSSScores()
:   VirtScores(),
    sssweight_(0.0f)
{
}

// -------------------------------------------------------------------------
// destructor:
//
SSSScores::~SSSScores()
{
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void SSSScores::ReadScores( const char* filename )
{
    FILE* fp = NULL;
    myruntime_error mre;
    mystring msg = "Reading ";
    msg += filename;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL )
        throw MYRUNTIME_ERROR( 
              mystring("SSSScores::ReadScores: Failed to open file ") + filename );

    try {
        message( msg.c_str());
        ReadLevels( fp );
        ReadScoreTables( fp );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );
    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// ReadLevels: read probability levels and those of effective number of 
// sequences
// TODO: reimplement large arrays to use heap (locbuffer)
//
void SSSScores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preamb = "SSSScores::ReadLevels: ";
    myruntime_error mre;
    const char*     p, *pp, *ppp;
    int             emsg;
    const char*     patstrssss = "ssss=";
    const char*     patstrpl = "+";
    float fpval;

    try {
        prblvs_.Clear();
        levels_.Clear();
        //try read levels of eff. no. sequences
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrssss )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrssss );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( pp = strstr( locbuffer, patstrpl )) != NULL ) {
            //read probability levels
            for( ; p < pp; )
            {
                if(( emsg = read_float( p, length - size_t( p - locbuffer ), &fpval, &rbts )) != 0 ) {
                    if( emsg == ERR_RD_NOVL )
                        break;
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                }
                if( fpval < 0.0f || 100.0f <= fpval )
                    throw MYRUNTIME_ERROR( preamb + "Invalid probability threshold.");

                prblvs_.Push( fpval/100.0f );
                for( p += rbts; p < pp && (*p == ' ' || *p == '\t'); p++);
            }
            p = pp + strlen( patstrpl );
        }

        ppp = locbuffer + length;
        for( ; p < ppp; p += rbts )
        {
            //read levels of eff. no. sequences
            if(( emsg = read_float( p, length - size_t( p - locbuffer ), &fpval, &rbts )) != 0 ) {
                if( emsg == ERR_RD_NOVL )
                    break;
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
            }
            levels_.Push( fpval );
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset()) {
        throw mre;
    }
}

// -------------------------------------------------------------------------
// ReadScoreTables: read score tables from files
// TODO: reimplement large arrays to use heap (locbuffer)
//
void SSSScores::ReadScoreTables( FILE* fp )
{
    size_t          length;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preamb = "SSSScores::ReadScoreTables: ";
    myruntime_error mre;
    const char*     p;
    int             emsg;
    char            levcmb[BUF_MAX] = {0};
    int szprblvs = SLC_MAX( 1, prblvs_.GetSize());
    int n, m, k1, k2, plvl, lvl, noplvs, notbls;

    if( fp == NULL )
        return;

    if( szprblvs < 1 || 10 < szprblvs )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of probability level values.");
    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid number of levels for effective number of sequences.");

    noplvs = ( szprblvs *( szprblvs + 1 )) >> 1;
    notbls = ( levels_.GetSize() *( levels_.GetSize() + 1 )) >> 1;
    NewScores( noplvs, notbls );

    plvl = 0;
    for( k1 = 0; k1 < szprblvs; k1++ ) {
        for( k2 = k1; k2 < szprblvs; k2++, plvl++ ) {
            lvl = 0;
            for( n = 0; n < levels_.GetSize(); n++ ) {
                for( m = n; m < levels_.GetSize(); m++ ) 
                {
                    if( prblvs_.GetSize())
                        sprintf( levcmb, "%02d%02d %d%d",
                                ( int )rintf( prblvs_.GetValueAt(k1) * 100.0f),
                                ( int )rintf( prblvs_.GetValueAt(k2) * 100.0f),
                                ( int )rintf( levels_.GetValueAt(n)),
                                ( int )rintf( levels_.GetValueAt(m)));
                    else
                        sprintf( levcmb, "0000 %d%d",
                                ( int )rintf( levels_.GetValueAt(n)),
                                ( int )rintf( levels_.GetValueAt(m)));

                    try {
                        //read next level values and then scores
                        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                        if( feof( fp ) || !length )
                            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

                        if(( strncmp( levcmb, p = locbuffer, strlen(levcmb))) != 0 )
                            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

                        p += strlen( levcmb );
                        for( ; p && (*p == ' '||*p == '\t'||*p == '\r'||*p == '\n'); p++);

                        if( size_t( p - locbuffer ) < length )
                            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

                    } catch( myexception const& ex ) {
                        mre = ex;
                    }

                    if( mre.isset())
                        throw mre;

                    ReadScoresHelper( fp, scores_[plvl][lvl++]);
                }
            }
        }
    }
}
