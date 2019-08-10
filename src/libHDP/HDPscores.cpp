/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "extsp/psl.h"
#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "liblib/myfiler.h"
#include "HDPscores.h"

using namespace extspsl;

// =========================================================================
//global object for application
HDPscores HDPSCORES;
// HDPscores HDPctSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPscores::HDPscores()
:   scores_( NULL ),
    noplvs_( 0 ),
    notbls_( 0 ),
    prblvs_(),
    levels_(),
    naval_( -9999.0f ),
    card_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
HDPscores::~HDPscores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void HDPscores::ReadScores( const char* filename )
{
    FILE* fp = NULL;
    myruntime_error mre;
    mystring msg = "Reading ";
    msg += filename;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL )
        throw MYRUNTIME_ERROR( mystring("HDPscores::ReadScores: Failed to open file ")+
                filename );
    try {
        message( msg.c_str());
        ReadLevels( fp );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );
    if( mre.isset())
        throw mre;

    ReadScoreTables( filename );
}

// -------------------------------------------------------------------------
// ReadScores: read probability levels and those of eff. no. sequences
//
void HDPscores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "HDPscores::ReadLevels: ";
    const char*     p, *pp, *ppp;
    int             emsg;
    const char*     patstrhdps = "hdps=";
    const char*     patstrpl = "+";
    float rdval;

    prblvs_.Clear();
    levels_.Clear();
    //read levels of eff. no. sequences
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));

    if( feof(fp) || !length )
        throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

    if(( p = strstr( locbuffer, patstrhdps )) == NULL )
        throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

    p += strlen( patstrhdps );

    if( length <= (size_t)( p - locbuffer ))
        throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

    if(( pp = strstr( locbuffer, patstrpl )) != NULL ) {
        //read probability levels
        for( ; p < pp; )
        {
            if(( emsg = read_float( p, length-(size_t)(p-locbuffer), &rdval, &rbts )) != 0 ) {
                if( emsg == ERR_RD_NOVL )
                    break;
                throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));
            }
            if( rdval < 0.0f || 100.0f <= rdval )
                throw MYRUNTIME_ERROR( preambl + "Invalid probability threshold.");

            prblvs_.Push( rdval * 0.01f );
            for( p += rbts; p < pp && (*p == ' ' || *p == '\t'); p++);
        }
        p = pp + strlen( patstrpl );
    }

    ppp = locbuffer + length;
    for( ; p < ppp; p += rbts )
    {
        //read levels of eff. no. sequences
        if(( emsg = read_float( p, length-(size_t)(p-locbuffer), &rdval, &rbts )) != 0 ) {
            if( emsg == ERR_RD_NOVL )
                break;
            throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));
        }
        levels_.Push( rdval );
    }
}

// -------------------------------------------------------------------------
// ReadScoreTables: read score tables from file(s)
//
void HDPscores::ReadScoreTables( const char* filename )
{
    FILE*       fp = NULL;
    const mystring preambl = "HDPscores::ReadScoreTables: ";
    myruntime_error mre;
    mystring    fname;
    char        locbuf[KBYTE] = {0};
    int szprblvs = SLC_MAX( 1, prblvs_.GetSize());
    int n, m, k1, k2, plvl, lvl, noplvs, notbls;

    if( !filename )
        return;

    if( szprblvs < 1 || 10 < szprblvs )
        throw MYRUNTIME_ERROR( preambl + 
                              "Invalid number of hdps probability level values.");
    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw MYRUNTIME_ERROR( preambl + 
                              "Invalid number of hdps level values.");

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
                        sprintf( locbuf, "%d%d%d%d",
                                ( int )rintf( prblvs_.GetValueAt(k1) * 100.0f),
                                ( int )rintf( prblvs_.GetValueAt(k2) * 100.0f),
                                ( int )rintf( levels_.GetValueAt(n)),
                                ( int )rintf( levels_.GetValueAt(m)));
                    else
                        sprintf( locbuf, "%d%d",( int )rintf( levels_.GetValueAt(n)),
                                ( int )rintf( levels_.GetValueAt(m)));

                    fname = filename + mystring( locbuf );
                    if(( fp = fopen( fname.c_str(), "r" )) == NULL )
                        throw MYRUNTIME_ERROR( preambl + mystring("Failed to open file ")+ fname );

                    try {
                        ReadScoresHelper( fp, scores_[plvl][lvl++]);
                    } catch( myexception const& ex ) {
                        mre = ex;
                    }

                    fclose( fp );
                    if( mre.isset())
                        throw mre;
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file descriptor
//
void HDPscores::ReadScoresHelper( FILE* fp, Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "HDPscores::ReadScoresHelper: ";
    myruntime_error mre;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = (int)strlen( patstrNA );
    int intval, noelems, card, n, c;
    float rdval;

    try {
        if( scos == NULL )
            throw MYRUNTIME_ERROR( preambl + "Memory access error.");

        //read cardinality
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrcard )) == NULL )
            throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

        p += strlen( patstrcard );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preambl + "Wrong file format.");

        if(( emsg = read_integer( p, length-(size_t)(p-locbuffer), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));

        if( intval < 1 || 4000 < intval )
            throw MYRUNTIME_ERROR( preambl + "Invalid cardinality value.");

        card = intval;
        if(!GetCardinality())
            SetCardinality( card );
        else if( card != GetCardinality())
            throw MYRUNTIME_ERROR( preambl + "Inconsistent score tables.");

        //allocate space for vector
        noelems = card + (( card * ( card-1 )) >> 1 );
        scos->Clear();
        scos->Allocate( noelems );

        //read scores
        for( n = 0; n < card; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));

            if( !buffer.length())
                throw MYRUNTIME_ERROR( preambl + "Short of scores.");

            p = buffer.c_str();
            pp = p + buffer.length();
            for( c = 0; c <= n; c++, p += rbts )
            {
                if( pp <= p )
                    throw MYRUNTIME_ERROR( preambl + "Short of scores.");

                //check for NA
                for( ; *p == ' ' || *p == '\t'; p++ );
                if( strncmp( p, patstrNA, lenpatstrNA ) == 0 ) {
                    //NA value
                    rbts = lenpatstrNA;
                    scos->Push( GetNAvalue());
                    continue;
                }

                if(( emsg = read_float( p, 
                      buffer.length()-(size_t)(p-buffer.c_str()), &rdval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preambl + TranslateReadError(emsg));

                scos->Push( rdval );
            }
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset()) {
        scos->Clear();
        throw mre;
    }
}

