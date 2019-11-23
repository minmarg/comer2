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

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "VirtScores.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
VirtScores::VirtScores()
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
VirtScores::~VirtScores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file descriptor
// TODO: reimplement large arrays to use heap (locbuffer)
//
void VirtScores::ReadScoresHelper( FILE* fp, extspsl::Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preamb = "VirtScores::ReadScoresHelper: ";
    myruntime_error mre;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = (int)strlen(patstrNA);
    int     intval, noelems, card, n, c;
    float   fpval;

    try {
        if( scos == NULL )
            throw MYRUNTIME_ERROR( preamb + "Memory access error.");

        //read cardinality
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrcard )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrcard );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( intval <= 0 || 4000 < intval )
            throw MYRUNTIME_ERROR( preamb + "Invalid cardinality.");

        card = intval;
        if(!GetCardinality())
            SetCardinality( card );
        else if( card != GetCardinality())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent score tables.");

        //allocate space for vector
        noelems = card + (( card * ( card-1 )) >> 1 );
        scos->Clear();
        scos->Allocate( noelems );

        //read scores
        for( n = 0; n < card; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( !buffer.length())
                throw MYRUNTIME_ERROR( preamb + "Short of scores.");

            p = buffer.c_str();
            pp = p + buffer.length();
            for( c = 0; c <= n; c++, p += rbts )
            {
                if( pp <= p )
                    throw MYRUNTIME_ERROR( preamb + "Short of scores.");

                //check for NA
                for( ; *p == ' ' || *p == '\t'; p++ );
                if( strncmp( p, patstrNA, lenpatstrNA ) == 0 ) {
                    //NA value
                    rbts = lenpatstrNA;
                    scos->Push( GetNAvalue());
                    continue;
                }

                if(( emsg = read_float( p, 
                            buffer.length() - size_t( p - buffer.c_str()), &fpval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                scos->Push( fpval );
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
