/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "extsp/ivector.h"
#include "DbProProbs.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
DbProProbs::DbProProbs()
:   probs_( NULL ),
    card_( 0 )
{
    float dv;
    int n, nv;
    int card;

    efflvs_.Allocate(10);
    lenlvs_.Allocate(10);
    midefflvs_.Allocate(10);
    midlenlvs_.Allocate(10);

    for( dv = 2.0f; dv <= 14.0f; dv += 2.0f ) efflvs_.Push(dv);
    for( nv = 50; nv <= 400; nv *= 2 ) lenlvs_.Push(nv);
    lenlvs_.Push(600);
    lenlvs_.Push(800);

    for( n = 0; n < efflvs_.GetSize()-1; n++ ) {
        dv = (efflvs_.GetValueAt(n) + efflvs_.GetValueAt(n+1)) * 0.5f;
        midefflvs_.Push(dv);
    }
    for( n = 0; n < lenlvs_.GetSize()-1; n++ ) {
        nv = (lenlvs_.GetValueAt(n) + lenlvs_.GetValueAt(n+1)) / 2;
        midlenlvs_.Push(nv);
    }

    card = efflvs_.GetSize() * lenlvs_.GetSize();

    margs_.Reserve( card );

    SetCardinality( card );
}

// -------------------------------------------------------------------------
// destructor:
//
DbProProbs::~DbProProbs()
{
    DestroyProbs();
}



// -------------------------------------------------------------------------
// ReadProbs: read probabilities from file
//
void DbProProbs::ReadProbs( const char* filename )
{
    FILE* fp = NULL;
    myruntime_error mre;
    mystring msg = "Reading probabilities...";
//     msg += filename;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL ) {
//         throw myruntime_error( mystring("Failed to open file ")+ filename );
        warning("Update of profile pair probabilities ignored.");
        return;
    }

    try {
        message( msg.c_str());
        NewProbs();
        ReadProbsHelper( fp, GetProbs());
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// ReadProbsHelper: helper method read probabilities from file descriptor
// TODO: reimplement local variables to use heap: locbuffer
//
void DbProProbs::ReadProbsHelper( FILE* fp, extspsl::Pslvector* probs )
{
    if( fp == NULL )
        return;

    const mystring  preamb = "DbProProbs::ReadProbsHelper: ";
    myruntime_error mre;
    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = (int)strlen(patstrNA);
    int     intval, noelems, card, n, c;
    float   fpval;

    try {
        if( probs == NULL )
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

        if( intval <= 0 || 100 < intval )
            throw MYRUNTIME_ERROR( preamb + "Invalid cardinality.");

        if( intval != efflvs_.GetSize() * lenlvs_.GetSize())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent cardinality.");

        card = intval;
        if(!GetCardinality())
            SetCardinality( card );
        else if( card != GetCardinality())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent size of probability matrix.");

        //allocate space for vector
        noelems = card + (( card * ( card-1 )) >> 1 );
        probs->Clear();
        probs->Allocate( noelems );

        //read probabilities
        for( n = 0; n < card; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( !buffer.length())
                throw MYRUNTIME_ERROR( preamb + "Missing probabilities.");

            p = buffer.c_str();
            pp = p + buffer.length();
            for( c = 0; c <= n; c++, p += rbts )
            {
                if( pp <= p )
                    throw MYRUNTIME_ERROR( preamb + "Missing probabilities.");

                //check for NA
                for( ; *p == ' ' || *p == '\t'; p++ );
                if( strncmp( p, patstrNA, lenpatstrNA ) == 0 ) {
                    //NA value
                    rbts = lenpatstrNA;
                    probs->Push(0.0f);
                    continue;
                }

                if(( emsg = read_float( p, buffer.length() - size_t( p - buffer.c_str()), &fpval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                probs->Push( fpval );
            }
        }

        VerifyProbs( probs, 1.e-3f );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset()) {
        probs->Clear();
        throw mre;
    }
}



// -------------------------------------------------------------------------
// UpdateMargs: update marginal counts of profile pairs
//
void DbProProbs::UpdateMargs( float E1, int L1 )
{
    const mystring  preamb = "DbProProbs::UpdateMargs: ";
    int e1, l1;//level indices
    int row;

    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Levels are not set.");

    for( e1 = 0; e1 < midefflvs_.GetSize(); e1++ )
        if( E1 < midefflvs_.GetValueAt(e1))
            break;
    for( l1 = 0; l1 < midlenlvs_.GetSize(); l1++ )
        if( L1 < midlenlvs_.GetValueAt(l1))
            break;

    row = e1 * lenlvs_.GetSize() + l1;
    if( margs_.GetSize() <= row )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    margs_.AddValueAt( row, 1 );
}

// -------------------------------------------------------------------------
// CalcProbs: calculate probabilities given marginal pair counts
//
void DbProProbs::CalcProbs()
{
    const mystring  preamb = "DbProProbs::CalcProbs: ";
    const int   card = GetCardinality();
    int         noelems, rowcnt, colcnt;
    int         count = 0;
    int         v, n, c;
    float       count2 = 0.0f;
    float       val;

    if( card < 1 )
        return;

    if( card != margs_.GetSize())
        throw MYRUNTIME_ERROR( preamb + "Inconsistent data sizes.");

    noelems = card + (( card * ( card-1 )) >> 1 );
    NewProbs();
    probs_->Allocate( noelems );

    for( n = 0; n < margs_.GetSize(); n++ )
        count += margs_.GetValueAt(n);

    if( count < 1 )
        throw MYRUNTIME_ERROR( preamb + "No marginal count data.");

    count2 = (float)count * (float)count;

    for( v = n = 0; n < card; n++ ) {
        rowcnt = margs_.GetValueAt(n);
        for( c = 0; c <= n; c++, v++ )
        {
            colcnt = margs_.GetValueAt(c);
            if( rowcnt < 1 || colcnt < 1 )
                val = 0.0f;
            else
                val = (float)rowcnt * (float)colcnt / count2;
            if( c < n )
                val *= 2.0f;
            probs_->Push(val);
        }
    }

    VerifyProbs( probs_ );
}



// -------------------------------------------------------------------------
// VerifyProbs: verify probabilities
//
void DbProProbs::VerifyProbs( extspsl::Pslvector* probs, float acc )
{
    const mystring  preamb = "DbProProbs::VerifyProbs: ";
    static char locbuf[BUF_MAX];
    float sum;

    if( probs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null probabilities.");

    sum = probs->Sum() - 1.0f;
    if( acc < fabsf(sum)) {
        sprintf( locbuf, "%s Probabilities not conserved: %g", preamb.c_str(), sum+1.0f );
        warning( locbuf );
    }
}



// -------------------------------------------------------------------------
// WriteProbs: write probabilities to file
//
void DbProProbs::WriteProbs( const char* filename )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if( GetProbs() == NULL || GetProbs()->GetSize() < 1 )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw MYRUNTIME_ERROR( 
        mystring("DbProProbs::WriteProbs: Failed to open file for writing: ") + filename );

    try {
        WriteProbsHelper( fp, GetProbs());
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// WriteProbsHelper: helper method to write probabilities to file descriptor
//
void DbProProbs::WriteProbsHelper( FILE* fp, extspsl::Pslvector* probs )
{
    if( fp == NULL )
        return;

    const mystring  preamb = "DbProProbs::WriteProbsHelper: ";
//     size_t          length, rbts;
//     const size_t    locsize = TIMES4( KBYTE );
//     char            locbuffer[locsize+1] = {0};
    mystring        buffer;
//     const char*     p, *pp;
//     int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
//     const int       lenpatstrNA = strlen( patstrNA );
    const int       card = GetCardinality();
    int noelems, v, n, c;
    float val;

    if( probs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    noelems = card + (( card * ( card-1 )) >> 1 );

    if( noelems != probs->GetSize())
        throw MYRUNTIME_ERROR( preamb + "Inconsistent probability vector size.");

    fprintf( fp, "%s %d%s", patstrcard, card, NL );

    //write probabilities
    for( v = n = 0; n < card; n++ )
    {
        for( c = 0; c <= n; c++, v++ )
        {
            val = probs->GetValueAt(v);
            if( val == 0.0f )
                fprintf( fp, " %8s", patstrNA );
            else {
                if( val < 5.e-7f )
                    fprintf( fp, " %8d", 0 );
                else
                    fprintf( fp, " %8.6f", val );
            }
        }
        fprintf( fp, "%s", NL );
    }
}
