/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "root.h"
#include "CtxtCoefficients.h"

// number of allocated positions by default
const size_t cAllocPoss = 100;

////////////////////////////////////////////////////////////////////////////
// CLASS CtxtCoefficients
//
// Default constructor
//

CtxtCoefficients::CtxtCoefficients( size_t size )
:   length_( 0 ),
    cweight_( 0.0f ),
    multip_( 0.0f ),
    coeffs_( NULL ),
    logcoeffs_( NULL )
{
    Init( size );
}

// Constructor
//

CtxtCoefficients::CtxtCoefficients( size_t size, float weight )
:   length_( 0 ),
    cweight_( 0.0f ),
    multip_( 0.0f ),
    coeffs_( NULL ),
    logcoeffs_( NULL )
{
    Init( size );
    SetCentralWeight( weight );
}

// Default constructor
//
CtxtCoefficients::CtxtCoefficients()
{
    throw MYRUNTIME_ERROR("CtxtCoefficients::CtxtCoefficients: Default construction forbidden.");
}

// -------------------------------------------------------------------------
// Destructor
//
CtxtCoefficients::~CtxtCoefficients()
{
    Destroy();
}

// -------------------------------------------------------------------------
// Init: allocate and initialize memory
//
void CtxtCoefficients::Init( size_t size )
{
    Destroy();

    coeffs_ = (float*)malloc( size * sizeof(float));
    logcoeffs_ = (float*)malloc( size * sizeof(float));

    if( !coeffs_ || !logcoeffs_ )
        throw MYRUNTIME_ERROR("CtxtCoefficients::Init: Not enough memory.");

    memset( coeffs_, 0, size * sizeof(float));
    memset( logcoeffs_, 0, size * sizeof(float));

    SetLength( size );
}

// -------------------------------------------------------------------------
// Destroy: destroy data
//
void CtxtCoefficients::Destroy()
{
    if( coeffs_ ) {
        free( coeffs_ );
        coeffs_ = NULL;
    }
    if( logcoeffs_ ) {
        free( logcoeffs_ );
        logcoeffs_ = NULL;
    }
}

// -------------------------------------------------------------------------
// fdffunction: function and its derivative to find the root of multiplier
//     f = 2w SUM x^i + w - 1
//    df = 2w SUM ix^(i-1), w is central weight
//
void CtxtCoefficients::fdfunction( float x, float* f, float* df, void* params )
{
    if( !params )
        return;

    float w = ((float*)params)[0];
    float w2 = w + w;
    float ldf = 0.0f;
    float lf = 1.0f;
    size_t length = (size_t)((float*)params)[1];
    size_t nd2 = length / 2;
    size_t n;

    for( n = 0; n < nd2; n++ ) {
        ldf = ldf * x + ( nd2 - n );
        lf  = lf  * x + 1.0f;
    }
    ldf = w2 * ldf;
    lf = w2 * ( lf - 1.0f ) + w - 1.0f;

    if( df ) *df = ldf;
    if( f )  *f  = lf;
    return;
}

// -------------------------------------------------------------------------
// ExpandCoefficients: expand coefficients for each position
//
float CtxtCoefficients::ExpandCoefficients()
{
    float cweight = GetCentralWeight();
    float multipl = GetMultiplier();
    float value = cweight;
    float sum = cweight;
    size_t length = GetLength();
    size_t nd2 = length / 2;
    size_t n;

    if(!( length % 2 ))
        return 0.0f;

    SetCoefficientAt( nd2, value );
    SetLogCoefficientAt( nd2 );

    for( n = 1; n <= nd2; n++ ) {
        value *= multipl;
        SetCoefficientAt( nd2 + n, value );
        SetCoefficientAt( nd2 - n, value );
        SetLogCoefficientAt( nd2 + n );//log
        SetLogCoefficientAt( nd2 - n );//log
        sum += value + value;
    }

    return sum;
}

// -------------------------------------------------------------------------
// FindCoefficients: find missing coefficients given central weight
//
void CtxtCoefficients::FindCoefficients()
{
    float limit = 0.001f;
    float accuracy = limit * 0.1f;
    float x1 = 0.0f + accuracy;
    float x2 = 5.0f;
    float lolim = 0.0f + limit;
    float uplim = 1.0f - limit;
    float cweight = GetCentralWeight();
    float params[2] = { cweight, (float)GetLength() };
    float multiplier = 0.0f;
    float consv = 0.0f;
    const char* emsg = NULL;
    char strbuf[BUF_MAX];
    mystring errmsg = "CtxtCoefficients::FindCoefficients: ";

    if(!( GetLength() % 2 ))
        throw MYRUNTIME_ERROR( errmsg + "Even number of coefficients.");

    if( cweight != 1.0f ) {
        if( cweight <= lolim || uplim <= cweight )
            throw MYRUNTIME_ERROR( errmsg + "Invalid central weight.");

        emsg = root_by_NR_and_bisection(
            &CtxtCoefficients::fdfunction,
            x1, x2, accuracy, MAXIT, params, &multiplier
        );
    }
    if( emsg != NULL )
        throw MYRUNTIME_ERROR( errmsg + emsg );

    if( multiplier < 0.0f )
        throw MYRUNTIME_ERROR( errmsg + "Failed to find the root of multiplier." );

    SetMultiplier( multiplier );

    consv = ExpandCoefficients();
    if( consv <= 1.0f - accuracy || 1.0f + accuracy <= consv ) {
        sprintf( strbuf, "Invalid coefficients: sum= %g.", consv );
        throw MYRUNTIME_ERROR( errmsg + strbuf );
    }
    return;
}

// =========================================================================
// Write: write coefficients to file
//
void CtxtCoefficients::Write( FILE* fp ) const
{
    size_t  n;
    if( fp == NULL )
        return;

    for( n = 0; n < GetLength(); n++ )
        fprintf( fp, " %8g", GetCoefficientAt( n ));
    fprintf( fp, "\n" );
}

// Write: write the logs of coefficients to file
//
void CtxtCoefficients::WriteLogs( FILE* fp ) const
{
    size_t  n;
    if( fp == NULL )
        return;

    for( n = 0; n < GetLength(); n++ )
        fprintf( fp, " %8g", GetLogCoefficientAt( n ));
    fprintf( fp, "\n" );
}
