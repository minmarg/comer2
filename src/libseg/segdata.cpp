/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdlib.h>
#include <math.h>

#include "extsp/psl.h"
#include "liblib/mybase.h"
#include "segdata.h"

#define MAXLOGFACT  10000

namespace SEG {

// Variable definitions

const size_t s_SZPRSQR =  101;
const size_t s_SZPRLOG =  200;
const size_t s_SZPRENT =  200;
const size_t s_SZPRLFT = 1000;

const _TSQUARES             PRESQUARES;
const _TLOGARITHMS          LOGARITHMS;
const _TPARTIAL_ENTROPIES   PRT_ENTROPIES;
_TLOG_FACT                  LOG_FACT;

// -------------------------------------------------------------------------
// _TPRECOMPUTED: precomputed values
//
_TPRECOMPUTED::_TPRECOMPUTED( size_t size )
:   no_vals_( size ),
    data_( NULL )
{
    data_ = (float*)malloc( sizeof(float) * no_vals_ );

    if( data_ == NULL )
        throw MYRUNTIME_ERROR( "_TPRECOMPUTED::_TPRECOMPUTED: Not enough memory." );
}

_TPRECOMPUTED::~_TPRECOMPUTED()
{
    if( data_ )
        free( data_ );
}

void _TPRECOMPUTED::Init()
{
    if( data_ == NULL )
        throw MYRUNTIME_ERROR( "_TPRECOMPUTED::Init: Memory access error." );

    for( size_t v = 0; v < no_vals_; v++ )
        data_[v] = Compute( v );
}

// -------------------------------------------------------------------------
// _TLOGARITHMS: precomputed squares
//
_TSQUARES::_TSQUARES( size_t size )
:   _TPRECOMPUTED( size )
{
    Init();
}

_TSQUARES::~_TSQUARES()
{
}

float _TSQUARES::Compute( size_t v ) const
{
    return SQUARE((float)v );
}

// -------------------------------------------------------------------------
// _TLOGARITHMS: precomputed logarithms
//
_TLOGARITHMS::_TLOGARITHMS( size_t size )
:   _TPRECOMPUTED( size )
{
    Init();
}

_TLOGARITHMS::~_TLOGARITHMS()
{
}

float _TLOGARITHMS::Compute( size_t v ) const
{
    if( v == 0 )
        return -9999.0f;

    return logf((float)v );
}

// -------------------------------------------------------------------------
// _TPARTIAL_ENTROPIES: precomputed entropy values
//
_TPARTIAL_ENTROPIES::_TPARTIAL_ENTROPIES( size_t size )
:   _TPRECOMPUTED( size )
{
    Init();
}

_TPARTIAL_ENTROPIES::~_TPARTIAL_ENTROPIES()
{
}

float _TPARTIAL_ENTROPIES::Compute( size_t v ) const
{
    if( v == 0 )
        return 0.0f;

    return -(float)v * logf((float)v) / SLC_LN2;
}

// -------------------------------------------------------------------------
// _TLOG_FACT: precomputed log-factorial values
//
_TLOG_FACT::_TLOG_FACT()
:   no_vals_( 0 ),
    data_( NULL )
{
    Precompute( s_SZPRLFT );
}

_TLOG_FACT::~_TLOG_FACT()
{
    if( data_ )
        free( data_ );
}

void _TLOG_FACT::Precompute( size_t newsize )
{
    size_t currentsize = GetSize();

    if( newsize <= currentsize )
        return;

    if( data_ == NULL )
        data_ = (float*)malloc( sizeof(float) * newsize );
    else
        data_ = (float*)realloc( data_, sizeof(float) * newsize );

    if( data_ == NULL )
        throw MYRUNTIME_ERROR( "_TLOG_FACT::Precompute: Not enough memory." );

    float   gammap1 = 0.0f;
    size_t  portion = 50;
    size_t  zerosps = 2;
    size_t  beginwith = ( currentsize < zerosps )? zerosps: currentsize;

    for( size_t v = currentsize; v < zerosps && v < no_vals_; v++ )
        data_[v] = 0.0f;

    for( size_t v = beginwith, pv = v; v < newsize; v = pv ) {
        gammap1 = 1.0f;

        //compute log-factorials in steps to avoid overflow
        for( pv = v; pv < v + portion && pv < newsize; pv++ ) {
            gammap1 *= (float)pv;
            data_[pv] = logf(gammap1) + data_[v-1];
        }
    }

    SetSize( newsize );
}

float _TLOG_FACT::GetValueOf( size_t value )
{
    if( GetSize() == 0 )
//         return 0.0f;
        throw MYRUNTIME_ERROR( "_TLOG_FACT::GetValueOf: Memory unallocated." );

    if( value < GetSize())
        return data_[value];

    if( value < GetSize() + GetSize() && value < MAXLOGFACT ) {
        Precompute( GetSize() + GetSize());
        return data_[value];
    }

    float lngammap1 = 0.0f;

    for( ; value >= GetSize(); value-- )
        lngammap1 += logf((float)value);

    return lngammap1 + data_[GetSize()-1];
}

}//namespace SEG
