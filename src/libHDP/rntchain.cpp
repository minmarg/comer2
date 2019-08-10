/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rntchain.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
RntChain::RntChain( int size )
:   values_( NULL ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: copy
//
RntChain::RntChain( const RntChain& right )
:   values_( NULL ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
RntChain::RntChain()
:   values_( NULL ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
RntChain::~RntChain()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
RntChain& RntChain::operator=( const RntChain& right )
{
    Clear();
    SetActualSize( 0 );

    ReserveRestaurants( right.GetSize());
    ReserveVacans( right.GetNoVacans());
    Copy( right );
    SetActualSize( right.GetActualSize());
    SetDestroy( false/*right.GetDestroy()*/);

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void RntChain::Realloc( int newcap )
{
    Restaurant** tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ <= 0 ) {
        tmp_values = ( Restaurant** )malloc( sizeof( void* ) * newcap );
    } else {
        tmp_values = ( Restaurant** )realloc( values_, sizeof( void* ) * newcap );
    }

    if( !tmp_values )
        throw MYRUNTIME_ERROR("RntChain::Realloc: Not enough memory.");

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void RntChain::ReallocVacans( int newcap )
{
    int* tmp_vacans = NULL;

    if( newcap <= capvac_ )
        return;

    if( capvac_ == 0 ) {
        tmp_vacans = ( int* )malloc( sizeof( int ) * newcap );
    } else {
        tmp_vacans = ( int* )realloc( vacancies_, sizeof( int ) * newcap );
    }

    if( !tmp_vacans )
        throw MYRUNTIME_ERROR("RntChain::ReallocVacans: Not enough memory.");

    vacancies_ = tmp_vacans;

    // fill uninitialized memory with -1
    memset( vacancies_ + capvac_, 0xff, sizeof( int ) * ( newcap - capvac_ ));
    capvac_ = newcap;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     where lengths of vectors are not equal
//
void RntChain::Copy( const RntChain& vector )
{
    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;
    if( GetCapVacans() < vector.GetNoVacans())
        return;

    if( GetSize() < vector.GetSize())
        SetSize( vector.GetSize());

    int noelms = GetSize();
    int n;

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( n = 0; n < noelms; n++ ) {
        SetRestaurantAt( n, vector.GetRestaurantAt( n ));
    }

    if( GetNoVacans() < vector.GetNoVacans())
        SetNoVacans( vector.GetNoVacans());

    noelms = GetNoVacans();

    if( vector.GetNoVacans() < noelms )
        noelms = vector.GetNoVacans();

    for( n = 0; n < noelms; n++ )
        SetVacantAt( n, vector.GetVacantAt( n ));
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void RntChain::Clear()
{
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        if( GetDestroy() && GetRestaurantAt( n ))
            delete GetRestaurantAt( n );
        SetRestaurantAt( n, NULL );
    }
    if( GetVacancies())
        for( n = 0; n < GetNoVacans(); n++ )
            SetVacantAt( n, -1 );
    SetSize( 0 );
    SetNoVacans( 0 );
}

// -------------------------------------------------------------------------
// Print: print vector to file
//
void RntChain::Print( FILE* fp ) const
{
    const int cvpl = 10;
    int n;
    if( fp == NULL )
        return;
    fprintf( fp, "Vector:%s", NL );
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " (%p)", GetRestaurantAt( n ));
        if(( n + 1 ) % cvpl == 0 ) {
            fprintf( fp, "%s", NL );
            if( n + 1 < GetSize())
                fprintf( fp, "    ");
        }
    }
    if( n % cvpl )
        fprintf( fp, "%s", NL );
    fprintf( fp, "Vacancies:%s", NL );
    for( n = 0; n < GetNoVacans(); n++ ) {
        fprintf( fp, " %d", GetVacantAt( n ));
        if(( n + 1 ) % cvpl == 0 ) {
            fprintf( fp, "%s", NL );
            if( n + 1 < GetNoVacans())
                fprintf( fp, "    ");
        }
    }
    if( n % cvpl )
        fprintf( fp, "%s", NL );
}

// =========================================================================
