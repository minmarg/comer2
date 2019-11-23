/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dish.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
Dish::Dish( int size )
:   time_(( time_t )0 ),
    basin_( NULL ),
    values_( NULL ),
    processed_( NULL ),
    dishproced_( false ),
    prosize_( -1 ),
    pronotbls_( -1 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 ),
    notables_( 0 ),
    tmpval_( 0.0f )
{
    Reserve( size );
    time( &time_ );
}

// -------------------------------------------------------------------------
// constructor: copy
//
Dish::Dish( const Dish& right )
:   time_(( time_t )0 ),
    basin_( NULL ),
    values_( NULL ),
    processed_( NULL ),
    dishproced_( false ),
    prosize_( -1 ),
    pronotbls_( -1 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 ),
    notables_( 0 ),
    tmpval_( 0.0f )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
Dish::Dish()
:   time_(( time_t )0 ),
    basin_( NULL ),
    values_( NULL ),
    processed_( NULL ),
    dishproced_( false ),
    prosize_( -1 ),
    pronotbls_( -1 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( false ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 ),
    notables_( 0 ),
    tmpval_( 0.0f )
{
    time( &time_ );
}

// -------------------------------------------------------------------------
// destructor:
//
Dish::~Dish()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Dish& Dish::operator=( const Dish& right )
{
    Copy( right );
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Dish::Realloc( int newcap )
{
    int*    tmp_values = NULL;
    bool*   tmp_processed = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ <= 0 ) {
        tmp_values = ( int* )malloc( sizeof( int ) * newcap );
        tmp_processed = ( bool* )malloc( sizeof( bool ) * newcap );
    } else {
        tmp_values = ( int* )realloc( values_, sizeof( int ) * newcap );
        tmp_processed = ( bool* )realloc( processed_, sizeof( bool ) * newcap );
    }

    if( !tmp_values || !tmp_processed )
        throw MYRUNTIME_ERROR("Dish::Realloc: Not enough memory.");

    values_ = tmp_values;
    processed_ = tmp_processed;

    // fill array elements with -1
    memset( values_ + capacity_, 0xff, sizeof( int ) * ( newcap - capacity_ ));
    memset( processed_ + capacity_, 0, sizeof( bool ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Dish::ReallocVacans( int newcap )
{
    int*    tmp_vacans = NULL;

    if( newcap <= capvac_ )
        return;

    if( capvac_ == 0 ) {
        tmp_vacans = ( int* )malloc( sizeof( int ) * newcap );
    } else {
        tmp_vacans = ( int* )realloc( vacancies_, sizeof( int ) * newcap );
    }

    if( !tmp_vacans )
        throw MYRUNTIME_ERROR("Dish::ReallocVacans: Not enough memory.");

    vacancies_ = tmp_vacans;

    // fill array elements with -1
    memset( vacancies_ + capvac_, 0xff, sizeof( int ) * ( newcap - capvac_ ));
    capvac_ = newcap;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     when lengths of vectors are not equal
//
void Dish::Copy( const Dish& vector )
{
    Clear();
    SetActualSize( 0 );

    ReserveValues( vector.GetSize());
    ReserveVacans( vector.GetNoVacans());

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
        SetVectorNIndAt( n, vector.GetVectorNIndAt( n ));
        SetProcessedAt( n, vector.GetProcessedAt( n ));
    }

    SetTime( vector.GetTime());
    SetBasin( vector.GetBasin());
    SetNoTables( vector.GetNoTables());
    SetTmpValue( vector.GetTmpValue());
    SetProcessed( vector.GetProcessed());
    SetReadSize( vector.GetReadSize());

    if( GetNoVacans() < vector.GetNoVacans())
        SetNoVacans( vector.GetNoVacans());

    noelms = GetNoVacans();

    if( vector.GetNoVacans() < noelms )
        noelms = vector.GetNoVacans();

    for( n = 0; n < noelms; n++ )
        SetVacantAt( n, vector.GetVacantAt( n ));

    SetActualSize( vector.GetActualSize());
    SetDestroy( true );
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void Dish::Clear()
{
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        SetVectorNIndAt( n, -1 );
        SetProcessedAt( n, false );
    }
    if( GetVacancies())
        for( n = 0; n < GetNoVacans(); n++ )
            SetVacantAt( n, -1 );
    SetSize( 0 );
    SetReadSize( -1 );
    SetActualSize( 0 );
    SetNoVacans( 0 );
    SetNoTables( 0 );
    SetProcessed( false );
}

// -------------------------------------------------------------------------
// Print: print vector to file
//
void Dish::Print( FILE* fp ) const
{
    const int   cvpl = 15;
    int         n;
    if( fp == NULL )
        return;
    fprintf( fp, "Vector:%s", NL );
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " (%d %d)", GetVectorNIndAt( n ), GetProcessedAt( n ));
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
