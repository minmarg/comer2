/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "restaurant.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
Restaurant::Restaurant( int size )
:   values_( NULL ),
    novecs_( 0 ),
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
Restaurant::Restaurant( const Restaurant& right )
:   values_( NULL ),
    novecs_( 0 ),
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
Restaurant::Restaurant()
:   values_( NULL ),
    novecs_( 0 ),
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
// -------------------------------------------------------------------------

Restaurant::~Restaurant()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Restaurant& Restaurant::operator=( const Restaurant& right )
{
    Copy( right );
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Restaurant::Realloc( int newcap )
{
    Table** tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ <= 0 ) {
        tmp_values = ( Table** )malloc( sizeof( void* ) * newcap );
    } else {
        tmp_values = ( Table** )realloc( values_, sizeof( void* ) * newcap );
    }

    if( !tmp_values )
        throw MYRUNTIME_ERROR("Restaurant::Realloc: Not enough memory.");

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Restaurant::ReallocVacans( int newcap )
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
        throw MYRUNTIME_ERROR("Restaurant::ReallocVacans: Not enough memory.");

    vacancies_ = tmp_vacans;

    // fill uninitialized memory with -1
    memset( vacancies_ + capvac_, 0xff, sizeof( int ) * ( newcap - capvac_ ));
    capvac_ = newcap;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     where lengths of vectors are not equal
//
void Restaurant::Copy( const Restaurant& vector )
{
    Clear();
    SetActualSize( 0 );

    ReserveTables( vector.GetSize());
    ReserveVacans( vector.GetNoVacans());

    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;
    if( GetCapVacans() < vector.GetNoVacans())
        return;

    if( GetSize() < vector.GetSize())
        SetSize( vector.GetSize());

    Table* tbl = NULL;
    int noelms = GetSize();
    int n;

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( n = 0; n < noelms; n++ ) {
        tbl = NULL;
        if( vector.GetTableAt( n ))
            tbl = new Table( *vector.GetTableAt( n ));
        SetTableAt( n, tbl );
    }

    SetNoVectors( vector.GetNoVectors());

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
void Restaurant::Clear()
{
    int n;
    if( GetTables())
        for( n = 0; n < GetSize(); n++ ) {
            if( GetDestroy() && GetTableAt( n ))
                delete GetTableAt( n );
            SetTableAt( n, NULL );
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
void Restaurant::Print( FILE* fp ) const
{
    const int   cvpl = 15;
    int         n;
    if( fp == NULL )
        return;
    fprintf( fp, "Vector:%s", NL );
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " %p", GetTableAt( n ));
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
