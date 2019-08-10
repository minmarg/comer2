/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basin.h"

using namespace extspsl;

// -------------------------------------------------------------------------
// constructor: initialization
//
Basin::Basin( int size )
:   values_( NULL ),
    processed_( NULL ),
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
Basin::Basin( const Basin& right )
:   values_( NULL ),
    processed_( NULL ),
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
Basin::Basin()
:   values_( NULL ),
    processed_( NULL ),
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
Basin::~Basin()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Basin& Basin::operator=( const Basin& right )
{
    Clear();
    SetActualSize( 0 );

    ReserveValues( right.GetSize());
    ReserveVacans( right.GetNoVacans());
    Copy( right );
    actlen_ = right.GetActualSize();
    destroy_ = false;//right.GetDestroy();

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Basin::Realloc( int newcap )
{
    Pslvector** tmp_values = NULL;
    bool*       tmp_processed = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ <= 0 ) {
        tmp_values = ( Pslvector** )malloc( sizeof( void* ) * newcap );
        tmp_processed = ( bool* )malloc( sizeof( bool ) * newcap );
    } else {
        tmp_values = ( Pslvector** )realloc( values_, sizeof( void* ) * newcap );
        tmp_processed = ( bool* )realloc( processed_, sizeof( bool ) * newcap );
    }

    if( !tmp_values || !tmp_processed )
        throw MYRUNTIME_ERROR("Basin::Realloc: Not enough memory.");

    values_ = tmp_values;
    processed_ = tmp_processed;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
    memset( processed_ + capacity_, 0, sizeof( bool ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Basin::ReallocVacans( int newcap )
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
        throw MYRUNTIME_ERROR("Basin::ReallocVacans: Not enough memory.");

    vacancies_ = tmp_vacans;

    // fill uninitialized memory with binary 1s
    memset( vacancies_ + capvac_, 0xff, sizeof( int ) * ( newcap - capvac_ ));
    capvac_ = newcap;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     when lengths of vectors are not equal
//
void Basin::Copy( const Basin& vector )
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
        SetValueAt( n, vector.GetValueAt( n ));
        SetProcessedAt( n, vector.GetProcessedAt( n ));
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
void Basin::Clear()
{
    int n;
    if( GetVector())
        for( n = 0; n < GetSize(); n++ ) {
            if( GetDestroy() && GetValueAt( n ))
                delete GetValueAt( n );
            SetValueAt( n, NULL );
            SetProcessedAt( n, false );
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
void Basin::Print( FILE* fp ) const
{
    const int   cvpl = 15;
    int         n;
    if( fp == NULL )
        return;
    fprintf( fp, "Vector:%s", NL );
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " %p", GetValueAt( n ));
//         fprintf( fp, " (%p %d)", GetValueAt( n ), GetProcessedAt( n ));
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
