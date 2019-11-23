/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pslerror.h"
#include "ivector.h"

namespace extspsl {

// -------------------------------------------------------------------------
// constructor: initialization
//
Ivector::Ivector( int size )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: copy
//
Ivector::Ivector( const Ivector& right )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
Ivector::Ivector()
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
Ivector::~Ivector()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Ivector& Ivector::operator=( const Ivector& right )
{
    Destroy();
    Reserve( right.GetSize());
    Copy( right );
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Ivector::Realloc( int newcap )
{
    int* tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = ( int* )malloc( sizeof( int ) * newcap );
    } else {
        tmp_values = ( int* )realloc( values_, sizeof( int ) * newcap );
    }

    if( !tmp_values )
        PRIVERROR("Ivector::Realloc: Not enough memory.");

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( int ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: insert value at the end of vector
//
void Ivector::Push( int value )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    values_[length_] = value;

    length_++;
}

// -------------------------------------------------------------------------
// InsertAt: insert value by shifting elements from position loc 
//     to the right
//
void Ivector::InsertAt( int loc, int value )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    if( loc < 0 || length_ < loc )
        PRIVERROR("Ivector::InsertAt: Unable to insert value.");

    for( int n = length_; n > loc; n-- )
        values_[n] = values_[n-1];

    values_[loc] = value;

    length_++;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     when lengths of vectors are not equal
//
void Ivector::Copy( const Ivector& vector )
{
    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;

    if( GetSize() < vector.GetSize()) {
        SetSize( vector.GetSize());
    }

    int noelms = GetSize();

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( int n = 0; n < noelms; n++ )
        SetValueAt( n, vector.GetValueAt( n ));
//     memcpy( values_, vector.GetVector(), sizeof( int ) * noelms );
}

// -------------------------------------------------------------------------
// Zero: assign all elements to zero
//
void Ivector::Zero()
{
    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, 0 );
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void Ivector::Clear()
{
    if( GetVector())
        for( int n = 0; n < GetSize(); n++ )
            SetValueAt( n, 0 );
//         memset( values_, 0, sizeof( int ) * capacity_ );
    SetSize( 0 );
}

// -------------------------------------------------------------------------
// Print: print vector to file
//
void Ivector::Print( FILE* fp ) const
{
    const int   cvpl = 10;
    int         n;
    if( fp == NULL )
        return;
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " %d", GetValueAt( n ));
        if(( n + 1 ) % cvpl == 0 ) {
            fprintf( fp, "\n");
            if( n + 1 < GetSize())
                fprintf( fp, "    ");
        }
    }
    if( n % cvpl )
        fprintf( fp, "\n");
}

// -------------------------------------------------------------------------
// Print: print vector to file
//
void Ivector::Print( FILE* fp, const char* format ) const
{
    if( fp == NULL )
        return;
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, format, GetValueAt( n ));
    }
    fprintf( fp, "\n");
}

// -------------------------------------------------------------------------
// SetAllToValue: set all elements to equal to `value'
//
void Ivector::SetAllToValue( int value )
{
    int n;
    for( n = 0; n < GetSize(); n++ )
        SetValueAt( n, value );
}

// =========================================================================
// Min: find minimum value
//
int Ivector::Min() const
{
    int min = 0;
    int val = 0;
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        val = GetValueAt( n );
        if( !n ) {
            min = val;
            continue;
        }
        if( val < min )
            min = val;
    }
    return min;
}

// Max: find maximum value
//
int Ivector::Max() const
{
    int max = 0;
    int val = 0;
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        val = GetValueAt( n );
        if( !n ) {
            max = val;
            continue;
        }
        if( max < val )
            max = val;
    }
    return max;
}

// -------------------------------------------------------------------------
// Sum: Sum all vector elements and return the result
//
int Ivector::Sum() const
{
    int sum = 0;
    int x;
    int n;

    for( n = 0; n < GetSize(); n++ ) {
        x = GetValueAt( n );

        if( x == 0 )
            continue;

        sum += x;
    }

    return sum;
}

// =========================================================================
// MultiplyBy: multibly vector by a scalar
//
int Ivector::MultiplyBy( int value )
{
    if( value == 1 )
        return 0;

    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, value * GetValueAt( n ));
    return 0;
}

}//namespace extspsl
