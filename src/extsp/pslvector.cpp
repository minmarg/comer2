/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rv/rvnorm.h"
#include "psl.h"
#include "pslerror.h"
#include "pslmatrix.h"
#include "pslvector.h"

namespace extspsl {

// -------------------------------------------------------------------------
// constructor: initialization
//
Pslvector::Pslvector( int size )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    stride_( 1 ),
    master_( true )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: copy
//
Pslvector::Pslvector( const Pslvector& right )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    stride_( 1 ),
    master_( true )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
Pslvector::Pslvector()
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    stride_( 1 ),
    master_( true )
{
}

// -------------------------------------------------------------------------
// destructor:
//
Pslvector::~Pslvector()
{
    Destroy();
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Pslvector& Pslvector::operator=( const Pslvector& right )
{
    Destroy();
    stride_ = 1;
    master_ = true;

    if( right.GetMaster()) {
        stride_ = right.GetStride();
        Reserve( right.GetSize());
        Copy( right );
    }
    else {
        master_ = false;
        values_ = right.values_;
        length_ = right.length_;
        stride_ = right.stride_;
    }
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Pslvector::Realloc( int newcap )
{
    if( !GetMaster())
        PRIVERROR("Pslvector::Realloc: Only master is allowed to manage memory.");

    float* tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = (float*)malloc( sizeof(float) * newcap );
    } else {
        tmp_values = (float*)realloc( values_, sizeof(float) * newcap );
    }

    if( !tmp_values )
        PRIVERROR("Pslvector::Realloc: Not enough memory");

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof(float) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: insert value at the end of vector
//
void Pslvector::Push( float value )
{
    if( !GetMaster())
        PRIVERROR("Pslvector::Push: Only master is allowed to push data.");

    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    values_[length_] = value;

    length_++;
}

// -------------------------------------------------------------------------
// InsertAt: insert value by shifting elements from the position
//     to the right
//
void Pslvector::InsertAt( int loc, float value )
{
    if( !GetMaster())
        PRIVERROR("Pslvector::InsertAt: Only master is allowed to insert data.");

    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    if( loc < 0 || length_ < loc )
        PRIVERROR("Pslvector::InsertAt: Unable to insert value");

    for( int n = length_; n > loc; n-- )
        values_[n] = values_[n-1];

    values_[loc] = value;

    length_++;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manage cases
//     where lengths of vectors are not equal
//
void Pslvector::Copy( const Pslvector& vector )
{
    if( vector.GetSize() <= 0 )
        return;

    if( GetMaster() && GetCapacity() < vector.GetSize())
        return;

    if( GetSize() < vector.GetSize()) {
        if( GetMaster())
            SetSize( vector.GetSize());
        else
            PRIVERROR("Pslvector::Copy: Slave has insufficient space.");
    }

    int noelms = GetSize();

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( int n = 0; n < noelms; n++ )
        SetValueAt( n, vector.GetValueAt( n ));
//     memcpy( values_, vector.GetVector(), sizeof(float) * noelms );
}

// -------------------------------------------------------------------------
// AssignTo: assign all elements to the given value
//
void Pslvector::AssignTo( float value )
{
    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, value );
}

// -------------------------------------------------------------------------
// Zero: assign all elements to zero
//
void Pslvector::Zero()
{
    AssignTo( 0.0f );
//     for( int n = 0; n < GetSize(); n++ )
//         SetValueAt( n, 0.0f );
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void Pslvector::Clear()
{
    if( !GetMaster())
        PRIVERROR("Pslvector::Clear: Only master is allowed to clear data.");

    if( GetVector())
        for( int n = 0; n < GetSize(); n++ )
            SetValueAt( n, 0.0f );
//         memset( values_, 0, sizeof(float) * capacity_ );
    SetSize( 0 );
}

// -------------------------------------------------------------------------
// Print: print vector to file
//
void Pslvector::Print( FILE* fp ) const
{
    const int cvpl = 10;
    int n;
    if( fp == NULL )
        return;
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " %g", GetValueAt( n ));
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
void Pslvector::Print( FILE* fp, const char* format ) const
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
// AddGNoise: add Gaussian noise with std. deviations in `stds'
//
int Pslvector::AddGNoise( const Pslvector& stds )
{
    if( GetSize() != stds.GetSize())
        return PSL_ERR_DIM;

    time_t  s_cnt;

    MTRng   rng;
    RVNorm  nrv( rng, RVNorm::TRVN_Ratio );
    float   val;
    int n, err;

    time( &s_cnt );
    rng.Set((unsigned int)(size_t)this + (unsigned int)s_cnt );
    for( n = 0; n < GetSize(); n++ ) {
        nrv.SetStd( stds.GetValueAt( n ));
        err = nrv.Gen( &val );
        if( err != 0 )
            return err;
        AddValueAt( n, val );
    }
    return 0;
}

// -------------------------------------------------------------------------
// SetAllToValue: set all elements to equal to `value'
//
void Pslvector::SetAllToValue( float value )
{
    int n;
    for( n = 0; n < GetSize(); n++ )
        SetValueAt( n, value );
}

// -------------------------------------------------------------------------
// SubVector: make a new vector to correspond to a subvector
//
const Pslvector Pslvector::SubVector( int offset, int n ) const
{
    Pslvector   sub;

    if( n < 1 || offset < 0 )
        PRIVERROR( "Pslvector::SubVector: invalid number of elements." );
    if( GetSize() < offset + n )
        PRIVERROR( "Pslvector::SubVector: invalid offset." );

    sub.SetVector( values_ + offset*stride_, n );
    sub.SetStride( stride_ );
    return sub;
}

// =========================================================================
// LINEAR ALGEBRA
//
// Min: find minimum value
//
float Pslvector::Min() const
{
    float min = 0.0f;
    float val = 0.0f;
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
float Pslvector::Max() const
{
    float max = 0.0f;
    float val = 0.0f;
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
float Pslvector::Sum() const
{
    float sum = 0.0f;
    float x;
    int n;

    for( n = 0; n < GetSize(); n++ ) {
        x = GetValueAt( n );

        if( x == 0.0f )
            continue;

        sum += x;
    }

    return sum;
}

// =========================================================================
// Norm2: euclidean norm
//
float Pslvector::Norm2() const
{
    float scale = 0.0f;
    float ssq = 1.0f;
    float x, ax;
    int n;

    if( GetSize() <= 0 )
        return 0.0f;

    if( GetSize() == 1 )
        return fabsf( GetValueAt( 0 ));

    for( n = 0; n < GetSize(); n++ ) {
        x = GetValueAt( n );

        if( x == 0.0f )
            continue;

        ax = fabsf( x );
        if( scale < ax ) {
            ssq = 1.0f + ssq * ( scale / ax ) * ( scale / ax );
            scale = ax;
        } else {
            ssq += ( ax / scale ) * ( ax / scale );
        }
    }

    return scale * sqrtf( ssq );
}

// =========================================================================
// DotProduct: calculate dot product
//
int Pslvector::DotProduct( const Pslvector& vect2, float* res ) const
{
    if( GetSize() != vect2.GetSize())
        return PSL_ERR_DIM;

    if( res == NULL )
        return PSL_ERR_ADDRESS;

    *res = 0.0f;
    for( int n = 0; n < GetSize(); n++ )
        *res += GetValueAt( n ) * vect2.GetValueAt( n );

    return 0;
}

// =========================================================================
// Superposition: compute vector . scalar product and add the result to this
//     vector
//
int Pslvector::Superposition( float alpha, const Pslvector& vect2 )
{
    if( GetSize() != vect2.GetSize())
        return PSL_ERR_DIM;

    if( alpha == 0.0f )
        return 0;

    int n;

    if( alpha == 1.0f )
        for( n = 0; n < GetSize(); n++ )
            AddValueAt( n, vect2.GetValueAt( n ));
    else if( alpha == -1.0f )
        for( n = 0; n < GetSize(); n++ )
            AddValueAt( n, -vect2.GetValueAt( n ));
    else
        for( n = 0; n < GetSize(); n++ )
            AddValueAt( n, alpha * vect2.GetValueAt( n ));
    return 0;
}

// =========================================================================
// MultiplyBy: multibly vector by a scalar
//
int Pslvector::MultiplyBy( float alpha )
{
    if( alpha == 1.0f )
        return 0;

    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, alpha * GetValueAt( n ));
    return 0;
}

// =========================================================================
// Mul: mutliply matrix by vector and save the result
//
int Pslvector::Mul( const Pslmatrix& mt, const Pslvector& v )
{
    int mtrows = mt.GetNoRows();
    int mtcols = mt.GetNoCols();
    int vsize = v.GetSize();
    int n, k;
    float sum;
    if( mtcols != vsize )
        return PSL_ERR_DIM;

    Reserve( mtrows );
    for( n = 0; n < mtrows; n++ ) {
        sum = 0.0f;
        for( k = 0; k < mtcols; k++ )
            sum += mt.GetValueAt( n, k ) * v.GetValueAt( k );
        SetValueAt( n, sum );
    }
    return PSL_SUCCESS;
}

// =========================================================================
// Transpose: apply the transpose operator to the vector
//
int Pslvector::Transpose( Pslmatrix& tr ) const
{
    int size = GetSize();
    int n;
    if( size < 1 )
        return PSL_SUCCESS;

    tr.Reserve( 1, size );
    for( n = 0; n < size; n++ )
        tr.SetValueAt( 0, n, GetValueAt( n ));
    return PSL_SUCCESS;
}

// =========================================================================
// Exp: exponentiate vector
//
int Pslvector::Exp()
{
    float value;
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        value = GetValueAt( n );
        if( !isfinite( value ))
            return PSL_ERR_DOMAIN;
        if( SLC_LOG_SP_MAX < value  )
            return PSL_ERR_DOMAIN;
        if( value < SLC_LOG_SP_MIN )
            value = 0.0f;
        else
            value = expf( value );
        SetValueAt( n, value );
    }
    return 0;
}

}//namespace extspsl
