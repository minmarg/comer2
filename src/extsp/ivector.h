/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_ivector__
#define __extspsl_ivector__

#include <stdio.h>
#include <stdlib.h>
#include "pslerror.h"

namespace extspsl {

class Ivector
{
public:
    Ivector( int size );
    Ivector( const Ivector& );
    explicit Ivector();
    ~Ivector();

    Ivector&  operator=( const Ivector& );

    int         GetSize() const { return length_; }
    int         GetCapacity() const { return capacity_; }

    void        SetAllToValue( int value );
    void        SetUnity() { SetAllToValue( 1 ); }

    int         GetLast() const;                    //get last value
    int         GetValueAt( int n ) const;          //get value at position n
    void        SetValueAt( int n, int value );     //set value at position n
    void        AddValueAt( int n, int value );     //add value at position n
    void        MulValueAt( int n, int value );     //multiply by value at position n

    void        InsertAt( int loc, int value );     //insert value at position loc
    void        Push( int value );                  //push value at the end

    void        Copy( const Ivector& vector );      //copy elements
    void        Zero();                             //assign all elements to zero
    void        Clear();                            //clear all elements
    void        Print( FILE* ) const;               //print vector
    void        Print( FILE*, const char* format ) const;//print vector

    int         Min() const;
    int         Max() const;
    int         Sum() const;                        //sum of all members
    int         MultiplyBy( int value );            //multiplication by a scalar

    void        Allocate( int cap );
    void        Reserve( int size );

    const int*  GetVector() const { return values_; }
    int*        GetVector() { return values_; }
    void        DecDim() { if( length_ ) length_--; }

protected:
    void        Realloc( int newcap );
    void        Destroy();

    void        SetSize( int value ) { length_ = value; }

private:
    int*        values_;        //values of vector
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the vector
};


// -------------------------------------------------------------------------
// Allocate: Allocate space for vector
//
inline
void Ivector::Allocate( int size )
{
    if( 0 < size )
        Realloc( size );
}

// -------------------------------------------------------------------------
// Reserve: Reserve space for vector
//
inline
void Ivector::Reserve( int size )
{
    if( 0 < size ) {
        Realloc( size );
        SetSize( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vector
//
inline
void Ivector::Destroy()
{
    if( values_ )
        free( values_ );
    values_ = NULL;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// GetLast: get last value
//
inline
int Ivector::GetLast() const
{
    if( !values_ || length_ < 1 )
        PRIVERROR("Ivector::GetLast: Memory access error.");
    return values_[length_-1];
}

// -------------------------------------------------------------------------
// GetValueAt: get value at position loc
//
inline
int Ivector::GetValueAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        PRIVERROR("Ivector::GetValueAt: Memory access error.");
#endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at position loc
//
inline
void Ivector::SetValueAt( int loc, int value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || capacity_ <= loc )
        PRIVERROR("Ivector::SetValueAt: Memory access error.");
//     if( length_ + 1 <= loc )
//         PRIVERROR("Ivector::SetValueAt: Memory access error.");
#endif
    if( length_ <= loc )
        length_ = loc + 1;
    values_[loc] = value;
}

// -------------------------------------------------------------------------
// AddValueAt: add value at position loc
//
inline
void Ivector::AddValueAt( int loc, int value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || capacity_ <= loc )
        PRIVERROR("Ivector::AddValueAt: Memory access error.");
//     if( length_ + 1 <= loc )
//         PRIVERROR("Ivector::AddValueAt: Memory access error.");
#endif
    if( length_ <= loc )
        length_ = loc + 1;
    values_[loc] += value;
}

// -------------------------------------------------------------------------
// MulValueAt: multiply by value at position loc
//
inline
void Ivector::MulValueAt( int loc, int value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || capacity_ <= loc )
        PRIVERROR("Ivector::MulValueAt: Memory access error.");
//     if( length_ + 1 <= loc )
//         PRIVERROR("Ivector::MulValueAt: Memory access error.");
#endif
    if( length_ <= loc )
        length_ = loc + 1;
    values_[loc] *= value;
}

}//namespace extspsl

#endif//__extspsl_ivector__
