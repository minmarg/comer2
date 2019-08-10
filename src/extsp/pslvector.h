/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_pslvector__
#define __extspsl_pslvector__

#include <stdio.h>
#include <stdlib.h>
#include "pslerror.h"

namespace extspsl {

class Pslmatrix;

class Pslvector
{
public:
    Pslvector( int size );
    Pslvector( const Pslvector& );
    explicit Pslvector();
    ~Pslvector();

    Pslvector&  operator=( const Pslvector& );

    int         GetSize() const { return length_; }
    int         GetCapacity() const { return capacity_; }

    int         AddGNoise( const Pslvector& stds );
    void        SetAllToValue( float value );
    void        SetUnity() { SetAllToValue( 1.0f ); }

    const Pslvector SubVector( int offset, int n ) const;

    float       GetValueAt( int n ) const;          //get value at position n
    void        SetValueAt( int n, float value );   //set value at position n
    void        AddValueAt( int n, float value );   //add value at position n
    void        MulValueAt( int n, float value );   //multiply by value at position n
    int         DivValueAt( int n, float value );   //divide by value at position n

    void        InsertAt( int loc, float value );   //insert value at position loc
    void        Push( float value );                //push value at the end

    void        Copy( const Pslvector& vector );    //copy elements
    void        AssignTo( float value );            //assign all elements to the given value
    void        Zero();                             //assign all elements to zero
    void        Clear();                            //clear all elements
    void        Print( FILE* ) const;               //print vector
    void        Print( FILE*, const char* format ) const;//print vector

    //LINEAR ALGEBRA
    float       Min() const;
    float       Max() const;
    float       Sum() const;//sum of all members
    float       Norm2() const;//norm
    int         DotProduct( const Pslvector& vect2, float* res ) const;//dot product of two vectors
    int         Superposition( float alpha, const Pslvector& vect2 );//addition of two vectors
    int         MultiplyBy( float alpha );//multiplication by scalar
    int         Mul( const Pslmatrix& mt, const Pslvector& v ); //multiplication of matrix by vector
    int         Transpose( Pslmatrix& tr ) const; //transpose of vector
    int         Exp();//exponentiate vector

    void        Allocate( int cap );
    void        Reserve( int size );

    const float*    GetVector() const { return values_; }
    float*          GetVector() { return values_; }
    void            DecDim() { if( !GetMaster()) return; if( length_ ) length_--; }

    void        SetVector( float* vect, int len );

protected:
    void        Realloc( int newcap );
    void        Destroy();

    bool        GetMaster() const { return master_; }
    void        SetMaster() { master_ = true; }
    void        SetSlave() { master_ = false; }

    int         GetStride() const { return stride_; }
    void        SetStride( int value ) { stride_ = value; }

    void        SetSize( int value ) { if( !GetMaster()) return; length_ = value; }

    friend class Pslmatrix;

private:
    float*      values_;        //single-precision values of vector
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the sequence
    int         stride_;        //stride of vector
    bool        master_;        //flag of mastering an object
};

// -------------------------------------------------------------------------
// Allocate: Allocate space for vector
//
inline
void Pslvector::Allocate( int size )
{
    if( !GetMaster())
        PRIVERROR("Pslvector::Allocate: Only master is allowed to allocate memory.");
    if( 0 < size ) {
        SetStride( 1 );
        Realloc( size );
    }
}

// -------------------------------------------------------------------------
// Reserve: Reserve space for vector
//
inline
void Pslvector::Reserve( int size )
{
    if( !GetMaster())
        PRIVERROR("Pslvector::Reserve: Only master is allowed to reserve memory.");
    if( 0 < size ) {
        SetStride( 1 );
        Realloc( size );
        SetSize( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vector
//
inline
void Pslvector::Destroy()
{
    if( master_ )
        if( values_ )
            free( values_ );
    values_ = NULL;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// SetVector: set vector and its size
//
inline
void Pslvector::SetVector( float* vect, int len )
{
    Destroy();
    SetSlave();
    length_ = len;
    values_ = vect;
}

// -------------------------------------------------------------------------
// GetValueAt: get value at position loc
//
inline
float Pslvector::GetValueAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
        PRIVERROR("Pslvector::GetValueAt: Memory access error");
#endif
    if( stride_ == 1 )
        return values_[loc];
    return values_[loc*stride_];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at position loc
//
inline
void Pslvector::SetValueAt( int loc, float value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::SetValueAt: Memory access error");
//         if( length_ + 1 <= loc )
//             PRIVERROR( "Pslvector::SetValueAt: Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::SetValueAt: Memory access error");
#endif
    }
    if( stride_ == 1 )
        values_[loc] = value;
    else
        values_[loc*stride_] = value;
}

// -------------------------------------------------------------------------
// AddValueAt: add value at the position
//
inline
void Pslvector::AddValueAt( int loc, float value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::AddValueAt: Memory access error");
//         if( length_ + 1 <= loc )
//             PRIVERROR( "Pslvector::AddValueAt: Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::AddValueAt: Memory access error");
#endif
    }
    if( stride_ == 1 )
        values_[loc] += value;
    else
        values_[loc*stride_] += value;
}

// -------------------------------------------------------------------------
// MulValueAt: multiply by value at position loc
//
inline
void Pslvector::MulValueAt( int loc, float value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::MulValueAt: Memory access error");
//         if( length_ + 1 <= loc )
//             PRIVERROR( "Pslvector::MulValueAt: Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::MulValueAt: Memory access error");
#endif
    }
    if( stride_ == 1 )
        values_[loc] *= value;
    else
        values_[loc*stride_] *= value;
}

// -------------------------------------------------------------------------
// DivValueAt: divide by value at position loc
//
inline
int Pslvector::DivValueAt( int loc, float value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::DivValueAt: Memory access error");
//         if( length_ + 1 <= loc )
//             PRIVERROR( "Pslvector::DivValueAt: Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            PRIVERROR("Pslvector::DivValueAt: Memory access error");
#endif
    }
    if( value == 0.0 )
        return PSL_ERR_ILLEGAL;

    if( stride_ == 1 )
        values_[loc] /= value;
    else
        values_[loc*stride_] /= value;

    return PSL_SUCCESS;
}

}//namespace extspsl

#endif//__extspsl_pslvector__
