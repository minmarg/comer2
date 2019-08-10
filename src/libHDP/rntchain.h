/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __rntchain__
#define __rntchain__

#include <stdio.h>
#include <stdlib.h>
#include "liblib/mybase.h"
#include "restaurant.h"

// -------------------------------------------------------------------------
// class RntChain: chain of restaurants 
// (the whole dataset; all groups/restaurants)
//
class RntChain
{
public:
    RntChain( int size );
    RntChain( const RntChain& );
    explicit RntChain();
    ~RntChain();

    RntChain&  operator=( const RntChain& );

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    void        CalcNoVectors();

    Restaurant* GetRestaurantAt( int n ) const;
    void        SetRestaurantAt( int n, Restaurant* value );
    int         NewRestaurant( Restaurant* value );//push restaurant
    void        RemRestaurantAt( int loc, Restaurant* value );//remove restaurant

    void        Copy( const RntChain& vector ); //copy elements
    void        Clear();                        //clear all elements
    void        Print( FILE* ) const;           //print vector

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    void        Reserve( int size );
    void        ReserveRestaurants( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyRestaurants();
    int         GetCapacity() const { return capacity_; }

    Restaurant** GetRestaurants() const { return values_; }

    void        ReallocVacans( int newcap );
    void        DestroyVacans();
    int         GetCapVacans() const { return capvac_; }

    const int*  GetVacancies() const { return vacancies_; }
    int         GetVacantAt( int n ) const;
    void        SetVacantAt( int n, int index );
    void        PushVacant( int index );

    int         GetNoVacans() const { return novacs_; }
    void        SetNoVacans( int value ) { novacs_ = value; }

    void        SetSize( int value ) { length_ = value; }
    void        SetActualSize( int value ) { actlen_ = value; }

private:
    Restaurant** values_;       //restaurants
    int         actlen_;        //actual number of values
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the vector
    bool        destroy_;       //flag of destroying object
    int*        vacancies_;     //vector of indices of vacancies
    int         novacs_;        //number of vacant elements
    int         capvac_;        //capcity of vacancies_
};

// -------------------------------------------------------------------------
// Reserve: Reserve space for vectors
//
inline
void RntChain::Reserve( int size )
{
    ReserveRestaurants( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveRestaurants: Reserve space for vector of restaurants
//
inline
void RntChain::ReserveRestaurants( int size )
{
    if( 0 < size ) {
        Realloc( size );
//         SetSize( size );
    }
}

// -------------------------------------------------------------------------
// ReserveVacans: Reserve space for vector of vacancies
//
inline
void RntChain::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void RntChain::Destroy()
{
    DestroyRestaurants();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyRestaurants: Destroy vector of restaurants
//
inline
void RntChain::DestroyRestaurants()
{
    int n;
    if( GetDestroy()) {
        if( values_ ) {
            for( n = 0; n < length_; n++ )
                if( values_[n] )
                    delete values_[n];
            free( values_ );
        }
    }
    values_ = NULL;
    actlen_ = 0;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// DestroyVacans: Destroy vector of vacancies
//
inline
void RntChain::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}

// -------------------------------------------------------------------------
// CalcNoVectors: calculate number of vectors in each of the restaurants
//
inline
void RntChain::CalcNoVectors()
{
    Restaurant* rest;
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        rest = GetRestaurantAt( n );
        if( rest == NULL )
            continue;
        rest->CalcNoVectors();
    }
}

// -------------------------------------------------------------------------
// GetRestaurantAt: get restaurant at position loc
//
inline
Restaurant* RntChain::GetRestaurantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("RntChain::GetRestaurantAt: Memory access error.");
#endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetRestaurantAt: set restaurant at position loc
//
inline
void RntChain::SetRestaurantAt( int loc, Restaurant* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("RntChain::SetRestaurantAt: Memory access error.");
#endif
    values_[loc] = value;
}

// -------------------------------------------------------------------------
// NewRestaurant: push restaurant
//
inline
int RntChain::NewRestaurant( Restaurant* value )
{
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw MYRUNTIME_ERROR("RntChain::NewRestaurant: Memory access error.");
        values_[loc] = value;
        novacs_--;
        actlen_++;
        return loc;
    }
    if( capacity_ <= ( loc = length_ ))
        Realloc( TIMES2( capacity_ + 1 ));
    values_[loc] = value;
    length_++;
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemRestaurantAt: remove restaurant at position loc
//
inline
void RntChain::RemRestaurantAt( int loc, Restaurant* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("RntChain::RemRestaurantAt: Memory access error.");
#endif
    if( values_[loc] != value )
        throw MYRUNTIME_ERROR("RntChain::RemRestaurantAt: Memory access error.");
    values_[loc] = NULL;
    PushVacant( loc );
    actlen_--;
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int RntChain::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("RntChain::GetVacantAt: Memory access error.");
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void RntChain::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("RntChain::SetVacantAt: Memory access error.");
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index into the vector of vacancies
//
inline
void RntChain::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__rntchain__
