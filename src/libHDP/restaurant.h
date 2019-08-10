/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __restaurant__
#define __restaurant__

#include <stdio.h>
#include <stdlib.h>
#include "liblib/mybase.h"
#include "table.h"

// -------------------------------------------------------------------------
// class Restaurant: reusable vector of tables (in one restaurant, or group!)
//
class Restaurant
{
public:
    Restaurant( int size );
    Restaurant( const Restaurant& );
    explicit Restaurant();
    ~Restaurant();

    Restaurant&  operator=( const Restaurant& );

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    int         GetNoVectors() const { return novecs_; }
    void        SetNoVectors( int value ) { novecs_ = value; }
    void        CalcNoVectors();

    Table*      GetTableAt( int n ) const;
    void        SetTableAt( int n, Table* value );
    int         NewTable( Table* value );//push table
    void        RemTableAt( int loc, Table* value );//remove table

    void        Copy( const Restaurant& vector );//copy elements
    void        Clear();                        //clear all elements
    void        Print( FILE* ) const;           //print vector

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    void        Reserve( int size );
    void        ReserveTables( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyTables();
    int         GetCapacity() const { return capacity_; }

    Table**     GetTables() const { return values_; }

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
    Table**     values_;        //tables
    int         novecs_;        //number of vectors
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
void Restaurant::Reserve( int size )
{
    ReserveTables( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveTables: Reserve space for vector of tables
//
inline
void Restaurant::ReserveTables( int size )
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
void Restaurant::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void Restaurant::Destroy()
{
    DestroyTables();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyTables: Destroy vector of tables
//
inline
void Restaurant::DestroyTables()
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
void Restaurant::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}

// -------------------------------------------------------------------------
// CalcNoVectors: calculate number of vectors in restaurant over all tables
//
inline
void Restaurant::CalcNoVectors()
{
    Table*  tbl;
    int     nvecs = 0;
    int     n;
    for( n = 0; n < GetSize(); n++ ) {
        tbl = GetTableAt( n );
        if( tbl == NULL )
            continue;
        nvecs += tbl->GetActualSize();
    }
    SetNoVectors( nvecs );
}

// -------------------------------------------------------------------------
// GetTableAt: get table at position loc
//
inline
Table* Restaurant::GetTableAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Restaurant::GetTableAt: Memory access error.");
#endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetTableAt: set table at position loc
//
inline
void Restaurant::SetTableAt( int loc, Table* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Restaurant::SetTableAt: Memory access error.");
#endif
    values_[loc] = value;
}

// -------------------------------------------------------------------------
// NewTable: push table
//
inline
int Restaurant::NewTable( Table* value )
{
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw MYRUNTIME_ERROR("Restaurant::NewTable: Memory access error.");
        values_[loc] = value;
        novacs_--;
    }
    else {
        loc = length_;
        if( capacity_ <= loc )
            Realloc( TIMES2( capacity_ + 1 ));
        values_[loc] = value;
        length_++;
    }
    actlen_++;
    if( value ) {
        if( value->GetDish() == NULL )
            throw MYRUNTIME_ERROR("Restaurant::NewTable: Null Table's dish.");
        value->GetDish()->IncNoTables();
    }
    return loc;
}

// -------------------------------------------------------------------------
// RemTableAt: remove the table attached at position loc
//
inline
void Restaurant::RemTableAt( int loc, Table* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Restaurant::RemTableAt: Memory access error.");
#endif
    if( values_[loc] != value )
        throw MYRUNTIME_ERROR("Restaurant::RemTableAt: Memory access error.");
//     if( GetDestroy())
//         delete values_[loc];
    values_[loc] = NULL;
    PushVacant( loc );
    actlen_--;
    if( value )
        value->GetDish()->DecNoTables();
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int Restaurant::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Restaurant::GetVacantAt: Memory access error.");
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void Restaurant::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Restaurant::SetVacantAt: Memory access error.");
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index into the vector of vacancies
//
inline
void Restaurant::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__restaurant__
