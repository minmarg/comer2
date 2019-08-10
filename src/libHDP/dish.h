/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __dish__
#define __dish__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "liblib/mybase.h"
#include "extsp/pslvector.h"
#include "basin.h"

// -------------------------------------------------------------------------
// class Dish: reusable vector of values (objects) from the same cluster
//
class Dish
{
public:
    Dish( int size );
    Dish( const Dish& );
    explicit Dish();
    ~Dish();

    Dish&  operator=( const Dish& );

    time_t      GetTime() const { return time_; }
    void        SetTime( time_t value ) { time_ = value; }

    int         GetReadSize() const { return prosize_; }
    void        SetReadSize( int value ) { prosize_ = value; }

    int         GetReadNoTables() const { return pronotbls_; }
    void        SetReadNoTables( int value ) { pronotbls_ = value; }

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    int         GetDishSize() const;

    int         GetVectorNIndAt( int loc ) const;
    extspsl::Pslvector* GetVectorNAt( int loc ) const;
    void        SetVectorNIndAt( int loc, int n );
    int         NewVectorNInd( int n );//push value
    void        RemValueAt( int loc, const extspsl::Pslvector* value );//remove value

    bool        GetProcessedAt( int loc ) const;
    void        SetProcessedAt( int loc, bool value );

    bool        GetProcessed() const { return dishproced_; }
    void        SetProcessed( bool value ) { dishproced_ = value; }

    void        Copy( const Dish& vector );     //copy elements
    void        Clear();                        //clear all elements
    void        Print( FILE* ) const;           //print vector

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    const Basin* GetBasin() const { return basin_; }
    void         SetBasin( const Basin* value ) { basin_ = value; }

    float       GetTmpValue() const { return tmpval_; }
    void        SetTmpValue( float value ) { tmpval_ = value; }

    int         GetNoTables() const { return notables_; }
    void        IncNoTables() { notables_++; }
    void        DecNoTables() { if( 0 < notables_ ) notables_--; }
protected:
    void        SetNoTables( int value ) { notables_ = value; }
public:
    void        Reserve( int size );
    void        ReserveValues( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyValues();
    int         GetCapacity() const { return capacity_; }

    int*        GetVector() const { return values_; }

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
    time_t      time_;          //time at which dish created
    const Basin* basin_;        //global basin of vectors
    int*        values_;        //values (indices of vectors) attributable to the dish
    bool*       processed_;     //flags of processed objects
    bool        dishproced_;    //`processed' flag of whole dish
    int         prosize_;       //proposal size of dish
    int         pronotbls_;     //proposal number of tables with this dish
    int         actlen_;        //actual number of values
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the vector
    bool        destroy_;       //flag of destroying object
    int*        vacancies_;     //vector of indices of vacancies
    int         novacs_;        //number of vacant elements
    int         capvac_;        //capacity of vacancies_
    int         notables_;      //number of tables associated with dish k
    float       tmpval_;        //temporary value used to store probability
};


// -------------------------------------------------------------------------
// GetDishSize: get dish size: either actual size of vectors or proposal
//  size if set
//
inline
int Dish::GetDishSize() const
{
    return ( 0 < GetActualSize())? GetActualSize(): GetReadSize();
}

// -------------------------------------------------------------------------
// Reserve: Reserve space for vectors
//
inline
void Dish::Reserve( int size )
{
    ReserveValues( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveValues: Reserve space for vector of values
//
inline
void Dish::ReserveValues( int size )
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
void Dish::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void Dish::Destroy()
{
    DestroyValues();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyValues: Destroy vector of values
//
inline
void Dish::DestroyValues()
{
    if( values_ )
        free( values_ );
    if( processed_ )
        free( processed_ );
    values_ = NULL;
    processed_ = NULL;
    actlen_ = 0;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// DestroyVacans: Destroy vector of vacancies
//
inline
void Dish::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}

// -------------------------------------------------------------------------
// GetVectorNIndAt: get the index of the vector from the basin, stored at 
//  `loc'
//
inline
int Dish::GetVectorNIndAt( int loc ) const
{
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Dish::GetVectorNIndAt: Memory access error.");
    return values_[loc];
}

// -------------------------------------------------------------------------
// GetVectorNAt: get vector from the basin with index stored at `loc'
//
inline
extspsl::Pslvector* Dish::GetVectorNAt( int loc ) const
{
    if( !basin_ || !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Dish::GetVectorNAt: Memory access error.");
    int ndx = values_[loc];

    if( ndx < 0 || basin_->GetSize() <= ndx )
        throw MYRUNTIME_ERROR("Dish::GetVectorNAt: Memory access error.");
    return basin_->GetValueAt( ndx );
}

// -------------------------------------------------------------------------
// SetVectorNIndAt: set vector's basin index at `loc'
//
inline
void Dish::SetVectorNIndAt( int loc, int n )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Dish::SetVectorNIndAt: Memory access error.");
#endif
    values_[loc] = n;
}

// -------------------------------------------------------------------------
// NewVectorNInd: push vector's basin index
//
inline
int Dish::NewVectorNInd( int n )
{
    if( !values_ )
        throw MYRUNTIME_ERROR("Dish::NewVectorNInd: Memory access error.");
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw MYRUNTIME_ERROR("Dish::NewVectorNInd: Memory access error.");
        values_[loc] = n;
        novacs_--;
    }
    else {
        if( capacity_ <= ( loc = length_ ))
            Realloc( TIMES2( capacity_ + 1 ));
        values_[loc] = n;
        length_++;
    }
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemValueAt: remove value at position loc
//
inline
void Dish::RemValueAt( int loc, const extspsl::Pslvector* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Dish::RemValueAt: Memory access error.");
#endif
    if( GetVectorNAt( loc ) != value )
        throw myruntime_error("Dish::RemValueAt: Memory access error.");
    values_[loc] = -1;
    processed_[loc] = false;
    PushVacant( loc );
    actlen_--;
}

// -------------------------------------------------------------------------
// GetProcessedAt: get a flag of processed object at position loc
//
inline
bool Dish::GetProcessedAt( int loc ) const
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR( "Dish::GetProcessedAt: Memory access error.");
#endif
    return processed_[loc];
}

// -------------------------------------------------------------------------
// SetProcessedAt: set a flag of processed object at position loc
//
inline
void Dish::SetProcessedAt( int loc, bool value )
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Dish::SetProcessedAt: Memory access error.");
#endif
    processed_[loc] = value;
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int Dish::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Dish::GetVacantAt: Memory access error.");
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void Dish::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Dish::SetVacantAt: Memory access error.");
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index in the vector of vacancies
//
inline
void Dish::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__dish__
