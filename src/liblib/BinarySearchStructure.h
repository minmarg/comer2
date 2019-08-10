/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __BinarySearchStructure__
#define __BinarySearchStructure__

#include "debug.h"
#include "mystring.h"
#include "myexception.h"

//Comparison function:
// functions are assumed to return zero if two keys are equal,
//  negative value if key1 is less than key2,
//  positive otherwise
typedef int ( *TComparator )( const void* key1, const void* key2 );
typedef int ( *TComparator1 )( const void* key1, const void* key2, void* pars );

// _________________________________________________________________________
// Class SimpleVector
//
class SimpleVector
{
public:
    SimpleVector( size_t size );
    virtual ~SimpleVector();

    size_t          GetSize() const { return length_; }

    const void*     GetValueAt( size_t loc ) const;                 //get value at position
    void            SetValueAt( size_t loc, const void* value );    //set value at position
    void            InsertValueAt( size_t loc, const void* value ); //insert value at position
    virtual bool    Push( const void* key, int* loc = NULL );       //save value in the structure
    bool            SimpleFind( const void* key, TComparator, size_t* loc = NULL ) const;
    const void*     RemoveLast();                                   //remove the last entry
    void            Clear();                                        //clear all elements

protected:
    explicit SimpleVector();

    void    Realloc( size_t newcap );
    size_t  GetCapacity() const { return capacity_; }

protected:
    const void**    values;     //table of values (void pointers); each value corresponds to a hash key
    size_t          length_;    //number of values present in the vector
    size_t          capacity_;  //capacity of the vector

};

// _________________________________________________________________________
// Class BinarySearchStructure
//
class BinarySearchStructure: public SimpleVector
{
public:
    BinarySearchStructure( TComparator, size_t size, bool keep = false );
    BinarySearchStructure( TComparator1, size_t size, bool keep = false, void* pars = NULL );
    virtual ~BinarySearchStructure();

    bool    KeepDuplicates() const  { return keep_duplicates; }

    virtual bool    Push( const void* key, int* loc = NULL );       //save value in the structure
    virtual bool    Find( const void* key, int* loc = NULL ) const; //find a key in the structure

    void*           GetParams() const { return params_; }
    void            SetParams( void* value ) { params_ = value; }

protected:
    explicit BinarySearchStructure();

private:
    TComparator     comparator;         //comparison function
    TComparator1    comparator1;        //comparison function
    bool            keep_duplicates;    //whether to keep duplicate keys
    void*           params_;            //private parameters
};

// INLINES ...

// -------------------------------------------------------------------------
// GetValueAt: get value at position
//
inline
const void* SimpleVector::GetValueAt( size_t loc ) const
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring("SimpleVector: Memory access error."), __EXCPOINT__ );
#endif
    return values[loc];
}

// -------------------------------------------------------------------------
// SetValueAt: Set value at position
//
inline
void SimpleVector::SetValueAt( size_t loc, const void* value )
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring("SimpleVector: Memory access error."), __EXCPOINT__ );
#endif
    values[loc] = value;
}

// -------------------------------------------------------------------------
// RemoveLast: remove the last entry
//
inline
const void* SimpleVector::RemoveLast()
{
    if( values == NULL || length_ < 1 )
        return NULL;
    length_--;
    return values[length_];
}

#endif//__BinarySearchStructure__
