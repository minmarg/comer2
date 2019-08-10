/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <string.h>
#include <stdlib.h>

#include "mystring.h"
#include "myexception.h"
#include "BinarySearchStructure.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
SimpleVector::SimpleVector( size_t size )
:   values( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Realloc( size );
}

// -------------------------------------------------------------------------
// constructor: default
//
SimpleVector::SimpleVector()
:   values( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    throw myruntime_error("SimpleVector: Default initialization prohibited.", __EXCPOINT__ );
}

// -------------------------------------------------------------------------
// destructor:
//
SimpleVector::~SimpleVector()
{
    if( values )
        free( values );
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void SimpleVector::Realloc( size_t newcap )
{
    const void** tmp_values;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = ( const void** )malloc( sizeof( void* ) * newcap );
    } else {
        tmp_values = ( const void** )realloc( values, sizeof( void* ) * newcap );
    }

    if( !tmp_values )
        throw myruntime_error("SimpleVector: Not enough memory.", __EXCPOINT__ );

    values = tmp_values;

    // fill uninitialized memory
    memset( values + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: push back value in the vector; 
//     return flag indicating whether or not the key has been inserted
//
bool SimpleVector::Push( const void* key, int* /*not used*/ )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    values[length_] = key;

    length_++;

    return true;
}

// -------------------------------------------------------------------------
// InsertValueAt: insert value at position by shifting elements at the
//     right to the right
//
void SimpleVector::InsertValueAt( size_t loc, const void* key )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    if( length_ < loc )
        throw myruntime_error("SimpleVector: Unable to insert value.", __EXCPOINT__ );

    for( size_t n = length_; n > loc; n-- )
        values[n] = values[n-1];

    values[loc] = key;

    length_++;
}

// -------------------------------------------------------------------------
// SimpleFind: simple find of a key using a given comparison function
//
bool SimpleVector::SimpleFind( const void* key, TComparator comp, size_t* loc ) const
{
    if( values == NULL )
        return false;
    for( size_t i = 0; i < GetSize(); i++ ) {
        if((*comp)(values[i], key) == 0 ) {
            if( loc )
                *loc = i;
            return true;
        }
    }
    return false;
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void SimpleVector::Clear()
{
    memset( values, 0, sizeof( void* ) * capacity_ );
    length_ = 0;
}


// /////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////
// CLASS BinarySearchStructure
//
// constructor: initialization
//
BinarySearchStructure::BinarySearchStructure( TComparator comp, size_t size, bool keep )
:   SimpleVector( size ),
    comparator( comp ),
    comparator1( NULL ),
    keep_duplicates( keep ),
    params_( NULL )
{
}

BinarySearchStructure::BinarySearchStructure( 
        TComparator1 comp, size_t size, bool keep, void* pars )
:   SimpleVector( size ),
    comparator( NULL ),
    comparator1( comp ),
    keep_duplicates( keep ),
    params_( pars )
{
}

// -------------------------------------------------------------------------
// constructor: default
//
BinarySearchStructure::BinarySearchStructure()
:   comparator( NULL ),
    comparator1( NULL ),
    keep_duplicates( false ),
    params_( NULL )
{
    throw myruntime_error("BinarySearchStructure: Default initialization prohibited.",
                          __EXCPOINT__ );
}

// -------------------------------------------------------------------------
// destructor:
//
BinarySearchStructure::~BinarySearchStructure()
{
}

// -------------------------------------------------------------------------
// Push: insert value in the structure so that the order is preserved;
//      return flag whether or not the key has been inserted; 
//      the key will not be inserted if a duplicate exists and the 
//      structure does not keep them;
//      loc will be set to point to the duplicate element
//
bool BinarySearchStructure::Push( const void* key, int* loc )
{
#ifdef __DEBUG__
    if( !comparator && !comparator1 )
        throw myruntime_error("BinarySearchStructure: Memory access error.", __EXCPOINT__ );
#endif

    int location;

    if( Find( key, &location )) {
        //key found 
        if( !KeepDuplicates()) {
            if( loc )
                *loc = location;
            return false; //nothing to do if not to keep duplicates
        }
    }

    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    for( int n = (int)length_; n > location; n-- )
        values[n] = values[n-1];

    values[location] = key;

    length_++;

    if( loc )
        *loc = location;

    return true;
}

// -------------------------------------------------------------------------
// Find: find a key in the structure and remember the location of the found 
//      key;
//      if no key is found, then the location will point to the nearest 
//      greater element
//
bool BinarySearchStructure::Find( const void* key, int* loc ) const
{
    int     left = 0;
    int     right = ( int )GetSize() - 1;
    int     middle = 0;
    int     comp = 1;

    while( left <= right )
    {
        middle = ( left + right ) >> 1;
        if( comparator )
            comp = ( *comparator )( values[middle], key );
        else if( comparator1 )
            comp = ( *comparator1 )( values[middle], key, GetParams());
        else
            throw myruntime_error("BinarySearchStructure: Memory access error.", __EXCPOINT__ );

        if( comp < 0 )      //key is greater than the middle element
            left  = middle + 1;
        else if( comp > 0 ) //key is less than the middle element
            right = middle - 1;
        else {
            if( loc )
                *loc = middle;
            return true;
        }
    }

    if( loc ) {
        if( comp < 0 )
            *loc = left;
        else if( comp > 0 )
            *loc = middle;
    }

    return false;
}

