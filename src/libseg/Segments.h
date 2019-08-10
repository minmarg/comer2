/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __segments_h__
#define __segments_h__

#include "liblib/mybase.h"

#include <stdio.h>

namespace SEG {

// _________________________________________________________________________
// Class Segments
//
class Segments
{
public:
    Segments( size_t size );
    virtual ~Segments();

    size_t  GetSize() const { return length_; }

    size_t  GetLeftAt( size_t loc ) const;//obtain left boundary value
    size_t  GetRightAt( size_t loc ) const;//obtain right boundary value

    bool    Push( size_t left, size_t right );//save segment boundaries
    void    Clear();//clear all elements

    void    Print( FILE* );//testing

protected:
    explicit Segments();

    void    Realloc( size_t newcap );
    size_t  GetCapacity() const { return capacity_; }

    bool    Find( size_t left, size_t right, int* loc = NULL ) const;//find segment

    void    SetLeftAt( size_t loc, size_t left );//set left boundary value
    void    SetRightAt( size_t loc, size_t right );//set right boundary value

    void    Remove( size_t );//remove segment at the given position
    void    RemoveAdjacentLeft( size_t );//remove intersecting segment at left
    void    RemoveAdjacentRight( size_t );//remove intersecting segment at right
    void    Expand( size_t loc, size_t left, size_t right );//expand segment at the given location if needed

    static int Compare( size_t oneleft, size_t oneright, size_t anoleft, size_t anoright );

protected:
    size_t*     segments_;//array of segments
    size_t      length_;//number of segments
    size_t      capacity_;//current capacity
};

// INLINES ...

// -------------------------------------------------------------------------
// GetLeftAt: get the left boundary value of segment at position loc
//
inline
size_t Segments::GetLeftAt( size_t loc ) const
{
#ifdef __DEBUG__
    if( !segments_ || length_ <= loc )
        throw MYRUNTIME_ERROR( "Segments::GetLeftAt: Memory access error." );
#endif
    return segments_[TIMES2(loc)];
}

// -------------------------------------------------------------------------
// GetRightAt: get the right boundary value of segment at position loc
//
inline
size_t Segments::GetRightAt( size_t loc ) const
{
#ifdef __DEBUG__
    if( !segments_ || length_ <= loc )
        throw MYRUNTIME_ERROR( "Segments::GetRightAt: Memory access error." );
#endif
    return segments_[TIMES2(loc)+1];
}

// -------------------------------------------------------------------------
// SetLeftAt: set the left boundary value of segment at position loc
//
inline
void Segments::SetLeftAt( size_t loc, size_t left )
{
#ifdef __DEBUG__
    if( !segments_ || length_ <= loc )
        throw MYRUNTIME_ERROR( "Segments::SetLeftAt: Memory access error." );
#endif
    segments_[TIMES2(loc)] = left;
}

// -------------------------------------------------------------------------
// SetRightAt: set the right boundary value of segment at position loc
//
inline
void Segments::SetRightAt( size_t loc, size_t right )
{
#ifdef __DEBUG__
    if( !segments_ || length_ <= loc )
        throw MYRUNTIME_ERROR( "Segments::SetRightAt: Memory access error." );
#endif
    segments_[TIMES2(loc)+1] = right;
}

// -------------------------------------------------------------------------
// Compare: compare segments; return 0 if two segments intersect, 
// 1 if one is greater than another, -1 otherwise;
// NOTE: assumes that boundaries are correct
//
inline
int Segments::Compare( size_t oneleft, size_t oneright, size_t anoleft, size_t anoright )
{
    if( oneright < anoleft )
        return -1;

    if( anoright < oneleft )
        return 1;

    return 0;
}

}//namespace SEG

#endif//__segments_h__
