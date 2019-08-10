/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <string.h>
#include <stdlib.h>

#include "Segments.h"

namespace SEG {

// -------------------------------------------------------------------------
// constructor: memory allocation
//
Segments::Segments( size_t size )
:   segments_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Realloc( size );
}

// -------------------------------------------------------------------------
// constructor: default
//
Segments::Segments()
:   segments_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    throw MYRUNTIME_ERROR( "Segments::Segments: Default construction prohibited." );
}

// -------------------------------------------------------------------------
// destructor:
//
Segments::~Segments()
{
    if( segments_ )
        free( segments_ );
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void Segments::Clear()
{
    memset( segments_, 0, sizeof(size_t) * TIMES2(capacity_));
    length_ = 0;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Segments::Realloc( size_t newcap )
{
    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        segments_ = ( size_t* )malloc( sizeof(size_t) * TIMES2(newcap));
    } else {
        segments_ = ( size_t* )realloc( segments_, sizeof(size_t) * TIMES2(newcap));
    }

    if( !segments_ )
        throw MYRUNTIME_ERROR( "Segments::Realloc: Not enough memory." );

    // fill uninitialized memory with zeros
    memset( segments_ + TIMES2(capacity_), 0, sizeof(size_t) * TIMES2(newcap-capacity_));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: push segment with boundaries left and right into their array
//
bool Segments::Push( size_t left, size_t right )
{
    int location = -1;

    if( right < left )
        throw MYRUNTIME_ERROR( "Segments::Push: Invalid segment." );

    if( Find( left, right, &location )) {
        //intersecting segment found, adjust by expanding it if needed
        Expand((size_t)location, left, right );

        RemoveAdjacentRight((size_t)location );
        RemoveAdjacentLeft((size_t)location );

        return false;
    }

    if( location < 0 )
        return false;

    if( capacity_ <= length_ ) {
        size_t newcap = TIMES2( capacity_ );
        if( newcap <= length_ )
            newcap = length_ + 1;
        Realloc( newcap );
    }

	length_++;

    for( int n = (int)(length_-1); n > location; n-- ) {
        SetLeftAt((size_t)n, GetLeftAt((size_t)( n-1 )));
        SetRightAt((size_t)n, GetRightAt((size_t)( n-1 )));
    }

    SetLeftAt((size_t)location, left );
    SetRightAt((size_t)location, right );

    return true;
}

// -------------------------------------------------------------------------
// Find: find a segment that intersects with the one given by arguments; 
// if not found, then set the location of the nearest greater segment
//
bool Segments::Find( size_t segleft, size_t segright, int* loc ) const
{
    int     left = 0;
    int     right = ( int )GetSize() - 1;
    int     middle = 0;
    int     comp = 1;
    size_t  midleft;        //left value of the middle segment
    size_t  midright;       //right value of the middle segment

    while( left <= right )
    {
        middle = ( left + right ) >> 1;
        midleft = GetLeftAt( middle );
        midright = GetRightAt( middle );
        comp = Compare( midleft, midright, segleft, segright );

        if( comp < 0 )      //if segment is greater than the middle element
            left = middle + 1;
        else if( comp > 0 ) //if segment is less than the middle element
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

// -------------------------------------------------------------------------
// RemoveAdjacentLeft: remove segment intersecting at left
//
void Segments::RemoveAdjacentLeft( size_t location )
{
    if( location == 0 || GetSize() <= location )
        return;

    size_t  prvleft = GetLeftAt( location - 1 );
    size_t  prvright = GetRightAt( location - 1 );
    size_t  locleft = GetLeftAt( location );
    size_t  locright = GetRightAt( location );

    if( Compare( prvleft, prvright, locleft, locright ) == 0 ) {
        Expand( location - 1, locleft, locright );
        Remove( location );
    }
}

// -------------------------------------------------------------------------
// RemoveAdjacentRight: remove segment intersecting at right
//
void Segments::RemoveAdjacentRight( size_t location )
{
    if( GetSize() <= location + 1 )
        return;

    size_t  locleft = GetLeftAt( location );
    size_t  locright = GetRightAt( location );
    size_t  nxtleft = GetLeftAt( location + 1 );
    size_t  nxtright = GetRightAt( location + 1 );

    if( Compare( locleft, locright, nxtleft, nxtright ) == 0 ) {
        Expand( location, nxtleft, nxtright );
        Remove( location + 1 );
    }
}

// -------------------------------------------------------------------------
// Expand: expand segment at the given location if needed
//
void Segments::Expand( size_t location, size_t left, size_t right )
{
    if( GetSize() <= location )
        return;

    size_t  fndleft = GetLeftAt( location );
    size_t  fndright = GetRightAt( location );

    if( left < fndleft )
        SetLeftAt( location, left );
    if( fndright < right )
        SetRightAt( location, right );
}

// -------------------------------------------------------------------------
// Remove: remove segment at the given position
//
void Segments::Remove( size_t loc )
{
    if( GetSize() <= loc )
        throw MYRUNTIME_ERROR( "Segments::Remove: Memory access error." );

    for( ; loc < GetSize() - 1; loc++ ) {
        SetLeftAt( loc, GetLeftAt( loc + 1 ));
        SetRightAt( loc, GetRightAt( loc + 1 ));
    }

    length_--;
}

// =========================================================================
// testing:
// Print: print array of segments
//
void Segments::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "%sSegments%s Segments:%s ", NL, NL, NL );
    for( size_t n = 0; n < GetSize(); n++ )
        fprintf( fp, " [%zu-%zu]", GetLeftAt( n ), GetRightAt( n ));

    fprintf( fp, "%s%s", NL, NL );
}

}//namespace SEG
