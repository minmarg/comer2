/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdlib.h>
#include "mysort.h"

// -------------------------------------------------------------------------
//
static inline 
void downheap( void* vector, const size_t N, size_t k, TSCompFunction cmfunc, TSSwapFunction swfunc )
{
    while( k <= ( N >> 1 )) {
        size_t j = k + k;

        if( j < N && ( *cmfunc )( vector, j, j + 1 ) < 0 )
            j++;

        if(( *cmfunc )( vector, k, j ) < 0 )
           ( *swfunc )( vector, j, k);
        else
            break;

        k = j;
    }
}

// -------------------------------------------------------------------------
// HeapSort: sort the array in ascending order; this is a true inplace
//  algorithm with N log N operations; worst case is ~20% slower
//  /GNU scientific library/
//
void HeapSort( void* vector, size_t size, TSCompFunction cmfunc, TSSwapFunction swfunc )
{
    size_t N, k;

    if( cmfunc == NULL || swfunc == NULL )
        return;
    if( size < 1 )
        return;

    //set N to the last element number
    N = size - 1;
    k = N >> 1;
    k++;//compensate the first use of k--
    do {
        k--;
        downheap( vector, N, k, cmfunc, swfunc );
    } while( 0 < k );

    while( 0 < N ) {
        //swap the elements
        ( *swfunc )( vector, 0, N );

        //then process the heap
        N--;
        downheap( vector, N, 0, cmfunc, swfunc );
    }
}

// =========================================================================
//
static inline
void downheapind( size_t* inds, const void* vector, const size_t N, size_t k, TSCompFunction cmfunc )
{
    const size_t indk = inds[k];

    while( k <= ( N >> 1 )) {
        size_t j = k + k;

        if( j < N && ( *cmfunc )( vector, inds[j], inds[j+1]) < 0 )
            j++;

        if( 0 <= ( *cmfunc )( vector, indk, inds[j] ))
            break;

        inds[k] = inds[j];
        k = j;
    }

    inds[k] = indk;
}

// -------------------------------------------------------------------------
// HeapSort: sort the array in ascending order placing indices in 
//  array `inds' of size `size'; NOTE: `inds' must be preallocated;
//  inplace algorithm with N log N operations; worst case is ~20% slower
//  /GNU scientific library/
//
void HeapSortInd( size_t* inds, const void* vector, size_t size, TSCompFunction cmfunc )
{
    size_t i, k, N;

    if( cmfunc == NULL )
        return;
    if( size <= 0 )
        return;

    for( i = 0; i < size; i++ )
        inds[i] = i;//set permutation to identity

    //set N to the last element number
    N = size - 1;
    k = N >> 1;
    k++;//compensate the first use of k--
    do {
        k--;
        downheapind( inds, vector, N, k, cmfunc );
    } while( 0 < k );

    while( 0 < N ) {
        //swap the elements
        size_t tmp = inds[0];
        inds[0] = inds[N];
        inds[N] = tmp;

        //then process the heap
        N--;
        downheapind( inds, vector, N, 0, cmfunc );
    }
}
