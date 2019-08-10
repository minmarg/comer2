/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __myalloc_h__
#define __myalloc_h__

#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "preproc.h"
#include "mystring.h"
#include "myexception.h"

// MYREALLOCBLK macro: reallocate memory for a block array (pointer to arrays) variable;
// NAME, name of variable;
// TYPE, type of array elements;
// BLKSIZE, block size of array;
// NREQSTD, requested number of elements (of size ELEMSIZE);
// NALLOCD, number of elements for which memory has been already allocated;
// MSG, message to throw on error;
#define MYREALLOCBLK( NAME, TYPE, BLKSIZE, NREQSTD, NALLOCD, MSG ) \
  if( NREQSTD > NALLOCD ) { \
    TYPE (*tmp_##NAME)[BLKSIZE] = NULL; \
    if( NALLOCD <= 0 ) \
        tmp_##NAME = ( TYPE(*)[BLKSIZE] )malloc( sizeof(TYPE) * BLKSIZE * NREQSTD ); \
    else \
        tmp_##NAME = ( TYPE(*)[BLKSIZE] )realloc( NAME, sizeof(TYPE) * BLKSIZE * NREQSTD ); \
    if( !tmp_##NAME ) \
        throw MYRUNTIME_ERROR( MSG ); \
    NAME = tmp_##NAME; \
    if( NALLOCD > 0 ) \
        tmp_##NAME += NALLOCD; \
    memset( tmp_##NAME, 0, sizeof(TYPE) * BLKSIZE * ( NREQSTD - NALLOCD )); \
  }

// myreallocblk: template version of reallocating memory for a block array variable;
// variable, address of variable to allocate memory for;
// T, type of array elements;
// S, type of size variables nreqstd and nallocd;
// blksize, block size of array;
// nreqstd, requested number of elements (of size ELEMSIZE);
// nallocd, number of elements for which memory has been already allocated;
template <typename T, typename S, int blksize>
void myreallocblk( T (**variable)[blksize], S nreqstd, S nallocd )
{
    if( nreqstd <= nallocd )
        return;
    if( variable == NULL )
        throw MYRUNTIME_ERROR("myreallocblk: Memory access error.");
    T (*tmp_variable)[blksize] = NULL;
    if( nallocd <= 0 )
        tmp_variable = ( T(*)[blksize] )malloc( sizeof(T) * blksize * nreqstd );
    else
        tmp_variable = ( T(*)[blksize] )realloc( *variable, sizeof(T) * blksize * nreqstd );
    if( !tmp_variable )
        throw MYRUNTIME_ERROR("myreallocblk: Not enough memory.");
    *variable = tmp_variable;
    if( nallocd > 0 )
        tmp_variable += nallocd;
    memset( tmp_variable, 0, sizeof(T) * blksize * ( nreqstd - nallocd ));
}

// MYREALLOC macro: reallocate memory for a (pointer) variable;
// NAME, name of variable;
// TYPE, type of elements;
// ELEMSIZE, size of one element a pointer points to;
// NREQSTD, requested number of elements (of size ELEMSIZE);
// NALLOCD, number of elements for which memory has been already allocated;
// MSG, message to throw on error;
#define MYREALLOC( NAME, TYPE, ELEMSIZE, NREQSTD, NALLOCD, MSG ) \
  if( NREQSTD > NALLOCD ) { \
    TYPE tmp_##NAME = NULL; \
    if( NALLOCD <= 0 ) \
        tmp_##NAME = ( TYPE )malloc( ELEMSIZE * NREQSTD ); \
    else \
        tmp_##NAME = ( TYPE )realloc( NAME, ELEMSIZE * NREQSTD ); \
    if( !tmp_##NAME ) \
        throw MYRUNTIME_ERROR( MSG ); \
    NAME = tmp_##NAME; \
    if( NALLOCD > 0 ) \
        tmp_##NAME += NALLOCD; \
    memset( tmp_##NAME, 0, ELEMSIZE * ( NREQSTD - NALLOCD )); \
  }

// myrealloc: template version of reallocating memory for a pointer variable;
// variable, address of variable to allocate memory for;
// T, type of elements;
// S, type of size variables nreqstd and nallocd;
// // // ELEMSIZE, size of one element a pointer points to;
// nreqstd, requested number of elements;
// nallocd, number of elements for which memory has been already allocated;
template <typename T, typename S>
void myrealloc( T** variable, S nreqstd, S nallocd )
{
    if( nreqstd <= nallocd )
        return;
    if( variable == NULL )
        throw MYRUNTIME_ERROR("myrealloc: Memory access error.");
    T* tmp_variable = NULL;
    if( nallocd <= 0 )
        tmp_variable = ( T* )malloc( sizeof(T) * nreqstd );
    else
        tmp_variable = ( T* )realloc( *variable, sizeof(T) * nreqstd );
    if( !tmp_variable )
        throw MYRUNTIME_ERROR("myrealloc: Not enough memory.");
    *variable = tmp_variable;
    if( nallocd > 0 )
        tmp_variable += nallocd;
    memset( tmp_variable, 0, sizeof(T) * ( nreqstd - nallocd ));
}

#endif//__myalloc_h__
