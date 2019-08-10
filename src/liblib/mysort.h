/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mysort_h__
#define __mysort_h__

#include <stdlib.h>

//definition of function types
typedef int ( *TSCompFunction )( const void*, size_t, size_t );
typedef int ( *TSSwapFunction )( void*, size_t, size_t );

//sorting of data with heap sort
void HeapSort( void* vector, size_t size, TSCompFunction cmfunc, TSSwapFunction swfunc );

//sorting of data with heap sort placing indices in array `inds' of size `size'
void HeapSortInd( size_t* inds, const void* vector, size_t size, TSCompFunction cmfunc );

#endif//__mysort_h__

