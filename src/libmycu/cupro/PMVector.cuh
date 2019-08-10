/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMVector_h__
#define __PMVector_h__

#include "libmycu/cupro/PMVectorFields.h"

// _________________________________________________________________________
// Class PMVector
// 
// Vector representing a position of profile model for scoring or 
// profile-specific data
//
class PMVector
{
public:
    __host__ __device__ PMVector( char* addr ): vector_(addr) {};
    __host__ __device__ explicit PMVector(): vector_(NULL) {};
    __host__ __device__ ~PMVector() {};

    __host__ __device__ PMVector& operator=( char* addr ) { vector_ = addr; return *this; }

    __host__ __device__ const char& operator[]( int idx ) const { return vector_[idx]; }

    __host__ __device__ char* GetField( int fldaddr ) { return vector_ + fldaddr; }

    __host__ __device__ void SetField( int fldaddr, FPTYPE value ) { 
            char* p = GetField( fldaddr );
            *(FPTYPE*)p = value; 
    }
#ifdef FPTYPEisINTYPE
#else
    __host__ __device__ void SetField( int fldaddr, INTYPE value ) { 
            char* p = GetField( fldaddr );
            *(INTYPE*)p = value; 
    }
#endif
    __host__ __device__ void SetField( int fldaddr, CHTYPE value ) { 
            char* p = GetField( fldaddr );
            *(CHTYPE*)p = value; 
    }

private:
    char* vector_;//complete data for one position
};

////////////////////////////////////////////////////////////////////////////
// PMVector inlines
//

#endif//__PMVector_h__
