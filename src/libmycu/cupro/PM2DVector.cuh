/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PM2DVector_h__
#define __PM2DVector_h__

#include "libmycu/cupro/PM2DVectorFields.h"

// _________________________________________________________________________
// Class PM2DVector
// 
// 2D Vector representing a profile model or a batch of them
//
class PM2DVector
{
public:
    __host__ __device__ PM2DVector( char** addr ): vector_(addr) {};
    __host__ __device__ explicit PM2DVector(): vector_(NULL) {};
    __host__ __device__ ~PM2DVector() {};

    __host__ __device__ PM2DVector& operator=( char** addr ) { vector_ = addr; return *this; }

    __host__ __device__ const char* operator[]( int idx ) const { return vector_[idx]; }

    __host__ __device__ char*& GetFieldVector( int vndx ) { return vector_[vndx]; }

    __host__ __device__ char* GetField( int vndx, int fldaddr ) { return GetFieldVector(vndx) + fldaddr; }


    //set a value for a field and move the pointer to the next position
    __host__ __device__ void SetFieldNext( int vndx, int ndx, FPTYPE value ) { 
            FPTYPE*& p = (FPTYPE*&)GetFieldVector(vndx);
            *(p++ + ndx) = value; 
    }
    __host__ __device__ void SetField( int vndx, int ndx, FPTYPE value ) { 
            FPTYPE* p = (FPTYPE*)GetFieldVector(vndx);
            p[ndx] = value; 
    }

#if !defined(LNTYPEisINTYPE) && !defined(FPTYPEisINTYPE)
    __host__ __device__ void SetFieldNext( int vndx, int ndx, LNTYPE value ) { 
            LNTYPE*& p = (LNTYPE*&)GetFieldVector(vndx);
            *(p++ + ndx) = value; 
    }
    __host__ __device__ void SetField( int vndx, int ndx, LNTYPE value ) { 
            LNTYPE* p = (LNTYPE*)GetFieldVector(vndx);
            p[ndx] = value; 
    }
#endif

#ifndef FPTYPEisINTYPE
    __host__ __device__ void SetFieldNext( int vndx, int ndx, INTYPE value ) { 
            INTYPE*& p = (INTYPE*&)GetFieldVector(vndx);
            *(p++ + ndx) = value; 
    }
    __host__ __device__ void SetField( int vndx, int ndx, INTYPE value ) { 
            INTYPE* p = (INTYPE*)GetFieldVector(vndx);
            p[ndx] = value; 
    }
#endif

#ifndef CHTYPEisINTYPE
    __host__ __device__ void SetFieldNext( int vndx, int ndx, CHTYPE value ) { 
            CHTYPE*& p = (CHTYPE*&)GetFieldVector(vndx);
            *(p++ + ndx) = value; 
    }
    __host__ __device__ void SetField( int vndx, int ndx, CHTYPE value ) { 
            CHTYPE* p = (CHTYPE*)GetFieldVector(vndx);
            p[ndx] = value; 
    }
#endif

private:
    char** vector_;//complete data
};

////////////////////////////////////////////////////////////////////////////
// PM2DVector inlines
//

#endif//__PM2DVector_h__
