/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CtxtCoefficients__
#define __CtxtCoefficients__

#include <stdio.h>
// #include <math.h>
#include <cmath>
#include "mybase.h"


// _________________________________________________________________________
// CLASS CtxtCoefficients
//
class CtxtCoefficients
{
public:
    CtxtCoefficients( size_t size );
    CtxtCoefficients( size_t size, float weight );

    ~CtxtCoefficients();

    size_t          GetLength() const { return length_; }

    float           GetCentralWeight() const { return cweight_; }
    void            SetCentralWeight( float value ) { cweight_ = value; }
    void            PutCentralWeight();

    float           GetMultiplier() const { return multip_; }

    const float*    GetCoefficients() const { return coeffs_; }
    float           GetCoefficientAt( size_t n ) const;
    void            SetCoefficientAt( size_t n, float value );

    const float*    GetLogCoefficients() const { return logcoeffs_; }
    float           GetLogCoefficientAt( size_t n ) const;
    void            SetLogCoefficientAt( size_t n );

    void            FindCoefficients();

//     bool            Read( FILE* );
    void            Write( FILE* ) const;
    void            WriteLogs( FILE* ) const;

protected:
    explicit CtxtCoefficients();

    static void     fdfunction( float x, float* f, float* df, void* );
    void            SetMultiplier( float value ) { multip_ = value; }
    float           ExpandCoefficients();

private:
    void            Init( size_t size );
    void            Destroy();
    void            SetLength( size_t value ) { length_ = value; }

private:
    size_t  length_;        //length of coefficient vector
    float   cweight_;       //central weight
    float   multip_;        //multiplier to be found
    float*  coeffs_;        //coefficients
    float*  logcoeffs_;     //precomputed logs of coefficients
};


////////////////////////////////////////////////////////////////////////////
// Class CtxtCoefficients inlines
//
// -------------------------------------------------------------------------
// GetCoefficientAt: get coefficient at a position
//
inline
float CtxtCoefficients::GetCoefficientAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !coeffs_ || GetLength() <= n )
        throw MYRUNTIME_ERROR("CtxtCoefficients::GetCoefficientAt: Memory access error.");
#endif
    return coeffs_[n];
}

// SetCoefficientAt: set coefficient value at a position
//
inline
void CtxtCoefficients::SetCoefficientAt( size_t n, float value )
{
#ifdef __DEBUG__
    if( !coeffs_ || GetLength() <= n )
        throw MYRUNTIME_ERROR("CtxtCoefficients::SetCoefficientAt: Memory access error.");
#endif
    coeffs_[n] = value;
}

// PutCoefficientAt: set coefficient value at the central position
//
inline
void CtxtCoefficients::PutCentralWeight()
{
    const float value = GetCentralWeight();
    const float precd = 0.0001f;
    const size_t cpost = GetLength() / 2;

    if( value <= 0.0f + precd || 1.0f - precd <= value )
        throw MYRUNTIME_ERROR("CtxtCoefficients::PutCentralWeight: Invalid central coefficient.");
    if( GetLength() < 1 || !( GetLength() % 2 ))
        throw MYRUNTIME_ERROR("CtxtCoefficients::PutCentralWeight: Even number of coefficients.");

    coeffs_[cpost] = value;
}

// -------------------------------------------------------------------------
// GetLogCoefficientAt: get the log of coefficient at a position
//
inline
float CtxtCoefficients::GetLogCoefficientAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !logcoeffs_ || GetLength() <= n )
        throw MYRUNTIME_ERROR("CtxtCoefficients::GetLogCoefficientAt: Memory access error.");
#endif
    return logcoeffs_[n];
}

// SetLogCoefficientAt: set the log of coefficient at a position
//
inline
void CtxtCoefficients::SetLogCoefficientAt( size_t n )
{
#ifdef __DEBUG__
    if( !logcoeffs_ || GetLength() <= n )
        throw MYRUNTIME_ERROR("CtxtCoefficients::SetLogCoefficientAt: Memory access error.");
#endif
    float value = GetCoefficientAt( n );
    if( 0.0f < value )
        logcoeffs_[n] = logf( value );
    else if( 0.0f == value )
        logcoeffs_[n] = -9999.0f;
    else
        throw MYRUNTIME_ERROR("CtxtCoefficients::SetLogCoefficientAt: Log of non-positive value.");
}

#endif//__CtxtCoefficients__
