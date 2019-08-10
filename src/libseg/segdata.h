/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __segdata_h__
#define __segdata_h__

#include <stdlib.h>

namespace SEG {

//CONST DATA
//
extern const size_t s_SZPRSQR;// =  101;
extern const size_t s_SZPRLOG;// =  200;
extern const size_t s_SZPRENT;// =  200;
extern const size_t s_SZPRLFT;// = 1000;

// _________________________________________________________________________
// Abstract structure of precomputed values
//
struct _TPRECOMPUTED
{
    _TPRECOMPUTED( size_t );
    virtual ~_TPRECOMPUTED();

    void    Init();
    float   GetValueOf( size_t v ) const
    {
        if( v < GetSize())
            return data_[v];
        return Compute( v );
    }

    size_t  GetSize() const  { return no_vals_; }

    virtual float Compute( size_t v ) const = 0;

private:
    const size_t    no_vals_;
    float*          data_;
};

// _________________________________________________________________________
// Precomputed squared values
//
struct _TSQUARES: public _TPRECOMPUTED
{
    _TSQUARES( size_t = s_SZPRSQR );
    virtual ~_TSQUARES();
    virtual float Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed logarithm values
//
struct _TLOGARITHMS: public _TPRECOMPUTED
{
    _TLOGARITHMS( size_t = s_SZPRLOG );
    virtual ~_TLOGARITHMS();
    virtual float Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed entropy values
//
struct _TPARTIAL_ENTROPIES: public _TPRECOMPUTED
{
    _TPARTIAL_ENTROPIES( size_t = s_SZPRENT );
    virtual ~_TPARTIAL_ENTROPIES();
    virtual float Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed log-factorial values
//
struct _TLOG_FACT
{
    _TLOG_FACT();
    ~_TLOG_FACT();

    float   GetValueOf( size_t v );
    size_t  GetSize() const  { return no_vals_; }

private:
    void    Precompute( size_t newsize );
    void    SetSize( size_t newsize ) { no_vals_ = newsize; }

    size_t  no_vals_;
    float*  data_;
};

// =========================================================================
// Declarations of extern constant variables
//
extern const struct _TSQUARES               PRESQUARES;
extern const struct _TLOGARITHMS            LOGARITHMS;
extern const struct _TPARTIAL_ENTROPIES     PRT_ENTROPIES;
extern struct _TLOG_FACT                    LOG_FACT;

}//namespace SEG

#endif//__segdata_h__
