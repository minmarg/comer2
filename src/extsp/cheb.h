/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_cheb__
#define __extspsl_cheb__

namespace extspsl {

// data type for a Chebyshev series over a given interval
//
typedef struct Tcheb_series_struct {
    float*      c;          // coefficients
    int         order;      // order of expansion
    float       a;          // lower interval point
    float       b;          // upper interval point
    int         order_sp;   // effective single precision order
} Tcheb_series;

int cheb_eval_e( const Tcheb_series* cs, const float x, float* result, float* err );

int psl_multiply_e( const float x, const float y, float* result, float* err );
int psl_multiply_err_e( const float x, const float dx, const float y, const float dy, float* result, float* err );

}//namespace extspsl

#endif//__extspsl_cheb__
