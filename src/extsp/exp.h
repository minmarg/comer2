/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_exp__
#define __extspsl_exp__

namespace extspsl {

// exp(x)
int psl_exp_e( const float x, float* result, float* err );

// exp(x) with add. error estimate
int psl_exp_err_e( const float x, const float dx, float* result, float* err );


// y exp(x)
int psl_exp_mult_e( const float x, const float y, float* result, float* err );

// y exp(x) with add. error estimate
int psl_exp_mult_err_e( const float x, const float dx, const float y, const float dy, float* result, float* err );


// exp(x)-1; accurate for small x
int psl_expm1_e( const float x, float* result, float* err );


// (exp(x)-1)/x; accurate for small x
int psl_exprel_e( const float x, float* result, float* err );

// (exp(x)-1)/x; using continued fraction representation
int psl_exprel_n_CF_e( const float N, const float x, float* result, float* err );

// 2(exp(x)-1-x)/x^2; accurate for small x
int psl_exprel_2_e( float x, float* result, float* err );

// N-relative exponential, n-th generalization of  exprel and exprel_2
int psl_exprel_n_e( const int N, const float x, float* result, float* err );

}//namespace extspsl

#endif//__extspsl_exp__
