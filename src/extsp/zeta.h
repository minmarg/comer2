/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_zeta__
#define __extspsl_zeta__

namespace extspsl {

// The Hurwitz zeta function, SUM (k+q)^(-s)
int psl_hzeta_e( const float s, const float q, float* result, float* err );


// The Riemann zeta function, SUM_k (k)^(-s) for s<>1
int psl_zeta_e( const float s, float* result, float* err );

// The Riemann zeta function, SUM_k (k)^(-n) for integer n<>1
int psl_zeta_int_e( const int n, float* result, float* err );


// The Riemann zeta function minus 1 for s<>1;
int psl_zetam1_e( const float s, float* result, float* err );

// The Riemann zeta function minus 1 for integer n<>1;
int psl_zetam1_int_e( const int n, float* result, float* err );


// The eta function, (1-2^(1-n))zeta(n) for integer n 
int psl_eta_int_e( int n, float* result, float* err );

// The eta function, (1-2^(1-s))zeta(s)
int psl_eta_e( const float s, float* result, float* err );

}//namespace extspsl

#endif//__extspsl_zeta__
