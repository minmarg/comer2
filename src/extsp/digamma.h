/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_digamma__
#define __extspsl_digamma__

namespace extspsl {

int psl_psi_int_e( const int n, float* result, float* err );
int psl_psi_e( const float x, float* result, float* err );

// Trigamma function, psi'(n), n>0
int psl_psi_1_int_e( const int n, float* result, float* err );

// Trigamma function, psi'(n) for general x
int psl_psi_1_e( const float x, float* result, float* err );

// Polygamma function, psi^(n)(x) for n>=0, x>0
int psl_psi_n_e( const int n, const float x, float* result, float* err );

}//namespace extspsl

#endif//__extspsl_digamma__
