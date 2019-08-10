/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_gamma__
#define __extspsl_gamma__

namespace extspsl {

// Log(Gamma(x)); for x<0 the real part of the result is returned
int psl_lngamma_e( float x, float* result, float* err );

// Log(Gamma(x)); value of Gamma can be reconstructed by sgn*exp(result)
int psl_lngamma_sgn_e( float x, float* result, float* err, float* sgn );

// Gamma(x); the maximum value of x is SLC_G_GAMMA_XMAX
int psl_gamma_e( const float x, float* result, float* err );

// Regulated Gamma function, Gamma*(x) for x>0
int psl_gammastar_e( const float x, float* result, float* err );

// Reciprocal of the Gamma function, 1/Gamma(x)
int psl_gammainv_e( const float x, float* result, float* err );


// Factorial, n! = Gamma(n+1); the maximum value of n is SLC_G_FACT_NMAX
int psl_fact_e( const unsigned int n, float* result, float* err );

// Double factorial, n!!; the maximum value of n is SLC_G_DOUBLEFACT_NMAX
int psl_doublefact_e( const unsigned int n, float* result, float* err );

// Log of factorial, log(n!) = log(Gamma(n+1))
int psl_lnfact_e( const unsigned int n, float* result, float* err );

// Log of double factorial, log(n!!)
int psl_lndoublefact_e( const unsigned int n, float* result, float* err );


// Combinatorial factor n choose m, n!/(m!-(n-m)!)
int psl_choose_e( unsigned int n, unsigned int m, float* result, float* err );

// Log of n choose m, log(n!) - log(m!) - log((n-m)!)
int psl_lnchoose_e( unsigned int n, unsigned int m, float* result, float* err );

// Taylor coefficient, x^n/n! for x>=0, n>=0
int psl_taylorcoeff_e( const int n, const float x, float* result, float* err );

}//namespace extspsl

#endif//__extspsl_gamma__
