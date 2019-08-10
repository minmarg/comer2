/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include "root.h"

// Greatest common divisor
// (the code adopted from the NCBI BLAST toolbox)
int mygcd( int n, int d )
{
    int t;

    if( d < 0 )
        d = -d;

    if( d > n ) {
        t = n;
        n = d;
        d = t;
    }

    while( d ) {
        t = n % d;
        n = d;
        d = t;
    }

    return n;
}

#define  IM1 2147483563
#define  IM2 2147483399
#define  AM (1.0f/IM1)
#define  IMM1 (IM1-1)
#define  IA1 40014
#define  IA2 40692
#define  IQ1 53668
#define  IQ2 52774
#define  IR1 12211
#define  IR2 3791
#define  NTAB 32
#define  NDIV (1+IMM1/NTAB)
#define  EPS 1.2e-7f
#define  RNMX (1.0f-EPS)

// Long period (> 2 x 10^18 ) random number generator of L'Ecuyer with
// Bays-Durham shuffle and added safeguards.
// Returns a uniform random deviate between 0.0 and 1.0 (exclusively).
// Call with idum a negative integer to initialize; thereafter, do not
// alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less
// than 1.
// This code is adopted from
// W.H.Press et al's. "Numerical Recipes in C".

float LecuyerRand( long* idum )
{
    int         j;
    long        k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if( *idum <= 0 ) {                          //Initialize
        if( -*idum < 1 )                        //Be sure to prevent idum = 0
             *idum = 1;
        else *idum = -*idum;
        idum2 = *idum;

        for( j = NTAB + 7; j >= 0; j-- ) {      //Load the shuffle table (after 8 warm-ups)
            k = *idum / IQ1;
            *idum = IA1 * ( *idum - k * IQ1 ) - k * IR1;
            if( *idum < 0 ) *idum += IM1;
            if( j < NTAB )  iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = *idum / IQ1;                            //Start here when not initializing
    *idum = IA1 * ( *idum - k * IQ1 ) - k * IR1;//Compute idum=(IA1*idum) % IM1 without...
                                                //      overflows by Schrage's method
    if( *idum < 0 ) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * ( idum2 - k * IQ2 ) - k * IR2;//Compute idum2=(IA2*idum) % IM2 likewise
    if( idum2 < 0 ) idum2 += IM2;
    j = iy / NDIV;                              //Will be in the range 0..NTAB-1
    iy = iv[j] - idum2;                         //Here idum is shuffled, idum and idum2 are
                                                //      combined to generate output
    iv[j] = *idum;
    if( iy < 1 ) iy += IMM1;
    if(( temp = AM * iy ) > RNMX )              //Because one doesn't expect endpoint values
        return RNMX;
    return temp;
}

// -------------------------------------------------------------------------
// find the root of a function bracketed between x1 and x2 using a
// combination of Newton-Raphson and bisection methods.
// A bisection step is taken whenever Newton-Raphson would take the solution
// out of bounds, or whenever it is not reducing the size of the brackets
// rapidly enough.
// The root will be refined until its accuracy is reached within ±xacc or
// maximum number of iteration, maxit, has been passed. The method must be
// supplied with a routine that returns both the function value and the
// first derivative of the function.
// The code is as found in
//     Numerical Recipes in C by W.H.Press et al.
//
const char* root_by_NR_and_bisection(
        TFdFfunc fdfunction,
        float x1,
        float x2,
        float xacc,
        int maxit,
        void* params,
        float* result )
{
    float df, dx, dxold, f, fh, fl;
    float temp, xh, xl, rts;
    int n;

    if( !result )
        return((const char*)0 );

    ( *fdfunction )( x1, &fl, &df, params );
    ( *fdfunction )( x2, &fh, &df, params );

    if(( fl > 0.0 && fh > 0.0 ) || ( fl < 0.0 && fh < 0.0 ))
        return("Root evaluation: Root must be bracketed in given interval.");

    if( fl < 0.0 ) {                //Orient the search so that f (xl) < 0.
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5f * ( x1 + x2 );		//Initialize the guess for root,
    dxold = fabsf( x2 - x1 );       //the "stepsize before last,"
    dx = dxold;                     //and the last step.

    (*fdfunction)( rts, &f, &df, params );

    for( n = 0; n < maxit; n++ )                            //Loop over allowed iterations.
    {                                                       //Bisect if Newton out of range,
        if( ((( rts - xh ) * df - f ) * (( rts - xl ) * df - f ) > 0.0 ) ||
               ( fabsf( 2.0f * f ) > fabsf( dxold * df ))) {//or not decreasing fast enough.
            dxold = dx;
            dx = 0.5f * ( xh - xl );
            rts = xl + dx;
            if( xl == rts ) break; //return rts;            //Change in root is negligible.
        } else {                                            //Newton step acceptable. Take it.
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if( temp == rts ) break; //return rts;
        }
        if( fabsf( dx ) < xacc ) break; //return rts;       //Convergence criterion.
        (*fdfunction)( rts, &f, &df, params );
        //The one new function evaluation per iteration Maintain the bracket on the root.
        if( f < 0.0 )
            xl = rts;
        else
            xh = rts;
    }

    if( maxit <= n )
        return( "Root evaluation: maximum number of iterations reached." );

    *result = rts;
    return((const char*)0 );
}
