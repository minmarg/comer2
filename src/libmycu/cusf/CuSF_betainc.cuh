/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 *                                                                         *
 * CUDA implementation of the LOGARITHM of the incomplete beta function    *
 * with positive parameters for a device. The code is based on the incbetf *
 * function from the cephes library as indicated below.                    *
 ***************************************************************************/

#ifndef __CuSF_betainc_h__
#define __CuSF_betainc_h__

#include "extsp/psl.h"
#include "libmycu/cucom/cucommon.h"
#include "CuSF_com.h"

#define CUSF_LBETAINC_FAILVAL 999.0f

/* Incomplete beta integral
 *
 *
 * SYNOPSIS:
 *
 * float a, b, x, y, incbetf();
 *
 * y = incbetf( a, b, x );
 *
 *
 * DESCRIPTION:
 *
 * Returns incomplete beta integral of the arguments, evaluated
 * from zero to x.  The function is defined as
 *
 *                  x
 *     -            -
 *    | (a+b)      | |  a-1     b-1
 *  -----------    |   t   (1-t)   dt.
 *   -     -     | |
 *  | (a) | (b)   -
 *                 0
 *
 * The domain of definition is 0 <= x <= 1.  In this
 * implementation a and b are restricted to positive values.
 * The integral from x to 1 may be obtained by the symmetry
 * relation
 *
 *    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
 *
 * The integral is evaluated by a continued fraction expansion.
 * If a < 1, the function calls itself recursively after a
 * transformation to increase a to a+1.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,x) with a and b in the indicated
 * interval and x between 0 and 1.
 *
 * arithmetic   domain     # trials      peak         rms
 * Relative error:
 *    IEEE       0,30       10000       3.7e-5      5.1e-6
 *    IEEE       0,100      10000       1.7e-4      2.5e-5
 * The useful domain for relative error is limited by underflow
 * of the single precision exponential function.
 * Absolute error:
 *    IEEE       0,30      100000       2.2e-5      9.6e-7
 *    IEEE       0,100      10000       6.5e-5      3.7e-6
 *
 * Larger errors may occur for extreme ratios of a and b.
 *
 * ERROR MESSAGES:
 *   message         condition      value returned
 * incbetf domain     x<0, x>1          0.0
 *
 *
 * Cephes Math Library, Release 2.2:  July, 1992
 * Copyright 1984, 1987, 1992 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

__device__ __forceinline__
float cusf_lbetaincf_helper( float a, float b, float x );

__device__ __forceinline__
float cusf_lbetainccff( float a, float b, float x );

__device__ __forceinline__
float cusf_lbetaincdf( float a, float b, float x );

__device__ __forceinline__
float cusf_lbetaincpsf( float a, float b, float x );

// -------------------------------------------------------------------------
// cusf_lbetaincf: main function to calculate the logarithm of the 
// incomplete beta function
//
__device__ __forceinline__
float cusf_lbetaincf( float aa, float bb, float xx )
{
    float ans, t;

    if((xx <= 0.0f) || (xx >= 1.0f)) {
        if( xx == 0.0f )
            return -999.0f;
        if( xx == 1.0f )
            return 0.0f;
        // domain error
        return CUSF_LBETAINC_FAILVAL;
    }

    if( aa <= 0.0f )
        return CUSF_LBETAINC_FAILVAL;

    /* transformation for small aa */
    if( aa <= 1.0f ) {
        ans = cusf_lbetaincf_helper( aa+1.0f, bb, xx );
        ans = ans > SLC_LOG_SP_MIN ? __expf(ans): 0.0f;
        t = aa * __logf(xx) + bb * __logf( 1.0f-xx )
            + lgammaf(aa+bb) - lgammaf(aa+1.0f) - lgammaf(bb);
        if( t > SLC_LOG_SP_MIN )
            ans += __expf(t);
        return ans? __logf(ans): -999.0f;
    }

    return cusf_lbetaincf_helper( aa, bb, xx );
}

// -------------------------------------------------------------------------
// cusf_lbetaincf: main function to calculate the logarithm of the 
// incomplete beta function
//
__device__ __forceinline__
float cusf_lbetaincf_helper( float a, float b, float x )
{
    float ans, t;
//     int flag;

    /* see if x is greater than the mean */
    // MM: do not apply the symmetry relation, 
    // always calculate the logarithm for the given parameters
//     if( x > __fdividef(a,a+b)) {
//         flag = 1;
//         myhdswap(a, b);
//         t = x;
//         x = 1.0f - x;
//     } else {
//         flag = 0;
        t = 1.0f - x;
//     }

    /* Choose expansion for optimal convergence */
    if( b > 10.0f ) {
        if( fabsf(__fdividef(b*x,a)) < 0.3f ) {
            t = cusf_lbetaincpsf( a, b, x );
            return t;
//             goto bdone;
        }
    }

    ans = x * __fdividef(a+b-2.0f, a-1.0f);
    if( ans < 1.0f ) {
        ans = cusf_lbetainccff( a, b, x );
        t = b * __logf(t);
    } else {
        ans = cusf_lbetaincdf( a, b, x );
        t = (b-1.0f) * __logf(t);
    }

    t += a*__logf(x) + lgammaf(a+b) - lgammaf(a) - lgammaf(b);
    t += __logf(__fdividef(ans,a));

// bdone:
//     if( flag )
//         t = 1.0f - t;

    return t;
}



// -------------------------------------------------------------------------
/* Continued fraction expansion #1
 * for incomplete beta integral
 */
__device__ __forceinline__
float cusf_lbetainccff( float a, float b, float x )
{
    float xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    float k1, k2, k3, k4, k5, k6, k7, k8;
    float r, t, ans;

    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1.0f;
    k5 = 1.0f;
    k6 = b - 1.0f;
    k7 = k4;
    k8 = a + 2.0f;

    pkm2 = 0.0f;
    qkm2 = 1.0f;
    pkm1 = 1.0f;
    qkm1 = 1.0f;
    ans = 1.0f;
    r = 0.0f;

    for( int n = 0; n < 100; n++ )
    {
        xk = -__fdividef( x * k1 * k2, k3 * k4 );
        pk = pkm1 +  pkm2 * xk;
        qk = qkm1 +  qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = __fdividef( x * k5 * k6, k7 * k8 );
        pk = pkm1 +  pkm2 * xk;
        qk = qkm1 +  qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if( qk != 0.0f )
            r = __fdividef( pk, qk );
        if( r != 0.0f ) {
            t = fabsf(__fdividef( ans - r, r ));
            ans = r;
        }
        else
            t = 1.0f;

        if( t < CUSF_MACHEPF )
            break;//goto cdone;

        k1 += 1.0f;
        k2 += 1.0f;
        k3 += 2.0f;
        k4 += 2.0f;
        k5 += 1.0f;
        k6 -= 1.0f;
        k7 += 2.0f;
        k8 += 2.0f;

        if((fabsf(qk) + fabsf(pk)) > CUSF_BIGSPNUM) {
            pkm2 *= CUSF_MACHEPF;
            pkm1 *= CUSF_MACHEPF;
            qkm2 *= CUSF_MACHEPF;
            qkm1 *= CUSF_MACHEPF;
        }
        if((fabsf(qk) < CUSF_MACHEPF) || (fabsf(pk) < CUSF_MACHEPF)) {
            pkm2 *= CUSF_BIGSPNUM;
            pkm1 *= CUSF_BIGSPNUM;
            qkm2 *= CUSF_BIGSPNUM;
            qkm1 *= CUSF_BIGSPNUM;
        }
    }

// cdone:
    return ans;
}

// -------------------------------------------------------------------------
/* Continued fraction expansion #2
 * for incomplete beta integral
 */
__device__ __forceinline__
float cusf_lbetaincdf( float a, float b, float x )
{
    float xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    float k1, k2, k3, k4, k5, k6, k7, k8;
    float r, t, ans, z;

    k1 = a;
    k2 = b - 1.0f;
    k3 = a;
    k4 = a + 1.0f;
    k5 = 1.0f;
    k6 = a + b;
    k7 = a + 1.0f;
    k8 = a + 2.0f;

    pkm2 = 0.0f;
    qkm2 = 1.0f;
    pkm1 = 1.0f;
    qkm1 = 1.0f;
    z = __fdividef(x, 1.0f-x);
    ans = 1.0f;
    r = 0.0f;

    for( int n = 0; n < 100; n++ )
    {
        xk = -__fdividef( z * k1 * k2, k3 * k4 );
        pk = pkm1 +  pkm2 * xk;
        qk = qkm1 +  qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = __fdividef( z * k5 * k6, k7 * k8 );
        pk = pkm1 +  pkm2 * xk;
        qk = qkm1 +  qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if( qk != 0.0f )
            r = __fdividef( pk, qk );
        if( r != 0.0f ) {
            t = fabsf(__fdividef( ans - r, r ));
            ans = r;
        }
        else
            t = 1.0f;

        if( t < CUSF_MACHEPF )
            break;//goto cdone;

        k1 += 1.0f;
        k2 -= 1.0f;
        k3 += 2.0f;
        k4 += 2.0f;
        k5 += 1.0f;
        k6 += 1.0f;
        k7 += 2.0f;
        k8 += 2.0f;

        if((fabsf(qk) + fabsf(pk)) > CUSF_BIGSPNUM) {
            pkm2 *= CUSF_MACHEPF;
            pkm1 *= CUSF_MACHEPF;
            qkm2 *= CUSF_MACHEPF;
            qkm1 *= CUSF_MACHEPF;
        }
        if((fabsf(qk) < CUSF_MACHEPF) || (fabsf(pk) < CUSF_MACHEPF)) {
            pkm2 *= CUSF_BIGSPNUM;
            pkm1 *= CUSF_BIGSPNUM;
            qkm2 *= CUSF_BIGSPNUM;
            qkm1 *= CUSF_BIGSPNUM;
        }
    }

// cdone:
    return ans;
}


// -------------------------------------------------------------------------
// MM: the logarithm of ...
/* power series */
__device__ __forceinline__
float cusf_lbetaincpsf( float a, float b, float x )
{
    float t, u, y, s;

    y = a * __logf(x) + (b-1.0f)*__logf(1.0f-x) - __logf(a);
    y -= lgammaf(a) + lgammaf(b);
    y += lgammaf(a+b);

    t = __fdividef(x, 1.0f-x);
    s = 0.0f;
    u = 1.0f;
    do {
        b -= 1.0f;
        if( b == 0.0f )
            break;
        a += 1.0f;
        u *= __fdividef( t*b, a);
        s += u;
    }
    while( fabsf(u) > CUSF_MACHEPF );

    s = y + __logf(1.0f+s);

    return s;
}

#endif//__CuSF_betainc_h__
