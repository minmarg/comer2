/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 *                                                                         *
 * CUDA implementation of the LOGARITHM of the incomplete gamma function   *
 * with positive parameters for a device. The code is based on the igamf   *
 * function from the cephes library as indicated below.                    *
 ***************************************************************************/

#ifndef __CuSF_gammainc_h__
#define __CuSF_gammainc_h__

#include "extsp/psl.h"
#include "libmycu/cucom/cucommon.h"
#include "CuSF_com.h"

#define CUSF_LGAMMAINC_FAILVAL 999.0f

/* Incomplete gamma integral
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamf();
 *
 * y = igamf( a, x );
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        20000       7.8e-6      5.9e-7
 *
 */
/* Complemented incomplete gamma integral
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamcf();
 *
 * y = igamcf( a, x );
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       7.8e-6      5.9e-7
 *
 *
 * Cephes Math Library Release 2.2: June, 1992
 * Copyright 1985, 1987, 1992 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

 __device__ __forceinline__
float cusf_lgammainc_Qf( float a, float x );

__device__ __forceinline__
float cusf_lgammainc_Pf( float a, float x );

// -------------------------------------------------------------------------
// cusf_lgammainc_Qf: main function to calculate the logarithm of the 
// normalized upper incomplete gamma function
//
__device__ __forceinline__
float cusf_lgammainc_Qf( float a, float x )
{
    float ans, c, yc, ax, y, z;
    float pk, pkm1, pkm2, qk, qkm1, qkm2;
    float r, t;

    if((x <= 0.0f) || (a <= 0.0f))
        return 0.0f;

    if((x < 1.0f) || (x < a)) {
        //NOTE: indirect recursion makes the compiler use stack
        ans = cusf_lgammainc_Pf(a,x);
        ans = ans > SLC_LOG_SP_MIN ? __expf(ans): 0.0f;
        return ans < 1.0f ? __logf(1.0f - ans): -999.0f;
    }

    ax = a * __logf(x) - x - lgammaf(a);
    /**/

    /* continued fraction */
    y = 1.0f - a;
    z = x + y + 1.0f;
    c = 0.0f;
    pkm2 = 1.0f;
    qkm2 = x;
    pkm1 = x + 1.0f;
    qkm1 = z * x;
    ans = __fdividef(pkm1, qkm1);

    do {
        c += 1.0f;
        y += 1.0f;
        z += 2.0f;
        yc = y * c;
        pk = pkm1 * z  -  pkm2 * yc;
        qk = qkm1 * z  -  qkm2 * yc;
        if( qk != 0.0f ) {
            r = __fdividef(pk, qk);
            t = fabsf(__fdividef( ans - r, r ));
            ans = r;
        }
        else
            t = 1.0f;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( fabsf(pk) > CUSF_BIGSPNUM ) {
            pkm2 *= CUSF_MACHEPF;
            pkm1 *= CUSF_MACHEPF;
            qkm2 *= CUSF_MACHEPF;
            qkm1 *= CUSF_MACHEPF;
        }
    }
    while( t > CUSF_MACHEPF );

    return __logf(ans) + ax;
}



// -------------------------------------------------------------------------
// cusf_lgammainc_Pf: main function to calculate the logarithm of the 
// normalized lower incomplete gamma function
//
/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */
__device__ __forceinline__
float cusf_lgammainc_Pf( float a, float x )
{
    float ans, ax, c, r;

    if((x <= 0.0f) || ( a <= 0.0f))
        return -999.0f;

    if((x > 1.0f) && (x > a)) {
        ans = cusf_lgammainc_Qf(a,x);
        ans = ans > SLC_LOG_SP_MIN ? __expf(ans): 0.0f;
        return ans < 1.0f ? __logf(1.0f - ans): -999.0f;
    }

    /* Compute  x**a * exp(-x) / gamma(a)  */
    ax = a * __logf(x) - x - lgammaf(a);
    /**/

    /* power series */
    r = a;
    c = 1.0f;
    ans = 1.0f;

    do {
        r += 1.0f;
        c *= __fdividef(x,r);
        ans += c;
    }
    while( __fdividef(c,ans) > CUSF_MACHEPF );

    // MM: ans>0
    return __logf(ans) + ax - __logf(a);
}

#endif//__CuSF_gammainc_h__
