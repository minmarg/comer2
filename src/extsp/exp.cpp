/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
// *** the code adopted from GSL ***

#include <math.h>
#include "psl.h"
#include "pslerror.h"
#include "gamma.h"
#include "exp.h"

namespace extspsl {

// -------------------------------------------------------------------------
// * Evaluate the continued fraction for exprel.
// * [Abramowitz+Stegun, 4.2.41]
//
static int exprel_n_CF( const float N, const float x, float* result, float* err )
{
    const float RECUR_BIG = SLC_SQRT_SP_MAX;
    const int maxiter = 5000;
    int   n = 1;
    float Anm2 = 1.0f;
    float Bnm2 = 0.0f;
    float Anm1 = 0.0f;
    float Bnm1 = 1.0f;
    float a1 = 1.0f;
    float b1 = 1.0f;
    float a2 = -x;
    float b2 = N + 1.0f;
    float an, bn;

    float fn;

    float An = b1 * Anm1 + a1 * Anm2;   // A1
    float Bn = b1 * Bnm1 + a1 * Bnm2;   // B1

    // * One explicit step, before we get to the main pattern
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    An = b2 * Anm1 + a2 * Anm2;   // A2
    Bn = b2 * Bnm1 + a2 * Bnm2;   // B2

    fn = An / Bn;

    while( n < maxiter ) {
        float old_fn;
        float del;
        n++;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        an = SLC_ODD( n )? (( n - 1 )/ 2 ) * x: -( N + ( n / 2 ) - 1 ) * x;
        bn = N + n - 1;
        An = bn * Anm1 + an * Anm2;
        Bn = bn * Bnm1 + an * Bnm2;

        if( RECUR_BIG < fabsf( An ) || RECUR_BIG < fabsf( Bn )) {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
            Bnm2 /= RECUR_BIG;
        }

        old_fn = fn;
        fn = An / Bn;
        del = old_fn / fn;

        if( fabsf( del - 1.0f ) < 2.0f * SLC_SP_EPSILON )
            break;
    }

    if( result ) {
        *result = fn;
        if( err )
            *err = 4.0f *( n + 1.0f )* SLC_SP_EPSILON * fabsf( fn );
    }
    if( n == maxiter )
        return PSL_MAXITERATS;
    return PSL_SUCCESS;
}


// --------------------- Functions with Error Codes ------------------------
// exp(x)
//
int psl_exp_e( const float x, float* result, float* err )
{
    if( SLC_LOG_SP_MAX < x ) {
        return PSL_ERR_OVERFLOW;
    }
    else if( x < SLC_LOG_SP_MIN ) {
        return PSL_ERR_UNDERFLOW;
    }
    else {
        if( result ) {
            *result = expf( x );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
}

// -------------------------------------------------------------------------
// y exp(x)
//
int psl_exp_mult_e( const float x, const float y, float* result, float* err )
{
    const float ay = fabsf( y );

    if( y == 0.0f ) {
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if(( x < 0.5f * SLC_LOG_SP_MAX   &&  0.5f * SLC_LOG_SP_MIN < x ) &&
           ( ay < 0.8f * SLC_SQRT_SP_MAX  &&  1.2f * SLC_SQRT_SP_MIN < ay )) {
        const float ex = expf( x );
        if( result ) {
            *result = y * ex;
            if( err )
                *err = ( 2.0f + fabsf( x )) * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        const float ly  = logf( ay );
        const float lnr = x + ly;

        if( SLC_LOG_SP_MAX - 0.01f < lnr ) {
            return PSL_ERR_OVERFLOW;
        }
        else if( lnr < SLC_LOG_SP_MIN + 0.01f ) {
            return PSL_ERR_UNDERFLOW;
        }
        else {
            const float sy   = (float)SLC_SIGN( y );
            const float M    = floorf( x );
            const float N    = floorf( ly );
            const float a    = x  - M;
            const float b    = ly - N;
            const float berr = 2.0f * SLC_SP_EPSILON * ( fabsf( ly ) + fabsf( N ));
            if( result ) {
                *result = sy * expf( M + N ) * expf( a + b );
                if( err ) {
                    *err = berr * fabsf( *result );
                    *err += 2.0f * SLC_SP_EPSILON * ( M + N + 1.0f ) * fabsf( *result );
                }
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// y exp(x) with add. error estimate
//
int psl_exp_mult_err_e( const float x, const float dx, const float y, const float dy, float* result, float* err )
{
    const float ay = fabsf(y);

    if( y == 0.0f ) {
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = fabsf( dy * expf( x ));
        }
        return PSL_SUCCESS;
    }
    else if(( x < 0.5f * SLC_LOG_SP_MAX   &&  0.5f * SLC_LOG_SP_MIN < x ) &&
           ( ay < 0.8f * SLC_SQRT_SP_MAX  &&  1.2f * SLC_SQRT_SP_MIN < ay )) {
        float ex = expf( x );
        if( result ) {
            *result = y * ex;
            if( err ) {
                *err = ex * ( fabsf( dy ) + fabsf( y * dx ));
                *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
        }
        return PSL_SUCCESS;
    }
    else {
        const float ly  = logf( ay );
        const float lnr = x + ly;

        if( SLC_LOG_SP_MAX - 0.01f < lnr ) {
            return PSL_ERR_OVERFLOW;
        }
        else if( lnr < SLC_LOG_SP_MIN + 0.01f ) {
            return PSL_ERR_UNDERFLOW;
        }
        else {
            const float sy  = (float)SLC_SIGN( y );
            const float M   = floorf( x );
            const float N   = floorf( ly );
            const float a   = x  - M;
            const float b   = ly - N;
            const float eMN = expf( M + N );
            const float eab = expf( a + b );
            if( result ) {
                *result = sy * eMN * eab;
                if( err ) {
                    *err = eMN * eab * 2.0f * SLC_SP_EPSILON;
                    *err += eMN * eab * fabsf( dy / y );
                    *err += eMN * eab * fabsf( dx );
                }
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// exp(x)-1; accurate for small x
//
int psl_expm1_e( const float x, float* result, float* err )
{
    const float cut = 0.002f;

    if( x < SLC_LOG_SP_MIN ) {
        if( result ) {
            *result = -1.0f;
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = expf( x ) - 1.0f;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = x *( 1.0f + 0.5f * x *( 1.0f + x / 3.0f *( 1.0f + 0.25f * x *( 1.0f + 0.2f * x ))));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    } 
    else if( x < SLC_LOG_SP_MAX ) {
        if( result ) {
            *result = expf( x ) - 1.0f;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// (exp(x)-1)/x; accurate for small x
//
int psl_exprel_e( const float x, float* result, float* err )
{
    const float cut = 0.002f;

    if( x < SLC_LOG_SP_MIN ) {
        if( result ) {
            *result = -1.0f / x;
            if( err )
                *err = SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = ( expf( x ) - 1.0f )/ x;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = ( 1.0f + 0.5f * x *( 1.0f + x / 3.0f *( 1.0f + 0.25f * x *( 1.0f + 0.2f * x ))));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    } 
    else if( x < SLC_LOG_SP_MAX ) {
        if( result ) {
            *result = ( expf( x ) - 1.0f )/ x;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// 2(exp(x)-1-x)/x^2; accurate for small x
//
int psl_exprel_2_e( float x, float* result, float* err )
{
    const float cut = 0.002f;

    if( x < SLC_LOG_SP_MIN ) {
        if( result ) {
            *result = -2.0f / x * ( 1.0f + 1.0f / x );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = 2.0f *( expf( x ) - 1.0f - x )/( x * x );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = ( 1.0f + 1.0f / 3.0f * x *( 1.0f + 0.25f * x *( 1.0f + 0.2f * x *( 1.0f + 1.0f / 6.0f * x ))));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < SLC_LOG_SP_MAX ) {
        if( result ) {
            *result = 2.0f *( expf( x ) - 1.0f - x )/( x * x );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// (exp(x)-1)/x; using continued fraction representation
//
int psl_exprel_n_CF_e( const float N, const float x, float* result, float* err )
{
    return exprel_n_CF( N, x, result, err );
}


// -------------------------------------------------------------------------
// N-relative exponential, n-th generalization of  exprel and exprel_2,
// N!/x^N (exp(x)-SUM(x^k/k!))
//
int psl_exprel_n_e( const int N, const float x, float* result, float* err )
{
    if( N < 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( x == 0.0f ) {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( fabsf( x ) < SLC_ROOT3_SP_EPSILON * N ) {
        if( result ) {
            *result = 1.0f + x /( N + 1 ) * ( 1.0f + x /( N + 2 ));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( N == 0 ) {
        return psl_exp_e( x, result, err );
    }
    else if( N == 1 ) {
        return psl_exprel_e( x, result, err );
    }
    else if( N == 2 ) {
        return psl_exprel_2_e( x, result, err );
    }
    else {
        if( N < x && ( -x + N *( 1.0f + logf( x / N )) < SLC_LOG_SP_EPSILON )) {
            // * x is much larger than N; ignore polynomial part, so
            // * exprel_N(x) ~= e^x N!/x^N
            float lnf_N, lnf_Nerr;
            float lnr_val;
            float lnr_err;
            float lnterm;
            psl_lnfact_e( N, &lnf_N, &lnf_Nerr );
            lnterm =  N * logf( x );
            lnr_val = x + lnf_N - lnterm;
            lnr_err = SLC_SP_EPSILON * ( fabsf( x ) + fabsf( lnf_N ) + fabsf( lnterm ));
            lnr_err += lnf_Nerr;
            return psl_exp_err_e( lnr_val, lnr_err, result, err );
        }
        else if( N < x ) {
            // * Write the identity
            // *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
            // * then use the asymptotic expansion
            // * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
            float ln_x = logf( x );
            float lnf_N, lnf_Nerr;
            float lg_N;
            float lnpre_val;
            float lnpre_err;
            psl_lnfact_e( N, &lnf_N, &lnf_Nerr ); // log(N!)
            lg_N = lnf_N - logf((float)N); // log(Gamma(N))
            lnpre_val  = x + lnf_N - N * ln_x;
            lnpre_err  = SLC_SP_EPSILON * ( fabsf( x ) + fabsf( lnf_N ) + fabsf( N * ln_x ));
            lnpre_err += lnf_Nerr;
            if( lnpre_val < SLC_LOG_SP_MAX - 5.0f ) {
                int eGstatus;
                float bigG_ratio, bigG_ratio_err;
                float pre = 0.0f, prerr = 0.0f;
                int exstatus = psl_exp_err_e( lnpre_val, lnpre_err, &pre, &prerr );
                float ln_bigG_ratio_pre = -x + ( N - 1 )* ln_x - lg_N;
                float bigGsum = 1.0f;
                float term = 1.0f;
                int k;
                for( k = 1; k < N; k++ ) {
                    term *= ( N - k )/ x;
                    bigGsum += term;
                }
                eGstatus = psl_exp_mult_e( ln_bigG_ratio_pre, bigGsum, &bigG_ratio, &bigG_ratio_err );
                if( eGstatus == PSL_SUCCESS ) {
                    if( result ) {
                        *result = pre * ( 1.0f - bigG_ratio );
                        if( err ) {
                            *err = pre * ( 2.0f * SLC_SP_EPSILON + bigG_ratio_err );
                            *err += prerr * fabsf( 1.0f - bigG_ratio );
                            *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
                        }
                    }
                    return exstatus;
                }
                else {
                    if( result ) {
                        *result = 0.0f;
                        if( err )
                            *err = 0.0f;
                    }
                    return eGstatus;
                }
            }
            else {
                return PSL_ERR_OVERFLOW;
            }
        }
        else if( -10.0f * N < x ) {
            return exprel_n_CF((float)N, x, result, err );
        }
        else {
            // * x -> -Inf asymptotic:
            // * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            // *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            float sum  = 1.0f;
            float term = 1.0f;
            int k;
            for( k = 1; k < N; k++ ) {
                term *= ( N - k )/ x;
                sum += term;
            }
            if( result ) {
                *result = -N / x * sum;
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// exp(x) with add. error estimate
//
int psl_exp_err_e( const float x, const float dx, float* result, float* err )
{
    const float adx = fabsf( dx );

    if( SLC_LOG_SP_MAX < x + adx ) {
        return PSL_ERR_OVERFLOW;
    }
    else if( x - adx < SLC_LOG_SP_MIN ) {
        return PSL_ERR_UNDERFLOW;
    }
    else {
        const float ex  = expf( x );
        const float edx = expf( adx );
        if( result ) {
            *result = ex;
            if( err ) {
                *err = ex * SLC_MAX( SLC_SP_EPSILON, edx - 1.0f / edx );
                *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
        }
        return PSL_SUCCESS;
    }
}

// -------------------------------------------------------------------------

}//namespace extspsl
