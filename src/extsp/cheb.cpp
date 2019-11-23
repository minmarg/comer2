/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

/* Author: G. Jungman */
// *** the code adopted from GSL ***

// #include <math.h>
#include <cmath>
#include "psl.h"
#include "pslerror.h"
#include "cheb.h"

namespace extspsl {

// -------------------------------------------------------------------------
//
int cheb_eval_e( const Tcheb_series* cs, const float x, float* result, float* err )
{
    int j;
    float d  = 0.0f;
    float dd = 0.0f;
    float y  = ( 2.0f * x - cs->a - cs->b )/( cs->b - cs->a );
    float y2 = 2.0f * y;
    float e = 0.0f;
    float temp;

    for( j = cs->order; j >= 1; j-- ) {
        temp = d;
        d = y2 * d - dd + cs->c[j];
        e += fabsf( y2 * temp ) + fabsf( dd ) + fabsf( cs->c[j] );
        dd = temp;
    }

        temp = d;
        d = y * d - dd + 0.5f * cs->c[0];
        e += fabsf( y * temp ) + fabsf( dd ) + 0.5f * fabsf( cs->c[0] );

    if( result )
        *result = d;
    if( err )
        *err = SLC_SP_EPSILON * e + fabsf( cs->c[cs->order] );

    return PSL_OK;
}

// -------------------------------------------------------------------------
//
int psl_multiply_e( const float x, const float y, float* result, float* err )
{
    const float ax = fabsf( x );
    const float ay = fabsf( y );

    if( x == 0.0f || y == 0.0f ) {
        // It is necessary to eliminate this immediately
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if(( ax <= 1.0f && 1.0f <= ay ) || ( ay <= 1.0f && 1.0f <= ax )) {
        // Straddling 1.0 is always safe
        if( result ) {
            *result = x * y;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        const float f = 1.0f - 2.0f * SLC_SP_EPSILON;
        const float min = SLC_MIN( fabsf( x ), fabsf( y ));
        const float max = SLC_MAX( fabsf( x ), fabsf( y ));
        if( max < 0.9f * SLC_SQRT_SP_MAX || min < ( f * SLC_SP_MAX )/ max) {
            if( result ) {
                *result = ( x * y );
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
                if( fabsf( *result ) < SLC_SP_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }
        else {
            return PSL_ERR_OVERFLOW;
        }
    }
}


// -------------------------------------------------------------------------
//
int psl_multiply_err_e( const float x, const float dx, const float y, const float dy, float* result, float* err )
{
    int status = psl_multiply_e( x, y, result, err );
    if( err )
        *err += fabsf( dx * y ) + fabsf( dy * x );
    return status;
}

}//namespace extspsl
