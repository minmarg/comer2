/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logitnormal.h"

// -------------------------------------------------------------------------
// LogitNormalErrCorr: correct error of logit-normal variable so that it 
//  sums exactly to 1;
// NOTE: result is written back to `f'
//
void LogitNormalErrCorr( double* f, size_t fsz, double tol )
{
    double          vv, vsum, corr;
    size_t          r, non;
    char            strbuf[1024];
    const double    accuracy = 1.0e-6;
    mystring errmsg = "LogitNormalErrCorr: ";

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw MYRUNTIME_ERROR( errmsg + "Invalid dimensionality.");

    vsum = 0.0;
    non = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0 )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal value.");
        if( vv != 0.0 )
            non++;
        vsum += vv;
    }

    corr = vsum - 1.0;
    if( corr == 0.0 )
        return;

    if( tol < corr ) {
        sprintf( strbuf, "Error tolerance exceeded: %g.", corr );
        throw MYRUNTIME_ERROR( errmsg + strbuf );
    }

    if( fabs( corr ) <= accuracy ) {
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            //correct small error using the first value met
            if(( 0.0 < corr && corr <= vv ) ||
               ( corr < 0.0 && (( non < fsz && vv == 0.0 ) || 
                                ( fsz <= non && vv - corr <= 1.0 ))))
            {
                vv -= corr;
                f[r] = vv;//write
                break;
            }
        }
    }
    else if( accuracy < fabs( corr ) && non ) {
        corr /= (double)non;
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            if( vv == 0.0 )
                continue;
            vv -= corr;
            non--;
            //correct overflow/undeflow error
            if( vv < 0.0 ) {
                corr -= vv / (double)non;
                vv = 0.0;
            }
            else if( 1.0 < vv ) {
                corr += ( 1.0 - vv )/(double)non;
                vv = 1.0;
            }
            f[r] = vv;//write
        }
    }
}

// -------------------------------------------------------------------------
// LogitNormal2Normal: tranform logit-normally distributed variable to 
//      normally distributed variable;
//  NOTE: tranformed variable is written back to `f' and 
//  NOTE: dimensionality is reduced by 1
//
void LogitNormal2Normal( double* f, size_t fsz, double tol, bool correct )
{
    double          vv, vlast, vsum, corr;
    size_t          r, noz, non;
    double          fake = 1.0e-4;//fake value
    const double    accuracy = tol;//1.0e-6;
    mystring errmsg = "LogitNormal2Normal: ";

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw MYRUNTIME_ERROR( errmsg + "Invalid dimensionality.");

    for( r = 0, vsum = 0.0; r < fsz; r++ )
        vsum += f[r];
    if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum ) {
        if( correct )
            LogitNormalErrCorr( f, fsz, tol );
        else
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable.");
    }

    vsum = 0.0;
    non = noz = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0 )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal value.");
        if( vv == 0.0 ) {
            noz++;
            continue;
        }
        non++;
        vsum += vv;
    }

    if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum )
        throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable.");
    if( non < 1 || noz + non != fsz )
        throw MYRUNTIME_ERROR( errmsg + "Logit-normal variable contains invalid values.");

    vlast = f[fsz-1];

    if( correct ) {
        //make sure there are no null frequencies
        corr = fake * (double)noz /(double)non;
        if( 0.0 < corr ) {
            vsum = 0.0;
            for( r = 0; r < fsz; r++ ) {
                vv = f[r];
                if( vv == 0.0 ) {
                    vsum += fake;
                    //null frequency, add fake value
                    f[r] = fake;//write
                    continue;
                }
                vv -= corr;
                non--;
                if( vv <= 0.0 ) {
                    vsum += fake;
                    //less than correction value; recalc. `corr'
                    f[r] = fake;//write
                    corr += ( fake - vv ) /(double)non;
                    continue;
                }
                vsum += vv;
                f[r] = vv;//write
            }
        }
        if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum || vlast <= 0.0 )
            throw MYRUNTIME_ERROR( errmsg + "Correction failed: Invalid logit-normal variable.");
    }
    else
        if( 0 < noz )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable: Zero values.");


    //apply transformation to normal variable
    for( r = 0; r < fsz - 1; r++ ) {
        vv = f[r];
        if( vv <= 0.0 )
            throw MYRUNTIME_ERROR( errmsg + "Transformation failed: Invalid logit-normal value.");
        vv = log( vv / vlast );
        f[r] = vv;//write
    }
    f[fsz-1] = 0.0;
}

// -------------------------------------------------------------------------
// Normal2LogitNormal: apply inverse logistic transformation to tranform 
//      normal random variable to logistic-normal random variable;
//  NOTE: result is written to lnv whose size lnvsz should be 1+ the 
//  size nvsz of the source nv
//
void Normal2LogitNormal( const double* nv, size_t nvsz, double* lnv, size_t lnvsz )
{
    if( nv == NULL || lnv == NULL )
        throw MYRUNTIME_ERROR("Normal2LogitNormal: Null arguments.");
    if( nvsz < 1 || nvsz+1 != lnvsz )
        throw MYRUNTIME_ERROR("Normal2LogitNormal: Invalid sizes.");

    size_t n;
    double val, sum, tsum;
    sum = 0.0;
    for( n = 0; n < nvsz; n++ ) {
        val = exp( nv[n]);
        lnv[n] = val;
        sum += val;
    }
    tsum = 0.0;
    for( n = 0; n < nvsz; n++ ) {
        lnv[n] /= 1.0 + sum;
        tsum += lnv[n];
    }
    lnv[n] = 1.0 - tsum;
}

// =========================================================================
// Float-type functional equivalents ---------------------------------------
// -------------------------------------------------------------------------
// LogitNormalErrCorr_f: correct error of logit-normal variable so that it 
//  sums exactly to 1;
// NOTE: result is written back to `f'
//
void LogitNormalErrCorr_f( float* f, size_t fsz, float tol )
{
    float vv, vsum, corr;
    size_t r, non;
    char strbuf[1024];
    const float accuracy = 1.0e-6f;
    mystring errmsg = "LogitNormalErrCorr_f: ";

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw MYRUNTIME_ERROR( errmsg + "Invalid dimensionality.");

    vsum = 0.0f;
    non = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0f )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal value.");
        if( vv != 0.0f )
            non++;
        vsum += vv;
    }

    corr = vsum - 1.0f;
    if( corr == 0.0f )
        return;

    if( tol < corr ) {
        sprintf( strbuf, "Error tolerance exceeded: %g.", corr );
        throw MYRUNTIME_ERROR( errmsg + strbuf );
    }

    if( fabsf( corr ) <= accuracy ) {
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            //correct small error using the first value met
            if(( 0.0 < corr && corr <= vv ) ||
               ( corr < 0.0f && (( non < fsz && vv == 0.0f ) || 
                                ( fsz <= non && vv - corr <= 1.0f ))))
            {
                vv -= corr;
                f[r] = vv;//write
                break;
            }
        }
    }
    else if( accuracy < fabsf( corr ) && non ) {
        corr /= (float)non;
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            if( vv == 0.0f )
                continue;
            vv -= corr;
            non--;
            //correct overflow/undeflow error
            if( vv < 0.0f ) {
                corr -= vv / (float)non;
                vv = 0.0f;
            }
            else if( 1.0f < vv ) {
                corr += ( 1.0f - vv )/(float)non;
                vv = 1.0f;
            }
            f[r] = vv;//write
        }
    }
}

// -------------------------------------------------------------------------
// LogitNormal2Normal_f: tranform logit-normally distributed variable to 
//      normally distributed variable;
//  NOTE: tranformed variable is written back to `f' and 
//  NOTE: dimensionality is reduced by 1
//
void LogitNormal2Normal_f( float* f, size_t fsz, float tol, bool correct )
{
    float vv, vlast, vsum, corr;
    size_t r, noz, non;
    float fake = 1.0e-4f;//fake value
    const float accuracy = tol;//1.0e-6f;
    mystring errmsg = "LogitNormal2Normal_f: ";

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw MYRUNTIME_ERROR( errmsg + "Invalid dimensionality.");

    for( r = 0, vsum = 0.0f; r < fsz; r++ )
        vsum += f[r];
    if( vsum < 1.0f - accuracy || 1.0f + accuracy < vsum ) {
        if( correct )
            LogitNormalErrCorr_f( f, fsz, tol );
        else
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable.");
    }

    vsum = 0.0f;
    non = noz = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0f )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal value.");
        if( vv == 0.0f ) {
            noz++;
            continue;
        }
        non++;
        vsum += vv;
    }

    if( vsum < 1.0f - accuracy || 1.0f + accuracy < vsum )
        throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable.");
    if( non < 1 || noz + non != fsz )
        throw MYRUNTIME_ERROR( errmsg + "Logit-normal variable contains invalid values.");

    vlast = f[fsz-1];

    if( correct ) {
        //make sure there are no null frequencies
        corr = fake * (float)noz /(float)non;
        if( 0.0f < corr ) {
            vsum = 0.0f;
            for( r = 0; r < fsz; r++ ) {
                vv = f[r];
                if( vv == 0.0f ) {
                    vsum += fake;
                    //null frequency, add fake value
                    f[r] = fake;//write
                    continue;
                }
                vv -= corr;
                non--;
                if( vv <= 0.0f ) {
                    vsum += fake;
                    //less than correction value; recalc. `corr'
                    f[r] = fake;//write
                    corr += ( fake - vv ) /(float)non;
                    continue;
                }
                vsum += vv;
                f[r] = vv;//write
            }
        }
        if( vsum < 1.0f - accuracy || 1.0f + accuracy < vsum || vlast <= 0.0f )
            throw MYRUNTIME_ERROR( errmsg + "Correction failed: Invalid logit-normal variable.");
    }
    else
        if( 0 < noz )
            throw MYRUNTIME_ERROR( errmsg + "Invalid logit-normal variable: Zero values.");


    //apply transformation to normal variable
    for( r = 0; r < fsz - 1; r++ ) {
        vv = f[r];
        if( vv <= 0.0f )
            throw MYRUNTIME_ERROR( errmsg + "Transformation failed: Invalid logit-normal value.");
        vv = logf( vv / vlast );
        f[r] = vv;//write
    }
    f[fsz-1] = 0.0f;
}

// -------------------------------------------------------------------------
// Normal2LogitNormal_f: apply inverse logistic transformation to tranform 
//      normal random variable to logistic-normal random variable;
//  NOTE: result is written to lnv whose size lnvsz should be 1+ the 
//  size nvsz of the source nv
//
void Normal2LogitNormal_f( const float* nv, size_t nvsz, float* lnv, size_t lnvsz )
{
    if( nv == NULL || lnv == NULL )
        throw MYRUNTIME_ERROR("Normal2LogitNormal_f: Null arguments.");
    if( nvsz < 1 || nvsz+1 != lnvsz )
        throw MYRUNTIME_ERROR("Normal2LogitNormal_f: Invalid sizes.");

    size_t n;
    float val, sum, tsum;
    sum = 0.0f;
    for( n = 0; n < nvsz; n++ ) {
        val = expf( nv[n]);
        lnv[n] = val;
        sum += val;
    }
    tsum = 0.0f;
    for( n = 0; n < nvsz; n++ ) {
        lnv[n] /= 1.0f + sum;
        tsum += lnv[n];
    }
    lnv[n] = 1.0f - tsum;
}
