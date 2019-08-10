/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __myfastapproxmath_h__
#define __myfastapproxmath_h__

// myfastapexpf: fast approximate calculation of the exponential function by
// Schraudolph's algorithm (Neural Computation 11, 1999)
//
__device__ __forceinline__ 
float myfastapexpf( float arg )
{
    union { float f; int i; } spr;
    //i = 2^23/ln(2) * arg + 2^23*127 - C,
    // where C is adjustment constant, 2^23*127=1065353216;
    // C=486411 is an optimal constant found by Schraudolph using the RMS criterion 
    // (C=60801.48) and multiplied by 2^3 to adjust for single precision exponent
    //spr.i = (int)(12102203.16f * arg + 1064866805.0f);
//     //NOTE: constant C=-1*2^3 for the max<e^y error to be 0%:
//     //spr.i = (int)(12102203.16f * arg + 1065353224.0f);
    //NOTE: constant C=90252.34*2^3 for the max>e^y error to be 0%
    spr.i = (int)(12102203.16f * arg + 1064631197.0f);
    return spr.f;
}

// myfastaplogf: inverse of Schraudolph's algorithm to rapidly calculate log
//
__device__ __forceinline__ 
float myfastaplogf( float arg )
{
    union { float f; int i; } spr = {arg};
    //multiply by 1/12102203.16
    //return (spr.i - 1064866805) * 8.2629582959339e-08f;
//     //NOTE: the case of constant C=-1*2^3 for the max<e^y error==0%
//     //return (spr.i - 1065353224) * 8.2629582959339e-08f;
    //NOTE: the case of constant C=90252.34*2^3 for the max>e^y error==0%
    return (spr.i - 1064631197) * 8.2629582959339e-08f;
}

// myNJfastapexpf: fast approximate calculation of the exponential function by
// Njuffa; more accurate but much slower than Schraudolph's algorithm 
//
__device__ __forceinline__ 
float myNJfastapexpf(float x)
 {
    union { float f; int i; } spr;
    //exp(x) = 2^i * 2^f; i = floor(log2(e) * x), 0<=f<=1
    float t = x * 1.442695041f;
    float fi = floorf(t);
    float f = t - fi;
    int i = (int)fi;
    //compute 2^f by Taylor's series:
    spr.f = (0.3371894346f * f + 0.657636276f) * f + 1.00172476f;
    spr.i += (i << 23);//scale by 2^i
    return spr.f;
 }

#endif//__myfastapproxmath_h__
