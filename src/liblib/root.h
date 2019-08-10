/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __root_h__
#define __root_h__

// typedefs
typedef void (*TFdFfunc)( float, float*, float*, void* );

// interface
//
int mygcd( int, int );//the greatest common divisor (PC_Gcd)
float LecuyerRand( long* );//L'Ecuyer's random number generator
const char* root_by_NR_and_bisection(//root finding by Newton-Raphson and bisection methods
        TFdFfunc fdfunction,
        float x1,
        float x2,
        float xacc,
        int maxit,
        void* params,
        float* result
);

#endif//__root_h__
