/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_psl__
#define __extspsl_psl__

#include <limits.h>

// single-precision library
namespace extspsl {

#define SLC_E           2.71828182845904523536028747135f

//pi, sqrt(pi), ln(pi), log(2*pi), log(sqrt(pi)), and ln(sqrt(2*pi))
#define SLC_PI          3.14159265358979323846264338328f
#define SLC_SQRTPI      1.77245385090551602729816748334f
#define SLC_LNPI        1.14472988584940017414342735135f
#define SLC_LN2PI       1.83787706640934548356065947281f
#define SLC_LNSQRTPI    0.57236494292470008706f
#define SLC_LNSQRT2PI   0.9189385332046727418f

// Euler constant
#define SLC_EULER       0.57721566490153286060651209008f

//sqrt(2), and ln(2)
#define SLC_SQRT2       1.41421356237309504880168872421f
#define SLC_LN2         0.69314718055994530941723212146f
#define SLC_LN1K        6.90775527898213681510242167860f

// machine dependent constants
#define SLC_SP_EPSILON          ( 2.3841858e-07f )
#define SLC_SQRT_SP_EPSILON     ( 0.00048828125f )
#define SLC_ROOT3_SP_EPSILON    ( 0.0062007854f )
#define SLC_ROOT4_SP_EPSILON    ( 0.022097087f )
#define SLC_ROOT5_SP_EPSILON    ( 0.047366143f )
#define SLC_ROOT6_SP_EPSILON    ( 0.078745066f )
#define SLC_LOG_SP_EPSILON      (-1.5249238e+01f )

#define SLC_SP_MIN              ( 1.1754944e-38f )
#define SLC_SQRT_SP_MIN         ( 1.0842022e-19f )
#define SLC_LOG_SP_MIN          (-8.7336545e+01f )

#define SLC_SP_MAX              ( 3.4028234e+38f )
#define SLC_SQRT_SP_MAX         ( 1.8446742e+19f )
#define SLC_LOG_SP_MAX          ( 8.8722839e+01f )

#define SLC_LOG_DP_MAX          ( 7.0978271289338397e+02f )
//

// functions
#define SLC_POW2( X )   (( X )*( X ))
#define SLC_POW3( X )   (( X )*( X )*( X ))
#define SLC_POW4( X )   ( SLC_POW2( X ) * SLC_POW2( X ))
#define SLC_POW5( X )   ( SLC_POW3( X ) * SLC_POW2( X ))
#define SLC_POW6( X )   ( SLC_POW3( X ) * SLC_POW3( X ))

#define SLC_ODD( X )    (( X ) & 1 )
#define SLC_EVEN( X )   (!SLC_ODD( X ))
#define SLC_SIGN( X )   (( 0 <= ( X ))? 1: -1 )

#define SLC_MIN( X, Y ) ((( X ) < ( Y ))? ( X ): ( Y ))
#define SLC_MAX( X, Y ) ((( X ) < ( Y ))? ( Y ): ( X ))
#define SLC_MAX3( X1, X2, X3 )  SLC_MAX( X1, SLC_MAX( X2, X3 ))
#define SLC_MAX4( X1, X2, X3, X4 )  SLC_MAX( X1, SLC_MAX3( X2, X3, X4 ))
#define SLC_MAX5( X1, X2, X3, X4, X5 )  SLC_MAX( X1, SLC_MAX4( X2, X3, X4, X5 ))

}//namespace extspsl

#endif//__extspsl_psl__
