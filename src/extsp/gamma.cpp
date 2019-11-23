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

// #include <math.h>
#include <cmath>
#include "psl.h"
#include "pslerror.h"
#include "cheb.h"
#include "digamma.h"
#include "exp.h"
#include "gamma.h"

// * The maximum x such that gamma(x) is not considered an overflow
#define SLC_G_GAMMA_XMAX        ( 34.0f )//     ( 171.0 )

//* The maximum n such that psl_fact(n) does not give an overflow
#define SLC_G_FACT_NMAX         ( 33 )//        ( 170 )

//* The maximum n such that psl_doublefact(n) does not give an overflow
#define SLC_G_DOUBLEFACT_NMAX   ( 55 )//        ( 297 )

namespace extspsl {

// -------------------------- Private Section ------------------------------
//
static float fact_table[SLC_G_FACT_NMAX+1] = {
//   f                                           //       n
// -------------------------------------------           ----
    1.0f,                                        //       0,
    1.0f,                                        //       1,
    2.0f,                                        //       2,
    6.0f,                                        //       3,
    24.0f,                                       //       4,
    120.0f,                                      //       5,
    720.0f,                                      //       6,
    5040.0f,                                     //       7,
    40320.0f,                                    //       8,

    362880.0f,                                   //       9,
    3628800.0f,                                  //       10,
    39916800.0f,                                 //       11,
    479001600.0f,                                //       12,

    6227020800.0f,                               //       13,
    87178291200.0f,                              //       14,
    1307674368000.0f,                            //       15,
    20922789888000.0f,                           //       16,
    355687428096000.0f,                          //       17,
    6402373705728000.0f,                         //       18,
    121645100408832000.0f,                       //       19,
    2432902008176640000.0f,                      //       20,
    51090942171709440000.0f,                     //       21,
    1124000727777607680000.0f,                   //       22,
    25852016738884976640000.0f,                  //       23,
    620448401733239439360000.0f,                 //       24,
    15511210043330985984000000.0f,               //       25,
    403291461126605635584000000.0f,              //       26,
    10888869450418352160768000000.0f,            //       27,
    304888344611713860501504000000.0f,           //       28,
    8841761993739701954543616000000.0f,          //       29,
    265252859812191058636308480000000.0f,        //       30,
    8222838654177922817725562880000000.0f,       //       31,
    263130836933693530167218012160000000.0f,     //       32,
    8683317618811886495518194401280000000.0f     //       33,
};

static float doub_fact_table[SLC_G_DOUBLEFACT_NMAX+1] = {
//   f                                           //       n
// -------------------------------------------           ----
    1.000000000000000000000000000f,              //       0,
    1.000000000000000000000000000f,              //       1,
    2.000000000000000000000000000f,              //       2,
    3.000000000000000000000000000f,              //       3,
    8.000000000000000000000000000f,              //       4,
    15.00000000000000000000000000f,              //       5,
    48.00000000000000000000000000f,              //       6,
    105.0000000000000000000000000f,              //       7,
    384.0000000000000000000000000f,              //       8,
    945.0000000000000000000000000f,              //       9,
    3840.000000000000000000000000f,              //       10,
    10395.00000000000000000000000f,              //       11,
    46080.00000000000000000000000f,              //       12,
    135135.0000000000000000000000f,              //       13,
    645120.00000000000000000000000f,             //       14,
    2.02702500000000000000000000000e6f,          //       15,
    1.03219200000000000000000000000e7f,          //       16,
    3.4459425000000000000000000000e7f,           //       17,
    1.85794560000000000000000000000e8f,          //       18,
    6.5472907500000000000000000000e8f,           //       19,
    3.7158912000000000000000000000e9f,           //       20,
    1.37493105750000000000000000000e10f,         //       21,
    8.1749606400000000000000000000e10f,          //       22,
    3.1623414322500000000000000000e11f,          //       23,
    1.96199055360000000000000000000e12f,         //       24,
    7.9058535806250000000000000000e12f,          //       25,
    5.1011754393600000000000000000e13f,          //       26,
    2.13458046676875000000000000000e14f,         //       27,
    1.42832912302080000000000000000e15f,         //       28,
    6.1902833536293750000000000000e15f,          //       29,
    4.2849873690624000000000000000e16f,          //       30,
    1.91898783962510625000000000000e17f,         //       31,
    1.37119595809996800000000000000e18f,         //       32,
    6.3326598707628506250000000000e18f,          //       33,
    4.6620662575398912000000000000e19f,          //       34,
    2.21643095476699771875000000000e20f,         //       35,
    1.67834385271436083200000000000e21f,         //       36,
    8.2007945326378915593750000000e21f,          //       37,
    6.3777066403145711616000000000e22f,          //       38,
    3.1983098677287777081562500000e23f,          //       39,
    2.55108265612582846464000000000e24f,         //       40,
    1.31130704576879886034406250000e25f,         //       41,
    1.07145471557284795514880000000e26f,         //       42,
    5.6386202968058350994794687500e26f,          //       43,
    4.7144007485205310026547200000e27f,          //       44,
    2.53737913356262579476576093750e28f,         //       45,
    2.16862434431944426122117120000e29f,         //       46,
    1.19256819277443412353990764062e30f,         //       47,
    1.04093968527333324538616217600e31f,         //       48,
    5.8435841445947272053455474391e31f,          //       49,
    5.2046984263666662269308108800e32f,          //       50,
    2.98022791374331087472622919392e33f,         //       51,
    2.70644318171066643800402165760e34f,         //       52,
    1.57952079428395476360490147278e35f,         //       53,
    1.46147931812375987652217169510e36f,         //       54,
    8.6873643685617511998269581003e36f           //       55,
};


// * Chebyshev coefficients for Gamma*(3/4(t+1)+1/2), -1<t<1
//
static float gstar_a_data[30] = {
    2.16786447866463034423060819465f,
   -0.05533249018745584258035832802f,
    0.01800392431460719960888319748f,
   -0.00580919269468937714480019814f,
    0.00186523689488400339978881560f,
   -0.00059746524113955531852595159f,
    0.00019125169907783353925426722f,
   -0.00006124996546944685735909697f,
    0.00001963889633130842586440945f,
   -6.3067741254637180272515795142e-06f,
    2.0288698405861392526872789863e-06f,
   -6.5384896660838465981983750582e-07f,
    2.1108698058908865476480734911e-07f,
   -6.8260714912274941677892994580e-08f,
    2.2108560875880560555583978510e-08f,
   -7.1710331930255456643627187187e-09f,
    2.3290892983985406754602564745e-09f,
   -7.5740371598505586754890405359e-10f,
    2.4658267222594334398525312084e-10f,
   -8.0362243171659883803428749516e-11f,
    2.6215616826341594653521346229e-11f,
   -8.5596155025948750540420068109e-12f,
    2.7970831499487963614315315444e-12f,
   -9.1471771211886202805502562414e-13f,
    2.9934720198063397094916415927e-13f,
   -9.8026575909753445931073620469e-14f,
    3.2116773667767153777571410671e-14f,
   -1.0518035333878147029650507254e-14f,
    3.4144405720185253938994854173e-15f,
   -1.0115153943081187052322643819e-15f
};
static Tcheb_series gstar_a_cs = {
    gstar_a_data,
    29,
    -1.f, 1.f,
    17
};


// * Chebyshev coefficients for
// * x^2(Gamma*(x) - 1 - 1/(12x)), x = 4(t+1)+2, -1 < t < 1
static float gstar_b_data[] = {
    0.0057502277273114339831606096782f,
    0.0004496689534965685038254147807f,
   -0.0001672763153188717308905047405f,
    0.0000615137014913154794776670946f,
   -0.0000223726551711525016380862195f,
    8.0507405356647954540694800545e-06f,
   -2.8671077107583395569766746448e-06f,
    1.0106727053742747568362254106e-06f,
   -3.5265558477595061262310873482e-07f,
    1.2179216046419401193247254591e-07f,
   -4.1619640180795366971160162267e-08f,
    1.4066283500795206892487241294e-08f,
   -4.6982570380537099016106141654e-09f,
    1.5491248664620612686423108936e-09f,
   -5.0340936319394885789686867772e-10f,
    1.6084448673736032249959475006e-10f,
   -5.0349733196835456497619787559e-11f,
    1.5357154939762136997591808461e-11f,
   -4.5233809655775649997667176224e-12f,
    1.2664429179254447281068538964e-12f,
   -3.2648287937449326771785041692e-13f,
    7.1528272726086133795579071407e-14f,
   -9.4831735252566034505739531258e-15f,
   -2.3124001991413207293120906691e-15f,
    2.8406613277170391482590129474e-15f,
   -1.7245370321618816421281770927e-15f,
    8.6507923128671112154695006592e-16f,
   -3.9506563665427555895391869919e-16f,
    1.6779342132074761078792361165e-16f,
   -6.0483153034414765129837716260e-17f
};
static Tcheb_series gstar_b_cs = {
    gstar_b_data,
    29,
    -1.f, 1.f,
    18
};


// * coefficients for gamma=7, kmax=8  Lanczos method
static float lanczos_7_c[9] = {
    0.99999999999980993227684700473478f,
    676.520368121885098567009190444019f,
   -1259.13921672240287047156078755283f,
    771.3234287776530788486528258894f,
   -176.61502916214059906584551354f,
    12.507343278686904814458936853f,
   -0.13857109526572011689554707f,
    9.984369578019570859563e-6f,
    1.50563273514931155834e-7f
};


// -------------------------------------------------------------------------
// * Lanczos method for real x > 0;
// * gamma=7, truncated at 1/(z+8) 
// * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
//
static int lngamma_lanczos( float x, float* result, float* err )
{
    int     k;
    float   Ag;
    float   term0, term1, term2;

    x -= 1.0f; // Lanczos writes z! instead of Gamma(z)

    Ag = lanczos_7_c[0];
    for( k = 1; k <= 8; k++ )
        Ag += lanczos_7_c[k]/( x + k );

//     // * (x+0.5)*log(x+7.5) - (x+7.5) + SLC_LNSQRT2PI + log(Ag(x))
//     term1 = ( x + 0.5f ) * logf(( x + 7.5f ) / SLC_E );
//     term2 = SLC_LNSQRT2PI + logf( Ag );

    term0 = logf( (x+7.5f)/SLC_E );
    if( SLC_LOG_SP_MAX < term0 + 5.0f )
        return PSL_ERR_OVERFLOW;

    term1 = ( x + 0.5f ) * term0;
    term2 = SLC_LNSQRT2PI + logf( Ag );

    if( result ) {
        *result = term1 + ( term2 - 7.0f );
        if( err ) {
            *err = 2.0f * SLC_SP_EPSILON * ( fabsf( term1 ) + fabsf( term2 ) + 7.0f );
            *err += SLC_SP_EPSILON * fabsf( *result );
        }
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// * x = eps near zero
// * gives single-precision for |eps| < 0.02
//
static int lngamma_sgn_0( float eps, float* result, float* err, float* sgn )
{
    // * calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps/2
    const float c1  = -0.07721566490153286061f;
    const float c2  = -0.01094400467202744461f;
    const float c3  =  0.09252092391911371098f;
    const float c4  = -0.01827191316559981266f;
    const float c5  =  0.01800493109685479790f;
    const float c6  = -0.00685088537872380685f;
    const float c7  =  0.00399823955756846603f;
    const float c8  = -0.00189430621687107802f;
    const float c9  =  0.00097473237804513221f;
    const float c10 = -0.00048434392722255893f;
    const float g6  = c6 + eps *( c7 + eps *( c8 + eps *( c9 + eps * c10 )));
    const float g   = eps *( c1 + eps *( c2 + eps *( c3 + eps *( c4 + eps *( c5 + eps * g6 )))));

    // * calculate Gamma(eps) eps, a positive quantity
    const float gee = g + 1.0f / ( 1.0f + eps ) + 0.5f * eps;

    if( result ) {
        *result = logf( gee / fabsf( eps ));
        if( err )
            *err = 4.0f * SLC_SP_EPSILON * fabsf( *result );
        if( sgn )
            *sgn = (float)SLC_SIGN( eps );
    }
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// * x near a negative integer
// * Calculates sign as well as log(|gamma(x)|).
// * x = -N + eps
// * assumes N >= 1
//
static int lngamma_sgn_sing( int N, float eps, float* result, float* err, float* sgn )
{
    if( eps == 0.0f ) {
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
            if( sgn )
                *sgn = 0.0f;
        }
        return PSL_ERR_DOMAIN;
    }
    else if( N == 1 ) {
        // * calculate series for
        // * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
        // * single-precision for |eps| < 0.02
        const float c0 =  0.07721566490153286061f;
        const float c1 =  0.08815966957356030521f;
        const float c2 = -0.00436125434555340577f;
        const float c3 =  0.01391065882004640689f;
        const float c4 = -0.00409427227680839100f;
        const float c5 =  0.00275661310191541584f;
        const float c6 = -0.00124162645565305019f;
        const float c7 =  0.00065267976121802783f;
        const float c8 = -0.00032205261682710437f;
        const float c9 =  0.00016229131039545456f;
        const float g5 = c5 + eps *( c6 + eps *( c7 + eps *( c8 + eps * c9 )));
        const float g  = eps *( c0 + eps *( c1 + eps *( c2 + eps *( c3 + eps *( c4 + eps * g5 )))));

        // * calculate eps gamma(-1+eps), a negative quantity
        const float gam_e = g - 1.0f - 0.5f * eps * ( 1.0f + 3.0f * eps ) / ( 1.0f - eps * eps );

        if( result ) {
            *result = logf( fabsf( gam_e ) / fabsf( eps ));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            if( sgn )
                *sgn = (( 0.0f < eps )? -1.0f : 1.0f );
        }
        return PSL_SUCCESS;
    }
    else {
        float g;

        // * series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
        // * single-precision for |eps| < 0.02
        const float cs1 = -1.6449340668482264365f;
        const float cs2 =  0.8117424252833536436f;
        const float cs3 = -0.1907518241220842137f;
        const float cs4 =  0.0261478478176548005f;
        const float cs5 = -0.0023460810354558236f;
        const float e2  = eps * eps;
        const float sin_ser = 1.0f + e2 *( cs1 + e2 *( cs2 + e2 *( cs3 + e2 *( cs4 + e2 * cs5 ))));

        // * calculate series for ln(gamma(1+N-eps))
        // * single-precision for |eps| < 0.02
        float aeps = fabsf( eps );
        float c1, c2, c3, c4, c5, c6, c7;
        float lng_ser;
        float c0, c0err;
        float psi_0, psi_0err;
        float psi_1, psi_1err;
        float psi_2 = 0.0f, psi_2err;
        float psi_3 = 0.0f, psi_3err;
        float psi_4 = 0.0f, psi_4err;
        float psi_5 = 0.0f, psi_5err;
        float psi_6 = 0.0f, psi_6err;

        psl_lnfact_e( N, &c0, &c0err );
        psl_psi_int_e( N+1, &psi_0, &psi_0err );
        psl_psi_1_int_e( N+1, &psi_1, &psi_1err );
        if( 0.00001f < aeps )   psl_psi_n_e( 2, N+1.0f, &psi_2, &psi_2err );
        if( 0.0002f  < aeps )   psl_psi_n_e( 3, N+1.0f, &psi_3, &psi_3err );
        if( 0.001f   < aeps )   psl_psi_n_e( 4, N+1.0f, &psi_4, &psi_4err );
        if( 0.005f   < aeps )   psl_psi_n_e( 5, N+1.0f, &psi_5, &psi_5err );
        if( 0.01f    < aeps )   psl_psi_n_e( 6, N+1.0f, &psi_6, &psi_6err );
        c1 = psi_0;
        c2 = psi_1 / 2.0f;
        c3 = psi_2 / 6.0f;
        c4 = psi_3 / 24.0f;
        c5 = psi_4 / 120.0f;
        c6 = psi_5 / 720.0f;
        c7 = psi_6 / 5040.0f;
        lng_ser = c0 - eps *( c1 - eps *( c2 - eps *( c3 - eps *( c4 - eps *( c5 - eps *( c6 - eps * c7 ))))));

        // * calculate
        // * g = ln(|eps gamma(-N+eps)|)
        // *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
        g = -lng_ser - logf( sin_ser );

        if( result ) {
            *result = g - logf( fabsf( eps ));
            if( err )
                *err = c0err + 2.0f * SLC_SP_EPSILON * ( fabsf( g ) + fabsf( *result ));
            if( sgn )
                *sgn = ( SLC_ODD( N )? -1.0f: 1.0f ) * (( 0.0f < eps )? 1.0f: -1.0f );
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
//
inline
static int lngamma_1_pade( const float eps, float* result, float* err )
{
    // * Use (2,2) Pade for Log[Gamma[1+eps]]/eps
    // * plus a correction series.
    const float n1 = -1.0017419282349508699871138440f;
    const float n2 =  1.7364839209922879823280541733f;
    const float d1 =  1.2433006018858751556055436011f;
    const float d2 =  5.0456274100274010152489597514f;
    const float num = ( eps + n1 ) * ( eps + n2 );
    const float den = ( eps + d1 ) * ( eps + d2 );
    const float pade = 2.0816265188662692474880210318f * num / den;
    const float c0 =  0.004785324257581753f;
    const float c1 = -0.01192457083645441f;
    const float c2 =  0.01931961413960498f;
    const float c3 = -0.02594027398725020f;
    const float c4 =  0.03141928755021455f;
    const float eps5 = eps * eps * eps * eps * eps;
    const float corr = eps5 *(c0 + eps *( c1 + eps *( c2 + eps *( c3 + c4 * eps ))));
    if( result ) {
        *result = eps * ( pade + corr );
        if( err )
            *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
//
inline
static int lngamma_2_pade( const float eps, float* result, float* err )
{
    // * Use (2,2) Pade for Log[Gamma[2+eps]]/eps
    // * plus a correction series.
    const float n1 = 1.000895834786669227164446568f;
    const float n2 = 4.209376735287755081642901277f;
    const float d1 = 2.618851904903217274682578255f;
    const float d2 = 10.85766559900983515322922936f;
    const float num = ( eps + n1 ) * ( eps + n2 );
    const float den = ( eps + d1 ) * ( eps + d2 );
    const float pade = 2.85337998765781918463568869f * num / den;
    const float c0 =  0.0001139406357036744f;
    const float c1 = -0.0001365435269792533f;
    const float c2 =  0.0001067287169183665f;
    const float c3 = -0.0000693271800931282f;
    const float c4 =  0.0000407220927867950f;
    const float eps5 = eps * eps * eps * eps * eps;
    const float corr = eps5 *( c0 + eps * (c1 + eps *( c2 + eps *( c3 + c4 * eps ))));
    if( result ) {
        *result = eps * ( pade + corr );
        if( err )
            *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
    }
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// * series for gammastar(x)
// * single-precision for x > 10.0
//
static int gammastar_ser( const float x, float* result, float* err )
{
    // * Use the Stirling series for the correction to Log(Gamma(x)),
    // * which is better behaved and easier to compute than the
    // * regular Stirling series for Gamma(x). 
    const float y = 1.0f / ( x * x );
    const float c0 =  1.0f / 12.0f;
    const float c1 = -1.0f / 360.0f;
    const float c2 =  1.0f / 1260.0f;
    const float c3 = -1.0f / 1680.0f;
    const float c4 =  1.0f / 1188.0f;
    const float c5 = -691.0f / 360360.0f;
    const float c6 =  1.0f / 156.0f;
    const float c7 = -3617.0f / 122400.0f;
    const float ser = c0 + y *( c1 + y *( c2 + y *( c3 + y *( c4 + y *( c5 + y *( c6 + y * c7 ))))));
    if( result ) {
        *result = expf( ser / x );
        if( err )
            *err = 2.0f * SLC_SP_EPSILON * *result * SLC_MAX( 1.0f, ser/x );
    }
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// * Chebyshev expansion for log(gamma(x)/gamma(8))
// * 5 < x < 10
// * -1 < t < 1
//
static float gamma_5_10_data[24] = {
   -1.5285594096661578881275075214f,
    4.8259152300595906319768555035f,
    0.2277712320977614992970601978f,
   -0.0138867665685617873604917300f,
    0.0012704876495201082588139723f,
   -0.0001393841240254993658962470f,
    0.0000169709242992322702260663f,
   -2.2108528820210580075775889168e-06f,
    3.0196602854202309805163918716e-07f,
   -4.2705675000079118380587357358e-08f,
    6.2026423818051402794663551945e-09f,
   -9.1993973208880910416311405656e-10f,
    1.3875551258028145778301211638e-10f,
   -2.1218861491906788718519522978e-11f,
    3.2821736040381439555133562600e-12f,
   -5.1260001009953791220611135264e-13f,
    8.0713532554874636696982146610e-14f,
   -1.2798522376569209083811628061e-14f,
    2.0417711600852502310258808643e-15f,
   -3.2745239502992355776882614137e-16f,
    5.2759418422036579482120897453e-17f,
   -8.5354147151695233960425725513e-18f,
    1.3858639703888078291599886143e-18f,
   -2.2574398807738626571560124396e-19f
};
static const Tcheb_series gamma_5_10_cs = {
    gamma_5_10_data,
    23,
    -1.f, 1.f,
    11
};


// -------------------------------------------------------------------------
// * gamma(x) for x >= 1/2
// * assumes x >= 1/2
//
static int gamma_xgthalf( const float x, float* result, float* err )
{
    if( x == 0.5f ) {
        if( result ) {
            *result = 1.77245385090551602729817f;
            if( err )
                *err = SLC_SP_EPSILON * *result;
        }
        return PSL_SUCCESS;
    } else if ( x <= ( SLC_G_FACT_NMAX + 1.0 ) && x == floorf( x )) {
        int n = ( int )floorf( x );
        if( result ) {
            *result = fact_table[n-1];
            if( err )
                *err = SLC_SP_EPSILON * *result;
        }
        return PSL_SUCCESS;
    }
    else if( fabsf( x - 1.0f ) < 0.01f ) {
        // * Use series for Gamma[1+eps] - 1/(1+eps).
        const float eps = x - 1.0f;
        const float c1 =  0.4227843350984671394f;
        const float c2 = -0.01094400467202744461f;
        const float c3 =  0.09252092391911371098f;
        const float c4 = -0.018271913165599812664f;
        const float c5 =  0.018004931096854797895f;
        const float c6 = -0.006850885378723806846f;
        const float c7 =  0.003998239557568466030f;
        if( result ) {
            *result = 1.0f / x + eps *( c1 + eps *( c2 + eps *( c3 + eps *( c4 + eps *( c5 + eps *( c6 + eps * c7 ))))));
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( fabsf( x - 2.0f ) < 0.01f ) {
        // * Use series for Gamma[1 + eps].
        const float eps = x - 2.0f;
        const float c1 =  0.4227843350984671394f;
        const float c2 =  0.4118403304264396948f;
        const float c3 =  0.08157691924708626638f;
        const float c4 =  0.07424901075351389832f;
        const float c5 = -0.00026698206874501476832f;
        const float c6 =  0.011154045718130991049f;
        const float c7 = -0.002852645821155340816f;
        const float c8 =  0.0021039333406973880085f;
        if( result ) {
            *result = 1.0f + 
              eps *( c1 + eps *( c2 + eps *( c3 + eps *( c4 + eps *( c5 + eps *( c6 + eps *( c7 + eps * c8 )))))));
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( x < 5.0f ) {
        // * Exponentiating the logarithm is fine, as
        // * long as the exponential is not so large
        // * that it greatly amplifies the error.
        float lg, lgerr;
        int status = lngamma_lanczos( x, &lg, &lgerr );
        if( status != PSL_SUCCESS )
            return status;
        if( result ) {
            *result = expf( lg );
            if( err )
                *err = *result * ( lgerr + 2.0f * SLC_SP_EPSILON );
        }
        return PSL_SUCCESS;
    }
    else if( x < 10.0f ) {
        // * This is a sticky area. The logarithm
        // * is too large and the gammastar series
        // * is not good.
        const float gamma_8 = 5040.0f;
        const float t = ( 2.0f * x - 15.0f ) / 5.0f;
        float c, cerr;
        cheb_eval_e( &gamma_5_10_cs, t, &c, &cerr );
        if( result ) {
            *result = expf( c ) * gamma_8;
            if( err ) {
                *err = *result * cerr;
                *err += 2.0f * SLC_SP_EPSILON * *result;
            }
        }
        return PSL_SUCCESS;
    }
    else if( x < SLC_G_GAMMA_XMAX ) {
        // * We do not want to exponentiate the logarithm
        // * if x is large because of the inevitable
        // * inflation of the error. So we carefully
        // * use pow() and exp() with exact quantities.
        float p = powf( x, 0.5f * x );
        float e = expf( -x );
        float q = ( p * e ) * p;
        float pre = SLC_SQRT2 * SLC_SQRTPI * q / sqrtf( x );
        float gstar, gstarerr;
        int gstatus = gammastar_ser( x, &gstar, &gstarerr );
        if( result ) {
            *result = pre * gstar;
            if( err )
                *err = ( x + 2.5f ) * SLC_SP_EPSILON * *result;
        }
        return gstatus;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// --------------------- Functions with Error Codes ------------------------
//
// Log(Gamma(x)), x not being a negative integer or zero; for x<0 the real 
// part of log(Gamma(x)) is returned; 
// computed using the real Lanczos method
//
int psl_lngamma_e( float x, float* result, float* err )
{
    int status;
    if( fabsf( x - 1.0f ) < 0.01f ) {
        // * Note that we must amplify the errors
        // * from the Pade evaluations because of
        // * the way we must pass the argument, i.e.
        // * writing (1-x) is a loss of precision
        // * when x is near 1.
        status = lngamma_1_pade( x - 1.0f, result, err );
        if( err )
            *err *= 1.0f / ( SLC_SP_EPSILON + fabsf( x - 1.0f ));
        return status;
    }
    else if( fabsf( x - 2.0f ) < 0.01f ) {
        status = lngamma_2_pade( x - 2.0f, result, err );
        if( err )
            *err *= 1.0f / ( SLC_SP_EPSILON + fabsf( x - 2.0f ));
        return status;
    }
    else if( 0.5f <= x ) {
        return lngamma_lanczos( x, result, err );
    }
    else if( x == 0.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else if( fabsf( x ) < 0.02f ) {
        float sgn;
        return lngamma_sgn_0( x, result, err, &sgn );
    }
    else if( -0.5f /( SLC_SP_EPSILON * SLC_PI ) < x ) {
        // * Try to extract a fractional part from x
        float z  = 1.0f - x;
        float s  = sinf( SLC_PI * z );
        float as = fabsf( s );
        if( s == 0.0f ) {
            return PSL_ERR_DOMAIN;
        }
        else if( as < SLC_PI * 0.015f ) {
            // * x is near a negative integer, -N
            if( x < INT_MIN + 2.0f ) {
                if( result ) {
                    *result = 0.0f;
                    if( err )
                        *err = 0.0f;
                }
                return PSL_ERR_ROUNDOFF;
            }
            else {
                int N = -( int )( x - 0.5f );
                float eps = x + N;
                float sgn;
                return lngamma_sgn_sing( N, eps, result, err, &sgn);
            }
        }
        else {
            float lg_z, lg_zerr;
            int status = lngamma_lanczos( z, &lg_z, &lg_zerr );
            if( status != PSL_SUCCESS )
                return status;
            if( result ) {
                *result = SLC_LNPI - ( logf( as ) + lg_z );
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result ) + lg_zerr;
            }
            return PSL_SUCCESS;
        }
    }
    else {
        // * |x| was too large to extract any fractional part
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_ERR_ROUNDOFF;
    }
}

// -------------------------------------------------------------------------
// Log(Gamma(x)) with sign of Gamma, x not being a negative integer or zero; 
// computed using the real Lanczos method;
// value of Gamma can be reconstructed by sgn*exp(result); 
//
int psl_lngamma_sgn_e( float x, float* result, float* err, float* sgn )
{
    int status;
    if( fabsf( x - 1.0f ) < 0.01f ) {
        status = lngamma_1_pade( x - 1.0f, result, err );
        if( err )
            *err *= 1.0f /( SLC_SP_EPSILON + fabsf( x - 1.0f ));
        if( sgn )
            *sgn = 1.0f;
        return status;
    }
    else if( fabsf( x - 2.0f ) < 0.01f ) {
        status = lngamma_2_pade( x - 2.0f, result, err );
        if( err )
            *err *= 1.0f /( SLC_SP_EPSILON + fabsf( x - 2.0f ));
        if( sgn )
            *sgn = 1.0f;
        return status;
    }
    else if( 0.5f <= x ) {
        if( sgn )
            *sgn = 1.0f;
        return lngamma_lanczos( x, result, err );
    }
    else if( x == 0.0f ) {
        if( sgn )
            *sgn = 0.0f;
        return PSL_ERR_DOMAIN;
    }
    else if( fabsf( x ) < 0.02f ) {
        return lngamma_sgn_0( x, result, err, sgn);
    }
    else if( -0.5f/( SLC_SP_EPSILON * SLC_PI ) < x ) {
        // * Try to extract a fractional part from x
        float z = 1.0f - x;
        float s = sinf( SLC_PI * x );///x vs. z
        float as = fabsf( s );
        if( s == 0.0f ) {
            if( sgn )
                *sgn = 0.0f;
            return PSL_ERR_DOMAIN;
        }
        else if( as < SLC_PI * 0.015f ) {
            // * x is near a negative integer, -N
            if( x < INT_MIN + 2.0f ) {
                if( result )
                    *result = 0.0f;
                if( err )
                    *err = 0.0f;
                if( sgn )
                    *sgn = 0.0f;
                return PSL_ERR_ROUNDOFF;
            }
            else {
                int N = -( int )( x - 0.5f );
                float eps = x + N;
                return lngamma_sgn_sing( N, eps, result, err, sgn );
            }
        }
        else {
            float lg_z, lg_zerr;
            int status = lngamma_lanczos( z, &lg_z, &lg_zerr );
            if( status != PSL_SUCCESS )
                return status;
            if( sgn )
                *sgn = (( 0.0f < s )? 1.0f: -1.0f );
            if( result ) {
                *result = SLC_LNPI - ( logf( as ) + lg_z );
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result ) + lg_zerr;
            }
            return PSL_SUCCESS;
        }
    }
    else {
        // * |x| was too large to extract any fractional part
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
            if( sgn )
                *sgn = 0.0f;
        }
        return PSL_ERR_ROUNDOFF;
    }
}


// -------------------------------------------------------------------------
// Gamma(x), x not being a negative integer or zero; 
// computed using the real Lanczos method;
// The maximum value of x is SLC_G_GAMMA_XMAX
//
int psl_gamma_e( const float x, float* result, float* err )
{
    if( x < 0.5f ) {
        int     rint_x = ( int )floorf( x + 0.5f );
        float f_x = x - rint_x;
        float sgn_gamma = ( SLC_EVEN( rint_x )? 1.0f : -1.0f );
        float sin_term = sgn_gamma * sinf( SLC_PI * f_x ) / SLC_PI;

        if( sin_term == 0.0f ) {
            return PSL_ERR_DOMAIN;
        }
        else if( -169.0f < x ) {
            float g, gerr;
            gamma_xgthalf( 1.0f - x, &g, &gerr );
            if( fabsf( sin_term) * g * SLC_SP_MIN < 1.0f ) {
                if( result ) {
                    *result = 1.0f /( sin_term * g );
                    if( err ) {
                        *err = fabsf( gerr / g ) * fabsf( *result );
                        *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
                    }
                }
                return PSL_SUCCESS;
            }
            else {
                return PSL_ERR_UNDERFLOW;
            }
        }
        else {
            // * It is hard to control it here.
            // * We can only exponentiate the
            // * logarithm and eat the loss of
            // * precision.
            float lng, lngerr;
            float sgn;
            int lngstatus = psl_lngamma_sgn_e( x, &lng, &lngerr, &sgn );
            int estatus   = psl_exp_mult_err_e( lng, lngerr, sgn, 0.0f, result, err );
            if( estatus != PSL_SUCCESS )
                return estatus;
            if( lngstatus != PSL_SUCCESS )
                return lngstatus;
            return PSL_SUCCESS;
        }
    }
    else {
        return gamma_xgthalf( x, result, err );
    }
}


// -------------------------------------------------------------------------
// Regulated Gamma function Gamma*(x) for x>0
//
int psl_gammastar_e( const float x, float* result, float* err )
{
    if( x <= 0.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else if( x < 0.5f ) {
        float lg = 0.0f, lgerr = 0.0f;
        const int lgstatus = psl_lngamma_e( x, &lg, &lgerr );
        const float lx = logf( x );
        const float c  = 0.5f *( SLC_LN2 + SLC_LNPI );
        const float lnr_val = lg - ( x - 0.5f ) * lx + x - c;
        const float lnr_err = lgerr + 2.0f * SLC_SP_EPSILON *(( x + 0.5f )* fabsf( lx ) + c );
        const int estatus  = psl_exp_err_e( lnr_val, lnr_err, result, err );
        if( lgstatus != PSL_SUCCESS )
            return lgstatus;
        if( estatus != PSL_SUCCESS )
            return estatus;
        return PSL_SUCCESS;
    }
    else if( x < 2.0f ) {
        const float t = 4.0f / 3.0f *( x - 0.5f ) - 1.0f;
        return cheb_eval_e( &gstar_a_cs, t, result, err );
    }
    else if( x < 10.0f ) {
        const float t = 0.25f *( x - 2.0f ) - 1.0f;
        float c, cerr;
        cheb_eval_e( &gstar_b_cs, t, &c, &cerr );
        if( result ) {
            *result = c /( x * x ) + 1.0f + 1.0f /( 12.0f * x );
            if( err ) {
                *err = cerr /( x * x );
                *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
        }
        return PSL_SUCCESS;
    }
    else if( x < 1.0f / SLC_ROOT4_SP_EPSILON ) {
        return gammastar_ser( x, result, err );
    }
    else if( x < 1.0f / SLC_SP_EPSILON ) {
        // * Use Stirling formula for Gamma(x)
        const float xi = 1.0f / x;
        if( result ) {
            *result = 1.0f + xi / 12.0f *( 1.0f + xi / 24.0f *( 1.0f - xi *( 139.0f / 180.0f + 571.0f / 8640.0f * xi )));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = 1.0f / x;
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Reciprocal of the Gamma function, 1/Gamma(x)
// computed using the real Lanczos method;
//
int psl_gammainv_e( const float x, float* result, float* err )
{
    if( x <= 0.0f && x == floor( x )) { // negative integer
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( x < 0.5f ) {
        float lng, lngerr;
        float sgn;
        int lngstatus = psl_lngamma_sgn_e( x, &lng, &lngerr, &sgn );
        if( lngstatus == PSL_ERR_DOMAIN ) {
            if( result ) {
                *result = 0.0f;
                if( err )
                    *err = 0.0f;
            }
            return PSL_SUCCESS;
        }
        else if( lngstatus != PSL_SUCCESS ) {
            if( result ) {
                *result = 0.0f;
                if( err )
                    *err = 0.0f;
            }
            return lngstatus;
        }
        else {
            return psl_exp_mult_err_e( -lng, lngerr, sgn, 0.0f, result, err );
        }
    }
    else {
        float g, gerr;
        int gstatus = gamma_xgthalf( x, &g, &gerr );
        if( gstatus == PSL_ERR_OVERFLOW ) {
            return PSL_ERR_UNDERFLOW;
        }
        else {
            if( result ) {
                *result = 1.0f / g;
                if( err ) {
                    *err = fabsf( gerr / g ) * fabsf( *result );
                    *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
                }
                if( fabsf( *result ) < SLC_SP_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// Taylor coefficient, x^n/n! for x>=0, n>=0
//
int psl_taylorcoeff_e( const int n, const float x, float* result, float* err )
{
    if( x < 0.0f || n < 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n == 0 ) {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( n == 1 ) {
        if( result ) {
            *result = x;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( x == 0.0f ) {
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else {
        const float log2pi = SLC_LNPI + SLC_LN2;
        const float ln_test = n *( logf( x ) + 1.0f ) + 1.0f - ( n + 0.5f )* logf( n + 1.0f ) + 0.5f * log2pi;

        if( ln_test < SLC_LOG_SP_MIN + 1.0f ) {
            return PSL_ERR_UNDERFLOW;
        }
        else if( SLC_LOG_SP_MAX - 1.0f < ln_test ) {
            return PSL_ERR_OVERFLOW;
        }
        else {
            float product = 1.0f;
            int k;
            for( k = 1; k <= n; k++ ) {
                product *= x / k;
            }
            if( result ) {
                *result = product;
                if( err )
                    *err = n * SLC_SP_EPSILON * product;
                if( fabsf( *result ) < SLC_SP_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }    
    }
}


// -------------------------------------------------------------------------
// Factorial, n! = Gamma(n+1);
// The maximum value of n is SLC_G_FACT_NMAX
//
int psl_fact_e( const unsigned int n, float* result, float* err )
{
    if( n < 12 ) {
        if( result ) {
            *result = fact_table[n];
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( n <= SLC_G_FACT_NMAX ) {
        if( result ) {
            *result = fact_table[n];
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
// Double factorial, n!! = n(n-2)(n-4)...;
// The maximum value of n is SLC_G_DOUBLEFACT_NMAX
//
int psl_doublefact_e( const unsigned int n, float* result, float* err )
{
    if( n < 18 ) {
        if( result ) {
            *result = doub_fact_table[n];
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( n <= SLC_G_DOUBLEFACT_NMAX ) {
        if( result ) {
            *result = doub_fact_table[n];
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
// Log of factorial, log(n!) = log(Gamma(n+1));
//
int psl_lnfact_e( const unsigned int n, float* result, float* err )
{
    if( n <= SLC_G_FACT_NMAX ) {
        if( result ) {
            *result = logf( fact_table[n] );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        psl_lngamma_e( n + 1.0f, result, err );
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Log of double factorial, log(n!!)
//
int psl_lndoublefact_e( const unsigned int n, float* result, float* err )
{
    if( n <= SLC_G_DOUBLEFACT_NMAX ) {
        if( result ) {
            *result = logf( doub_fact_table[n] );
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( SLC_ODD( n )) {
        float lg, lgerr;
        int status = psl_lngamma_e( 0.5f *( n + 2.0f ), &lg, &lgerr );
        if( status != PSL_SUCCESS )
            return status;
        if( result ) {
            *result = 0.5f *( n + 1.0f ) * SLC_LN2 - 0.5f * SLC_LNPI + lg;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result) + lgerr;
        }
        return PSL_SUCCESS;
    }
    else {
        float lg, lgerr;
        int status = psl_lngamma_e( 0.5f * n + 1.0f, &lg, &lgerr );
        if( status != PSL_SUCCESS )
            return status;
        if( result ) {
            *result = 0.5f * n * SLC_LN2 + lg;
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result ) + lgerr;
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Log of n choose m, log(n!) - log(m!) - log((n-m)!)
//
int psl_lnchoose_e( unsigned int n, unsigned int m, float* result, float* err )
{
    if( n < m ) {
        return PSL_ERR_DOMAIN;
    }
    else if( m == n || m == 0 ) {
        if( result ) {
            *result = 0.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else {
        float nf, nferr;
        float mf, mferr;
        float nmmf, nmmferr;
        if( n < m * 2 ) 
            m = n - m;
        psl_lnfact_e( n, &nf, &nferr );
        psl_lnfact_e( m, &mf, &mferr );
        psl_lnfact_e( n - m, &nmmf, &nmmferr );
        if( result ) {
            *result = nf - mf - nmmf;
            if( err ) {
                *err = nferr + mferr + nmmferr;
                *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Combinatorial factor n choose m, n!/(m!-(n-m)!)
//
int psl_choose_e( unsigned int n, unsigned int m, float* result, float* err )
{
    if( n < m ) {
        return PSL_ERR_DOMAIN;
    }
    else if( m == n || m == 0 ) {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = 0.0f;
        }
        return PSL_SUCCESS;
    }
    else if( n <= SLC_G_FACT_NMAX ) {
        if( result ) {
            *result = ( fact_table[n] / fact_table[m] ) / fact_table[n-m];
            if( err )
                *err = 6.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        if( m * 2 < n )
            m = n - m;

        if( n - m < 64 ) { // compute product for a manageable number of terms
            float prod = 1.0f;
            unsigned int k;

            for( k = n; m + 1 <= k; k-- ) {
                float tk = ( float )k / ( float )( k - m );
                if( SLC_SP_MAX / prod < tk ) {
                    return PSL_ERR_OVERFLOW;
                }
                prod *= tk;
            }
            if( result ) {
                *result = prod;
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * prod * fabsf((float)(n-m));
            }
            return PSL_SUCCESS;
        }
        else {
            float lc = 0.0f, lcerr = 0.0f;
            const int lcstatus = psl_lnchoose_e( n, m, &lc, &lcerr );
            const int estatus  = psl_exp_err_e( lc, lcerr, result, err );
            if( lcstatus != PSL_SUCCESS )
                return lcstatus;
            if( estatus != PSL_SUCCESS )
                return estatus;
            return PSL_SUCCESS;
        }
    }
}

// -------------------------------------------------------------------------

}//namespace extspsl
