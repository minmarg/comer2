/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rvnorm.h"

namespace extspsl {

// =========================================================================
// RVNorm
//
RVNorm::RVNorm( Rng& rng, TRVN met )
:   RVar( rng ),
    met_( met ),
    std_( 1.0f ),
    mean_( 0.0f )
{
}

RVNorm::~RVNorm()
{
}

// -------------------------------------------------------------------------
// GenPolar: generate random value of normally distributed variable;
//  Polar form of the Box-Mueller method (GNU SL)
//  (Knuth TAOCP v2 p.122)
//
int RVNorm::GenPolar( float* rv1, float* rv2 )
{
    const int maxn = 100;
    float r1, r2;
    float x, y, R, f;
    bool in = false;
    int n;

    if( rv1 == NULL )
        return PSL_ERR_ADDRESS;

    //1.27 times on average (Knuth TAOCP v2)
    for( n = 0; !in && n < maxn; n++ ) {
        //get x,y in uniform square: (-1,1)
        r1 = GetRng().GetSingle();
        r2 = GetRng().GetSingle();
        x = -1.0f + r1 + r1;
        y = -1.0f + r2 + r2;
        R = x * x + y * y;
        //check if in the unit circle
        in = R <= 1.0f && R != 0.0f;
    }

    if( !in )
        r1 = 0.01f;
    else
        r1 = logf(R) / R;
    //Box-Muller transform
    f = sqrtf( -r1 - r1 );
    *rv1 = x * f;
    if( rv2 )
        *rv2 = y * f;
    if( !in )
        return PSL_MAXITERATS;
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenRatio: generate random value of normally distributed variable;
//  Ratio method (Kinderman-Monahan), Leva implementation (adopted from GSL)
//  Leva JL. (1992) ACM Trans Math Software 18, 449-455
//  (Knuth TAOCP v2 p.130)
//
int RVNorm::GenRatio( float* rv )
{
    const int maxn = 100;
    float u, v, x, y, Q;
    const float s = 0.449871f;//Leva constants
    const float t = -0.386595f;
    const float a = 0.19600f;
    const float b = 0.25472f;
    const float r1 = 0.27597f;
    const float r2 = 0.27846f;
    int n = 0;

    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    //1.369 times on average (GSL)
    do {
        if( maxn <= n++ )
            break;
        //generate point P=(u,v) uniform in rectangle
        //enclosing K&M region v^2 <= - 4 u^2 log(u)

        //u in (0,1] to avoid singularity at u=0
        u = 1.0f - GetRng().GetSingle1();
        if( !u )
            return PSL_ERR_DOMAIN;

        //v is in asymmetric interval [-0.5,0.5);  however v=-0.5 is
        //rejected in the last part of the while clause;
        //resulting normal variate is strictly symmetric about 0
        //(provided that v is symmetric once v = -0.5 is excluded)
        v = GetRng().GetSingle() - 0.5f;

        //constant 1.7156 > sqrt(8/e) (for accuracy);
        //but not by too much (for efficiency)
        v *= 1.7156f;

        //compute Leva's quadratic form Q
        x = u - s;
        y = fabsf(v) - t;
        Q = x*x + y*( a*y - b*x );

        //accept P if Q < r1; /Leva/
        //reject P if Q > r2; /Leva/
        //accept if v^2 <= -4 u^2 log(u); /K&M/
        //this final test is executed 0.012 times on average (GSL)
    } while( r1 <= Q && ( r2 < Q || -4.0f*u*u*logf(u) < v*v ));

    //save slope
    if( rv )
        *rv = v / u;
    if( maxn < n )
        return PSL_MAXITERATS;
    return PSL_OK;
}

// =========================================================================
// DATA FOR ZIGGURAT METHOD BY MARSAGLIA & TSANG; GSL implementation
//
//position of right-most step
static const float grvn_ZM_R = 3.44428647676f;//R
static const float grvn_1oZM_R = 0.29033589591f;//1/R

//tabulated values for the heigt of the Ziggurat levels
static const float grvn_ytab[128] = {
  1.0f, 0.963598623011f, 0.936280813353f, 0.913041104253f,
  0.892278506696f, 0.873239356919f, 0.855496407634f, 0.838778928349f,
  0.822902083699f, 0.807732738234f, 0.793171045519f, 0.779139726505f,
  0.765577436082f, 0.752434456248f, 0.739669787677f, 0.727249120285f,
  0.715143377413f, 0.703327646455f, 0.691780377035f, 0.68048276891f,
  0.669418297233f, 0.65857233912f, 0.647931876189f, 0.637485254896f,
  0.62722199145f, 0.617132611532f, 0.607208517467f, 0.597441877296f,
  0.587825531465f, 0.578352913803f, 0.569017984198f, 0.559815170911f,
  0.550739320877f, 0.541785656682f, 0.532949739145f, 0.524227434628f,
  0.515614886373f, 0.507108489253f, 0.498704867478f, 0.490400854812f,
  0.482193476986f, 0.47407993601f, 0.466057596125f, 0.458123971214f,
  0.450276713467f, 0.442513603171f, 0.434832539473f, 0.427231532022f,
  0.419708693379f, 0.41226223212f, 0.404890446548f, 0.397591718955f,
  0.390364510382f, 0.383207355816f, 0.376118859788f, 0.369097692334f,
  0.362142585282f, 0.355252328834f, 0.348425768415f, 0.341661801776f,
  0.334959376311f, 0.328317486588f, 0.321735172063f, 0.31521151497f,
  0.308745638367f, 0.302336704338f, 0.29598391232f, 0.289686497571f,
  0.283443729739f, 0.27725491156f, 0.271119377649f, 0.265036493387f,
  0.259005653912f, 0.253026283183f, 0.247097833139f, 0.241219782932f,
  0.235391638239f, 0.229612930649f, 0.223883217122f, 0.218202079518f,
  0.212569124201f, 0.206983981709f, 0.201446306496f, 0.195955776745f,
  0.190512094256f, 0.185114984406f, 0.179764196185f, 0.174459502324f,
  0.169200699492f, 0.1639876086f, 0.158820075195f, 0.153697969964f,
  0.148621189348f, 0.143589656295f, 0.138603321143f, 0.133662162669f,
  0.128766189309f, 0.123915440582f, 0.119109988745f, 0.114349940703f,
  0.10963544023f, 0.104966670533f, 0.100343857232f, 0.0957672718266f,
  0.0912372357329f, 0.0867541250127f, 0.082318375932f, 0.0779304915295f,
  0.0735910494266f, 0.0693007111742f, 0.065060233529f, 0.0608704821745f,
  0.056732448584f, 0.05264727098f, 0.0486162607163f, 0.0446409359769f,
  0.0407230655415f, 0.0368647267386f, 0.0330683839378f, 0.0293369977411f,
  0.0256741818288f, 0.0220844372634f, 0.0185735200577f, 0.0151490552854f,
  0.0118216532614f, 0.00860719483079f, 0.00553245272614f, 0.00265435214565f
};

//tabulated values for 2^24 times x[i]/x[i+1],
//used to accept for U*x[i+1]<=x[i] without any floating point operations
static const unsigned int grvn_ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

//tabulated values of 2^{-24}*x[i]
static const float grvn_wtab[128] = {
  1.62318314817e-08f, 2.16291505214e-08f, 2.54246305087e-08f, 2.84579525938e-08f,
  3.10340022482e-08f, 3.33011726243e-08f, 3.53439060345e-08f, 3.72152672658e-08f,
  3.8950989572e-08f, 4.05763964764e-08f, 4.21101548915e-08f, 4.35664624904e-08f,
  4.49563968336e-08f, 4.62887864029e-08f, 4.75707945735e-08f, 4.88083237257e-08f,
  5.00063025384e-08f, 5.11688950428e-08f, 5.22996558616e-08f, 5.34016475624e-08f,
  5.44775307871e-08f, 5.55296344581e-08f, 5.65600111659e-08f, 5.75704813695e-08f,
  5.85626690412e-08f, 5.95380306862e-08f, 6.04978791776e-08f, 6.14434034901e-08f,
  6.23756851626e-08f, 6.32957121259e-08f, 6.42043903937e-08f, 6.51025540077e-08f,
  6.59909735447e-08f, 6.68703634341e-08f, 6.77413882848e-08f, 6.8604668381e-08f,
  6.94607844804e-08f, 7.03102820203e-08f, 7.11536748229e-08f, 7.1991448372e-08f,
  7.2824062723e-08f, 7.36519550992e-08f, 7.44755422158e-08f, 7.52952223703e-08f,
  7.61113773308e-08f, 7.69243740467e-08f, 7.77345662086e-08f, 7.85422956743e-08f,
  7.93478937793e-08f, 8.01516825471e-08f, 8.09539758128e-08f, 8.17550802699e-08f,
  8.25552964535e-08f, 8.33549196661e-08f, 8.41542408569e-08f, 8.49535474601e-08f,
  8.57531242006e-08f, 8.65532538723e-08f, 8.73542180955e-08f, 8.8156298059e-08f,
  8.89597752521e-08f, 8.97649321908e-08f, 9.05720531451e-08f, 9.138142487e-08f,
  9.21933373471e-08f, 9.30080845407e-08f, 9.38259651738e-08f, 9.46472835298e-08f,
  9.54723502847e-08f, 9.63014833769e-08f, 9.71350089201e-08f, 9.79732621669e-08f,
  9.88165885297e-08f, 9.96653446693e-08f, 1.00519899658e-07f, 1.0138063623e-07f,
  1.02247952126e-07f, 1.03122261554e-07f, 1.04003996769e-07f, 1.04893609795e-07f,
  1.05791574313e-07f, 1.06698387725e-07f, 1.07614573423e-07f, 1.08540683296e-07f,
  1.09477300508e-07f, 1.1042504257e-07f, 1.11384564771e-07f, 1.12356564007e-07f,
  1.13341783071e-07f, 1.14341015475e-07f, 1.15355110887e-07f, 1.16384981291e-07f,
  1.17431607977e-07f, 1.18496049514e-07f, 1.19579450872e-07f, 1.20683053909e-07f,
  1.21808209468e-07f, 1.2295639141e-07f, 1.24129212952e-07f, 1.25328445797e-07f,
  1.26556042658e-07f, 1.27814163916e-07f, 1.29105209375e-07f, 1.30431856341e-07f,
  1.31797105598e-07f, 1.3320433736e-07f, 1.34657379914e-07f, 1.36160594606e-07f,
  1.37718982103e-07f, 1.39338316679e-07f, 1.41025317971e-07f, 1.42787873535e-07f,
  1.44635331499e-07f, 1.4657889173e-07f, 1.48632138436e-07f, 1.50811780719e-07f,
  1.53138707402e-07f, 1.55639532047e-07f, 1.58348931426e-07f, 1.61313325908e-07f,
  1.64596952856e-07f, 1.68292495203e-07f, 1.72541128694e-07f, 1.77574279496e-07f,
  1.83813550477e-07f, 1.92166040885e-07f, 2.05295471952e-07f, 2.22600839893e-07f
};

// -------------------------------------------------------------------------
// GenZigMsg00: generate random normal variates by implementation of the
//  ziggurat algorithm from
//  Marsaglia G, Tsang WW. (2000) Journal of Statistical Software 5(8), 1-7.
//  (Usually the Zigurrat algorithm is the fastest one)
// 
// GNU SL modified the algorithm:
// 
// 1) use 128 steps instead of 256 to decrease the amount of static
// data necessary.  
// 
// 2) use an acceptance sampling from an exponential wedge
// exp(-R*(x-R/2)) for the tail of the base strip to simplify the
// implementation.  The area of exponential wedge is used in
// calculating 'v' and the coefficients in ziggurat table, so the
// coefficients differ slightly from those in the Marsaglia and Tsang
// paper.
// 
// For r.n.g. issues see
// Leong et al. (2005) Journal of Statistical Software 5(7), 1-4.
// For issues of using r.n. in the rejection method see
// Doornik et al. (2005) Tech report, University of Oxford.
// 
int RVNorm::GenZigMsg00( float* rv )
{
    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    const int maxn = 100;
    unsigned int i, j;
    float x, y;
    int n, sign;

    const unsigned int range = GetRng().GetMax() - GetRng().GetMin();
    const unsigned int offset = GetRng().GetMin();

    for( n = 0; n < maxn; n++ ) {
        if( 0xFFFFFFFF <= range ) {
            unsigned int k = GetRng().Get() - offset;
            i = k & 0xFF;
            j = ( k >> 8 ) & 0xFFFFFF;
        }
        else if( 0x00FFFFFF <= range ) {
            unsigned int k1 = GetRng().Get() - offset;
            unsigned int k2 = GetRng().Get() - offset;
            i = (k1 & 0xFF);
            j = (k2 & 0x00FFFFFF);
        }
        else {
            i = GetRng().GetI( 256 );//choose the step
            j = GetRng().GetI( 16777216 );//sample from 2^24
        }

        sign = ( i & 0x80 )? 1: -1;
        i &= 0x7f;

        x = j * grvn_wtab[i];

        if( j < grvn_ktab[i])
            break;

        if( i < 127 ) {
            float y0, y1, U1;
            y0 = grvn_ytab[i];
            y1 = grvn_ytab[i+1];
            U1 = GetRng().GetSingle();
            y = y1 + ( y0 - y1 )* U1;
        }
        else {
            float U1, U2;
            U1 = 1.0f - GetRng().GetSingle0();
            U2 = GetRng().GetSingle();
            x = grvn_ZM_R - logf(U1) * grvn_1oZM_R;
            y = expf( -grvn_ZM_R *( x - 0.5f * grvn_ZM_R )) * U2;
        }

        if( y < expf( -0.5f * x * x ))
            break;
    }
    if( sign < 0 )
        *rv = -x;
    else
        *rv = x;
    if( maxn <= n )
        return PSL_MAXITERATS;
    return PSL_OK;
}

}//namespace extspsl
