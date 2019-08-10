/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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
#include "cheb.h"
#include "gamma.h"
#include "exp.h"
#include "zeta.h"

namespace extspsl {

// ------------------------- Private Section -------------------------------

// * chebyshev fit for (s(t)-1)Zeta[s(t)]
// * s(t)= (t+1)/2
// * -1 <= t <= 1
static float zeta_xlt1_data[14] = {
    1.48018677156931561235192914649f,
    0.25012062539889426471999938167f,
    0.00991137502135360774243761467f,
   -0.00012084759656676410329833091f,
   -4.7585866367662556504652535281e-06f,
    2.2229946694466391855561441361e-07f,
   -2.2237496498030257121309056582e-09f,
   -1.0173226513229028319420799028e-10f,
    4.3756643450424558284466248449e-12f,
   -6.2229632593100551465504090814e-14f,
   -6.6116201003272207115277520305e-16f,
    4.9477279533373912324518463830e-17f,
   -1.0429819093456189719660003522e-18f,
    6.9925216166580021051464412040e-21f
};
static Tcheb_series zeta_xlt1_cs = {
    zeta_xlt1_data,
    13,
    -1.f, 1.f,
    8
};

// * chebyshev fit for (s(t)-1)Zeta[s(t)]
// * s(t)= (19t+21)/2
// * -1 <= t <= 1
static float zeta_xgt1_data[30] = {
    19.3918515726724119415911269006f,
    9.1525329692510756181581271500f,
    0.2427897658867379985365270155f,
    -0.1339000688262027338316641329f,
    0.0577827064065028595578410202f,
    -0.0187625983754002298566409700f,
    0.0039403014258320354840823803f,
    -0.0000581508273158127963598882f,
    -0.0003756148907214820704594549f,
    0.0001892530548109214349092999f,
    -0.0000549032199695513496115090f,
    8.7086484008939038610413331863e-6f,
    6.4609477924811889068410083425e-7f,
    -9.6749773915059089205835337136e-7f,
    3.6585400766767257736982342461e-7f,
    -8.4592516427275164351876072573e-8f,
    9.9956786144497936572288988883e-9f,
    1.4260036420951118112457144842e-9f,
    -1.1761968823382879195380320948e-9f,
    3.7114575899785204664648987295e-10f,
    -7.4756855194210961661210215325e-11f,
    7.8536934209183700456512982968e-12f,
    9.9827182259685539619810406271e-13f,
    -7.5276687030192221587850302453e-13f,
    2.1955026393964279988917878654e-13f,
    -4.1934859852834647427576319246e-14f,
    4.6341149635933550715779074274e-15f,
    2.3742488509048340106830309402e-16f,
    -2.7276516388124786119323824391e-16f,
    7.8473570134636044722154797225e-17f
};
static Tcheb_series zeta_xgt1_cs = {
    zeta_xgt1_data,
    29,
    -1.f, 1.f,
    17
};


// * chebyshev fit for Ln[Zeta[s(t)] - 1 - 2^(-s(t))]
// * s(t)= 10 + 5t
// * -1 <= t <= 1; 5 <= s <= 15
static float zetam1_inter_data[24] = {
    -21.7509435653088483422022339374f,
    -5.63036877698121782876372020472f,
    0.0528041358684229425504861579635f,
    -0.0156381809179670789342700883562f,
    0.00408218474372355881195080781927f,
    -0.0010264867349474874045036628282f,
    0.000260469880409886900143834962387f,
    -0.0000676175847209968878098566819447f,
    0.0000179284472587833525426660171124f,
    -4.83238651318556188834107605116e-6f,
    1.31913788964999288471371329447e-6f,
    -3.63760500656329972578222188542e-7f,
    1.01146847513194744989748396574e-7f,
    -2.83215225141806501619105289509e-8f,
    7.97733710252021423361012829496e-9f,
    -2.25850168553956886676250696891e-9f,
    6.42269392950164306086395744145e-10f,
    -1.83363861846127284505060843614e-10f,
    5.25309763895283179960368072104e-11f,
    -1.50958687042589821074710575446e-11f,
    4.34997545516049244697776942981e-12f,
    -1.25597782748190416118082322061e-12f,
    3.61280740072222650030134104162e-13f,
    -9.66437239205745207188920348801e-14f
};
static Tcheb_series zetam1_inter_cs = {
    zetam1_inter_data,
    22,
    -1.f, 1.f,
    12
};



// -------------------------------------------------------------------------
// assumes s >= 0 and s != 1.0 
//
inline
static int riemann_zeta_sgt0( float s, float* result, float* err )
{
    if( s < 1.0f ) {
        float c, cerr;
        cheb_eval_e( &zeta_xlt1_cs, 2.0f * s - 1.0f, &c, &cerr );
        if( result ) {
            *result = c /( s - 1.0f );
            if( err )
                *err = cerr / fabsf( s - 1.0f ) + SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else if( s <= 20.0f ) {
        float x = ( 2.0f * s - 21.0f )/ 19.0f;
        float c, cerr;
        cheb_eval_e( &zeta_xgt1_cs, x, &c, &cerr );
        if( result ) {
            *result = c /( s - 1.0f );
            if( err )
                *err = cerr /( s - 1.0f ) + SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        float f2 = 1.0f - powf( 2.0f, -s );
        float f3 = 1.0f - powf( 3.0f, -s );
        float f5 = 1.0f - powf( 5.0f, -s );
        float f7 = 1.0f - powf( 7.0f, -s );
        if( result ) {
            *result = 1.0f/( f2 * f3 * f5 * f7 );
            if( err )
                *err = 3.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
//
inline
static int riemann_zeta1ms_slt0( float s, float* result, float* err )
{
    if( -19.0f < s ) {
        float x = ( -19.0f - 2.0f * s )/ 19.0f;
        float c, cerr;
        cheb_eval_e( &zeta_xgt1_cs, x, &c, &cerr );
        if( result ) {
            *result = c /( -s );
            if( err )
                *err = cerr /( -s ) + SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        float f2 = 1.0f - powf( 2.0f, -( 1.0f - s ));
        float f3 = 1.0f - powf( 3.0f, -( 1.0f - s ));
        float f5 = 1.0f - powf( 5.0f, -( 1.0f - s ));
        float f7 = 1.0f - powf( 7.0f, -( 1.0f - s ));
        if( result ) {
            *result = 1.0f /( f2 * f3 * f5 * f7 );
            if( err )
                *err = 3.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// works for 5 < s < 15
//
static int riemann_zeta_minus_1_intermediate_s( float s, float* result, float* err )
{
    float t = ( s - 10.0f )/ 5.0f;
    float c, cerr;
    cheb_eval_e( &zetam1_inter_cs, t, &c, &cerr );
    if( result ) {
        *result = expf( c ) + powf( 2.0f, -s );
        if( err )
            *err = ( cerr + 2.0f * SLC_SP_EPSILON ) * *result;
    }
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// assumes s is large and positive;
// write: zeta(s) - 1 = zeta(s) * (1 - 1/zeta(s))
// and expand a few terms of the product formula to evaluate 1 - 1/zeta(s);
//
// works well for s > 15
//
static int riemann_zeta_minus1_large_s( float s, float* result, float* err )
{
    float a = powf( 2.0f, -s );
    float b = powf( 3.0f, -s );
    float c = powf( 5.0f, -s );
    float d = powf( 7.0f, -s );
    float e = powf(11.0f, -s );
    float f = powf(13.0f, -s );
    float t1 = a + b + c + d + e + f;
    float t2 = a *( b + c + d + e + f ) + b *( c + d + e + f ) + 
               c *( d + e + f ) + d *( e + f ) + e * f;
    /*
    float t3 = a*(b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f) +
                b*(c*(d+e+f) + d*(e+f) + e*f) +
                c*(d*(e+f) + e*f) +
                d*e*f;
    float t4 = a*(b*(c*(d + e + f) + d*(e + f) + e*f) + c*(d*(e+f) + e*f) + d*e*f) +
                b*(c*(d*(e+f) + e*f) + d*e*f) +
                c*d*e*f;
    float t5 = b*c*d*e*f + a*c*d*e*f+ a*b*d*e*f+ a*b*c*e*f+ a*b*c*d*f+ a*b*c*d*e;
    float t6 = a*b*c*d*e*f;
    */
    float numt = t1 - t2 /* + t3 - t4 + t5 - t6 */;
    float zeta = 1.0f /(( 1.0f - a )*( 1.0f - b )*( 1.0f - c )*( 1.0f - d )*( 1.0f - e )*( 1.0f - f ));
    if( result ) {
        *result = numt * zeta;
        if( err )
            *err = ( 15.0f / s + 1.0f ) * 6.0f * SLC_SP_EPSILON * *result;
    }
    return PSL_SUCCESS;
}


// zeta(n) - 1
#define ZETA_POS_TABLE_NMAX     100
static float zetam1_pos_int_table[ZETA_POS_TABLE_NMAX+1] = {
   -1.5f,                               // zeta(0) 
    0.0f,       // FIXME: Infinity      // zeta(1) - 1 
    0.644934066848226436472415166646f,  // zeta(2) - 1 
    0.202056903159594285399738161511f,
    0.082323233711138191516003696541f,
    0.036927755143369926331365486457f,
    0.017343061984449139714517929790f,
    0.008349277381922826839797549849f,
    0.004077356197944339378685238508f,
    0.002008392826082214417852769232f,
    0.000994575127818085337145958900f,
    0.000494188604119464558702282526f,
    0.000246086553308048298637998047f,
    0.000122713347578489146751836526f,
    0.000061248135058704829258545105f,
    0.000030588236307020493551728510f,
    0.000015282259408651871732571487f,
    7.6371976378997622736002935630e-6f,
    3.8172932649998398564616446219e-6f,
    1.9082127165539389256569577951e-6f,
    9.5396203387279611315203868344e-7f,
    4.7693298678780646311671960437e-7f,
    2.3845050272773299000364818675e-7f,
    1.1921992596531107306778871888e-7f,
    5.9608189051259479612440207935e-8f,
    2.9803503514652280186063705069e-8f,
    1.4901554828365041234658506630e-8f,
    7.4507117898354294919810041706e-9f,
    3.7253340247884570548192040184e-9f,
    1.8626597235130490064039099454e-9f,
    9.3132743241966818287176473502e-10f,
    4.6566290650337840729892332512e-10f,
    2.3283118336765054920014559759e-10f,
    1.1641550172700519775929738354e-10f,
    5.8207720879027008892436859891e-11f,
    2.9103850444970996869294252278e-11f,
    1.4551921891041984235929632245e-11f,
    7.2759598350574810145208690123e-12f,
    3.6379795473786511902372363558e-12f,
    1.8189896503070659475848321007e-12f,
    9.0949478402638892825331183869e-13f,
    4.5474737830421540267991120294e-13f,
    2.2737368458246525152268215779e-13f,
    1.1368684076802278493491048380e-13f,
    5.6843419876275856092771829675e-14f,
    2.8421709768893018554550737049e-14f,
    1.4210854828031606769834307141e-14f,
    7.1054273952108527128773544799e-15f,
    3.5527136913371136732984695340e-15f,
    1.7763568435791203274733490144e-15f,
    8.8817842109308159030960913863e-16f,
    4.4408921031438133641977709402e-16f,
    2.2204460507980419839993200942e-16f,
    1.1102230251410661337205445699e-16f,
    5.5511151248454812437237365905e-17f,
    2.7755575621361241725816324538e-17f,
    1.3877787809725232762839094906e-17f,
    6.9388939045441536974460853262e-18f,
    3.4694469521659226247442714961e-18f,
    1.7347234760475765720489729699e-18f,
    8.6736173801199337283420550673e-19f,
    4.3368086900206504874970235659e-19f,
    2.1684043449972197850139101683e-19f,
    1.0842021724942414063012711165e-19f,
    5.4210108624566454109187004043e-20f,
    2.7105054312234688319546213119e-20f,
    1.3552527156101164581485233996e-20f,
    6.7762635780451890979952987415e-21f,
    3.3881317890207968180857031004e-21f,
    1.6940658945097991654064927471e-21f,
    8.4703294725469983482469926091e-22f,
    4.2351647362728333478622704833e-22f,
    2.1175823681361947318442094398e-22f,
    1.0587911840680233852265001539e-22f,
    5.2939559203398703238139123029e-23f,
    2.6469779601698529611341166842e-23f,
    1.3234889800848990803094510250e-23f,
    6.6174449004244040673552453323e-24f,
    3.3087224502121715889469563843e-24f,
    1.6543612251060756462299236771e-24f,
    8.2718061255303444036711056167e-25f,
    4.1359030627651609260093824555e-25f,
    2.0679515313825767043959679193e-25f,
    1.0339757656912870993284095591e-25f,
    5.1698788284564313204101332166e-26f,
    2.5849394142282142681277617708e-26f,
    1.2924697071141066700381126118e-26f,
    6.4623485355705318034380021611e-27f,
    3.2311742677852653861348141180e-27f,
    1.6155871338926325212060114057e-27f,
    8.0779356694631620331587381863e-28f,
    4.0389678347315808256222628129e-28f,
    2.0194839173657903491587626465e-28f,
    1.0097419586828951533619250700e-28f,
    5.0487097934144756960847711725e-29f,
    2.5243548967072378244674341938e-29f,
    1.2621774483536189043753999660e-29f,
    6.3108872417680944956826093943e-30f,
    3.1554436208840472391098412184e-30f,
    1.5777218104420236166444327830e-30f,
    7.8886090522101180735205378276e-31f
};


#define ZETA_NEG_TABLE_NMAX     63  //99
#define ZETA_NEG_TABLE_SIZE     32  //50
static float zeta_neg_int_table[ZETA_NEG_TABLE_SIZE] = {
   -0.083333333333333333333333333333f,      // zeta(-1)
    0.008333333333333333333333333333f,      // zeta(-3)
   -0.003968253968253968253968253968f,      // ... 
    0.004166666666666666666666666667f,
   -0.007575757575757575757575757576f,
    0.021092796092796092796092796093f,
   -0.083333333333333333333333333333f,
    0.44325980392156862745098039216f,
   -3.05395433027011974380395433027f,
    26.4562121212121212121212121212f,
   -281.460144927536231884057971014f,
    3607.5105463980463980463980464f,
   -54827.583333333333333333333333f,
    974936.82385057471264367816092f,
   -2.0052695796688078946143462272e+07f,
    4.7238486772162990196078431373e+08f,
   -1.2635724795916666666666666667e+10f,
    3.8087931125245368811553022079e+11f,
   -1.2850850499305083333333333333e+13f,
    4.8241448354850170371581670362e+14f,
   -2.0040310656516252738108421663e+16f,
    9.1677436031953307756992753623e+17f,
   -4.5979888343656503490437943262e+19f,
    2.5180471921451095697089023320e+21f,
   -1.5001733492153928733711440151e+23f,
    9.6899578874635940656497942895e+24f,
   -6.7645882379292820990945242302e+26f,
    5.0890659468662289689766332916e+28f,
   -4.1147288792557978697665486068e+30f,
    3.5666582095375556109684574609e+32f,    // ...
   -3.3066089876577576725680214670e+34f,    // zeta(-61)
    3.2715634236478716264211227016e+36f     // zeta(-63)
};


// coefficients for Maclaurin summation in hzeta()
// B_{2j}/(2j)!
static float hzeta_c[15] = {
    1.00000000000000000000000000000f,
    0.083333333333333333333333333333f,
   -0.00138888888888888888888888888889f,
    0.000033068783068783068783068783069f,
   -8.2671957671957671957671957672e-07f,
    2.0876756987868098979210090321e-08f,
   -5.2841901386874931848476822022e-10f,
    1.3382536530684678832826980975e-11f,
   -3.3896802963225828668301953912e-13f,
    8.5860620562778445641359054504e-15f,
   -2.1748686985580618730415164239e-16f,
    5.5090028283602295152026526089e-18f,
   -1.3954464685812523340707686264e-19f,
    3.5347070396294674716932299778e-21f,
   -8.9535174270375468504026113181e-23f
};

#define ETA_POS_TABLE_NMAX  24  //100
static float eta_pos_int_table[ETA_POS_TABLE_NMAX+1] = {
    0.50000000000000000000000000000f,   // eta(0)
    SLC_LN2,                            // eta(1)
    0.82246703342411321823620758332f,   // ... 
    0.90154267736969571404980362113f,
    0.94703282949724591757650323447f,
    0.97211977044690930593565514355f,
    0.98555109129743510409843924448f,
    0.99259381992283028267042571313f,
    0.99623300185264789922728926008f,
    0.99809429754160533076778303185f,
    0.99903950759827156563922184570f,
    0.99951714349806075414409417483f,
    0.99975768514385819085317967871f,
    0.99987854276326511549217499282f,
    0.99993917034597971817095419226f,
    0.99996955121309923808263293263f,
    0.99998476421490610644168277496f,
    0.99999237829204101197693787224f,
    0.99999618786961011347968922641f,
    0.99999809350817167510685649297f,
    0.99999904661158152211505084256f,
    0.99999952325821554281631666433f,
    0.99999976161323082254789720494f,
    0.99999988080131843950322382485f,
    0.99999994039889239462836140314f,
};


#define ETA_NEG_TABLE_NMAX  47  //99
#define ETA_NEG_TABLE_SIZE  24  //50
static float eta_neg_int_table[ETA_NEG_TABLE_SIZE] = {
    0.25000000000000000000000000000f,   // eta(-1)
   -0.12500000000000000000000000000f,   // eta(-3)
    0.25000000000000000000000000000f,   // ... 
   -1.06250000000000000000000000000f,
    7.75000000000000000000000000000f,
   -86.3750000000000000000000000000f,
    1365.25000000000000000000000000f,
   -29049.0312500000000000000000000f,
    800572.750000000000000000000000f,
   -2.7741322625000000000000000000e+7f,
    1.1805291302500000000000000000e+9f,
   -6.0523980051687500000000000000e+10f,
    3.6794167785377500000000000000e+12f,
   -2.6170760990658387500000000000e+14f,
    2.1531418140800295250000000000e+16f,
   -2.0288775575173015930156250000e+18f,
    2.1708009902623770590275000000e+20f,
   -2.6173826968455814932120125000e+22f,
    3.5324148876863877826668602500e+24f,
   -5.3042033406864906641493838981e+26f,
    8.8138218364311576767253114668e+28f,
   -1.6128065107490778547354654864e+31f,
    3.2355470001722734208527794569e+33f,// ...
   -7.0876727476537493198506645215e+35f,// eta(-47)
//     1.6890450341293965779175629389e+38,
};


// --------------------- Functions with Error Codes ------------------------
//
// The Hurwitz zeta function, SUM (k+q)^(-s)
//
int psl_hzeta_e( const float s, const float q, float* result, float* err )
{
    if( s <= 1.0f || q <= 0.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else {
        const float max_bits = 54.0f;
        const float ln_term0 = -s * logf( q );  

        if( ln_term0 < SLC_LOG_SP_MIN + 1.0f ) {
            return PSL_ERR_UNDERFLOW;
        }
        else if( SLC_LOG_SP_MAX - 1.0 < ln_term0 ) {
            return PSL_ERR_OVERFLOW;
        }
        else if(( max_bits < s && q < 1.0f ) || ( 0.5f * max_bits < s && q < 0.25f )) {
            if( result ) {
                *result = powf( q, -s );
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
        else if( 0.5f * max_bits < s && q < 1.0f ) {
            const float p1 = powf( q, -s );
            const float p2 = powf( q /( 1.0f + q ), s );
            const float p3 = powf( q /( 2.0f + q ), s );
            if( result ) {
                *result = p1 * ( 1.0f + p2 + p3 );
                if( err )
                    *err = SLC_SP_EPSILON * ( 0.5f * s + 2.0f ) * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
        else {
            // * Euler-Maclaurin summation formula 
            // * [Moshier, p. 400, with several typo corrections]
            const int jmax = 12;
            const int kmax = 10;
            int j, k;
            const float pmax = powf( kmax + q, -s );
            float scp = s;
            float pcp = pmax /( kmax + q );
            float ans = pmax *(( kmax + q )/( s - 1.0f ) + 0.5f );

            for( k = 0; k < kmax; k++ ) {
                ans += powf( k + q, -s );
            }

            for( j = 0; j <= jmax; j++ ) {
                float delta = hzeta_c[j+1] * scp * pcp;
                ans += delta;
                if( fabsf( delta / ans ) < 0.5f * SLC_SP_EPSILON )
                    break;
                scp *= ( s + 2*j + 1 )*( s + 2*j + 2);
                pcp /= ( kmax + q )*( kmax + q );
            }
            if( result ) {
                *result = ans;
                if( err )
                    *err = 2.0f * ( jmax + 1.0f ) * SLC_SP_EPSILON * fabsf( ans );
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// The Riemann zeta function, SUM_k (k)^(-s) for s<>1
//
int psl_zeta_e( const float s, float* result, float* err )
{
    if( s == 1.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else if( 0.0f <= s ) {
        return riemann_zeta_sgt0( s, result, err );
    }
    else {
        // reflection formula, [Abramowitz+Stegun, 23.2.5]
        float zeta_one_minus_s, zeta_one_minus_s_err;
        const int zomstatus = riemann_zeta1ms_slt0( s, &zeta_one_minus_s, &zeta_one_minus_s_err );
        const float sin_term = 
          ( fmodf( s, 2.0f ) == 0.0f )? 0.0f: sinf( 0.5f * SLC_PI * fmodf( s, 4.0f ))/ SLC_PI;

        if( sin_term == 0.0f ) {
            if( result ) {
                *result = 0.0f;
                if( err )
                    *err = 0.0f;
            }
            return PSL_SUCCESS;
        }
        else if( -40.0f < s ) {   //( -170 < s ) {
            // * We have to be careful about losing digits
            // * in calculating pow(2 Pi, s). The gamma
            // * function is fine because we were careful
            // * with that implementation.
            // * We keep an array of (2 Pi)^(10 n).
            const float twopi_pow[5] = { 1.0f,
                                      9.589560061550901348e+007f,
                                      9.195966217409212684e+015f,
                                      8.818527036583869903e+023f,
                                      8.456579467173150313e+031f
            };
            const int n = (int)floorf(( -s )/ 10.0f );
            const float fs = s + 10.0f * n;
            const float p = powf( 2.0f * SLC_PI, fs ) / twopi_pow[n];

            float g, gerr;
            const int gstatus = psl_gamma_e( 1.0f - s, &g, &gerr );
            if( result ) {
                *result = p * g * sin_term * zeta_one_minus_s;
                if( err ) {
                    *err = fabsf( p * g * sin_term ) * zeta_one_minus_s_err;
                    *err += fabsf( p * sin_term * zeta_one_minus_s ) * gerr;
                    *err += SLC_SP_EPSILON *( fabsf( s ) + 2.0f ) * fabsf( *result );
                }
            }
            if( gstatus != PSL_SUCCESS )
                return gstatus;
            if( zomstatus != PSL_SUCCESS )
                return zomstatus;
            return PSL_SUCCESS;
        }
        else {
            // * The actual zeta function may or may not
            // * overflow here. But we have no easy way
            // * to calculate it when the prefactor(s)
            // * overflow. Trying to use log's and exp
            // * is no good because we loose a couple
            // * digits to the exp error amplification.
            // * When we gather a little more patience,
            // * we can implement something here. Until
            // * then just give up.
            return PSL_ERR_OVERFLOW;
        }
    }
}


// -------------------------------------------------------------------------
// The Riemann zeta function, SUM_k (k)^(-n) for integer n<>1
//
int psl_zeta_int_e( const int n, float* result, float* err )
{
    if( n < 0 ) {
        if( !SLC_ODD( n )) {
            if( result ) {
                *result = 0.0f; // exactly zero at even negative integers
                if( err )
                    *err = 0.0f;
            }
            return PSL_SUCCESS;
        }
        else if( -ZETA_NEG_TABLE_NMAX < n ) {
            if( result ) {
                *result = zeta_neg_int_table[-(n+1)/2];
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
        else {
            return psl_zeta_e( (float)n, result, err );
        }
    }
    else if( n == 1 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= ZETA_POS_TABLE_NMAX ) {
        if( result ) {
            *result = 1.0f + zetam1_pos_int_table[n];
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// The Riemann zeta function minus 1 for s<>1;
// for large positive argument, Riemann zeta -> 1; thus the 
// fractional part is evaluated explicitly
//
int psl_zetam1_e( const float s, float* result, float* err )
{
    if( s <= 5.0f ) {
        int status = psl_zeta_e( s, result, err );
        if( result )
            *result -= 1.0f;
        return status;
    }
    else if( s < 15.0f ) {
        return riemann_zeta_minus_1_intermediate_s( s, result, err );
    }
    else {
        return riemann_zeta_minus1_large_s( s, result, err );
    }
}


// -------------------------------------------------------------------------
// The Riemann zeta function minus 1 for integer n<>1;
// for large positive argument, Riemann zeta -> 1; thus the 
// fractional part is evaluated explicitly
//
int psl_zetam1_int_e( const int n, float* result, float* err )
{
    if( n < 0 ) {
        if( !SLC_ODD( n )) {
            if( result ) {
                *result = -1.0f; //at even negative integers zetam1 == -1 since zeta is exactly zero
                if( err )
                    *err = 0.0f;
            }
            return PSL_SUCCESS;
        }
        else if( -ZETA_NEG_TABLE_NMAX < n ) {
            if( result ) {
                *result = zeta_neg_int_table[-(n+1)/2] - 1.0f;
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
        else {
            // could use gsl_sf_zetam1_e here but subtracting 1 makes no difference
            // for such large values, so go straight to the result
            return psl_zeta_e( (float)n, result, err );
        }
    }
    else if( n == 1 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= ZETA_POS_TABLE_NMAX ) {
        if( result ) {
            *result = zetam1_pos_int_table[n];
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return psl_zetam1_e( (float)n, result, err );
    }
}


// -------------------------------------------------------------------------
// The eta function, (1-2^(1-n))zeta(n) for integer n 
//
int psl_eta_int_e( int n, float* result, float* err )
{
    if( ETA_POS_TABLE_NMAX < n ) {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( 0 <= n ) {
        if( result ) {
            *result = eta_pos_int_table[n];
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        // n < 0
        if( !SLC_ODD( n )) {
            // exactly zero at even negative integers
            if( result ) {
                *result = 0.0f;
                if( err )
                    *err = 0.0f;
            }
            return PSL_SUCCESS;
        }
        else if( -ETA_NEG_TABLE_NMAX < n ) {
            if( result ) {
                *result = eta_neg_int_table[-(n+1)/2];
                if( err )
                    *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            return PSL_SUCCESS;
        }
        else {
            float z, zerr;
            float p, perr;
            int zstatus = psl_zeta_int_e( n, &z, &zerr );
            int pstatus = psl_exp_e(( 1.0f - n ) * SLC_LN2, &p, &perr );
            int mstatus = psl_multiply_e( -p, z, result, err );
            if( err ) {
                *err = fabsf( perr *( SLC_LN2 *( 1.0f - n )) * z ) + zerr * fabsf( p );
                *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
            }
            if( mstatus != PSL_SUCCESS )
                return mstatus;
            if( pstatus != PSL_SUCCESS )
                return pstatus;
            if( zstatus != PSL_SUCCESS )
                return zstatus;
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// The eta function, (1-2^(1-s))zeta(s)
//
int psl_eta_e( const float s, float* result, float* err )
{
    if( 100.0f < s ) {
        if( result ) {
            *result = 1.0f;
            if( err )
                *err = SLC_SP_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( fabsf( s - 1.0f ) < 10.0f * SLC_ROOT5_SP_EPSILON ) {
        float del = s - 1.0f;
        float c0  = SLC_LN2;
        float c1  = SLC_LN2 * ( SLC_EULER - 0.5f * SLC_LN2 );
        float c2  = -0.0326862962794492996f;
        float c3  =  0.0015689917054155150f;
        float c4  =  0.00074987242112047532f;
        if( result ) {
            *result = c0 + del *( c1 + del *( c2 + del *( c3 + del * c4 )));
            if( err )
                *err = 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        float z = 0.0f, zerr = 0.0f;
        float p = 0.0f, perr = 0.0f;
        int zstatus = psl_zeta_e( s, &z, &zerr );
        int pstatus = psl_exp_e(( 1.0f - s ) * SLC_LN2, &p, &perr );
        int mstatus = psl_multiply_e( 1.0f - p, z, result, err );
        if( err ) {
            *err = fabsf( perr *( SLC_LN2 *( 1.0f - s )) * z ) + zerr * fabsf( p );
            *err += 2.0f * SLC_SP_EPSILON * fabsf( *result );
        }
        if( mstatus != PSL_SUCCESS )
            return mstatus;
        if( pstatus != PSL_SUCCESS )
            return pstatus;
        if( zstatus != PSL_SUCCESS )
            return zstatus;
        return PSL_SUCCESS;
    }
}

// -------------------------------------------------------------------------

}//namespace extspsl
