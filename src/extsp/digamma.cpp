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
#include "exp.h"
#include "gamma.h"
#include "zeta.h"
#include "digamma.h"

namespace extspsl {

// Private-Section----------------------------------------------------------
//

/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
                                       with weighted error   2.03E-17
                                        log weighted error  16.69
                              significant figures required  16.39
                                   decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
                                       with weighted error   5.54E-17
                                        log weighted error  16.26
                              significant figures required  14.42
                                   decimal places required  16.86

*/

static float psics_data[23] = {
  -.038057080835217922f,
   .491415393029387130f,
  -.056815747821244730f,
   .008357821225914313f,
  -.001333232857994342f,
   .000220313287069308f,
  -.000037040238178456f,
   .000006283793654854f,
  -.000001071263908506f,
   .000000183128394654f,
  -.000000031353509361f,
   .000000005372808776f,
  -.000000000921168141f,
   .000000000157981265f,
  -.000000000027098646f,
   .000000000004648722f,
  -.000000000000797527f,
   .000000000000136827f,
  -.000000000000023475f,
   .000000000000004027f,
  -.000000000000000691f,
   .000000000000000118f,
  -.000000000000000020f
};

static Tcheb_series psi_cs = {
    psics_data,
    22,
    -1.f, 1.f,
    17
};

static float apsics_data[16] = {    
  -.0204749044678185f,
  -.0101801271534859f,
   .0000559718725387f,
  -.0000012917176570f,
   .0000000572858606f,
  -.0000000038213539f,
   .0000000003397434f,
  -.0000000000374838f,
   .0000000000048990f,
  -.0000000000007344f,
   .0000000000001233f,
  -.0000000000000228f,
   .0000000000000045f,
  -.0000000000000009f,
   .0000000000000002f,
  -.0000000000000000f
};    

static Tcheb_series apsi_cs = {
    apsics_data,
    15,
    -1.f, 1.f,
    9
};

#define PSI_TABLE_NMAX 100
static float psi_table[PSI_TABLE_NMAX+1] = {
  0.0f,  // Infinity                // psi(0)
 -SLC_EULER,                        // psi(1)
  0.42278433509846713939348790992f, // ...
  0.92278433509846713939348790992f,
  1.25611766843180047272682124325f,
  1.50611766843180047272682124325f,
  1.70611766843180047272682124325f,
  1.87278433509846713939348790992f,
  2.01564147795560999653634505277f,
  2.14064147795560999653634505277f,
  2.25175258906672110764745616389f,
  2.35175258906672110764745616389f,
  2.44266167997581201673836525479f,
  2.52599501330914535007169858813f,
  2.60291809023222227314862166505f,
  2.67434666166079370172005023648f,
  2.74101332832746036838671690315f,
  2.80351332832746036838671690315f,
  2.86233685773922507426906984432f,
  2.91789241329478062982462539988f,
  2.97052399224214905087725697883f,
  3.02052399224214905087725697883f,
  3.06814303986119666992487602645f,
  3.11359758531574212447033057190f,
  3.15707584618530734186163491973f,
  3.1987425128519740085283015864f,
  3.2387425128519740085283015864f,
  3.2772040513135124700667631249f,
  3.3142410883505495071038001619f,
  3.3499553740648352213895144476f,
  3.3844381326855248765619282407f,
  3.4177714660188582098952615740f,
  3.4500295305349872421533260902f,
  3.4812795305349872421533260902f,
  3.5115825608380175451836291205f,
  3.5409943255438998981248055911f,
  3.5695657541153284695533770196f,
  3.5973435318931062473311547974f,
  3.6243705589201332743581818244f,
  3.6506863483938174848844976139f,
  3.6763273740348431259101386396f,
  3.7013273740348431259101386396f,
  3.7257176179372821503003825420f,
  3.7495271417468059598241920658f,
  3.7727829557002943319172153216f,
  3.7955102284275670591899425943f,
  3.8177324506497892814121648166f,
  3.8394715810845718901078169905f,
  3.8607481768292527411716467777f,
  3.8815815101625860745049801110f,
  3.9019896734278921969539597029f,
  3.9219896734278921969539597029f,
  3.9415975165651470989147440166f,
  3.9608282857959163296839747858f,
  3.9796962103242182164764276160f,
  3.9982147288427367349949461345f,
  4.0163965470245549168131279527f,
  4.0342536898816977739559850956f,
  4.0517975495308205809735289552f,
  4.0690389288411654085597358518f,
  4.0859880813835382899156680552f,
  4.1026547480502049565823347218f,
  4.1190481906731557762544658694f,
  4.1351772229312202923834981274f,
  4.1510502388042361653993711433f,
  4.1666752388042361653993711433f,
  4.1820598541888515500147557587f,
  4.1972113693403667015299072739f,
  4.2121367424746950597388624977f,
  4.2268426248276362362094507330f,
  4.2413353784508246420065521823f,
  4.2556210927365389277208378966f,
  4.2697055997787924488475984600f,
  4.2835944886676813377364873489f,
  4.2972931188046676391063503626f,
  4.3108066323181811526198638761f,
  4.3241399656515144859531972094f,
  4.3372978603883565912163551041f,
  4.3502848733753695782293421171f,
  4.3631053861958823987421626300f,
  4.3757636140439836645649474401f,
  4.3882636140439836645649474401f,
  4.4006092930563293435772931191f,
  4.4128044150075488557724150703f,
  4.4248526077786331931218126607f,
  4.4367573696833950978837174226f,
  4.4485220755657480390601880108f,
  4.4601499825424922251066996387f,
  4.4716442354160554434975042364f,
  4.4830078717796918071338678728f,
  4.4942438268358715824147667492f,
  4.5053549379469826935258778603f,
  4.5163439489359936825368668713f,
  4.5272135141533849868846929582f,
  4.5379662023254279976373811303f,
  4.5486045001977684231692960239f,
  4.5591308159872421073798223397f,
  4.5695474826539087740464890064f,
  4.5798567610044242379640147796f,
  4.5900608426370772991885045755f,
  4.6001618527380874001986055856f
};

#define PSI_1_TABLE_NMAX 100
static float psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0f,  // Infinity                // psi(1,0)
  SLC_PI * SLC_PI / 6.0f,           // psi(1,1)
  0.644934066848226436472415f,      // ... 
  0.394934066848226436472415f,
  0.2838229557371153253613041f,
  0.2213229557371153253613041f,
  0.1813229557371153253613041f,
  0.1535451779593375475835263f,
  0.1331370146940314251345467f,
  0.1175120146940314251345467f,
  0.1051663356816857461222010f,
  0.0951663356816857461222010f,
  0.0869018728717683907503002f,
  0.0799574284273239463058557f,
  0.0740402686640103368384001f,
  0.0689382278476838062261552f,
  0.0644937834032393617817108f,
  0.0605875334032393617817108f,
  0.0571273257907826143768665f,
  0.0540409060376961946237801f,
  0.0512708229352031198315363f,
  0.0487708229352031198315363f,
  0.0465032492390579951149830f,
  0.0444371335365786562720078f,
  0.0425467743683366902984728f,
  0.0408106632572255791873617f,
  0.0392106632572255791873617f,
  0.0377313733163971768204978f,
  0.0363596312039143235969038f,
  0.0350841209998326909438426f,
  0.0338950603577399442137594f,
  0.0327839492466288331026483f,
  0.0317433665203020901265817f,
  0.03076680402030209012658168f,
  0.02984853037475571730748159f,
  0.02898347847164153045627052f,
  0.02816715194102928555831133f,
  0.02739554700275768062003973f,
  0.02666508681283803124093089f,
  0.02597256603721476254286995f,
  0.02531510384129102815759710f,
  0.02469010384129102815759710f,
  0.02409521984367056414807896f,
  0.02352832641963428296894063f,
  0.02298749353699501850166102f,
  0.02247096461137518379091722f,
  0.02197713745088135663042339f,
  0.02150454765882086513703965f,
  0.02105185413233829383780923f,
  0.02061782635456051606003145f,
  0.02020133322669712580597065f,
  0.01980133322669712580597065f,
  0.01941686571420193164987683f,
  0.01904704322899483105816086f,
  0.01869104465298913508094477f,
  0.01834810912486842177504628f,
  0.01801753061247172756017024f,
  0.01769865306145131939690494f,
  0.01739086605006319997554452f,
  0.01709360088954001329302371f,
  0.01680632711763538818529605f,
  0.01652854933985761040751827f,
  0.01625980437882562975715546f,
  0.01599965869724394401313881f,
  0.01574770606433893015574400f,
  0.01550356543933893015574400f,
  0.01526687904880638577704578f,
  0.01503731063741979257227076f,
  0.01481454387422086185273411f,
  0.01459828089844231513993134f,
  0.01438824099085987447620523f,
  0.01418415935820681325171544f,
  0.01398578601958352422176106f,
  0.01379288478501562298719316f,
  0.01360523231738567365335942f,
  0.01342261726990576130858221f,
  0.01324483949212798353080444f,
  0.01307170929822216635628920f,
  0.01290304679189732236910755f,
  0.01273868124291638877278934f,
  0.01257845051066194236996928f,
  0.01242220051066194236996928f,
  0.01226978472038606978956995f,
  0.01212106372098095378719041f,
  0.01197590477193174490346273f,
  0.01183418141592267460867815f,
  0.01169577311142440471248438f,
  0.01156056489076458859566448f,
  0.01142844704164317229232189f,
  0.01129931481023821361463594f,
  0.01117306812421372175754719f,
  0.01104961133409026496742374f,
  0.01092885297157366069257770f,
  0.01081070552355853781923177f,
  0.01069508522063334415522437f,
  0.01058191183901270133041676f,
  0.01047110851491297833872701f,
  0.01036260157046853389428257f,
  0.01025632035036012704977199f,    // ... 
  0.01015219706839427948625679f,    // psi(1,99) 
  0.01005016666333357139524567f     // psi(1,100)
};

//
//
//

// -------------------------------------------------------------------------
// * digamma for x both positive and negative; we do both
// * cases here because of the way we use even/odd parts
// * of the function
//
static int psi_x( const float x, float* result, float* err )
{
    const float y = fabsf( x );

    if( x == 0.0f || x == -1.0f || x == -2.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else if( y >= 2.0f ) {
        const float t = 8.0f/( y * y ) - 1.0f;
        float pres, perr;
        cheb_eval_e( &apsi_cs, t, &pres, &perr );
        if( x < 0.0f ) {
            const float s = sinf( SLC_PI * x );
            const float c = cosf( SLC_PI * x );
            if( fabsf( s ) < 2.0f * SLC_SQRT_SP_MIN ) {
                return PSL_ERR_DOMAIN;
            }
            else {
                if( result ) {
                    *result = logf( y ) - 0.5f / x + pres - SLC_PI * c / s;
                    if( err ) {
                        *err  = SLC_PI * fabsf( x ) * SLC_SP_EPSILON /( s * s );
                        *err += perr;
                        *err += SLC_SP_EPSILON * fabsf( *result );
                    }
                }
                return PSL_OK;
            }
        }
        else {
            if( result ) {
                *result = logf( y ) - 0.5f / x + pres;
                if( err ) {
                    *err  = perr;
                    *err += SLC_SP_EPSILON * fabsf( *result );
                }
            }
            return PSL_OK;
        }
    }
    else { // -2 < x < 2
        float pres, perr;

        if( x < -1.0f ) { // x = -2 + v
            const float v  = x + 2.0f;
            const float t1 = 1.0f / x;
            const float t2 = 1.0f /( x + 1.0f );
            const float t3 = 1.0f / v;
            cheb_eval_e( &psi_cs, 2.0f * v - 1.0f, &pres, &perr );

            if( result ) {
                *result = -( t1 + t2 + t3 ) + pres;
                if( err ) {
                    *err  = SLC_SP_EPSILON * ( fabsf(t1) + fabsf(x/(t2*t2)) + fabsf(x/(t3*t3)) );
                    *err += perr;
                    *err += SLC_SP_EPSILON * fabsf( *result );
                }
            }
            return PSL_OK;
        }
        else if( x < 0.0f ) { // x = -1 + v
            const float v  = x + 1.0f;
            const float t1 = 1.0f / x;
            const float t2 = 1.0f / v;
            cheb_eval_e( &psi_cs, 2.0f * v - 1.0f, &pres, &perr );

            if( result ) {
                *result = -( t1 + t2 ) + pres;
                if( err ) {
                    *err  = SLC_SP_EPSILON * ( fabsf(t1) + fabsf(x/(t2*t2)) );
                    *err += perr;
                    *err += SLC_SP_EPSILON * fabsf( *result );
                }
            }
            return PSL_OK;
        }
        else if( x < 1.0f ) { // x = v
            const float t1 = 1.0f / x;
            cheb_eval_e( &psi_cs, 2.0f * x - 1.0f, &pres, &perr );

            if( result ) {
                *result = -t1 + pres;
                if( err ) {
                    *err  = SLC_SP_EPSILON * t1;
                    *err += perr;
                    *err += SLC_SP_EPSILON * fabsf( *result );
                }
            }
            return PSL_OK;
        }
        else { // x = 1 + v
            const float v = x - 1.0f;
            return cheb_eval_e( &psi_cs, 2.0f * v - 1.0f, result, err );
        }
    }
}


// -------------------------------------------------------------------------
// generic polygamma, (d/dx)^n [psi(x)];
// assumes n >= 0 and x > 0
//
static int psi_n_xg0( const int n, const float x, float* result, float* err )
{
    if( n == 0 ) {
        return psl_psi_e( x, result, err );
    }
    else {
        // Abramowitz + Stegun 6.4.10 
        float ln_nf, ln_nferr;
        float hzeta, hzetaerr;
        int hzstatus = psl_hzeta_e( n + 1.0f, x, &hzeta, &hzetaerr );
        int nfstatus = psl_lnfact_e((unsigned int) n, &ln_nf, &ln_nferr );
        int estatus  = psl_exp_mult_err_e( ln_nf, ln_nferr, hzeta, hzetaerr, result, err );
        if( result )
            if( SLC_EVEN( n ))
                *result = -*result;
        if( estatus != PSL_SUCCESS )
            return estatus;
        if( nfstatus != PSL_SUCCESS )
            return nfstatus;
        if( hzstatus != PSL_SUCCESS )
            return hzstatus;
        return PSL_SUCCESS;
    }
}

// =========================================================================
// Psi of natural numbers
//
int psl_psi_int_e( const int n, float* result, float* err )
{
    if( n <= 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= PSI_TABLE_NMAX ) {
        if( result ) {
            *result = psi_table[n];
            if( err )
                *err = SLC_SP_EPSILON * fabsf( *result );
        }
        return PSL_OK;
    }
    else {
        // Abramowitz+Stegun 6.3.18
        const float c2 = -1.0f / 12.0f;
        const float c3 =  1.0f / 120.0f;
        const float c4 = -1.0f / 252.0f;
        const float c5 =  1.0f / 240.0f;
        const float ni2 = ( 1.0f / n ) * ( 1.0f / n );
        const float ser = ni2 * ( c2 + ni2 * ( c3 + ni2 * ( c4 + ni2 * c5 )));
        if( result ) {
            *result = logf((float)n) - 0.5f / n + ser;
            if( err ) {
                *err  = SLC_SP_EPSILON * ( fabsf(logf((float)n)) + fabsf(0.5f/n) + fabsf(ser) );
                *err += SLC_SP_EPSILON * fabsf( *result );
            }
        }
        return PSL_OK;
    }
}

// -------------------------------------------------------------------------
// Psi of real numbers
//
int psl_psi_e( const float x, float* result, float* err )
{
    return psi_x( x, result, err );
}


// -------------------------------------------------------------------------
// Trigamma function, psi'(n), n>0
//
int psl_psi_1_int_e( const int n, float* result, float* err )
{
    if( n <= 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= PSI_1_TABLE_NMAX ) {
        if( result ) {
            *result = psi_1_table[n];
            if( err )
                *err = SLC_SP_EPSILON * *result;
        }
        return PSL_SUCCESS;
    }
    else {
        // Abramowitz+Stegun 6.4.12; single-precision for n > 100
        const float c0 = -1.0f / 30.0f;
        const float c1 =  1.0f / 42.0f;
        const float c2 = -1.0f / 30.0f;
        const float ni2 = ( 1.0f / n ) * ( 1.0f / n );
        const float ser =  ni2 * ni2 *( c0 + ni2 *( c1 + c2 * ni2 ));
        if( result ) {
            *result = ( 1.0f + 0.5f / n + 1.0f /( 6.0f * n * n ) + ser ) / n;
            if( err )
                *err = SLC_SP_EPSILON * *result;
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Trigamma function, psi'(n) for general x
//
int psl_psi_1_e( const float x, float* result, float* err )
{
    if( x == 0.0f || x == -1.0f || x == -2.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else if( 0.0f < x ) {
        return psi_n_xg0( 1, x, result, err );
    }
    else if( -5.0f < x ) {
        // Abramowitz + Stegun 6.4.6
        int M = -(int)floorf( x );
        float fx = x + M;
        float sum = 0.0;
        int m;

        if( fx == 0.0f )
            return PSL_ERR_DOMAIN;

        for( m = 0; m < M; m++ )
            sum += 1.0f /(( x + m )*( x + m ));

        int psistatus = psi_n_xg0( 1, fx, result, err );
        if( result ) {
            *result += sum;
            if( err )
                *err += M * SLC_SP_EPSILON * sum;
        }
        return psistatus;
    }
    else {
        // Abramowitz + Stegun 6.4.7
        const float sin_px = sinf( SLC_PI * x );
        const float d = SLC_PI * SLC_PI /( sin_px * sin_px );
        float r, rerr;
        int psistatus = psi_n_xg0( 1, 1.0f-x, &r, &rerr );
        if( result ) {
            *result = d - r;
            if( err )
                *err = rerr + 2.0f * SLC_SP_EPSILON * d;
        }
        return psistatus;
    }
}


// -------------------------------------------------------------------------
// Polygamma function, psi^(n)(x) for n>=0, x>0
//
int psl_psi_n_e( const int n, const float x, float* result, float* err )
{
    if( n == 0 ) {
        return psl_psi_e( x, result, err );
    }
    else if( n == 1 ) {
        return psl_psi_1_e( x, result, err );
    }
    else if( n < 0 || x <= 0.0f ) {
        return PSL_ERR_DOMAIN;
    }
    else {
        float ln_nf, ln_nferr;
        float hzeta, hzetaerr;
        int hzstatus = psl_hzeta_e( n+1.0f, x, &hzeta, &hzetaerr );
        int nfstatus = psl_lnfact_e((unsigned int)n, &ln_nf, &ln_nferr );
        int estatus  = psl_exp_mult_err_e( ln_nf, ln_nferr, hzeta, hzetaerr, result, err );
        if( result )
            if( SLC_EVEN( n ))
                *result = -*result;
        if( estatus != PSL_SUCCESS )
            return estatus;
        if( nfstatus != PSL_SUCCESS )
            return nfstatus;
        if( hzstatus != PSL_SUCCESS )
            return hzstatus;
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Functions w/ Natural Prototypes
//

// float psl_psi_int( const int n )
// {
//     float result, err;
//     if( psl_psi_int_e( n, &result, &err ) != PSL_OK )
//         ;
//     return result;
// }
// 
// float psl_psi( const float x )
// {
//     float result, err;
//     if( psl_psi_e( x, &result, &err ) != PSL_OK )
//         ;
//     return result;
// }

}//namespace extspsl
