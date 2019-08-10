/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/psl.h"
#include "extsp/gamma.h"
#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "liblib/logitnormal.h"
#include "liblib/BinarySearchStructure.h"
#include "HDPbase.h"

using namespace extspsl;

// =========================================================================
//global object for application only
HDPbase  HDPBASE;
// HDPbase  HDPctBASE;
//
void SetHDPBASEhlp( HDPbase& obj, const char* filename )
{
    mystring msg = "Reading ";
    msg += filename;
    message( msg.c_str());
    obj.ReadParameters( filename );
    obj.SetPriorProbs();
}
void SetHDPBASE( const char* filename )
{
    SetHDPBASEhlp( HDPBASE, filename );
}
// void SetHDPctBASE( const char* filename )
// {
//     SetHDPBASEhlp( HDPctBASE, filename );
// }

// -------------------------------------------------------------------------

float HDPbase::s_defdpmtau_ = 1.0f;//determines the rate of new tables
float HDPbase::s_defdpmgamma_ = 1.0f;//determines the rate of new dishes

//default dimensionality of samples
int HDPbase::s_defsampledim_ = NUMAA;
//default dish size
int HDPbase::s_defdishsize_ = 50;
//default table size
int HDPbase::s_deftablesize_ = 3;
//default basin size
int HDPbase::s_defbasinsize_ = 1000;

float HDPbase::s_defkappa_pp_a_ = 1.0f;//prior param. a for kappa
float HDPbase::s_defkappa_pp_b_ = 1.0f;//prior param. b for kappa
//default value of the NIW kappa parameter (no. measurements)
float HDPbase::s_defkappa_ = 1.0f;
//default value of the NIW nu parameter (deg. of freedom)
float HDPbase::s_defnu_ = 1.0f;

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPbase::HDPbase()
:   uninfprior_( true ),
    clst4each( false ),
    basin_( NULL ),
    menu_( NULL ),
    chain_( NULL ),
    ctxtsize_( 0 ),
    rdnorsts_( 0 ),
    totalnos_( 0 ),
    restype_( TRT_novalues ),
    mhupdate_( 0 ),
    gibbsit_( 0 ),
    restarted_( false ),
    mixweight_( 1.0f ),
    nosupclsts_( -1 ),
    adjweight_( 0.1f ),
    scores_( NULL ),
    tau_( 1.0f ),
    gamma_( 1.0f ),
    iterread_( 0 ),
    mlpd_( -1.e9f ),
    lastlpd_( -1.e9 ),
    S0scalefac_( 1.0f ),
    adjdegf_( 0.0f )
{
}

// -------------------------------------------------------------------------
// destructor:
//
HDPbase::~HDPbase()
{
    DestroyBasin();
    DestroyMenu();
    DestroyChain();
    ctxtsize_ = 0;
    totalnos_ = 0;
    adjdegf_ = 0.0f;
}

// =========================================================================
// Table2Mtx: copy data from table to matrix
//
void HDPbase::Table2Mtx( const Table* table, Pslmatrix* mtx ) const
{
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR("HDPbase::Table2Mtx: Null menu.");

    if( table == NULL || mtx == NULL )
        throw MYRUNTIME_ERROR("HDPbase::Table2Mtx: Memory access error.");

    const int   dim = GetMenu()->GetDim();
    const int   nov = table->GetActualSize();//number of vectors at the table
    const Pslvector* vec = NULL;
    int n, i, j;
    mtx->Reserve( dim, nov );

    for( n = 0, j = 0; n < table->GetSize(); n++ ){
        if( table->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = table->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        for( i = 0; i < vec->GetSize(); i++ )
            mtx->SetValueAt( i, j, vec->GetValueAt( i ));
        j++;
    }
}





// -------------------------------------------------------------------------
// SetPriorProbs: calculate and set prior probabilities for each dish
//
void HDPbase::SetPriorProbs()
{
    mystring        preamb = "HDPbase::SetPriorProbs: ";
    Menu*           menu = GetMenu();
    const float     gamma = GetDPMGamma();

    if( menu == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP menu size.");
    if( GetDPMGamma() <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid value of the HDP parameter gamma.");

    const Dish* dish = NULL;
    float sum;
    int nodshs = menu->GetSize();
    int k, nk;

    //normalize prior probabilities
    sum = gamma;
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish size.");
        sum += (float)nk;
    }
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        nk = dish->GetDishSize();
        menu->SetPriorProbAt( k, (float)nk / sum );
    }
    menu->SetPriorProbNewDish( gamma / sum );
}

// -------------------------------------------------------------------------
// PPrCompare: compare posterior probabilities
static inline
int PPrCompare( const void* key1, const void* key2, void* pars )
{
    void**  params = (void**)pars;
    if( params == NULL )
        throw MYRUNTIME_ERROR("PPrCompare: Null parameters.");
    const Menu* menu = (const Menu*)params[0];
    const float* lprior = (const float*)params[1];
    if( menu == NULL || lprior == NULL )
        throw MYRUNTIME_ERROR("PPrCompare: Null parameters.");
    Dish*   dish = NULL;
    int     nodshs = menu->GetSize();
    int     k1 = (int)(ssize_t)key1;
    int     k2 = (int)(ssize_t)key2;
    float   pr1 = *lprior;
    float   pr2 = *lprior;
    float   diff;
    if( 0 < k1 ) {
        if( nodshs <= k1 || ( dish = menu->GetDishAt(k1)) == NULL )
            throw MYRUNTIME_ERROR("PPrCompare: Null dish.");
        pr1 = dish->GetTmpValue();
    }
    if( 0 < k2 ) {
        if( nodshs <= k2 || ( dish = menu->GetDishAt(k2)) == NULL )
            throw MYRUNTIME_ERROR("PPrCompare: Null dish.");
        pr2 = dish->GetTmpValue();
    }
    //key2-key1 -- to sort in descending order
    diff = pr2 - pr1;
    return ( diff < 0.0f )? -1: (( 0.0f < diff )? 1: 0 );
}

// =========================================================================
// CalcPPProbs: calculate posterior predictive probabilities of central 
//  vector of logit-normal variable;
//  lnvar -- input logit-normal variable
//  bppr -- background posterior probability of vector `lnvar'
//  ppprobs -- output posterior probabilities
//  cindcs -- indices of clusters 
//  NOTE: input value `lnvar' will be changed on exit
//
void HDPbase::CalcPPProbs( Pslvector& lnvar, float* bppr, Pslvector* ppprobs, Ivector* cindcs,
                           float lpfact, bool usepriors, bool tonormal ) const
{
    mystring errstr;
    mystring preamb = "HDPbase::CalcPPProbs: ";
    if( GetDPMGamma() <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid value of the HDP parameter gamma.");

    const Menu*     menu = GetMenu();
    const float     tau0 = GetDPMTau();
    const float     gamma = GetDPMGamma();
    int             nospcls = GetNoSupClusters();
//     const float     cdLPFACT = 0.02f;//log-space factor in calculating probability threshold; TRAINING value
    const float     cdLPFACT = lpfact;//0.02f;//log-space factor in calculating probability threshold; WORKING value
//     const float     cdMINPRB = 0.1f;//minimum probability threshold
    const bool      cbPRIORS = usepriors;//if false, posteriors are approximated to likelihoods

    if( menu == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP menu size.");

    int     nogrps = GetReadNoGroups();
    int     nodshs = menu->GetSize();
    int     dim = menu->GetDim();
    int     adim = dim + 1;
    int     ctxtsz = menu->GetCtx();

    Dish*       dish = NULL;
    Pslvector   tvec( dim*ctxtsz );//normal transform
    float*  pv, *ptv;
    float   thlprob;//threshold log probability
    float   lprob, lprior, maxlp;
    float   pnorm;
    bool    pset;
    size_t  c;
    int     n, d, k, nk, mk, mtot;

    void*   params[] = { (void*)menu, (void*)&lprior };
    BinarySearchStructure dshnxs( PPrCompare, (size_t)nodshs+1, true/*keep*/, (void*)params );

    if( dim < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP dimensions.");
    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile context (fragment) size.");
    if( nogrps < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of groups.");
    if( nospcls == -1 )
        nospcls = nodshs + 1;//1 for prior
    if( nospcls < 2 || nodshs+1 < nospcls )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of HDP support clusters.");

    if( bppr == NULL || ppprobs == NULL || cindcs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null output data.");

    if( tonormal ) {
        if( lnvar.GetSize() != adim * ctxtsz )
            throw MYRUNTIME_ERROR( preamb + "Invalid size of input vector.");
        //set normal vector
        pv = lnvar.GetVector();
        ptv = tvec.GetVector();
        for( n = 0; n < ctxtsz; n++, pv += adim, ptv += dim ) {
            //TRANSFORM to multivariate normal
            ::LogitNormal2Normal_f( pv, (size_t)adim, 1.e-1f, false );
            //reduce dimenions of resulting vector
            for( d = 0; d < dim; d++ )
                ptv[d] = pv[d];
        }
    }
    else {
        if( lnvar.GetSize() != dim * ctxtsz )
            throw MYRUNTIME_ERROR( preamb + "Invalid size of input vector.");
        tvec = lnvar;
    }

    ppprobs->Allocate( nospcls );
    cindcs->Allocate( nospcls );
    ppprobs->Clear();
    cindcs->Clear();
    *bppr = 0.0f;

    mtot = 0; //total number of tables (local clusters)
    maxlp = -LOC_DBL_MAX;
    //calculate and save log prior probability to lprior
    PriorProbVec( &tvec, &lprior );
    if( !isfinite( lprior ))
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP prior probability value.");
    if( maxlp < lprior )
        maxlp = lprior;
    //calculate max exponent in log posterior predictive probabilities over all dishes
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        mk = dish->GetReadNoTables();
        if( nk < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish size.");
        if( mk < 1 || nk < mk )
            throw MYRUNTIME_ERROR( preamb + "Invalid number of tables associated with HDP dish.");
        mtot += mk;
        ProbVecOfDish( &tvec, k, &lprob );
        //density values can be >1
        if( !isfinite( lprob ))
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP posterior predictive probability value.");
        if( maxlp < lprob )
            maxlp = lprob;
        dish->SetTmpValue( lprob );
    }
    //prior probability and set threshold probability
    pnorm = 0.0f;
    lprior -= maxlp;
    if( lprior < SLC_LOG_SP_MIN || !isfinite( lprior ))
        lprior = 0.0f;
    else {
        lprior = expf( lprior );
        if( cbPRIORS )///
            ///lprior *= gamma;//apprx.
            lprior *= ((float)nogrps)*tau0*gamma;
        if( lprior )
            pnorm += lprior;
    }
    //exponentiate log probabilities
    for( k = 0; k < nodshs; k++ ) {
        pset = true;
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        mk = dish->GetReadNoTables();
        if( nk < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish size.");
        if( mk < 1 || nk < mk || mtot < mk )
            throw MYRUNTIME_ERROR( preamb + "Invalid number of tables associated with HDP dish.");
        lprob = dish->GetTmpValue();
        lprob -= maxlp;
        if( lprob < SLC_LOG_SP_MIN || !isfinite( lprob ))
            lprob = 0.0f;
        else {
            lprob = expf( lprob );
            if( cbPRIORS )///
                ///lprob *= (float)nk;//apprx.
                lprob *= (float)nk *((float)mtot+gamma ) + (float)(nogrps)*(float)mk*tau0;
        }
        if( pset )
            dish->SetTmpValue( lprob );
        if( lprob )
            pnorm += lprob;
    }

    if( pnorm <= 0.0f )
        return;

    //normalize probabilities
    if( pnorm <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Failed to normalize HDP posterior probabilities.");
    if( lprior )
        lprior /= pnorm;
    dshnxs.Push((const void*)-1 );
    thlprob = powf( lprior, cdLPFACT );
//     thlprob = cdMINPRB;

    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        if( dish->GetTmpValue()) {
            lprob = dish->GetTmpValue() / pnorm;
//             if( lprob < thlprob ) {
//                 dish->SetTmpValue( 0.0f );
//                 continue;
//             }
            dish->SetTmpValue( lprob );
            dshnxs.Push((const void*)(ssize_t)k );
        }
    }

    if( dshnxs.GetSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "No clusters to calculate probabilities for.");

    //write top probabilities and calculate background posterior probability
    *bppr = 0.0f;
    for( c = 0; c < dshnxs.GetSize() && c < (size_t)(nospcls-1); c++ ) {
        k = (int)(ssize_t)dshnxs.GetValueAt(c);
        if( k < 0 )
            continue;//prior will be placed at the end
        if( nodshs < k || ( dish = menu->GetDishAt(k)) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish.");
        if( c == 0 && dish->GetTmpValue() < thlprob )
            break;//the greatest probability is less than the threshold
        cindcs->Push( k );
        ppprobs->Push( dish->GetTmpValue());
        *bppr += dish->GetTmpValue() * ( cbPRIORS? 1.0f: menu->GetPriorProbAt(k));//already mixed
    }
    cindcs->Push( -1 );
    ppprobs->Push( lprior );
    *bppr += lprior * ( cbPRIORS? 1.0f: menu->GetPriorProbNewDish());//already mixed
}

// =========================================================================
// MixLNVar: mix central vector of logit-normal variable using clusters 
//  encapsulated in this class object;
//  lnvar -- input logit-normal variable
//  lnmixed -- output mixed logit-normal variable
//  NOTE: input value `lnvar' will be changed on exit
  //
void HDPbase::MixCLNVar( Pslvector& lnvar, Pslvector* lnmixed ) const
{
    mystring preamb = "HDPbase::MixCLNVar: ";
    if( GetDPMGamma() <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid value of the HDP parameter gamma.");

    const Menu*     menu = GetMenu();
    const float     gamma = GetDPMGamma();
    const float     mixwgt = GetMixWeight();
    const float     lpfact = 0.1f;//log-space factor in calculating probability threshold
    const bool      sngdsh = true;//use the single dish of the highest probability in mixing

    if( menu == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP menu size.");

    if( mixwgt <= 0.0f )
        return;
    if( 1.0f < mixwgt )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP mixing weight.");

    int     nodshs = menu->GetSize();
    int     dim = menu->GetDim();
    int     adim = dim + 1;
    int     ctxtsz = menu->GetCtx();
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    const Pslvector* mu = NULL;
    Dish*       dish = NULL;
    Pslvector   tvec( dim*ctxtsz );//normal transform
    Pslvector   subv;
    Pslvector   mxdv;
    float*  pv, *ptv;
    float   thlprob;//threshold log probability
    float   lprob, lprior, maxlp;
    float   pnorm, sum, lnv;
    int     n, d, k, kmaxp, nk;
    int     err;

    if( dim < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP dimensions.");
    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile context (fragment) size.");

    if( lnvar.GetSize() != adim * ctxtsz )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of input vector.");
    if( lnmixed == NULL || lnmixed->GetSize() != adim )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of output vector.");

    //copy input central vector to output vector
    subv = lnvar.SubVector( mid*adim, adim );
    lnmixed->Copy( subv );

    //set normal vector
    pv = lnvar.GetVector();
    ptv = tvec.GetVector();
    for( n = 0; n < ctxtsz; n++, pv += adim, ptv += dim ) {
        //TRANSFORM to multivariate normal
        ::LogitNormal2Normal_f( pv, adim, 1.e-1f, false );
        //reduce dimenions of resulting vector
        for( d = 0; d < dim; d++ )
            ptv[d] = pv[d];
    }

    maxlp = -LOC_DBL_MAX;
    //calculate and save log prior probability to lprior
    PriorProbVec( &tvec, &lprior );
    if( !isfinite( lprior ))
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP prior probability value.");
    if( maxlp < lprior )
        maxlp = lprior;
    //get threshold probability
    thlprob = lprior * lpfact;
    //calculate log posterior predictive probabilities for each dish
    for( k = 0, kmaxp = -1; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish size.");
        ProbVecOfDish( &tvec, k, &lprob );
        //density values can be >1
        if( !isfinite( lprob ))
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP posterior predictive probability value.");
        if( maxlp < lprob && thlprob <= lprob ) {
            maxlp = lprob;
            kmaxp = k;
        }
        dish->SetTmpValue( lprob );
    }
    //exponentiate log probabilities
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        lprob = dish->GetTmpValue();
        if( lprob < thlprob ) {
            dish->SetTmpValue(0.0f);
            continue;
        }
        if( sngdsh && kmaxp != k ) {
            dish->SetTmpValue(0.0f);
            continue;
        }
        lprob -= maxlp;
        if( lprob < SLC_LOG_SP_MIN || !isfinite( lprob ))
            lprob = 0.0f;
        else
            lprob = expf( lprob );
        dish->SetTmpValue( lprob );
    }
    //check threshold of prior probability
    if( lprior < thlprob || ( sngdsh && 0 < kmaxp ))
        lprior = 0.0f;
    else {
        lprior -= maxlp;
        if( lprior < SLC_LOG_SP_MIN || !isfinite( lprior ))
            lprior = 0.0f;
        else
            lprior = expf( lprior );
    }
    //normalize probabilities
    pnorm = 0.0f;
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid HDP dish size.");
        if( dish->GetTmpValue()) {
            dish->SetTmpValue((float)nk * dish->GetTmpValue());
            pnorm += dish->GetTmpValue();
        }
    }
    if( lprior ) {
        lprior *= gamma;
        pnorm += lprior;
    }

    if( pnorm <= 0.0f )
        return;

    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish.");
        if( dish->GetTmpValue())
            dish->SetTmpValue( dish->GetTmpValue() / pnorm );
    }
    if( lprior )
        lprior /= pnorm;

    //mix central vector with mean vectors of dishes
    mxdv = lnmixed->SubVector( 0, dim );
    mxdv.Zero();
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        mu = menu->GetMuVectorAt(k);
        if( mu == NULL || dish == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null HDP dish mean vector.");
        subv = mu->SubVector( mid*dim, dim );
        if(( err = mxdv.Superposition( dish->GetTmpValue(), subv )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    }
    subv = tvec.SubVector( mid*dim, dim );
    if(( err = mxdv.Superposition( lprior, subv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    //apply weight of mixing
    if( mixwgt < 1.0f ) {
        if(( err = mxdv.MultiplyBy( mixwgt )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        if(( err = mxdv.Superposition( 1.0f - mixwgt, subv )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    }

    //back-transform to logit-normal space
    if(( err = mxdv.Exp()) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    pnorm = 1.0f + mxdv.Sum();
    if( !isfinite( pnorm ))
        throw MYRUNTIME_ERROR( preamb + "Failed to back-transform HDP-mixed vector.");
    sum = 0.0f;
    for( d = 0; d < dim; d++ ) {
        lnv = mxdv.GetValueAt(d) / pnorm;
        lnmixed->SetValueAt( d, lnv );
        sum += lnv;
    }
    lnmixed->SetValueAt( dim, lnv = 1.0f - sum );
}





// =========================================================================
// CalcPriorParams: calculate prior parameters of NIW distribution;
//  return flag of adding noise
//
bool HDPbase::CalcPriorParams( bool addnoise, float stdfact )
{
    SPDmatrix*  S0 = NULL;
    Pslvector   stds;
    mystring    preamb = "HDPbase::CalcPriorParams: ";
    float val;
    int novs;
    int n, err;

    if( GetBasin() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null basin.");
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null menu.");

    if(( err = CalcPriorParams( stds )) == 0 )
        return false;

    if( !addnoise )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    novs = GetBasin()->GetSize();
    S0 = GetMenu()->GetS0();
    if( S0 == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null scatter matrix.");
    if( novs < 1 )
        throw MYRUNTIME_ERROR( preamb + "Basin contains no vectors.");
    if( stdfact < 0 )
        throw MYRUNTIME_ERROR( preamb + "Negative standard deviation factor.");
    //set standard deviations for each dimension
    stds.Reserve( S0->GetNoRows());
    for( n = 0; n < stds.GetSize(); n++ ) {
        val = S0->GetValueAt( n, n ) / (float)novs;
        if( val < 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Negative standard deviation.");
        stds.SetValueAt( n, stdfact * sqrtf(val));
    }

    if(( err = CalcPriorParams( stds )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    return true;
}

// =========================================================================
// CalcPriorParams: calculate prior parameters of NIW distribution
//
int HDPbase::CalcPriorParams( const Pslvector& stds )
{
    if( GetBasin() == NULL )
        return PSL_ERR_ADDRESS;
    if( GetMenu() == NULL )
        return PSL_ERR_ADDRESS;
    if( GetS0ScaleFac() <= 0.0f )
        throw MYRUNTIME_ERROR("HDPbase::CalcPriorParams: Invalid scale factor.");

    const Pslvector* frq = NULL;
    Pslvector*  mu = NULL;
    SPDmatrix*  S0 = NULL;
    SPDmatrix*  S0inv = NULL;
    mystring    preamb = "HDPbase::CalcPriorParams: ";
    myruntime_error mre;
    float ldet = 0.0f;
    float inos = 0.0f;
    int err, dim = 0;
    int n, nos;

    try {
        for( n = 0, nos = 0; n < GetBasin()->GetSize(); n++ ) {
            frq = GetBasin()->GetValueAt( n );
            if( frq == NULL )
                continue;
            if( mu == NULL ) {
                mu = new Pslvector( dim = frq->GetSize());
                if( mu == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
                if( stds.GetSize() && frq->GetSize() != stds.GetSize())
                    throw MYRUNTIME_ERROR( preamb + "Inconsistent vector sizes.");
            }
            if( frq->GetSize() != dim )
                throw MYRUNTIME_ERROR( preamb + "Inconsistent vector lengths." );
            if(( err = mu->Superposition( 1.0f, *frq )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            nos++;
        }
        if( mu == NULL || nos == 0 || dim == 0 )
            return 0;

        inos = 1.0f / (float)nos;

        //mean vector
        if(( err = mu->MultiplyBy( inos )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        Pslvector   dev( dim );
        Pslvector   adev( dim );
        Pslmatrix   Sdev( dim, dim );

        S0 = new SPDmatrix( dim );
        if( S0 == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        S0inv = new SPDmatrix( dim );
        if( S0inv == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        for( n = 0; n < GetBasin()->GetSize(); n++ ) {
            frq = GetBasin()->GetValueAt( n );
            if( frq == NULL )
                continue;

            //deviation of sample from the mean vector
            dev = *frq;
            if( stds.GetSize())
                dev.AddGNoise( stds );
            if(( err = dev.Superposition( -1.0f, *mu )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            //sum of deviations
            if(( err = adev.Superposition( 1.0f, dev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = S0->Add( Sdev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }

        //square of sum of deviations
        if(( err = Sdev.Mul( adev, adev )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        //divide by the number of observations
        if(( err = Sdev.Scale( inos )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        //correct the round off error of covariance matrix S0
        if(( err = S0->Sub( Sdev )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        //S0 is now covariance matrix multiplied by the number of observations;
        //multiplication below is valid only if S0 is cov. matrix;
        ////scale covariance matrix S0 to obtain valid NIW dist. parameter of scale matrix
        ////if(( err = S0->Scale( GetDefKappaParam())) != 0 )
        ////    throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            //NOTE:divide scale matrix by the number of observations
            if(( err = S0->Scale( inos )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        if(( err = S0->Scale( GetS0ScaleFac())) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        //SAVE scatter matrix here
        GetMenu()->SetS0ScaleFac( GetS0ScaleFac());
        GetMenu()->SetS0( S0 );
        //solve for inverse matrix and calculate determinant
        *S0inv = *S0;
        if(( err = S0inv->CholeskyDecompose()) == 0 ) {
            if(( err = S0inv->CDedLogDet( &ldet )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = S0inv->CDedInvert()) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }
        else {
            return err;//return for adding noise
            SPDmatrix   S0invhlp( dim );//for LU-decomposition
            S0invhlp = *S0;
            if(( err = S0invhlp.LUDecompose()) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = S0invhlp.LUDedLogDet( &ldet )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = S0invhlp.LUDedInvert( *S0inv )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }


        //SAVE
        /*GetMenu()->SetS0( S0 );*/ S0 = NULL;//saved above
        GetMenu()->SetInvS0( S0inv ); S0inv = NULL;
        GetMenu()->SetLDetS0( ldet );
        GetMenu()->SetMu0( mu ); mu = NULL;
        GetMenu()->SetKappa0_pp_a( GetDefKappa_pp_a());
        GetMenu()->SetKappa0_pp_b( GetDefKappa_pp_b());
        if( !GetMenu()->GetKappa0())
            GetMenu()->SetKappa0( inos / GetS0ScaleFac()/*1*//*dim*//*GetDefKappaParam()*/);
        if( !GetMenu()->GetNu0())
            GetMenu()->SetNu0((float)GetBasin()->GetSize()/*1*//*dim*//*GetDefNuParam()*/);
        GetMenu()->SetDim( dim );
        GetMenu()->SetCtx( GetCtxtSize());

        //probability factors
        CalcPriorProbFact();

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mu ) { delete mu; mu = NULL; }
    if( S0 ) { delete S0; S0 = NULL; }
    if( S0inv ) { delete S0inv; S0inv = NULL; }

    if( mre.isset())
        throw mre;

    return 0;
}

// -------------------------------------------------------------------------
// SetUninfPriorParamsS0I: set uninformative prior parameters of NIW 
//     distribution; scale matrix S0 is proportional to I
//
int HDPbase::SetUninfPriorParamsS0I( int dim, int ctx, const Pslvector& mvec, 
                                     float kp0preset, float nu0preset )
{
    mystring preamb = "HDPbase::SetUninfPriorParamsS0I: ";
    myruntime_error mre;

    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null menu.");
    if( GetS0ScaleFac() <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale factor.");
    if( dim < 1 || ctx < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( dim != mvec.GetSize())
        throw MYRUNTIME_ERROR( preamb + "Inconsistent dimensionality.");

    Pslvector* mu0 = NULL;
    SPDmatrix* S0 = NULL;
    SPDmatrix* S0inv = NULL;
    float scale = GetS0ScaleFac();
    float ldet = 0.0f;
    int err;

    try {
        mu0 = new Pslvector( mvec );
        S0 = new SPDmatrix( dim );
        S0inv = new SPDmatrix( dim );
        if( mu0 == NULL || S0 == NULL || S0inv == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        S0->SetIdentity();
        S0inv->SetIdentity();

        if(( err = S0->Scale( scale )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        if(( err = S0inv->Scale( 1.0f / scale )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        ldet = (float)dim * logf(scale);

        GetMenu()->SetS0ScaleFac( GetS0ScaleFac());
        GetMenu()->SetS0( S0 ); S0 = NULL;
        GetMenu()->SetInvS0( S0inv ); S0inv = NULL;
        GetMenu()->SetLDetS0( ldet );
        GetMenu()->SetMu0( mu0 ); mu0 = NULL;
        GetMenu()->SetKappa0_pp_a( GetDefKappa_pp_a());
        GetMenu()->SetKappa0_pp_b( GetDefKappa_pp_b());
        if( 0.0f < kp0preset )
                GetMenu()->SetKappa0( kp0preset );
        else    GetMenu()->SetKappa0( 1.0f / GetS0ScaleFac());
        if( 0.0f < nu0preset )
                GetMenu()->SetNu0( nu0preset );
        else    GetMenu()->SetNu0( GetS0ScaleFac());
        GetMenu()->SetDim( dim );
        GetMenu()->SetCtx( ctx );

        //probability factors
        CalcPriorProbFact();

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mu0 ) { delete mu0; mu0 = NULL; }
    if( S0 ) { delete S0; S0 = NULL; }
    if( S0inv ) { delete S0inv; S0inv = NULL; }

    if( mre.isset())
        throw mre;

    return 0;
}

// =========================================================================
// CalcPriorProbFact: recalculate probability factors
//
void HDPbase::CalcPriorProbFact()
{
    mystring preamb = "HDPbase::CalcPriorProbFact: ";

    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null menu.");
    if( GetMenu()->GetKappa0() <= 0.0f || GetMenu()->GetNu0() <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Invalid degrees of freedom.");

    float pfact = 0.0f;
    float term1, term2;
    float operr;
    int err;
    const float do2 = (float)GetMenu()->GetDim() * 0.5f;
    const float nu0p1o2 = (GetMenu()->GetNu0()+1.0f) * 0.5f;
    const float kp0 = GetMenu()->GetKappa0();
    const float kp0p1 = kp0 + 1.0f;

    const float nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    pfact = term1 - term2;
    pfact -= do2 * SLC_LNPI;
    pfact += do2 * ( logf(kp0) - logf(kp0p1));

    GetMenu()->SetLogPriorProbFact( pfact );
}

// =========================================================================
// RecalcMenuParams: recalculate parameters for each dish of the menu
//
void HDPbase::RecalcMenuParams()
{
    if( GetMenu() == NULL )
        return;

    Dish* dish = NULL;
    int k;

    for( k = 0; k < GetMenu()->GetSize(); k++ )
    {
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        RecalcDishParams( k );
    }
}

// -------------------------------------------------------------------------
// RecalcDishParams: recalculate parameters for the given dish of the menu
//
void HDPbase::RecalcDishParams( int k )
{
    mystring preamb = "HDPbase::RecalcDishParams: ";
    myruntime_error mre;

    if( GetBasin() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null basin.");
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null menu.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();

    if( 0.0f < kp0 && ( mu0 == NULL || L0 == NULL ))
        throw MYRUNTIME_ERROR( preamb + "Null prior parameters.");
    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");

    const Pslvector*    val = NULL;
    Pslvector*          mu = NULL;
    SPDmatrix*          L = NULL;
    SPDmatrix*          Linv = NULL;
    Dish*       dish = NULL;
    char        bufstr[KBYTE];
    float       ldet = 0.0f;
    float       inos = 0.0f;
    int         err;
    int         n, nos;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid dish index.");

    dish = GetMenu()->GetDishAt( k );
    if( dish == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null dish k.");

    try {
        mu = new Pslvector( dim );
        if( mu == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
        L = new SPDmatrix( dim );
        if( L == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
#ifdef PROBMTXINVSM
        Linv = new SPDmatrix( dim );
        if( Linv == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
#endif
        for( n = 0, nos = 0; n < dish->GetSize(); n++ ) {
            if( dish->GetVectorNIndAt( n ) < 0 )
                continue;
            val = dish->GetVectorNAt( n );
            if( val == NULL )
                continue;
            if( val->GetSize() != dim )
                throw MYRUNTIME_ERROR( preamb + "Inconsistent vector dimensionality.");
            if(( err = mu->Superposition( 1.0f, *val )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            nos++;
        }
        if( dish->GetActualSize() != nos ) {
            sprintf( bufstr, "Inconsistent size of dish %d.", k );
            throw MYRUNTIME_ERROR( preamb + bufstr );
        }
        if( nos != 0 &&
            nos != dish->GetDishSize()) {
            //when proposal size is defined
            if(( err = mu->MultiplyBy( (float)dish->GetDishSize()/(float)nos )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            nos = dish->GetDishSize();
        }
        if( nos == 0 ) {
            sprintf( bufstr, "Dish %d has no members.", k );
            throw MYRUNTIME_ERROR( preamb + bufstr );
        }
        inos = 1.0f /(float)nos;

        //mean vector
        if(( err = mu->MultiplyBy( inos )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

        Pslvector   eta( *mu );//unweighted mean vector
        float       tot = kp0 + (float)nos;
        float       tpr = kp0 * (float)nos;

        if( 0.0f < kp0 ) {
            Pslvector   m0w( *mu0 );
            //weigh prior mean vector
            if(( err = m0w.MultiplyBy( kp0 / tot )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            //weigh dish mean vector
            if(( err = mu->MultiplyBy( (float)nos / tot )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            //weighted mean vector
            if(( err = mu->Superposition( 1.0f, m0w )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }

        //proceed with SCALE matrix...
        Pslvector   dev( dim );
        Pslmatrix   Sdev( dim, dim );
        SPDmatrix   S( dim );

        //1.calculate S
        for( n = 0; n < dish->GetSize(); n++ ) {
            if( dish->GetVectorNIndAt( n ) < 0 )
                continue;
            val = dish->GetVectorNAt( n );
            if( val == NULL )
                continue;

            //deviation of sample from the mean vector
            dev = *val;
            if(( err = dev.Superposition( -1.0f, eta )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = S.Add( Sdev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }

        *L = S;
        if( 0.0f < kp0 ) {
            //2.term of deviation of means
            dev = eta;
            if(( err = dev.Superposition( -1.0f, *mu0 )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            if(( err = Sdev.Scale( tpr / tot )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

            //3.final parameter scale matrix
            if(( err = L->Add( *L0 )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = L->Add( Sdev )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        }

        //4.solve for inverse matrix and calculate determinant
#ifndef PROBMTXINVSM//if NOT defined
        Linv = &S;
#endif
        *Linv = *L;
        if(( err = Linv->CholeskyDecompose()) == 0 ) {
            if(( err = Linv->CDedLogDet( &ldet )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
#ifdef PROBMTXINVSM
            if(( err = Linv->CDedInvert()) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
#endif
        }
        else {
            SPDmatrix   Linvhlp( dim );//for LU-decomposition
            Linvhlp = *L;
            if(( err = Linvhlp.LUDecompose()) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = Linvhlp.LUDedLogDet( &ldet )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
#ifdef PROBMTXINVSM
            if(( err = Linvhlp.LUDedInvert( *Linv )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
#endif
        }

        //SAVE
        GetMenu()->SetSMatrixAt( k, L ); L = NULL;
#ifdef PROBMTXINVSM
        GetMenu()->SetInvSMatrixAt( k, Linv ); Linv = NULL;
#endif
        GetMenu()->SetLDetSMAt( k, ldet );
        GetMenu()->SetMuVectorAt( k, mu ); mu = NULL;

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mu ) { delete mu; mu = NULL; }
    if( L ) { delete L; L = NULL; }
#ifdef PROBMTXINVSM
    if( Linv ) { delete Linv; Linv = NULL; }
#endif
    if( mre.isset())
        throw mre;
}





// =========================================================================
// CalcProbOfData: calculate log probability of all data given current 
//  configuration of clusters
//
void HDPbase::CalcProbOfData( float* lprob )
{
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR("HDPbase::CalcProbOfData: Null menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR("HDPbase::CalcProbOfData: Memory access error.");

    Dish* dish = NULL;
    float lp;
    int k;

    *lprob = 0.0f;
    for( k = 0; k < GetMenu()->GetSize(); k++ ) {
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        CalcProbOfDish( k, &lp );
        GetMenu()->SetProbAt( k, lp );
        *lprob += lp;
    }
    SetLastLProbData( *lprob );
}

// =========================================================================
// CalcProbOfDish: calculate log probability of dish k
//
void HDPbase::CalcProbOfDish( int k, float* lprob )
{
    mystring preamb = "HDPbase::CalcProbOfDish: ";
    int err;

    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null menu.");

    *lprob = -1.e6f;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid dish index.");

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const float         ldet = GetMenu()->GetLDetSMAt( k );
    const float         L0ldet = GetMenu()->GetLDetS0();

    if( dishk == NULL || mu == NULL || L == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null cluster data.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetActualSize();
    const float         kp = kp0 + (float)nk;
    const float         nu = nu0 + (float)nk;

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( nk <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Dish has no members.");
    if( kp0 <= 0.0f || nu0 <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Invalid prior parameters kappa and nu.");

    float lnsum = 0.0f;
    float term1, term2;
    float operr;
    float j;
    const float do2 = (float)dim * 0.5f;
    const float nuo2 = nu * 0.5f;
    const float nko2 = (float)nk * 0.5f;

    const float nu_t_o2 = nuo2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const float dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5f;//(deg.of freedom+1 + dim.) over 2
    const float dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - nko2;//(deg.of freedom+1 + dim.- L) / 2
    const float dofdom_pdim_mv_o2 = dofdom_pdim_o2 - nko2;//(deg.of freedom + dim.- L) over 2

    //calculate probability next
    if( nk == 1 ) {
        if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
        lnsum += term1 - term2;
    }
    else
        for( j = 0.5f; j <= do2; j += 0.5f )
        {
            if(( err = psl_lngamma_e( dofdom_pdim1_o2 - j, &term1, &operr )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            if(( err = psl_lngamma_e( dofdom_pdim1_mv_o2 - j, &term2, &operr )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
            lnsum += term1 - term2;
        }
    lnsum -= (float)nk * do2 * SLC_LNPI;
    lnsum += do2 * ( logf(kp0) - logf(kp));
    lnsum += dofdom_pdim_mv_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
}





// =========================================================================
// TestMenu: test menu structure
//
bool HDPbase::TestMenu() const
{
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    const Restaurant*   rest = NULL;
    const Table*        tbl = NULL;
    const Dish*         dsh = NULL;
    const Pslvector*    frq = NULL;
    mystring    preamb = "HDPbase::TestMenu: ";
    mystring    merror;
    const int   nods = GetMenu()->GetSize();
    int*    dishs = NULL;
    int*    vecs = NULL;
    int     notbls = 0, novecs = 0;
    int     nts = 0, nvs = 0;
    int     ddx;
    int     r, t, k, v;

    try {
        dishs = ( int* )malloc( nods * sizeof( int ));
        vecs = ( int* )malloc( nods * sizeof( int ));
        if( dishs == NULL || vecs == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        memset( dishs, 0, nods * sizeof( int ));
        memset( vecs, 0, nods * sizeof( int ));

        for( r = 0; r < GetChain()->GetSize(); r++ ) {
            rest = GetChain()->GetRestaurantAt( r );
            if( rest == NULL )
                throw MYRUNTIME_ERROR( preamb + "Null restaurant.");
            for( t = 0; t < rest->GetSize(); t++ ) {
                tbl = rest->GetTableAt( t );
                if( tbl == NULL )
                    continue;
                k = tbl->GetDishIndex();
                dsh = tbl->GetDish();
                if( k < 0 || nods <= k || dsh == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Invalid table's dish index.");
                dishs[k]++;
                vecs[k] += tbl->GetActualSize();
                notbls++;
                novecs += tbl->GetActualSize();
                for( v = 0; v < tbl->GetSize(); v++ ) {
                    if( tbl->GetVectorNIndAt( v ) < 0 )
                        continue;
                    frq = tbl->GetVectorNAt( v );
                    ddx = tbl->GetVectorNDishIndAt( v );
                    if( ddx < 0 || dsh->GetSize() <= ddx )
                        throw MYRUNTIME_ERROR( preamb + "Invalid dish index saved in table vector.");
                    if( dsh->GetVectorNAt( ddx ) != frq )
                        throw MYRUNTIME_ERROR( preamb + "Table vector not found in dish.");
                }
            }
        }
        for( k = 0; k < nods; k++ ) {
            dsh = GetMenu()->GetDishAt( k );
            if( dsh == NULL )
                continue;
            if( dsh->GetNoTables() < 1 || dishs[k] != dsh->GetNoTables())
                throw MYRUNTIME_ERROR( preamb + "Dish has invalid number of tables.");
            if( dsh->GetActualSize() < 1 || vecs[k] != dsh->GetActualSize())
                throw MYRUNTIME_ERROR( preamb + "Dish has invalid number of tables.");
            nts += dsh->GetNoTables();
            nvs += dsh->GetActualSize();
        }
        if( notbls != nts )
            throw MYRUNTIME_ERROR( preamb + "Inconsistent number of tables.");
        if( novecs != nvs )
            throw MYRUNTIME_ERROR( preamb + "Inconsistent number of vectors.");

    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    if( dishs ) {
        free( dishs );
        dishs = NULL;
    }
    if( vecs ) {
        free( vecs );
        vecs = NULL;
    }
    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    return true;
}

