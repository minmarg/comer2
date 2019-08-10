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
#include "liblib/alpha.h"
#include "HDPbase.h"

using namespace extspsl;

// =========================================================================
// PriorProbVec: compute prior probability of t-distributed vector
//
void HDPbase::PriorProbVec( const Pslvector* vec, float* lprob ) const
{
    PriorProbVecHlpObs( vec, lprob );
//     PriorProbVecHlp( vec, lprob );
}

// -------------------------------------------------------------------------
// PriorProbVecHlpObs: compute prior probability of t-distributed vector;
//  a version requiring the inverse of scale matrix
//
void HDPbase::PriorProbVecHlpObs( const Pslvector* vec, float* lprob ) const
{
    mystring    preamb = "HDPbase::PriorProbVecHlpObs: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const SPDmatrix*    L0inv = GetMenu()->GetInvS0();
    const float         L0ldet = GetMenu()->GetLDetS0();
    const float         lpfact = GetMenu()->GetLogPriorProbFact();

    if( mu0 == NULL || L0 == NULL || L0inv == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( kp0 <= 0.0f || nu0 <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Too small number of degrees of freedom.");

    float lnsum = 0.0f;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu0 * 0.5f;
    const float nu0p1o2 = no2 + 0.5f;//(nu+1)/2
    const float kp01 = kp0 + 1.0f;

    const float nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;

    dev = *vec;
    if(( err = dev.Superposition( -1.0f, *mu0 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = dev.Transpose( trv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = pdt.Mul( trv, *L0inv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if( scl.GetValueAt( 0 ) <= -kp01 / kp0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid determinant." );

    lnsum = lpfact;
    lnsum -= 0.5f * L0ldet;
    lnsum -= dofdom_pdim_o2 * logf( 1.0f + kp0 / kp01  * scl.GetValueAt( 0 ));

    *lprob = lnsum;
//     *prob = expf( lnsum );
}

// -------------------------------------------------------------------------
// PriorProbVecHlp: compute prior probability of t-distributed vector;
//
void HDPbase::PriorProbVecHlp( const Pslvector* vec, float* lprob ) const
{
    mystring    preamb = "HDPbase::PriorProbVecHlp: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const float         L0ldet = GetMenu()->GetLDetS0();
    const float         lpfact = GetMenu()->GetLogPriorProbFact();

    if( mu0 == NULL || L0 == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( kp0 <= 0.0f || nu0 <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Too small number of degrees of freedom.");

    float lnsum = 0.0f;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu0 * 0.5f;
    const float nu0p1o2 = no2 + 0.5f;//(nu+1)/2
    const float kp01 = kp0 + 1.0f;

    const float nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const float dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5f;//(deg.of freedom-1 + dim.) over 2

    Pslvector   dev;
    SPDmatrix   MCM( dim );
    float       ldetMCM = 0.0f;

    dev = *vec;
    if(( err = dev.Superposition( -1.0f, *mu0 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = MCM.Mul( dev, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = MCM.Scale( kp0 / kp01 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = MCM.Add( *L0 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    //factorize by Cholesky decomp. and calculate determinant
    if(( err = MCM.CholeskyDecompose()) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = MCM.CDedLogDet( &ldetMCM )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));


    lnsum = lpfact;
    lnsum += dofdom_pdim_m1_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = expf( lnsum );
}





// =========================================================================
// ProbVecOfDish: compute probability of t-distributed vector to get dish k
//  (to belong to global cluster k)
//
void HDPbase::ProbVecOfDish( const Pslvector* vec, int k, float* lprob ) const
{
#ifdef PROBMTXINVSM
    ProbVecOfDishHlpObs( vec, k, lprob );
#else
    ProbVecOfDishHlp( vec, k, lprob );
#endif
}

#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbVecOfDishHlpObs: compute probability of t-distributed vector to get 
//     dish k; requires the inverse of scale matrix
//
void HDPbase::ProbVecOfDishHlpObs( const Pslvector* vec, int k, float* lprob ) const
{
    mystring    preamb = "HDPbase::ProbVecOfDishHlpObs: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid cluster index.");

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const SPDmatrix*    Linv = GetMenu()->GetInvSMatrixAt( k );
    const float         ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL || Linv == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const float         kp = kp0 + nk;
    const float         nu = nu0 + nk;

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( nk <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Cluster has no members.");
    if( kp <= 0.0f || nu <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Cluster has too few members.");

    float lnsum = 0.0f;
    float term1, term2;
    float operr;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu * 0.5f;
    const float nup1o2 = no2 + 0.5f;//(nu+1)/2
    const float kp1 = kp + 1.0f;

    const float nu_t_o2 = nup1o2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;

    dev = *vec;
    if(( err = dev.Superposition( -1.0f, *mu )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = dev.Transpose( trv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = pdt.Mul( trv, *Linv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError(err));

    if( scl.GetValueAt( 0 ) <= -kp1 / kp )
        throw MYRUNTIME_ERROR( preamb + "Invalid determinant.");

    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( logf( kp ) - logf( kp1 ));
    lnsum -= 0.5f * ldet;
    lnsum -= dofdom_pdim_o2 * logf( 1.0f + kp / kp1  * scl.GetValueAt( 0 ));

    *lprob = lnsum;
//     *prob = expf( lnsum );
}

#endif

// -------------------------------------------------------------------------
// ProbVecOfDishHlp: compute probability of t-distributed vector to get 
//     dish k (to belong to global cluster k);
//
void HDPbase::ProbVecOfDishHlp( const Pslvector* vec, int k, float* lprob ) const
{
    mystring    preamb = "HDPbase::ProbVecOfDishHlp: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin.");
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid cluster index.");

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const float         ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const float         kp = kp0 + nk;
    const float         nu = nu0 + nk;

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
//     if( nk <= 0 )
//         throw MYRUNTIME_ERROR( preamb + "Cluster has no members.");
    if( kp <= 0.0f || nu <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Cluster has too few members.");

    float lnsum = 0.0f;
    float term1, term2;
    float operr;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu * 0.5f;
    const float nup1o2 = no2 + 0.5f;//(nu+1)/2
    const float kp1 = kp + 1.0f;

    const float nu_t_o2 = nup1o2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const float dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5f;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    Pslvector   dev;
    SPDmatrix   MCM( dim );
    float       ldetMCM = 0.0f;

    dev = *vec;
    if(( err = dev.Superposition( -1.0f, *mu )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = MCM.Mul( dev, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = MCM.Scale( kp / kp1 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = MCM.Add( *L )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //factorize by Cholesky decomp. and calculate determinant
    if(( err = MCM.CholeskyDecompose()) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = MCM.CDedLogDet( &ldetMCM )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));


    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( logf( kp ) - logf( kp1 ));
    lnsum += dofdom_pdim_m1_o2 * ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = expf( lnsum );
}





// =========================================================================
// ProbVecOfDishExc: compute probability of t-distributed vector to get 
//  dish k, excluding that vector
//
void HDPbase::ProbVecOfDishExc( const Pslvector* vec, int k, float* lprob ) const
{
#ifdef PROBMTXINVSM
    ProbVecOfDishExcHlpObs( vec, k, lprob );
#else
    ProbVecOfDishExcHlp( vec, k, lprob );
#endif
}

#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbVecOfDishExcHlpObs: compute probability of t-distributed member 
//  vector to get dish k, excluding the given member vector;
//  requires the inverse of scale matrix;
// NOTE: member vector n is supposed to belong to dish k already
//
void HDPbase::ProbVecOfDishExcHlpObs( const Pslvector* vec, int k, float* lprob ) const
{
    mystring    preamb = "HDPbase::ProbVecOfDishExcHlpObs: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid cluster index.");

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const SPDmatrix*    Linv = GetMenu()->GetInvSMatrixAt( k );
    const float         ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL || Linv == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const float         kp = kp0 + nk;
    const float         nu = nu0 + nk;

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( nk <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Cluster has no members.");
    if( kp <= 0.0f || nu <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Cluster has too few members.");

    if( nk == 1 ) {
        //dish has only one member; probability equals the prior
        PriorProbVec( vec, lprob );
        return;
    }

    float lnsum = 0.0f;
    float term1, term2;
    float operr;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu * 0.5f;
    const float num1o2 = no2 - 0.5f;//(nu-1)/2
    const float km1 = kp - 1.0f;

    const float nu_t_o2 = no2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const float dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5f;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    Pslvector   muexc;//mean vector of dish k, excluding vec
    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;//scale matrix of dish k, excluding vec

    //mu_{n_k-1}
    muexc = *mu;
    if(( err = muexc.MultiplyBy( kp / km1 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = muexc.Superposition( -1.0f / km1, *vec )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //vec-mu_{n_k-1}
    dev = *vec;
    if(( err = dev.Superposition( -1.0f, muexc )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //|L_{n_k-1}|
    if(( err = dev.Transpose( trv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = pdt.Mul( trv, *Linv )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if( scl.GetValueAt( 0 ) <= kp / km1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid determinant.");

    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( logf( km1 ) - logf( kp ));
    lnsum -= dofdom_pdim_o2 * logf( 1.0f + km1 / kp * scl.GetValueAt( 0 ));
    lnsum -= 0.5f * ldet;

    *lprob = lnsum;
//     *prob = expf( lnsum );
}

#endif

// -------------------------------------------------------------------------
// ProbVecOfDishExcHlp: compute probability of t-distributed member 
//  vector to get dish k, excluding the given vector;
// NOTE: member vector n is supposed to belong to dish k already
//
void HDPbase::ProbVecOfDishExcHlp( const Pslvector* vec, int k, float* lprob ) const
{
    mystring    preamb = "HDPbase::ProbVecOfDishExcHlp: ";
    mystring    errstr;
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null Menu.");
    if( lprob == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    if( vec == NULL )
        throw MYRUNTIME_ERROR( preamb + "Invalid vector address.");

    *lprob = -1.e6f;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw MYRUNTIME_ERROR( preamb + "Invalid cluster index.");

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const float         ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    const int           dim = GetMenu()->GetDim();
    const float         kp0 = GetMenu()->GetKappa0();
    const float         nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const float         kp = kp0 + nk;
    const float         nu = nu0 + nk;

    if( dim <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality.");
    if( nk <= 0 )
        throw MYRUNTIME_ERROR( preamb + "Cluster has no members.");
    if( kp <= 0.0f || nu <= 0.0f/*(float)dim*/ )
        throw MYRUNTIME_ERROR( preamb + "Cluster has too few members.");

    if( nk == 1 ) {
        //dish has only one member; probability equals the prior
        PriorProbVec( vec, lprob );
        return;
    }

    float lnsum = 0.0;
    float term1, term2;
    float operr;
    const float do2 = (float)dim * 0.5f;
    const float no2 = nu * 0.5f;
//     const float num1o2 = no2 - 0.5f;//(nu-1)/2
    const float km1 = kp - 1.0f;

    const float nu_t_o2 = no2;//student's t nu over 2

//     const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_o2 = NU_t_o2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const float dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5f;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    Pslvector   muexc;//mean vector of dish k excluded vec
    Pslvector   dev;
    SPDmatrix   Lexc( dim );//scale matrix of dish k excluded vec
    float       Lexcldet = 0.0f;

    //mu_{n_k-1}
    muexc = *mu;
    if(( err = muexc.MultiplyBy( kp / km1 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = muexc.Superposition( -1.0f / km1, *vec )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //vec-mu_{n_k-1}
    dev = *vec;
    if(( err = dev.Superposition( -1.0f, muexc )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //L_{n_k-1}
    if(( err = Lexc.Mul( dev, dev )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = Lexc.Scale( -km1 / kp )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = Lexc.Add( *L )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    //factor by Cholesky decomp. and calculate determinant
    if(( err = Lexc.CholeskyDecompose()) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));

    if(( err = Lexc.CDedLogDet( &Lexcldet )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslatePSLError( err ));


    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( logf( km1 ) - logf( kp ));
    lnsum += dofdom_pdim_m1_o2 * Lexcldet;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
//     *prob = expf( lnsum );
}
