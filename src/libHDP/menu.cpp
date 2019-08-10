/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "menu.h"

using namespace extspsl;

// -------------------------------------------------------------------------
// constructor: initialization
//
Menu::Menu( int size )
:   values_( NULL ),
    probs_( NULL ),
    pprobs_( NULL ),
    npprob_( 0.0f ),
    mu_( NULL ),
    lambda_( NULL ),
#ifdef PROBMTXINVSM
    lambdainv_( NULL ),
#endif
    ldetlambda_( NULL ),
    lmvrprobfact_( 0.0f ),
    lmtxprobfact_( 0.0f ),
    S0scalefac_( 1.0f ),
    S0_( NULL ),
    S0inv_( NULL ),
    ldetS0_( 0.0f ),
    m0_( NULL ),
    a_k0_( 0.0f ),
    b_k0_( 0.0f ),
    k0_( 0.0f ),
    v0_( 0.0f ),
    dim_( 0 ),
    ctx_( 0 ),
    notbls_( 0 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( true ),
    destroypriors_( true ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: copy
//
Menu::Menu( const Menu& right )
:   values_( NULL ),
    probs_( NULL ),
    pprobs_( NULL ),
    npprob_( 0.0f ),
    mu_( NULL ),
    lambda_( NULL ),
#ifdef PROBMTXINVSM
    lambdainv_( NULL ),
#endif
    ldetlambda_( NULL ),
    lmvrprobfact_( 0.0f ),
    lmtxprobfact_( 0.0f ),
    S0scalefac_( 1.0f ),
    S0_( NULL ),
    S0inv_( NULL ),
    ldetS0_( 0.0f ),
    m0_( NULL ),
    a_k0_( 0.0f ),
    b_k0_( 0.0f ),
    k0_( 0.0f ),
    v0_( 0.0f ),
    dim_( 0 ),
    ctx_( 0 ),
    notbls_( 0 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( true ),
    destroypriors_( true ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
Menu::Menu()
:   values_( NULL ),
    probs_( NULL ),
    pprobs_( NULL ),
    npprob_( 0.0f ),
    mu_( NULL ),
    lambda_( NULL ),
#ifdef PROBMTXINVSM
    lambdainv_( NULL ),
#endif
    ldetlambda_( NULL ),
    lmvrprobfact_( 0.0f ),
    lmtxprobfact_( 0.0f ),
    S0scalefac_( 1.0f ),
    S0_( NULL ),
    S0inv_( NULL ),
    ldetS0_( 0.0f ),
    m0_( NULL ),
    a_k0_( 0.0f ),
    b_k0_( 0.0f ),
    k0_( 0.0f ),
    v0_( 0.0f ),
    dim_( 0 ),
    ctx_( 0 ),
    notbls_( 0 ),
    actlen_( 0 ),
    length_( 0 ),
    capacity_( 0 ),
    destroy_( true ),
    destroypriors_( true ),
    vacancies_( NULL ),
    novacs_( 0 ),
    capvac_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
Menu::~Menu()
{
    Destroy();
    if( S0_ && GetDestroyPriors()) 
        delete S0_; 
    S0_ = NULL;
    if( S0inv_ && GetDestroyPriors()) 
        delete S0inv_; 
    S0inv_ = NULL;
    if( m0_ && GetDestroyPriors()) 
        delete m0_; 
    m0_ = NULL;
}

// -------------------------------------------------------------------------
// operator=: assignment
//
Menu& Menu::operator=( const Menu& right )
{
    HardCopy( right );
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Menu::Realloc( int newcap )
{
    Dish**      tmp_values = NULL;
    float*      tmp_probs = NULL;
    float*      tmp_pprobs = NULL;
    Pslvector** tmp_mu = NULL;
    SPDmatrix** tmp_lambda = NULL;
#ifdef PROBMTXINVSM
    SPDmatrix** tmp_lambdainv = NULL;
#endif
    float*      tmp_ldetlambda = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ <= 0 ) {
        tmp_values = ( Dish** )malloc( sizeof( void* ) * newcap );
        tmp_probs = (float*)malloc( sizeof( float ) * newcap );
        tmp_pprobs = (float*)malloc( sizeof( float ) * newcap );
        tmp_mu = ( Pslvector** )malloc( sizeof( void* ) * newcap );
        tmp_lambda = ( SPDmatrix** )malloc( sizeof( void* ) * newcap );
#ifdef PROBMTXINVSM
        tmp_lambdainv = ( SPDmatrix** )malloc( sizeof( void* ) * newcap );
#endif
        tmp_ldetlambda = ( float* )malloc( sizeof( float ) * newcap );
    } else {
        tmp_values = ( Dish** )realloc( values_, sizeof( void* ) * newcap );
        tmp_probs = ( float* )realloc( probs_, sizeof( float ) * newcap );
        tmp_pprobs = ( float*)realloc( pprobs_, sizeof( float ) * newcap );
        tmp_mu = ( Pslvector** )realloc( mu_, sizeof( void* ) * newcap );
        tmp_lambda = ( SPDmatrix** )realloc( lambda_, sizeof( void* ) * newcap );
#ifdef PROBMTXINVSM
        tmp_lambdainv = ( SPDmatrix** )realloc( lambdainv_, sizeof( void* ) * newcap );
#endif
        tmp_ldetlambda = ( float* )realloc( ldetlambda_, sizeof( float ) * newcap );
    }

    if( !tmp_values || !tmp_probs || !tmp_pprobs ||
        !tmp_mu || !tmp_lambda || !tmp_ldetlambda )
        throw MYRUNTIME_ERROR("Menu::Realloc: Not enough memory.");
#ifdef PROBMTXINVSM
    if( !tmp_lambdainv )
        throw MYRUNTIME_ERROR("Menu::Realloc: Not enough memory.");
#endif

    values_ = tmp_values;
    probs_ = tmp_probs;
    pprobs_ = tmp_pprobs;
    mu_ = tmp_mu;
    lambda_ = tmp_lambda;
#ifdef PROBMTXINVSM
    lambdainv_ = tmp_lambdainv;
#endif
    ldetlambda_ = tmp_ldetlambda;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
    memset( probs_ + capacity_, 0, sizeof( float ) * ( newcap - capacity_ ));
    memset( pprobs_ + capacity_, 0, sizeof( float ) * ( newcap - capacity_ ));
    memset( mu_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
    memset( lambda_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
#ifdef PROBMTXINVSM
    memset( lambdainv_ + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));
#endif
    memset( ldetlambda_ + capacity_, 0, sizeof( float ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
//
void Menu::ReallocVacans( int newcap )
{
    int*    tmp_vacans = NULL;

    if( newcap <= capvac_ )
        return;

    if( capvac_ == 0 ) {
        tmp_vacans = ( int* )malloc( sizeof( int ) * newcap );
    } else {
        tmp_vacans = ( int* )realloc( vacancies_, sizeof( int ) * newcap );
    }

    if( !tmp_vacans )
        throw MYRUNTIME_ERROR("Menu::ReallocVacans: Not enough memory.");

    vacancies_ = tmp_vacans;

    // fill uninitialized memory with zeros
    memset( vacancies_ + capvac_, 0xff, sizeof( int ) * ( newcap - capvac_ ));
    capvac_ = newcap;
}

// -------------------------------------------------------------------------
// SoftCopy: copy addresses of argument vector to this vector; manage cases
//     when lengths of vectors are not equal
//
void Menu::SoftCopy( const Menu& vector )
{
    Clear();
    SetDestroy( false );

    SetActualSize( 0 );
    ReserveDishes( vector.GetSize());
    ReserveVacans( vector.GetNoVacans());

    SoftCopyPriors( vector );

    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;
    if( GetCapVacans() < vector.GetNoVacans())
        return;

    if( GetSize() < vector.GetSize())
        SetSize( vector.GetSize());

    int noelms = GetSize();
    int n;

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( n = 0; n < noelms; n++ ) {
        SetDishAt( n, vector.GetDishAt( n ));
        SetProbAt( n, vector.GetProbAt( n ));
        SetPriorProbAt( n, vector.GetPriorProbAt( n ));
        SetMuVectorAt( n, vector.GetMuVectorAt( n ));
        SetSMatrixAt( n, vector.GetSMatrixAt( n ));
#ifdef PROBMTXINVSM
        SetInvSMatrixAt( n, vector.GetInvSMatrixAt( n ));
#endif
        SetLDetSMAt( n, vector.GetLDetSMAt( n ));
    }

    SetNoTables( vector.GetNoTables());

    if( GetNoVacans() < vector.GetNoVacans())
        SetNoVacans( vector.GetNoVacans());

    noelms = GetNoVacans();

    if( vector.GetNoVacans() < noelms )
        noelms = vector.GetNoVacans();

    for( n = 0; n < noelms; n++ )
        SetVacantAt( n, vector.GetVacantAt( n ));

    SetActualSize( vector.GetActualSize());
}

// -------------------------------------------------------------------------
// SoftCopyPriors: copy addresses of prior parameters
//
void Menu::SoftCopyPriors( const Menu& vector )
{
    SetLogPriorProbFact( vector.GetLogPriorProbFact());
    SetLogMtxPriorProbFact( vector.GetLogMtxPriorProbFact());
    SetPriorProbNewDish( vector.GetPriorProbNewDish());
    SetS0ScaleFac( vector.GetS0ScaleFac());
    SetS0( vector.GetS0());
    SetInvS0( vector.GetInvS0());
    SetLDetS0( vector.GetLDetS0());
    SetMu0( vector.GetMu0());
    SetKappa0_pp_a( vector.GetKappa0_pp_a());
    SetKappa0_pp_b( vector.GetKappa0_pp_b());
    SetKappa0( vector.GetKappa0());
    SetNu0( vector.GetNu0());
    SetDim( vector.GetDim());
    SetCtx( vector.GetCtx());
    SetDestroyPriors( false );
}

// -------------------------------------------------------------------------
// HardCopy: copy elements of argument vector to this vector; manage cases
//     when lengths of vectors are not equal
//
void Menu::HardCopy( const Menu& vector )
{
    Dish*       dk = NULL;
    Pslvector*  mu = NULL;
    SPDmatrix*  sm = NULL;
#ifdef PROBMTXINVSM
    SPDmatrix*  smi = NULL;
#endif

    Clear();
    SetDestroy( true );

    SetActualSize( 0 );
    ReserveDishes( vector.GetSize());
    ReserveVacans( vector.GetNoVacans());

    HardCopyPriors( vector );

    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;
    if( GetCapVacans() < vector.GetNoVacans())
        return;

    if( GetSize() < vector.GetSize())
        SetSize( vector.GetSize());

    int noelms = GetSize();
    int n;

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    for( n = 0; n < noelms; n++ ) {
        dk = NULL;
        mu = NULL;
        sm = NULL;
#ifdef PROBMTXINVSM
        smi = NULL;
#endif
        if( vector.GetDishAt( n ))
            dk = new Dish( *vector.GetDishAt( n ));
        SetDishAt( n, dk );
        SetProbAt( n, vector.GetProbAt( n ));
        SetPriorProbAt( n, vector.GetPriorProbAt( n ));
        if( vector.GetMuVectorAt( n ))
            mu = new Pslvector( *vector.GetMuVectorAt( n ));
        SetMuVectorAt( n, mu );
        if( vector.GetSMatrixAt( n ))
            sm = new SPDmatrix( *vector.GetSMatrixAt( n ));
        SetSMatrixAt( n, sm );
#ifdef PROBMTXINVSM
        if( vector.GetInvSMatrixAt( n ))
            smi = new SPDmatrix( *vector.GetInvSMatrixAt( n ));
        SetInvSMatrixAt( n, smi );
#endif
        SetLDetSMAt( n, vector.GetLDetSMAt( n ));
    }

    SetNoTables( vector.GetNoTables());

    if( GetNoVacans() < vector.GetNoVacans())
        SetNoVacans( vector.GetNoVacans());

    noelms = GetNoVacans();

    if( vector.GetNoVacans() < noelms )
        noelms = vector.GetNoVacans();

    for( n = 0; n < noelms; n++ )
        SetVacantAt( n, vector.GetVacantAt( n ));

    SetActualSize( vector.GetActualSize());
}

// -------------------------------------------------------------------------
// HardCopyPriors: copy elements of prior parameters
//
void Menu::HardCopyPriors( const Menu& vector )
{
    Pslvector*  m = NULL;
    SPDmatrix*  s = NULL;
    SPDmatrix*  si = NULL;
    SetLogPriorProbFact( vector.GetLogPriorProbFact());
    SetLogMtxPriorProbFact( vector.GetLogMtxPriorProbFact());
    SetPriorProbNewDish( vector.GetPriorProbNewDish());
    SetS0ScaleFac( vector.GetS0ScaleFac());
    if( vector.GetS0()) s = new SPDmatrix( *vector.GetS0());
    SetS0( s );
    if( vector.GetInvS0()) si = new SPDmatrix( *vector.GetInvS0());
    SetInvS0( si );
    SetLDetS0( vector.GetLDetS0());
    if( vector.GetMu0()) m = new Pslvector( *vector.GetMu0());
    SetMu0( m );
    SetKappa0_pp_a( vector.GetKappa0_pp_a());
    SetKappa0_pp_b( vector.GetKappa0_pp_b());
    SetKappa0( vector.GetKappa0());
    SetNu0( vector.GetNu0());
    SetDim( vector.GetDim());
    SetCtx( vector.GetCtx());
    SetDestroyPriors( true );
}

// -------------------------------------------------------------------------
// Clear: clear all elements
//
void Menu::Clear()
{
    int n;
    for( n = 0; n < GetSize(); n++ ) {
        if( GetDestroy() && GetDishAt( n ))
            delete GetDishAt( n );
        SetDishAt( n, NULL );
        SetProbAt( n, 0.0f );
        SetPriorProbAt( n, 0.0f );
        if( GetDestroy() && GetMuVectorAt( n ))
            delete GetMuVectorAt( n );
        SetMuVectorAt( n, NULL );
        if( GetDestroy() && GetSMatrixAt( n ))
            delete GetSMatrixAt( n );
        SetSMatrixAt( n, NULL );
#ifdef PROBMTXINVSM
        if( GetDestroy() && GetInvSMatrixAt( n ))
            delete GetInvSMatrixAt( n );
        SetInvSMatrixAt( n, NULL );
#endif
        SetLDetSMAt( n, 0.0f );
    }

    SetLogPriorProbFact( 0.0f );
    SetLogMtxPriorProbFact( 0.0f );
    SetPriorProbNewDish( 0.0f );
    SetS0ScaleFac( 1.0f );
    SetS0( NULL );
    SetInvS0( NULL );
    SetLDetS0( 0.0f );
    SetMu0( NULL );
    SetKappa0_pp_a( 0.0f );
    SetKappa0_pp_b( 0.0f );
    SetKappa0( 0.0f );
    SetNu0( 0.0f );

    SetNoTables( 0 );

    if( GetVacancies())
        for( n = 0; n < GetNoVacans(); n++ )
            SetVacantAt( n, -1 );
    SetSize( 0 );
    SetNoVacans( 0 );
    SetDestroy( true );
    SetDestroyPriors( true );
}
