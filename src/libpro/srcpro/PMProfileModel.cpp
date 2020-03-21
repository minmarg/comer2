/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <math.h>
#include <cmath>

#include <fstream>

#include "extsp/psl.h"
#include "extsp/gamma.h"
#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "liblib/logitnormal.h"
#include "libpro/srcpro/MOptions.h"
#include "liblib/CtxtCoefficients.h"
#include "libHDP/HDPbase.h"
#include "SUBSTABLE.h"
#include "TCTXVECT.h"
#include "PMTransModel.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "PMProfileModelBase.h"
#include "PMProfileModel.h"

// #include "Serializer.h"

namespace pmodel {

// COMER profile format version
const char* PMProfileModel::dataversion = "2.0";
const char* PMProfileModel::bindataversion = "2.2";

////////////////////////////////////////////////////////////////////////////
// CLASS PMProfileModel
//
// Constructor
// 
PMProfileModel::PMProfileModel()
{
    private_buffer_ = NULL;
    private_buffer_ = ( char* )malloc( DescBufferSize * sizeof(char));
    if( private_buffer_ == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::PMProfileModel: Not enough memory.");
    init();
}

// Destructor
// 
PMProfileModel::~PMProfileModel()
{
    if( private_buffer_ )
        free( private_buffer_ );
    destroy();
}

// -------------------------------------------------------------------------
// init: initialize members
//
void PMProfileModel::init()
{
    PMProfileModelBase::init();
    //
// //     trgfrqs_ = NULL;
    obsfrqs_ = NULL;
    freqweights_ = NULL;
    information_ = NULL;
    expMIDs_ = NULL;
    //
    bppprob_ = NULL;
    ppprobs_ = NULL;
    pppndxs_ = NULL;
    noppps_ = NULL;
    //
    ctpsset_ = false;
    ctbppprob_ = NULL;
    ctppprobs_ = NULL;
    ctpppndxs_ = NULL;
    noctppps_ = NULL;
    //
    ctxvecset_ = false;
    lppfact_ = 0.0f;
    ctxvecnorm2_ = NULL;
    ctxveclpprb_ = NULL;
    ctxvecs_ = NULL;
    ctxvecsize_ = 0;
    //
    sssset_ = false;
    sssp3_ = false;
    sss_ = NULL;
    sssprob_ = NULL;
    //
    name_ = NULL;
    description_ = NULL;
    szname_ = 0;
    szdescription_ = 0;

    nosequences_ = 0;
    effnosequences_ = 0.0f;

    memset( backprobs_, 0, sizeof(float) * PVDIM );
    memset( postprobs_, 0, sizeof(float) * PVDIM );

    referenceLambda_ = -1.0f;
    referenceK_ = -1.0f;
    lambda_ = -1.0f;
    entropy_ = -1.0f;
    parameterK_ = -1.0f;
    expscore_ = -1.0f;

//     CalcLogPProbFactor();
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
//
void PMProfileModel::destroy()
{
    int n;
// //     if( trgfrqs_ )      { free( trgfrqs_ ); trgfrqs_ = NULL; }
    if( obsfrqs_ )      { free( obsfrqs_ ); obsfrqs_ = NULL; }
    if( freqweights_ )  { free( freqweights_ ); freqweights_ = NULL; }
    if( information_ )  { free( information_ ); information_ = NULL; }
    if( expMIDs_ )      { free( expMIDs_ ); expMIDs_ = NULL; }
    if( bppprob_ )      { free( bppprob_ ); bppprob_ = NULL; }
    if( noppps_ )       { free( noppps_ ); noppps_ = NULL; }
    if( ctbppprob_ )    { free( ctbppprob_ ); ctbppprob_ = NULL; }
    if( noctppps_ )     { free( noctppps_ ); noctppps_ = NULL; }
    if( ctxvecnorm2_ )  { free( ctxvecnorm2_ ); ctxvecnorm2_ = NULL; }
    if( ctxveclpprb_ )  { free( ctxveclpprb_ ); ctxveclpprb_ = NULL; }
    sssset_ = false;
    sssp3_ = false;
    if( sss_ )          { free( sss_ ); sss_ = NULL; }
    if( sssprob_ )      { free( sssprob_ ); sssprob_ = NULL; }
    szname_ = 0; 
    szdescription_ = 0; 
    if( name_ )         { free( name_ ); name_ = NULL; }
    if( description_ )  { free( description_ ); description_ = NULL; }
    //
    if( ppprobs_ ) {
        for( n = 0; n < GetSize(); n++ )
            if( ppprobs_[n])
                free( ppprobs_[n]);
        free( ppprobs_ );
        ppprobs_ = NULL;
    }
    if( pppndxs_ ) {
        for( n = 0; n < GetSize(); n++ )
            if( pppndxs_[n])
                free( pppndxs_[n]);
        free( pppndxs_ );
        pppndxs_ = NULL;
    }
    //
    ctpsset_ = false;
    if( ctppprobs_ ) {
        for( n = 0; n < GetSize(); n++ )
            if( ctppprobs_[n])
                free( ctppprobs_[n]);
        free( ctppprobs_ );
        ctppprobs_ = NULL;
    }
    if( ctpppndxs_ ) {
        for( n = 0; n < GetSize(); n++ )
            if( ctpppndxs_[n])
                free( ctpppndxs_[n]);
        free( ctpppndxs_ );
        ctpppndxs_ = NULL;
    }
    //
    ctxvecset_ = false;
    lppfact_ = 0.0f;
    ctxvecsize_ = 0;
    if( ctxvecs_ ) {
        for( n = 0; n < GetSize(); n++ )
            if( ctxvecs_[n])
                free( ctxvecs_[n]);
        free( ctxvecs_ );
        ctxvecs_ = NULL;
    }
    //
    nosequences_ = 0;
    effnosequences_ = 0.0f;

    memset( backprobs_, 0, sizeof(float) * PVDIM );
    memset( postprobs_, 0, sizeof(float) * PVDIM );

    referenceLambda_ = -1.0f;
    referenceK_ = -1.0f;
    lambda_ = -1.0f;
    entropy_ = -1.0f;
    parameterK_ = -1.0f;
    expscore_ = -1.0f;
    //
    PMProfileModelBase::destroy();
}

// -------------------------------------------------------------------------
// Clear: erase all information but leave the space allocated
//
void PMProfileModel::Clear()
{
    int n;
    if( allocated_ ) {
        for( n = 0; n < GetSize(); n++ ) {
            if( ppprobs_ && ppprobs_[n]) free( ppprobs_[n]);
            if( pppndxs_ && pppndxs_[n]) free( pppndxs_[n]);
            if( ctppprobs_ && ctppprobs_[n]) free( ctppprobs_[n]);
            if( ctpppndxs_ && ctpppndxs_[n]) free( ctpppndxs_[n]);
            if( ctxvecs_ && ctxvecs_[n]) free( ctxvecs_[n]);
        }
// //         memset( trgfrqs_, 0, sizeof(float) * NUMAA * allocated_ );
        memset( obsfrqs_, 0, sizeof(float) * PVDIM * allocated_ );
        memset( freqweights_, 0, sizeof(float) * allocated_ );
        memset( information_, 0, sizeof(float) * allocated_ );
        memset( expMIDs_, 0, sizeof(float) * PS_NSTATES * ( allocated_ + 1 ));
        memset( bppprob_, 0, sizeof(float) * allocated_ );
        memset( ppprobs_, 0, sizeof(float*) * allocated_ );
        memset( pppndxs_, 0, sizeof(int*) * allocated_ );
        memset( noppps_,  0, sizeof(size_t) * allocated_ );
        ctpsset_ = false;
        memset( ctbppprob_, 0, sizeof(float) * allocated_ );
        memset( ctppprobs_, 0, sizeof(float*) * allocated_ );
        memset( ctpppndxs_, 0, sizeof(int*) * allocated_ );
        memset( noctppps_,  0, sizeof(size_t) * allocated_ );
        ctxvecset_ = false;
        lppfact_ = 0.0f;
        ctxvecsize_ = 0;
        memset( ctxvecnorm2_,  0, sizeof(float) * allocated_ );
        memset( ctxveclpprb_,  0, sizeof(float) * allocated_ );
        memset( ctxvecs_,  0, sizeof(float*) * allocated_ );
        sssset_ = false;
        sssp3_ = false;
        memset( sss_, 0, sizeof(char) * allocated_ );
        memset( sssprob_, 0, sizeof(float) * SS_NSTATES * allocated_ );
    }

    if( name_ )         { free( name_ ); name_ = NULL; szname_ = 0; }
    if( description_ )  { free( description_ ); description_ = NULL; szdescription_ = 0; }

    nosequences_ = 0;
    effnosequences_ = 0.0f;

    //do not erase probabilities!
//     memset( backprobs_, 0, sizeof(float) * PVDIM );
//     memset( postprobs_, 0, sizeof(float) * PVDIM );

    referenceLambda_ = -1.0f;
    referenceK_ = -1.0f;

    lambda_ = -1.0f;
    entropy_ = -1.0f;
    parameterK_ = -1.0f;
    expscore_ = -1.0f;

    PMProfileModelBase::Clear();
}

// -------------------------------------------------------------------------
// reallocate: (re)allocate memory for class members
//
void PMProfileModel::reallocate( int size )
{
// //     myreallocblk<float,int,NUMAA>( &trgfrqs_, size, allocated_ );
    myreallocblk<float,int,PVDIM>( &obsfrqs_, size, allocated_ );
    myrealloc<float,int>( &freqweights_, size, allocated_ );
    myrealloc<float,int>( &information_, size, allocated_ );
    myreallocblk<float,int,PS_NSTATES>( &expMIDs_, size+1, allocated_+1 );
    myrealloc<float,int>( &bppprob_, size, allocated_ );
    myrealloc<float*,int>( &ppprobs_, size, allocated_ );
    myrealloc<int*,int>( &pppndxs_, size, allocated_ );
    myrealloc<size_t,int>( &noppps_, size, allocated_ );
    myrealloc<float,int>( &ctbppprob_, size, allocated_ );
    myrealloc<float*,int>( &ctppprobs_, size, allocated_ );
    myrealloc<int*,int>( &ctpppndxs_, size, allocated_ );
    myrealloc<size_t,int>( &noctppps_, size, allocated_ );
    myrealloc<float,int>( &ctxvecnorm2_, size, allocated_ );
    myrealloc<float,int>( &ctxveclpprb_, size, allocated_ );
    myrealloc<float*,int>( &ctxvecs_, size, allocated_ );
    myrealloc<char,int>( &sss_, size, allocated_ );
    myreallocblk<float,int,SS_NSTATES>( &sssprob_, size, allocated_ );
    if( allocated_ <= 0 )
        memset( expMIDs_, 0, sizeof(float) * PS_NSTATES );
    //
    PMProfileModelBase::reallocate( size );
}

// // // -------------------------------------------------------------------------
// // // reallocate: (re)allocate memory for class members
// // //
// // void PMProfileModel::reallocate( int size )
// // {
// //     if( size <= allocated_ )
// //         return;
// // 
// // // //     float   (*tmp_trgfrqs)[NUMAA] = NULL;
// //     float   (*tmp_obsfrqs)[PVDIM] = NULL;
// //     float*  tmp_freqweights = NULL;
// //     float*  tmp_information = NULL;
// //     float   (*tmp_expMIDs)[PS_NSTATES] = NULL;
// //     float*  tmp_bppprob = NULL;
// //     float** tmp_ppprobs = NULL;
// //     int**   tmp_pppndxs = NULL;
// //     size_t* tmp_noppps = NULL;
// //     float*  tmp_ctbppprob = NULL;
// //     float** tmp_ctppprobs = NULL;
// //     int**   tmp_ctpppndxs = NULL;
// //     size_t* tmp_noctppps = NULL;
// //     float*  tmp_ctxvecnorm2 = NULL;
// //     float*  tmp_ctxveclpprb = NULL;
// //     float** tmp_ctxvecs = NULL;
// //     char*   tmp_sss = NULL;
// //     float   (*tmp_sssprob)[SS_NSTATES] = NULL;
// // 
// //     if( allocated_ <= 0 ) {
// // // //         tmp_trgfrqs = ( float(*)[NUMAA] )malloc( sizeof(float) * NUMAA * size );
// //         tmp_obsfrqs = ( float(*)[PVDIM] )malloc( sizeof(float) * PVDIM * size );
// //         tmp_freqweights = ( float* )malloc( sizeof(float) * size );
// //         tmp_information = ( float* )malloc( sizeof(float) * size );
// //         tmp_expMIDs = ( float(*)[PS_NSTATES])malloc( sizeof(float)*PS_NSTATES*(size+1));
// //         tmp_bppprob = ( float* )malloc( sizeof(float) * size );
// //         tmp_ppprobs = ( float** )malloc( sizeof(float*) * size );
// //         tmp_pppndxs = ( int** )malloc( sizeof(int*) * size );
// //         tmp_noppps = ( size_t* )malloc( sizeof(size_t) * size );
// //         tmp_ctbppprob = ( float* )malloc( sizeof(float) * size );
// //         tmp_ctppprobs = ( float** )malloc( sizeof(float*) * size );
// //         tmp_ctpppndxs = ( int** )malloc( sizeof(int*) * size );
// //         tmp_noctppps = ( size_t* )malloc( sizeof(size_t) * size );
// //         tmp_ctxvecnorm2 = ( float* )malloc( sizeof(float) * size );
// //         tmp_ctxveclpprb = ( float* )malloc( sizeof(float) * size );
// //         tmp_ctxvecs = ( float** )malloc( sizeof(float*) * size );
// //         tmp_sss = ( char* )malloc( sizeof(char) * size );
// //         tmp_sssprob = ( float(*)[SS_NSTATES])malloc( sizeof(float)*SS_NSTATES*size );
// //     } else {
// // // //         tmp_trgfrq = ( float(*)[NUMAA] )realloc( trgfrqs_, sizeof(float) * NUMAA * size );
// //         tmp_obsfrqs = ( float(*)[PVDIM] )realloc( obsfrqs_, sizeof(float) * PVDIM * size );
// //         tmp_freqweights = ( float* )realloc( freqweights_, sizeof(float) * size );
// //         tmp_information = ( float* )realloc( information_, sizeof(float) * size );
// //         tmp_expMIDs = ( float(*)[PS_NSTATES])realloc( expMIDs_, sizeof(float)*PS_NSTATES*(size+1));
// //         tmp_bppprob = ( float* )realloc( bppprob_, sizeof(float) * size );
// //         tmp_ppprobs = ( float** )realloc( ppprobs_, sizeof(float*) * size );
// //         tmp_pppndxs = ( int** )realloc( pppndxs_, sizeof(int*) * size );
// //         tmp_noppps = ( size_t* )realloc( noppps_, sizeof(size_t) * size );
// //         tmp_ctbppprob = ( float* )realloc( ctbppprob_, sizeof(float) * size );
// //         tmp_ctppprobs = ( float** )realloc( ctppprobs_, sizeof(float*) * size );
// //         tmp_ctpppndxs = ( int** )realloc( ctpppndxs_, sizeof(int*) * size );
// //         tmp_noctppps = ( size_t* )realloc( noctppps_, sizeof(size_t) * size );
// //         tmp_ctxvecnorm2 = ( float* )realloc( ctxvecnorm2_, sizeof(float) * size );
// //         tmp_ctxveclpprb = ( float* )realloc( ctxveclpprb_, sizeof(float) * size );
// //         tmp_ctxvecs = ( float** )realloc( ctxvecs_, sizeof(float*) * size );
// //         tmp_sss = ( char* )realloc( sss_, sizeof(char) * size );
// //         tmp_sssprob = ( float(*)[SS_NSTATES])realloc( sssprob_, sizeof(float)*SS_NSTATES*size );
// //     }
// // 
// //     if( /*!tmp_trgfrqs || */!tmp_obsfrqs || !tmp_freqweights || !tmp_information || !tmp_expMIDs ||
// //         !tmp_bppprob || !tmp_ppprobs || !tmp_pppndxs || !tmp_noppps ||
// //         !tmp_ctbppprob || !tmp_ctppprobs || !tmp_ctpppndxs || !tmp_noctppps ||
// //         !tmp_ctxvecnorm2 || !tmp_ctxveclpprb || !tmp_ctxvecs ||
// //         !tmp_sss || !tmp_sssprob ) 
// //     {
// // // //         if( tmp_trgfrqs ) { free( tmp_trgfrqs ); tmp_trgfrqs = NULL; }
// //         if( tmp_obsfrqs ) { free( tmp_obsfrqs ); tmp_obsfrqs = NULL; }
// //         if( tmp_freqweights ) { free( tmp_freqweights ); tmp_freqweights = NULL; }
// //         if( tmp_information ) { free( tmp_information ); tmp_information = NULL; }
// //         if( tmp_expMIDs ) { free( tmp_expMIDs ); tmp_expMIDs = NULL; }
// //         if( tmp_bppprob ) { free( tmp_bppprob ); tmp_bppprob = NULL; }
// //         if( tmp_ppprobs ) { free( tmp_ppprobs ); tmp_ppprobs = NULL; }
// //         if( tmp_pppndxs ) { free( tmp_pppndxs ); tmp_pppndxs = NULL; }
// //         if( tmp_noppps ) { free( tmp_noppps ); tmp_noppps = NULL; }
// //         if( tmp_ctbppprob ) { free( tmp_ctbppprob ); tmp_ctbppprob = NULL; }
// //         if( tmp_ctppprobs ) { free( tmp_ctppprobs ); tmp_ctppprobs = NULL; }
// //         if( tmp_ctpppndxs ) { free( tmp_ctpppndxs ); tmp_ctpppndxs = NULL; }
// //         if( tmp_noctppps ) { free( tmp_noctppps ); tmp_noctppps = NULL; }
// //         if( tmp_ctxvecnorm2 ) { free( tmp_ctxvecnorm2 ); tmp_ctxvecnorm2 = NULL; }
// //         if( tmp_ctxveclpprb ) { free( tmp_ctxveclpprb ); tmp_ctxveclpprb = NULL; }
// //         if( tmp_ctxvecs ) { free( tmp_ctxvecs ); tmp_ctxvecs = NULL; }
// //         if( tmp_sss ) { free( tmp_sss ); tmp_sss = NULL; }
// //         if( tmp_sssprob ) { free( tmp_sssprob ); tmp_sssprob = NULL; }
// //         throw MYRUNTIME_ERROR("PMProfileModel::reallocate: Not enough memory.");
// //     }
// // 
// // // //     trgfrqs_ = tmp_trgfrqs;
// //     obsfrqs_ = tmp_obsfrqs;
// //     freqweights_ = tmp_freqweights;
// //     information_ = tmp_information;
// //     expMIDs_ = tmp_expMIDs; tmp_expMIDs++;
// //     bppprob_ = tmp_bppprob;
// //     ppprobs_ = tmp_ppprobs;
// //     pppndxs_ = tmp_pppndxs;
// //     noppps_ = tmp_noppps;
// //     ctbppprob_ = tmp_ctbppprob;
// //     ctppprobs_ = tmp_ctppprobs;
// //     ctpppndxs_ = tmp_ctpppndxs;
// //     noctppps_ = tmp_noctppps;
// //     ctxvecnorm2_ = tmp_ctxvecnorm2;
// //     ctxveclpprb_ = tmp_ctxveclpprb;
// //     ctxvecs_ = tmp_ctxvecs;
// //     sss_ = tmp_sss;
// //     sssprob_ = tmp_sssprob;
// // 
// //     if( 0 < allocated_ ) {
// // // //         tmp_trgfrqs += allocated_;
// //         tmp_obsfrqs += allocated_;
// //         tmp_freqweights += allocated_;
// //         tmp_information += allocated_;
// //         tmp_expMIDs += allocated_;
// //         tmp_bppprob += allocated_;
// //         tmp_ppprobs += allocated_;
// //         tmp_pppndxs += allocated_;
// //         tmp_noppps += allocated_;
// //         tmp_ctbppprob += allocated_;
// //         tmp_ctppprobs += allocated_;
// //         tmp_ctpppndxs += allocated_;
// //         tmp_noctppps += allocated_;
// //         tmp_ctxvecnorm2 += allocated_;
// //         tmp_ctxveclpprb += allocated_;
// //         tmp_ctxvecs += allocated_;
// //         tmp_sss += allocated_;
// //         tmp_sssprob += allocated_;
// //     } else {
// //         memset( expMIDs_, 0, sizeof(float) * PS_NSTATES );
// //     }
// // 
// //     int difference = size - allocated_;
// //     // fill uninitialized memory with zeros
// // // //     memset( tmp_trgfrqs, 0, sizeof(float) * NUMAA * difference );
// //     memset( tmp_obsfrqs, 0, sizeof(float) * PVDIM * difference );
// //     memset( tmp_freqweights, 0, sizeof(float) * difference );
// //     memset( tmp_information, 0, sizeof(float) * difference );
// //     memset( tmp_expMIDs, 0, sizeof(float) * PS_NSTATES * difference );
// //     memset( tmp_bppprob, 0, sizeof(float) * difference );
// //     memset( tmp_ppprobs, 0, sizeof(float*) * difference );
// //     memset( tmp_pppndxs, 0, sizeof(int*) * difference );
// //     memset( tmp_noppps, 0, sizeof(size_t) * difference );
// //     memset( tmp_ctbppprob, 0, sizeof(float) * difference );
// //     memset( tmp_ctppprobs, 0, sizeof(float*) * difference );
// //     memset( tmp_ctpppndxs, 0, sizeof(int*) * difference );
// //     memset( tmp_noctppps, 0, sizeof(size_t) * difference );
// //     memset( tmp_ctxvecnorm2, 0, sizeof(float) * difference );
// //     memset( tmp_ctxveclpprb, 0, sizeof(float) * difference );
// //     memset( tmp_ctxvecs, 0, sizeof(float*) * difference );
// //     memset( tmp_sss, 0, sizeof(char) * difference );
// //     memset( tmp_sssprob, 0, sizeof(float) * SS_NSTATES * difference );
// //     //
// //     PMProfileModelBase::reallocate( size );
// // }

// -------------------------------------------------------------------------
// Push: push information relevant to an additional position;
//
void PMProfileModel::Push( const float vals[PVDIM], const float frqs[PVDIM], 
                           char rs, float weight, float info, float expnobs[PS_NSTATES] )
{
    if( allocated_ <= GetSize()) {
        int newcap = TIMES2( allocated_ );
        if( newcap <= GetSize())
            newcap = GetSize() + 1;
        reallocate( newcap );
    }

    PushAt( vals, frqs, rs, weight, info, expnobs, GetSize());
}

// -------------------------------------------------------------------------
// PushAt: push information relevant to position pos
//
void PMProfileModel::PushAt( const float vals[PVDIM], const float frqs[PVDIM], 
    char rs, float weight, float info, float expnobs[PS_NSTATES], int pos )
{
    if( allocated_ <= pos )
        throw MYRUNTIME_ERROR("PMProfileModel::PushAt: Memory access error.");

    PMProfileModelBase::PushAt( vals, rs, pos );
    SetObsFreqsAt( pos, frqs );
    freqweights_[pos] = weight;
    information_[pos] = info;
    SetMIDExpNoObservationsAt( pos, expnobs );
}

// -------------------------------------------------------------------------
// PushAt: push partial information relevant to position pos
//
void PMProfileModel::PushAt( const float vals[PVDIM], const float frqs[PVDIM], char rs, int pos )
{
    if( allocated_ <= pos )
        throw MYRUNTIME_ERROR("PMProfileModel::PushAt: Memory access error.");
    PMProfileModelBase::PushAt( vals, rs, pos );
    SetObsFreqsAt( pos, frqs );
}

// // // -------------------------------------------------------------------------
// // // Finalize: finalize data; calculate target frequencies from scores
// // //
// // void PMProfileModel::Finalize()
// // {
// //     ConvertToTrgFreqs();
// // }





// -------------------------------------------------------------------------
// CalcLogPProbFactor: calculate log of prior probability factor used in 
// processing context vectors
//
void PMProfileModel::CalcLogPProbFactor()
{
    float lpfact = 0.0f;
    float term1, term2;
    float operr;
    int err;
    const float do2 = (float)CVS.DIM * 0.5f;
    const float nu0p1o2 = (CVS.NU0+1.0f) * 0.5f;
    const float kp0 = CVS.KAPPA0;
    const float kp0p1 = kp0+1.0f;

    const float nu_t_o2 = nu0p1o2;//nu over 2

    const float dofdom_o2 = nu_t_o2 - do2;
    const float dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = extspsl::psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw MYRUNTIME_ERROR( extspsl::TranslatePSLError( err ));
    if(( err = extspsl::psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw MYRUNTIME_ERROR( extspsl::TranslatePSLError( err ));

    lpfact = term1 - term2;
//     lpfact -= do2 * SLC_LNPI;//NOTE:do not include pi^{-1/2A}
    lpfact += do2 * (logf(kp0)-logf(kp0p1));

    lppfact_ = lpfact;
}





// -------------------------------------------------------------------------
// ConvertToTrgFreqs: calculate target frequencies from scores and 
// overwrite scores
//
void PMProfileModel::ConvertToTrgFreqs()
{
    int     effnoress = NUMAA;
    float (*trgfs )[PVDIM] = GetTrgFreqs();
    char    res;
    int     m, r;

    static const float lsclambda = STABLE.StatisParam( Ungapped, Lambda );

    if( trgfs == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::ConvertToTrgFreqs: Null target frequencies.");

    for( m = 0; m < GetSize(); m++, trgfs++ ) {
        res = GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;
        for( r = 0; r < effnoress; r++ )
            (*trgfs)[r] = expf((*trgfs)[r] * lsclambda ) * STABLE.PROBABility(r);
            //(*trgfs)[r] = expf( GetValueAt(m,r)*lsclambda ) * STABLE.PROBABility(r);
    }
}

// -------------------------------------------------------------------------
// ConvertToScores: calculate scores from target frequencies and overwrite 
// target frequencies 
//
void PMProfileModel::ConvertToScores()
{
    const mystring preamb = "PMProfileModel::ConvertToScores: ";
    //int     noress = NUMALPH;
    int     effnoress = NUMAA;
    float (*trgfs)[PVDIM] = GetTrgFreqs();
    float (*scos)[PVDIM] = trgfs;
    float   bprob/*, ratv*/;
    char    res;
    int     m, r;

    static const float lsclambda = STABLE.StatisParam( Ungapped, Lambda );

    if( scos == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");
    if( lsclambda <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid lambda value.");

    for( m = 0; m < GetSize(); m++, scos++ ) {
        res = GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;
        for( r = 0; r < effnoress; r++ ) {
            bprob = STABLE.PROBABility(r);
            if((*scos)[r] <= 0.0f || 1.0f <= (*scos)[r])
                throw MYRUNTIME_ERROR( preamb + "Invalid target frequencies.");
            if( bprob <= 0.0f || 1.0f <= bprob )
                throw MYRUNTIME_ERROR( preamb + "Invalid background probabilities.");
            (*scos)[r] = logf((*scos)[r] / bprob ) / lsclambda;
        }
//         for( ; r < noress; r++ ) {
//             ratv = STABLE.FreqRatio( res, r );
//             if( ratv <= 0.0f )
//                 (*scos)[r] = SCORE_MIN;
//             else
//                 (*scos)[r] = STABLE.PrecomputedEntry( res, r );
//         }
    }
}

// =========================================================================
// MixTrgFrequencies: mix target frequencies using the HDP mixture model
//
void PMProfileModel::MixTrgFrequencies( const HDPbase* hdpbase )
{
    mystring preamb = "PMProfileModel::MixTrgFrequencies: ";

    const int   noeffress = NUMAA;
    float   pme = 0.0f;//estimated probability
    float*  infvv = GetInformation();
    float  (*trgfs)[PVDIM] = GetTrgFreqs();

    if( infvv == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null information content.");
    if( trgfs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");
    if( hdpbase == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDPbase object.");

    int     ctxtsz = hdpbase->GetCtxtSize();
    int     prolen = GetSize();
    int     p;
    int     left, right, r;
    int     hlf = ctxtsz >> 1;
//     int     parity = ( ctxtsz & 1 ) ^ 1;
//     int     mid = hlf - parity;

    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz ) {
        warning(( preamb + "Profile length < context size.").c_str());
        return;
    }
    extspsl::Pslmatrix  promtx( prolen, noeffress );//profile matrix
    extspsl::Pslmatrix  ctxmtx;//context matrix
    extspsl::Pslvector  ctxsck;//context stack
    extspsl::Pslvector  ctxnrm( noeffress*ctxtsz );//normal transform
    extspsl::Pslvector  mixed( noeffress );//mixed vector
    float infrm;

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw MYRUNTIME_ERROR( preamb + "Invalid residue.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0f || 1.0f <= (*trgfs)[r])
                throw MYRUNTIME_ERROR( preamb + "Invalid target frequencies.");
            promtx.SetValueAt( p, r, (*trgfs)[r]);
        }
    }

    //iterate over all profile positions
    trgfs = GetTrgFreqs();
    for( p = 0; p < prolen; p++, trgfs++, infvv++ )
    {
        right = SLC_MIN( prolen-1, p+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //mix central vector of context;
        //input vector ctxnrm will be transformed;
        //mixed will be returned in logit-normal space
        hdpbase->MixCLNVar( ctxnrm, &mixed );

        //write PME back to target frequencies
        infrm = 0.0f;
        for( r = 0; r < noeffress; r++ ) {
            pme = mixed.GetValueAt(r);
            (*trgfs)[r] = pme;
            //calculate relative entropy
            if( STABLE.PROBABility(r) <= 0.0f )
                continue;
            if( 0.0f < pme )
                infrm += pme * logf( pme / STABLE.PROBABility(r));
        }
        //save relative entropy
        infrm /= SLC_LN2;
        infrm = SLC_MAX( 0.0f, infrm );
        *infvv = infrm;
    }
// //     //recalculate scores
// //     CalcScores();
}

// =========================================================================
// CalcTFPostPredictives: calculate posterior predictive probabilities of 
// target frequencies using the HDP mixture model
//
void PMProfileModel::CalcTFPostPredictives( const HDPbase* hdpbase )
{
    mystring preamb = "PMProfileModel::CalcTFPostPredictives: ";

    const int   noeffress = NUMAA;
    float   ppr = 0.0f;//posterior probability
    float  (*trgfs)[PVDIM] = GetTrgFreqs();

    if( trgfs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");
    if( hdpbase == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDPbase object.");

    int     nosupcls = hdpbase->GetNoSupClusters();
    int     ctxtsz = hdpbase->GetCtxtSize();
    int     prolen = GetSize();
    int     p;
    int     left, right, r;
    int     hlf = ctxtsz >> 1;
//     int     parity = ( ctxtsz & 1 ) ^ 1;
//     int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw MYRUNTIME_ERROR( preamb + "Wrong number of HDP support clusters.");
    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz ) {
        warning(( preamb + "Profile length < context size.").c_str());
        return;
    }
    extspsl::Pslmatrix  promtx( prolen, noeffress );//profile matrix
    extspsl::Pslmatrix  ctxmtx;//context matrix
    extspsl::Pslvector  ctxsck;//context stack
    extspsl::Pslvector  ctxnrm( noeffress*ctxtsz );//normal transform
    extspsl::Pslvector  ppprobs;//vector of posterior probabilities
    extspsl::Ivector    cindcs;//indices of clusters

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw MYRUNTIME_ERROR( preamb + "Invalid residue.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0f || 1.0f <= (*trgfs)[r])
                throw MYRUNTIME_ERROR( preamb + "Invalid target frequencies.");
            promtx.SetValueAt( p, r, (*trgfs)[r]);
        }
    }

    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = SLC_MIN( prolen-1, p+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //calculate posterior predictive probabilities;
        //input vector ctxnrm will be transformed;
        hdpbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs, 0.02f/*lpfact*/);

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw MYRUNTIME_ERROR( preamb + 
            "Inconsistent number of posterior predictive probabilities.");

        //save posteriors
        SetBckPPProbAt( ppr, p );
        SetPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}

// -------------------------------------------------------------------------
// Calc_ctTFPostPredictives: calculate posterior predictive probabilities of 
// target frequencies using the HDP mixture model; 
// (inferred distribution of context vectors)
//
void PMProfileModel::Calc_ctTFPostPredictives( const HDPbase* hdpctbase )
{
    mystring preamb = "PMProfileModel::Calc_ctTFPostPredictives: ";
    const int   noeffress = NUMAA;
    float   ppr = 0.0f;//posterior probability
    float   cfc;
    float  (*trgfs)[PVDIM] = GetTrgFreqs();
    const size_t cnCTXLEN = 21;//TODO: read ctx length from parameter file
    const float  cdCTCWGT = 0.5f;//TODO: read ctx central weight from parameter file
    CtxtCoefficients coeffs( cnCTXLEN, cdCTCWGT );
    coeffs.FindCoefficients();

    if( trgfs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");
    if( hdpctbase == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDPbase object.");

    int     nosupcls = hdpctbase->GetNoSupClusters();
    int     ctxtsz = hdpctbase->GetCtxtSize();
    int     prolen = GetSize();
    int     p, pn, c;
    int     left, right, r;
    int     hlf = cnCTXLEN >> 1;
    int     parity = ( cnCTXLEN & 1 ) ^ 1;
//     int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw MYRUNTIME_ERROR( preamb + "Wrong number of HDP support clusters.");
    if( ctxtsz < 1 || 1 < ctxtsz )
        throw MYRUNTIME_ERROR( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz || prolen < (int)(((cnCTXLEN-1)>>1) + 1) ) {
        warning(( preamb + "Profile length < context size.").c_str());
        return;
    }
    extspsl::Pslmatrix  promtx( prolen, noeffress-1 );//transformed profile matrix
    extspsl::Pslvector  tfrqs( noeffress );//target frequencies
    extspsl::Pslvector  ctxnrm( noeffress-1 );//positional representative
    extspsl::Pslvector  ppprobs;//vector of posterior probabilities
    extspsl::Ivector    cindcs;//indices of clusters

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw MYRUNTIME_ERROR( preamb + "Invalid residue.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0f || 1.0f <= (*trgfs)[r])
                throw MYRUNTIME_ERROR( preamb + "Invalid target frequencies.");
            tfrqs.SetValueAt( r, (*trgfs)[r]);
        }
        //transform to normal
        ::LogitNormal2Normal_f( tfrqs.GetVector(), noeffress, 1.e-1f, false );
        for( r = 0; r < noeffress-1; r++ ) {
            promtx.SetValueAt( p, r, tfrqs.GetValueAt(r));
        }
    }

    SetctPsSet( true );
    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = p + hlf;//SLC_MIN( prolen-1, p+hlf );
        left = p - hlf + parity;//SLC_MAX( 0, right-ctxtsz+1 );

        ctxnrm.Zero();
        for( r = left, c = 0; r <= right && c < (int)coeffs.GetLength(); r++, c++ ) {
            pn = r;
            if( prolen <= pn || pn < 0 )
                pn = p + p - pn;
            cfc = coeffs.GetCoefficientAt(c);
            ctxnrm.Superposition( cfc, promtx.RowVector(pn));
        }

        //calculate posterior predictive probabilities;
        hdpctbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs,
                                0.02f/*lpfact*/, true/*usepriors*/, false/*tonormal*/ );

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw MYRUNTIME_ERROR( preamb + 
            "Inconsistent number of posterior predictive probabilities.");

        //save posteriors
        SetctBckPPProbAt( ppr, p );
        SetctPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}





// =========================================================================
// CalcCtxVector: calculate context vector, its prior probability, and its 
// squared norm for each position
//
void PMProfileModel::CalcCtxVector()
{
    mystring preamb = "PMProfileModel::CalcCtxVector: ";
    const int   noeffress = NUMAA;
    float   cfc;
    float  (*trgfs)[PVDIM] = GetTrgFreqs();
    CtxtCoefficients coeffs( CVS.CTXLEN, CVS.CWGT );
    coeffs.FindCoefficients();

    if( CVS.CTXLEN < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid context length.");
    if( CVS.CWGT < 0.0f || 1.0f < CVS.CWGT )
        throw MYRUNTIME_ERROR( preamb + "Invalid context central weight.");

    CalcLogPProbFactor();

    if( trgfs == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");

    int     prolen = GetSize();
    int     p, pn, c, cc, d;
    int     left, right, r, err;
    int     parity = ( CVS.CTXLEN & 1 ) ^ 1;
    int     hlf = CVS.CTXLEN >> 1;
    int     mid = hlf - parity;

    if( prolen < CVS.CTXLEN || prolen < (int)(((CVS.CTXLEN-1)>>1) + 1) ) {
        warning(( preamb + "Profile length < context size.").c_str());
        return;
    }
    extspsl::Pslmatrix  promtx( prolen, noeffress-1 );//transformed profile matrix
    extspsl::Pslvector  tfrqs( noeffress );//target frequencies
    extspsl::Pslvector  ctxnrm( CVS.DIM );//positional representative
    float norm2;//squared norm of vector
    float pprob;//prior probability of vector

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw MYRUNTIME_ERROR( preamb + "Invalid residue.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0f || 1.0f <= (*trgfs)[r])
                throw MYRUNTIME_ERROR( preamb + "Invalid target frequencies.");
            tfrqs.SetValueAt( r, (*trgfs)[r]);
        }
        //transform to normal
        ::LogitNormal2Normal_f( tfrqs.GetVector(), noeffress, 1.e-1f, false );
        for( r = 0; r < noeffress-1; r++ ) {
            promtx.SetValueAt( p, r, tfrqs.GetValueAt(r));
        }
    }

    SetCtxVecSet( true );
    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = p + hlf;
        left = p - hlf + parity;

        ctxnrm.Zero(); cc = 0;
        for( r = left, c = 0; r <= right && c < (int)coeffs.GetLength(); r++, c++ ) {
            pn = r;
            if( prolen <= pn || pn < 0 )
                pn = p + p - pn;
            cfc = coeffs.GetCoefficientAt(c);
            if( CVS.MIX ) {
                if( CVS.AUG ) {
                    cc = 0;
                    if( mid == c ) cc = promtx.GetNoCols();
                    else if( mid < c ) cc = TIMES2( promtx.GetNoCols());
                    //
                    for( d = 0; d < promtx.GetNoCols(); d++ )
                        ctxnrm.AddValueAt( cc+d, cfc*(promtx.GetValueAt(pn,d)-CVS.MU0.GetValueAt(d)) );
                }
                else
                    if(( err = ctxnrm.Superposition( cfc, promtx.RowVector(pn))) != 0 )
                        throw MYRUNTIME_ERROR( preamb + extspsl::TranslatePSLError( err ));
            }
            else {
                for( d = 0; d < promtx.GetNoCols(); d++ )
                    ctxnrm.SetValueAt( cc++, cfc*(promtx.GetValueAt(pn,d)-CVS.MU0.GetValueAt(d)) );
            }
        }

        if( CVS.MIX && !CVS.AUG && CVS.MEAN ) {
            if(( err = ctxnrm.Superposition( -1.0f, CVS.MU0 )) != 0 )
                throw MYRUNTIME_ERROR( preamb + extspsl::TranslatePSLError( err ));
        }

        //vector's squared norm
        norm2 = ctxnrm.Norm2();
        norm2 *= norm2;

        //vector's prior probability scaled by pi^{1/2dim}, where dim=noeffress-1
        pprob = GetLogPProbFactor();
        pprob -= 0.5f*(CVS.NU0+1.0f) * logf(1.0f+CVS.KAPPA0/(CVS.KAPPA0+1.0f)*norm2);

        //save vector and calculated magnitudes
        SetCtxVecPlusAt( p, norm2, pprob, ctxnrm.GetVector(), ctxnrm.GetSize());
    }
}





// =========================================================================
// SetName: set profile name
//
void PMProfileModel::SetName( const char* newname )
{
    if( name_ ) {
        free( name_ );
        name_ = NULL;
        szname_ = 0;
    }

    if( !newname || !*newname )
        return;

    size_t newlength = strlen( newname );

    if( !newlength )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetName: Invalid argument." );

    name_ = ( char* )malloc( newlength+1 );
    if( !name_ )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetName: Not enough memory." );

    strncpy( name_, newname, newlength+1 ); // include the terminating null symbol
    szname_ = newlength;
}

// -------------------------------------------------------------------------
// SetDescription: set profile description
//
void PMProfileModel::SetDescription( const char* newdesc )
{
    if( description_ ) {
        free( description_ );
        description_ = NULL;
        szdescription_ = 0;
    }

    if( !newdesc )
        return;

    size_t newlength = strlen( newdesc );

    if( !newlength )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetDescription: Invalid argument." );

    description_ = ( char* )malloc( newlength+1 );
    if( !description_ )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetDescription: Not enough memory." );

    strncpy( description_, newdesc, newlength+1 ); // include the terminating null symbol
    szdescription_ = newlength;
}

// -------------------------------------------------------------------------
// Print: print information to file
//
void PMProfileModel::Print( const char* filename ) const
{
    int nores = PVDIM;//effective number of residues
    //
    myruntime_error mre;
    FILE* fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "PMProfileModel::Print: Failed to open file for writing." );

    try {
        fprintf( fp, "Profile (no. sequences %zu, eff. no. sequences %.3f)", 
                 GetNoSequences(), GetEffNoSequences());

        if( GetName())
            fprintf( fp, ": %s", GetName());

        fprintf( fp, "%s", NL );

        if( GetDescription())
            fprintf( fp, "Description: %s%s%s", GetDescription(), NL, NL );

        fprintf( fp, "%28c         Target frequencies       %22c Freq weights %c "
                     "Information %c Exp.Obs.%s", 32, 32, 32, 32, NL );
        fprintf( fp, "%9c", 32 );

        for( char r = 0; r < nores; r++ )
            fprintf( fp, "%3c", DehashCode( r ) );

        for( int p = 0; p < GetSize(); p++ ) {
            fprintf( fp, "%s%5d %c   ", NL, p+1, DehashCode( GetResidueAt(p)));

            for( char r = 0; r < nores; r++ )
                fprintf( fp, "%3d", (int)rintf( GetValueAt(p,r)));

            fprintf( fp, " %11d", (int)rintf( GetFrequencyWeightAt( p )));
            fprintf( fp, " %13.2f", GetInformationAt( p ));
            fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_M ));
            fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_I ));
            fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_D ));
        }
        fprintf( fp, "%s", NL );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( fp != stdout )
        fclose( fp );

    if( mre.isset())
        throw mre;
}

// // // -------------------------------------------------------------------------
// // // OutputProfile: output profile information in the text format
// // //
// // void OutputProfile( const char* filename,
// //         const FrequencyMatrix& freq,
// //         const LogOddsMatrix& odds,
// //         const GapScheme& gaps )
// // {
// //     FILE*           fp = stdout;
// //     size_t          l = 0;
// //     const size_t    res_count = NUMALPH - 1; //exclude gap symbol
// // 
// //     if( filename )
// //         fp = fopen( filename, "w" );
// // 
// //     if( fp == NULL )
// //         throw myruntime_error(
// //             mystring( "Failed to open file for writing." ));
// // 
// //     fprintf( fp, "Multiple alignment (no. sequences %d, eff. no. sequences %.3f)",
// //             odds.GetNoSequences(), odds.GetEffNoSequences());
// // 
// //     if( odds.GetName())
// //         fprintf( fp, ": %s", odds.GetName());
// // 
// //     fprintf( fp, "\n" );
// // 
// //     if( odds.GetDescription())
// //         fprintf( fp, "First sequence description: %s\n\n", odds.GetDescription());
// // 
// // 
// //     fprintf( fp,"%28c Position-specific scoring matrix "
// //                 "%53c Weighted observed frequencies %30c Gap weights %c Freq weights %c Information\n",
// //                 32, 32, 32, 32, 32 );
// // 
// //     fprintf( fp, "%9c", 32 );
// // 
// //     for( unsigned char r = 0; r < res_count; r++ )
// //         fprintf( fp, "%3c", DehashCode( r ) );
// // 
// //     for( unsigned char r = 0; r < res_count; r++ )
// //         fprintf( fp, "%4c", DehashCode( r ) );
// // 
// //     for( int p = 0; p < odds.GetColumns(); p++ ) {
// //         // omit unused positions and gaps in query
// //         if( odds.GetResidueAt( p ) == GAP )
// //             continue;
// // 
// //         fprintf( fp, "\n%5d %c   ", ++l, DehashCode( odds.GetResidueAt( p )));
// // 
// //         for( unsigned char r = 0; r < res_count; r++ )
// //             fprintf( fp, "%2d ", ( int )rint( odds.GetValueAt( p, r )));
// // 
// //         for( unsigned char r = 0; r < res_count; r++ )
// //             fprintf( fp, "%4d", ( int )rint( 100 * freq.GetValueAt( p, r )));
// // 
// //         fprintf( fp, " %6d", ( int )rint( 100 * gaps.GetWeightsAt( p )));
// //         fprintf( fp, " %14d", ( int )rint( odds.GetFrequencyWeightAt( p )));
// // 
// //         fprintf( fp, " %13.2f", odds.GetInformationAt( p ));
// //     }
// // 
// //     fprintf( fp, "\n\n" );
// // 
// // // #ifdef SCALE_PROFILES
// //     odds.PrintParameters( fp );
// // // #endif
// // 
// //     if( fp != stdout )
// //         fclose( fp );
// // }





// =========================================================================
// GetMaxAnnotationWidth: get maximum width for annotation
//
size_t PMProfileModel::GetMaxAnnotationWidth() const
{
    static size_t width = 2 * OUTPUTINDENT + OUTPUTWIDTH + 2;
    return width;
}

// -------------------------------------------------------------------------
// PrintAnnotation: print annotation of name and description to string 
// stream; 
// NOTE: space of stream must be pre-allocated!
//
void PMProfileModel::PrintAnnotation( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//to ensure printing starting at the beginning of the stream
    PrintAnnotation( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: print annotation to file
//
void PMProfileModel::PrintAnnotation( FILE* fp ) const
{
    PrintAnnotation( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: format and print annotation of name and description
//
void PMProfileModel::PrintAnnotation( TPrintFunction print_func, void* vpn ) const
{
    static size_t preamble = 0;
    static size_t textwidth = OUTPUTINDENT + OUTPUTWIDTH + 2;
    static size_t width = textwidth + preamble;
    static size_t max_length = width;
    static size_t max_rows = 1;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// GetMinimumRequiredSizeForDescription: get minimum size required to
// contain the whole description of the profile
//
size_t PMProfileModel::GetMinimumRequiredSizeForDescription() const
{
    static size_t indent = OUTPUTINDENT;//indentation
    static size_t textwidth = OUTPUTWIDTH;//output width

#if 1//def __DEBUG__
    if( !textwidth )
        return 0;
#endif

    size_t name_size = GetNameSize();//actual name size 
    size_t desc_size = GetDescriptionSize();//actual description size 
    size_t max_desc_size = GetPrivateBufferSize();//maximum size for description

    size_t title_size  = ( name_size + desc_size ) * ( 1 + ( indent + 2 ) / textwidth ) + 1;

    if( title_size >= max_desc_size )
        title_size  = max_desc_size;

    return title_size;
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description to file
//
void PMProfileModel::PrintDescriptionFirst( FILE* fp ) const
{
    PrintDescriptionFirst( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description
//
void PMProfileModel::PrintDescriptionFirst( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = 0;
    static size_t   textwidth = OUTPUTINDENT + OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to string stream;
// NOTE: memory for stream must be pre-allocated!
//
void PMProfileModel::PrintDescription( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//to print at the beginning of the stream
    PrintDescription( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to file
//
void PMProfileModel::PrintDescription( FILE* fp ) const
{
    PrintDescription( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description
//
void PMProfileModel::PrintDescription( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = OUTPUTINDENT;
    static size_t   textwidth = OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, false );
}



// -------------------------------------------------------------------------
// PrintDescriptionHelper: helper method for formating and printing the 
// description
//
void PMProfileModel::PrintDescriptionHelper(
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width, size_t max_rows, size_t max_length,
        bool annotation ) const
{
    if( vpn == NULL )
        return;

    bool            printname = false;

    char*           buffer = GetPrivateBuffer();
    const char*     name = printname? GetName(): NULL;
    const char*     description = GetDescription();
    const char*     ending = "...";
    const size_t    sznl = strlen(NL);

    size_t  sz_name = printname? GetNameSize(): 0;
    size_t  sz_description = GetDescriptionSize();
    size_t  szending = strlen( ending );
    size_t  size = 0;
    size_t  pos = 0;    //current position in line
    size_t  tot = 0;    //position in the buffer
    bool    term = false;

    size = sz_name + sz_description + ( name ? 3: 0 );//+parentheses and space
    size += ( size / textwidth + 1 ) * ( preamble*sznl + 1 );//#lines times preamble

    if( buffer == NULL || !max_length || GetPrivateBufferSize() < max_length )
        return;

    if( size >= max_length - max_rows * ( preamble*sznl + 1 )) {
        size = max_length - max_rows * ( preamble*sznl + 1 ) - szending - 1;
        term = true;
    }

    if( !annotation ) {
        *buffer++ = '>'; pos++; tot++;
    }
    if( name ) {
        *buffer++ = '('; pos++; tot++;
        FormatBuffer( buffer, name, tot, pos, size, preamble, width );
        if( tot < size ) { *buffer++ = ')'; pos++; tot++; }
        if( tot < size ) { *buffer++ = 32;  pos++; tot++; }
    }
    if( description ) {
        FormatBuffer( buffer, description, tot, pos, size, preamble, width );
    }

    if( term ) {
        if( width <= pos + szending ) {
// //             *buffer++ = '\n';
            sprintf( buffer, "%s", NL );
            buffer += sznl;
            for( size_t k = 0; k < preamble; k++ ) *buffer++ = 32;
        }
        for( const char* p = ending; *p; *buffer++ = *p++ );
    }

    if( !annotation )
// //         if( *( buffer - 1 ) != '\n' )
// //             *buffer++ = '\n';
        if( strncmp( buffer-sznl, NL, sznl ) != 0 ) {
            sprintf( buffer, "%s", NL );
            buffer += sznl;
        }

    *buffer++ = 0;

    print_func( vpn, "%s", GetPrivateBuffer());
}

// -------------------------------------------------------------------------
// FormatBuffer: auxiliary method to format character buffer; note that
// space required for the first argument must be pre-allocated
//
void PMProfileModel::FormatBuffer( char*& format, const char* desc,
        size_t& tot, size_t& pos,
        const size_t size,
        const size_t indent, const size_t width ) const
{
    size_t  loc = 0;
    size_t  lto = 0;
    const size_t sznl = strlen(NL);
    const char* beg = ( tot <= 2 )? NULL: desc;//if the very beginning, don't check for line feed

    while( *desc && tot < size ) {
        if( *desc == 32 || *desc == 9 || desc == beg ) {
            loc = pos + 1;
            lto = tot + 1;
            for( const char* p = desc + 1;
                    *p && *p != 32 && *p != 9 && loc < width && lto < size;
                     p++, loc++, lto++ );
            if( width <= loc && indent < pos ) {
// //                 *format++ = '\n';
                sprintf( format, "%s", NL );
                format += sznl;
                for( size_t k = 0; k < indent; k++ ) *format++ = 32;
                tot += indent;
                pos  = indent;
                if( *desc == 32 || *desc == 9 )
                    desc++;
                if( beg )
                    beg = NULL;
                continue;
            }
        }
        *format++ = *desc++;
        pos++; tot++;
    }
}





// =========================================================================

extern "C" {
/*static */const char* patstrDATBINVER = "COMER bin profile v";
/*static */const char* patstrDATVER = "COMER profile v";
/*static */const char* patstrDESC = "DESC:";
/*static */const char* patstrFILE = "FILE:";
/*static */const char* patstrCMD = "CMD:";
/*static */const char* patstrLEN = "LEN:";
/*static */const char* patstrNOS = "NOS:";
/*static */const char* patstrEFFNOS = "EFFNOS:";
/*static */const char* patstrSCALE = "SCALE:";
/*static */const char* patstrNULL = "NULL";
/*static */const char* patstrPOST = "POST";
/*static */const char* patstrCT = "CT:";
/*static */const char* patstrCV = "CV:";
/*static */const char* patstrSS = "SS:";
/*static */const char* patstrSPCOMP = "Computed  ungapped,";
/*static */const char* patstrSPREFR = "Reference ungapped,";
/*static */const char* patstrENTROPY = "Entropy,";
/*static */const char* patstrEXPECTED = "Expected,";
/*static */const char* patstrEXPNN = "Expected score per pos. is non-negative,";
/*static */const char* patstrEND = "*";
/*static */extern const size_t lenstrEND = strlen( patstrEND );
/*static */const int   INTSCALE = 1000;
}

// =========================================================================
// CalculateAndWriteAddressTable_v2_2: calculate and write address 
// table to a binary file of profile database;
// fp, output file stream;
// szpreamble, size reserved for preamble in the beginning;
// addrfields, file addresses of the beginnings of profile fields, to be 
// calculated and written to file;
// addrfname_end_addrs, address of end addresses of filenames to be 
// calculated;
// addrdesc_end_addrs, address of end addresses of descripptions to be 
// calculated;
// addrfname, address of profile filenames to be calculated;
// addrdesc, address of profile descriptions to be calculated;
// nprofiles, number of profiles to comprise a database;
// totalposs, total number of positions over all profiles in a database;
//
void CalculateAndWriteAddressTable_v2_2(
    std::ofstream& fp,
    const size_t szpreamble,
    size_t* addrfields,
    size_t& addrdesc_end_addrs, size_t& addrdesc,
    const size_t nprofiles, const size_t totalposs )
{
    if( addrfields == NULL )
        throw MYRUNTIME_ERROR( "CalculateAndWriteAddressTable_v2_2: Null file addresses." );

    if( MAXFILESIZE <= szpreamble )
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid reserved size.");

    const int noress = NUMAA;//PMProfileModelBase::PVDIM;
    const int vsize = NUMAA-1;//context vector length
    size_t table[BUF_MAX];//NOTE:BUF_MAX assumed to be greater than v2_2_NSECTIONS!
    size_t addr;
    int r, t, tt, n = 0;

    addr = szpreamble + PMProfileModel::v2_2_NSECTIONS * sizeof(table[0]);//start of the profile data

    table[n++] = addrfields[pps2DLen] = addr;//lengths
    addr += nprofiles * sizeof(int);//end of lengths
    table[n++] = addrfields[pps2DENO] = addr;//ENOs
    addr += nprofiles * sizeof(float);//end of ENOs
    for( r = 0; r < noress; r++ ) {
        table[n++] = addrfields[pps2DBkgPrbs+r] = addr;//bkg probabilities for each r
        addr += nprofiles * sizeof(float);//end
    }
    for( t = 0, tt = -1; t < P_NSTATES; t++ ) {
        if( t == P_ID || t == P_DI )
            continue;
        tt++;
        table[n++] = addrfields[ptr2DTrnPrbs+tt] = addr;//trn probabilities for each tt
        addr += (nprofiles + totalposs) * sizeof(float);//end
    }
    for( r = 0; r < noress; r++ ) {
        table[n++] = addrfields[pmv2DTrgFrqs+r] = addr;//target probabilities for each r
        addr += totalposs * sizeof(float);//end
    }
    for( t = 0; t < vsize; t++ ) {
        table[n++] = addrfields[pmv2DCVentrs+t] = addr;//context vector entries (for each t)
        addr += totalposs * sizeof(float);//end
    }
    table[n++] = addrfields[pmv2DCVprior] = addr;//context vector probabilities
    addr += totalposs * sizeof(float);//end
    table[n++] = addrfields[pmv2DCVnorm2] = addr;//squared norms of context vector
    addr += totalposs * sizeof(float);//end
    for( t = 0; t < SS_NSTATES; t++ ) {
        table[n++] = addrfields[pmv2DSSsprbs+t] = addr;//secondary structure state probabilities for each t
        addr += totalposs * sizeof(float);//end
    }
    table[n++] = addrfields[pmv2DHDP1prb] = addr;//HDP1 cluster membership probabilities
    addr += totalposs * sizeof(float);//end
    table[n++] = addrfields[pmv2DHDP1ind] = addr;//HDP1 cluster indices
    addr += totalposs * sizeof(int);//end
    table[n++] = addrfields[pmv2Daa] = addr;//amino acid sequences
    addr += totalposs * sizeof(char);//end
    table[n++] = addrfields[pmv2DSSstate] = addr;//secondary structure state sequences
    addr += totalposs * sizeof(char);//end
    //
//     addrfname_end_addrs = addr;//end addresses of filenames are not included!
//     addr += nprofiles * sizeof(size_t);//end
    table[n++] = addrdesc_end_addrs = addr;//end addresses of profile descriptions
    addr += nprofiles * sizeof(size_t);//end
    //
    table[n++] = addrdesc = addr;//profile descriptions; the end is indefinite!

    //write the address table
    if( n != PMProfileModel::v2_2_NSECTIONS )
        throw MYRUNTIME_ERROR("CalculateAndWriteAddressTable_v2_2: Invalid number of address table entries.");

    fp.seekp(std::streamoff(szpreamble));
    if(fp.fail())
        throw MYRUNTIME_ERROR("CalculateAndWriteAddressTable_v2_2: Failed to set file position.");
    fp.write(reinterpret_cast<const char*>(table), n * sizeof(table[0]));
    if(fp.bad())
        throw MYRUNTIME_ERROR("CalculateAndWriteAddressTable_v2_2: Write to file failed.");
}

// =========================================================================
// BinaryWriteProfile_v2_2: sparse write of profile data to binary file;
// fp, output file stream;
// addrfields, file addresses of profile fields, being updated on each write call;
// addrfname_end_addrs, file address of end addresses of filenames;
// addrdesc_end_addrs, file address of end addresses of descripptions;
// addrfname, file address of profile filename, being updated;
// addrdesc, file address of profile description, being updated;
// pssm, profile model;
// gaps, model of transition probabilities;
//
void BinaryWriteProfile_v2_2(
    std::ofstream& fp,
    size_t* addrfields,
    size_t& addrdesc_end_addrs, size_t& addrdesc,
    const PMProfileModel& pssm, const PMTransModel& gaps )
{
    if( pssm.GetSize() != gaps.GetSize())
        throw MYRUNTIME_ERROR( "BinaryWriteProfile_v2_2: Inconsistent Profile data." );

    if( pssm.GetSize() < 1 )
        throw MYRUNTIME_ERROR( "BinaryWriteProfile_v2_2: Invalid Profile length." );

    if( pssm.GetCtxVecSet() && pssm.GetCtxVecSize() != pmv2DNoCVEls )
        throw MYRUNTIME_ERROR( "BinaryWriteProfile_v2_2: Invalid context vector length." );

    if( addrfields == NULL )
        throw MYRUNTIME_ERROR( "BinaryWriteProfile_v2_2: Null file addresses." );

    const int noress = NUMAA;//PMProfileModelBase::PVDIM;

    mystring        sequence;
    extspsl::Pslvector spvals;
    extspsl::Ivector intvals;
    extspsl::Pslvector uninfctxvec(pmv2DNoCVEls);//uninformative context vector
    const int Xuninf = MOptions::GetX_UNINF();
    constexpr bool mildmasking = true;//use mild masking
    const float*    ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    int             noppps;//number of p.p.probability values
    //float           norm2;//context vector's squared norm
    float           lpprb;//log of context vector's prior probability
    int             vsize = pmv2DNoCVEls;//pssm.GetCtxVecSize();//size of vector
    const char ch0 = 0;
    char strbuf[BUF_MAX];
    size_t bytes;
    char res;
    int m, r, t, tt;

    spvals.Allocate(pssm.GetSize()+1);
    intvals.Allocate(pssm.GetSize());

    //init each entry to ensure CVS2S score==0 when CVS information is not to be in use
    uninfctxvec.AssignTo(0.17302947f);

    //length:
    if( MAXFILESIZE <= addrfields[pps2DLen])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid profile length address.");
    fp.seekp(std::streamoff(addrfields[pps2DLen]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    bytes = sizeof(pssm.length_);
    fp.write(reinterpret_cast<const char*>(&pssm.length_),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pps2DLen] += bytes;

    //ENO:
    if( MAXFILESIZE <= addrfields[pps2DENO])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid profile ENO address.");
    fp.seekp(std::streamoff(addrfields[pps2DENO]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    bytes = sizeof(pssm.effnosequences_);
    fp.write(reinterpret_cast<const char*>(&pssm.effnosequences_),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pps2DENO] += bytes;


    //background probabilities:
    for( r = 0; r < noress; r++ ) {
        if( MAXFILESIZE <= addrfields[pps2DBkgPrbs+r])
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of background probabilities.");
        fp.seekp(std::streamoff(addrfields[pps2DBkgPrbs+r]));
        if(fp.fail())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
        bytes = sizeof(pssm.backprobs_[0]);
        fp.write(reinterpret_cast<const char*>(pssm.backprobs_+r),bytes);
        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
        addrfields[pps2DBkgPrbs+r] += bytes;
    }

    //NOTE: posterior probabilities are not written
    //for( r = 0; r < noress; r++ )
    //    fp.write(reinterpret_cast<const char*>(pssm.postprobs_+r),sizeof(pssm.postprobs_[0]));


    //transition probabilities:
    for( t = 0, tt = -1; t < P_NSTATES; t++ ) {
        if( t == P_ID || t == P_DI )
            continue;
        tt++;
        spvals.Clear();
        bytes = 0;
        if( MAXFILESIZE <= addrfields[ptr2DTrnPrbs+tt])
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of transition probabilities.");
        fp.seekp(std::streamoff(addrfields[ptr2DTrnPrbs+tt]));
        if(fp.fail())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
        //beginning transition probabilities
        lpprb = gaps.GetOrgTrProbsAt(t,-1);
        //NOTE: calculate and write logarithms of transitions
        lpprb = lpprb? (P_MD<t? logf(lpprb)*MAIN_TRNPRB_LOGFCT: logf(lpprb)): -32768.0f;
        spvals.Push(lpprb);
        bytes += sizeof(lpprb);
        //positional transition probabilities
        for( m = 0; m < pssm.GetSize(); m++ ) {
            res = pssm.GetResidueAt(m);
            if( res == GAP )
                continue;//omit unused positions and gaps
            lpprb = gaps.GetOrgTrProbsAt(t,m);
            //NOTE: calculate and write logarithms of transitions
            lpprb = lpprb? (P_MD<t? logf(lpprb)*MAIN_TRNPRB_LOGFCT: logf(lpprb)): -32768.0f;
            spvals.Push(lpprb);
            bytes += sizeof(lpprb);
        }
        fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
        addrfields[ptr2DTrnPrbs+tt] += bytes;
    }


    //target probabilities:
    for( r = 0; r < noress; r++ ) {
        spvals.Clear();
        bytes = 0;
        if( MAXFILESIZE <= addrfields[pmv2DTrgFrqs+r])
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of target probabilities.");
        fp.seekp(std::streamoff(addrfields[pmv2DTrgFrqs+r]));
        if(fp.fail())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
        for( m = 0; m < pssm.GetSize(); m++ ) {
            res = pssm.GetResidueAt(m);
            if( res == GAP )
                continue;//omit unused positions and gaps
            lpprb = (Xuninf && res == X && !mildmasking)? pssm.backprobs_[r]: pssm.values_[m][r];
            spvals.Push(lpprb);
            bytes += sizeof(lpprb);
        }
        fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
        addrfields[pmv2DTrgFrqs+r] += bytes;
    }


    //{{context vectors:
    for( t = 0; t < vsize; t++ ) {
        spvals.Clear();
        bytes = 0;
        if( MAXFILESIZE <= addrfields[pmv2DCVentrs+t])
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of context vectors.");
        fp.seekp(std::streamoff(addrfields[pmv2DCVentrs+t]));
        if(fp.fail())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
        for( m = 0; m < pssm.GetSize(); m++ ) {
            res = pssm.GetResidueAt(m);
            if( res == GAP )
                continue;//omit unused positions and gaps
            if( pssm.GetCtxVecSet() && !(Xuninf && res == X))
                lpprb = pssm.GetCtxVecAt(m)? 
                    pssm.GetCtxVecAt(m)[t]: uninfctxvec.GetValueAt(t);//HUGE_VALF;
            else
                lpprb = uninfctxvec.GetValueAt(t);//HUGE_VALF;
            spvals.Push(lpprb);
            bytes += sizeof(lpprb);
        }
        fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
        addrfields[pmv2DCVentrs+t] += bytes;
    }

    //probabilities of context vectors
    spvals.Clear();
    bytes = 0;
    if( MAXFILESIZE <= addrfields[pmv2DCVprior])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of probabilities of context vectors.");
    fp.seekp(std::streamoff(addrfields[pmv2DCVprior]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt(m);
        if( res == GAP )
            continue;//omit unused positions and gaps
        if( pssm.GetCtxVecSet() && !(Xuninf && res == X))
            lpprb = pssm.GetCtxVecLpprobAt(m);
        else
            lpprb = 0.0f;//HUGE_VALF;
        spvals.Push(lpprb);
        bytes += sizeof(lpprb);
    }
    fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2DCVprior] += bytes;

    //squared norms of context vectors
    spvals.Clear();
    bytes = 0;
    if( MAXFILESIZE <= addrfields[pmv2DCVnorm2])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of norms of context vectors.");
    fp.seekp(std::streamoff(addrfields[pmv2DCVnorm2]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt(m);
        if( res == GAP )
            continue;//omit unused positions and gaps
        if( pssm.GetCtxVecSet() && !(Xuninf && res == X))
            lpprb = pssm.GetCtxVecNorm2At(m);
        else
            lpprb = 0.0f;//HUGE_VALF;
        spvals.Push(lpprb);
        bytes += sizeof(lpprb);
    }
    fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2DCVnorm2] += bytes;
    //}}//CV


    //{{secondary structure state probabilities:
    for( t = 0; t < SS_NSTATES; t++ ) {
        spvals.Clear();
        bytes = 0;
        if( MAXFILESIZE <= addrfields[pmv2DSSsprbs+t])
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of sec. strcture probabilities.");
        fp.seekp(std::streamoff(addrfields[pmv2DSSsprbs+t]));
        if(fp.fail())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
        for( m = 0; m < pssm.GetSize(); m++ ) {
            res = pssm.GetResidueAt(m);
            if( res == GAP )
                continue;//omit unused positions and gaps
            if( pssm.GetSSSSet() && !(Xuninf && res == X))
                lpprb = pssm.sssprob_[m]? pssm.sssprob_[m][t]: 0.0f;//HUGE_VALF;
            else
                lpprb = 0.0f;//so that SSS inf vector is uninformative//(HUGE_VALF;)
            spvals.Push(lpprb);
            bytes += sizeof(lpprb);
        }
        fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
        addrfields[pmv2DSSsprbs+t] += bytes;
    }
    //}}//SS


    //{{HDP1
    //cluster membership probabilities:
    spvals.Clear();
    bytes = 0;
    if( MAXFILESIZE <= addrfields[pmv2DHDP1prb])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of HDP1 cluster membership probabilities.");
    fp.seekp(std::streamoff(addrfields[pmv2DHDP1prb]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt(m);
        if( res == GAP )
            continue;//omit unused positions and gaps
        noppps = (int)pssm.GetNoPPProbsAt(m);
        if( noppps < 1 ) {
            sprintf(strbuf, "BinaryWriteProfile_v2_2: Invalid # HDP1 support clusters at %d.", m);
            throw MYRUNTIME_ERROR(strbuf);
        }
        ppprobs = pssm.GetPPProbsAt(m);
        if( ppprobs == NULL ) {
            sprintf(strbuf, "BinaryWriteProfile_v2_2: Null HDP1 cluster membership probabilities at %d.", m);
            throw MYRUNTIME_ERROR(strbuf);
        }
        if(!(Xuninf && res == X && !mildmasking))
            spvals.Push(ppprobs[0]);
        else
            spvals.Push(0.0f);//make uninformative
        bytes += sizeof(ppprobs[0]);
    }
    fp.write(reinterpret_cast<const char*>(spvals.GetVector()),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2DHDP1prb] += bytes;

    //cluster indices:
    intvals.Clear();
    bytes = 0;
    if( MAXFILESIZE <= addrfields[pmv2DHDP1ind])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of HDP1 cluster indices.");
    fp.seekp(std::streamoff(addrfields[pmv2DHDP1ind]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt(m);
        if( res == GAP )
            continue;//omit unused positions and gaps
        pppndxs = pssm.GetPPPIndsAt(m);
        if( pppndxs == NULL ) {
            sprintf(strbuf, "BinaryWriteProfile_v2_2: Null HDP1 cluster indices at %d.", m);
            throw MYRUNTIME_ERROR(strbuf);
        }
        if(!(Xuninf && res == X && !mildmasking))
            intvals.Push(pppndxs[0]);
        else
            intvals.Push(-1);//make uninformative
        bytes += sizeof(pppndxs[0]);
    }
    fp.write(reinterpret_cast<const char*>(intvals.GetVector()),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2DHDP1ind] += bytes;
    //}}//HDP1


    //profile amino acid sequence:
    sequence.reserve(pssm.GetSize()+4);
    for( m = 0; m < pssm.GetSize(); m++ )
        sequence += DehashCode(pssm.GetResidueAt(m));
    if( MAXFILESIZE <= addrfields[pmv2Daa])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid profile sequence address.");
    fp.seekp(std::streamoff(addrfields[pmv2Daa]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    bytes = pssm.GetSize();
    fp.write(sequence.c_str(),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2Daa] += bytes;

    //profile secondary structure state sequence:
    sequence.erase();
    for( m = 0; m < pssm.GetSize(); m++ )
        sequence += pssm.GetSSSSet()? DehashSSCodeLc(pssm.GetSSStateAt(m)): ' ';
    if( MAXFILESIZE <= addrfields[pmv2DSSstate])
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid sec. structure state sequence address.");
    fp.seekp(std::streamoff(addrfields[pmv2DSSstate]));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    bytes = pssm.GetSize();
    fp.write(sequence.c_str(),bytes);
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrfields[pmv2DSSstate] += bytes;


    //addresses of the ends of filenames are not included!
//     if( MAXFILESIZE <= addrfname_end_addrs)
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of filename end address.");
//     fp.seekp(std::streamoff(addrfname_end_addrs));
//     if(fp.fail())
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
//     bytes = addrfname + (pssm.GetName()? strlen(pssm.GetName())+1: 1);//including the terminal 0
//     fp.write(reinterpret_cast<const char*>(&bytes),sizeof(bytes));
//     if(fp.bad())
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
//     addrfname_end_addrs += sizeof(bytes);

    //address of the end of description:
    if( MAXFILESIZE <= addrdesc_end_addrs)
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid address of description end address.");
    fp.seekp(std::streamoff(addrdesc_end_addrs));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    bytes = addrdesc + (pssm.GetDescription()? strlen(pssm.GetDescription())+1: 1);//including the terminal 0
    fp.write(reinterpret_cast<const char*>(&bytes),sizeof(bytes));
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrdesc_end_addrs += sizeof(bytes);


    //profile filenames are not included!
//     if( MAXFILESIZE <= addrfname)
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid profile filename address.");
//     fp.seekp(std::streamoff(addrfname));
//     if(fp.fail())
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
//     if( pssm.GetName()) {
//         bytes = strlen(pssm.GetName())+1;//including the terminal 0
//         fp.write(pssm.GetName(),bytes);
//     } else {
//         bytes = 1;
//         fp.write(&ch0,bytes);
//     }
//     if(fp.bad())
//         throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
//     addrfname += bytes;

    //profile description:
    if( MAXFILESIZE <= addrdesc)
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Invalid profile description address.");
    fp.seekp(std::streamoff(addrdesc));
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Failed to set file position.");
    if( pssm.GetDescription()) {
        bytes = strlen(pssm.GetDescription())+1;//including the terminal 0
        fp.write(pssm.GetDescription(),bytes);
    } else {
        bytes = 1;
        fp.write(&ch0,bytes);
    }
    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile_v2_2: Write to file failed.");
    addrdesc += bytes;
}

// =========================================================================
// BinaryWriteProfile: write profile data to binary file
//
void BinaryWriteProfile( 
    std::ofstream& fp,
    const PMProfileModel& pssm, const PMTransModel& gaps )
{
    if( pssm.GetSize() != gaps.GetSize())
        throw MYRUNTIME_ERROR( "BinaryWriteProfile: Inconsistent Profile data." );

    const int noress = NUMAA;//PMProfileModelBase::PVDIM;

    extspsl::Pslvector uninfctxvec(pmv2DNoCVEls);//uninformative context vector
    const int Xuninf = MOptions::GetX_UNINF();
    constexpr bool mildmasking = true;//use mild masking
    mystring        sequence;
    float           bppprob;//background posterior predicitive
    const float*    ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    int             noppps;//number of p.p.probability values
    float           norm2;//context vector's squared norm
    float           lpprb;//log of context vector's prior probability
    const float*    ctxvec;//context vector
    int             vsize = pssm.GetCtxVecSize();//size of vector
    std::streamoff fppos = -1, fplst = -1;
    const float prob0 = 0.0f;
    const int ndxneg = -1;
    const char ch0 = 0;
    char res;
    int m, r, t;

    //init each entry to ensure CVS2S score==0 when CVS are not used
    uninfctxvec.AssignTo(0.17302947f);

    if((fppos = fp.tellp()) == (std::streamoff)-1 || fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Failed to determine file position.");

    //first, write total size of the profile (actual size to be written finally)
    int totprosize = 0;
    fp.write(reinterpret_cast<const char*>(&totprosize),sizeof(totprosize));

    //then, preamble
    fp.write(patstrDATBINVER,strlen(patstrDATBINVER));
    fp.write(PMProfileModel::dataversion,strlen(PMProfileModel::dataversion));
    //description and file:
    if( pssm.GetDescription())
        fp.write(pssm.GetDescription(),strlen(pssm.GetDescription())+1);//add 0
    else
        fp.write(&ch0,1);
    if( pssm.GetName())
        fp.write(pssm.GetName(),strlen(pssm.GetName())+1);//add 0
    else
        fp.write(&ch0,1);

    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

    //length:
    fp.write(reinterpret_cast<const char*>(&pssm.length_),sizeof(pssm.length_));
    //ENO:
    fp.write(reinterpret_cast<const char*>(&pssm.effnosequences_),sizeof(pssm.effnosequences_));

    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

    //background probabilities:
    for( r = 0; r < noress; r++ )
        fp.write(reinterpret_cast<const char*>(pssm.backprobs_+r),sizeof(pssm.backprobs_[0]));

    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

    //posterior probabilities:
    for( r = 0; r < noress; r++ )
        fp.write(reinterpret_cast<const char*>(pssm.postprobs_+r),sizeof(pssm.postprobs_[0]));

    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

    if( 0 < pssm.GetSize()) {
        //whole sequence at once
        sequence.reserve(pssm.GetSize()+4);
        for( m = 0; m < pssm.GetSize(); m++ )
            sequence += DehashCode(pssm.GetResidueAt(m));
        //fp.write(pssm.residues_,pssm.GetSize());
        fp.write(sequence.c_str(),pssm.GetSize());

        //beginning transition probabilities:
        for( t = 0; t < P_NSTATES; t++ ) {
            lpprb = gaps.GetOrgTrProbsAt(t,-1);
            //NOTE: calculate and write logarithms of transitions
            lpprb = lpprb? (P_MD<t? logf(lpprb)*MAIN_TRNPRB_LOGFCT: logf(lpprb)): -32768.0f;
            fp.write(reinterpret_cast<const char*>(&lpprb),sizeof(lpprb));
        }

        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");
    }

    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;

        //target probabilities:
        if(!(Xuninf && res == X && !mildmasking)) {
            for( r = 0; r < noress; r++ )
                fp.write(reinterpret_cast<const char*>(pssm.values_[m]+r),sizeof(pssm.values_[m][0]));
        } else {
            for( r = 0; r < noress; r++ )
                fp.write(reinterpret_cast<const char*>(pssm.backprobs_+r),sizeof(pssm.backprobs_[0]));
        }

        //transition probabilities:
        for( t = 0; t < P_NSTATES; t++ ) {
            lpprb = gaps.GetOrgTrProbsAt(t,m);
            //NOTE: calculate and write logarithms of transitions
            lpprb = lpprb? (P_MD<t? logf(lpprb)*MAIN_TRNPRB_LOGFCT: logf(lpprb)): -32768.0f;
            fp.write(reinterpret_cast<const char*>(&lpprb),sizeof(lpprb));
        }

        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

        //{{HDP1
        if(!(Xuninf && res == X && !mildmasking)) {
            bppprob = pssm.GetBckPPProbAt(m);
            ppprobs = pssm.GetPPProbsAt(m);
            pppndxs = pssm.GetPPPIndsAt(m);
            noppps = (int)pssm.GetNoPPProbsAt(m);
        } else {
            bppprob = 0.0f;
            ppprobs = &prob0;
            pppndxs = &ndxneg;
            noppps = 1;
        }

        //HDP1: background probability:
        fp.write(reinterpret_cast<const char*>(&bppprob),sizeof(bppprob));
        //HDP1: number of components:
        fp.write(reinterpret_cast<const char*>(&noppps),sizeof(noppps));

        if(fp.bad())
            throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

        if( 0 < noppps ) {
            if( ppprobs == NULL || pppndxs == NULL )
                throw MYRUNTIME_ERROR("BinaryWriteProfile: Null posterior probabilities.");

            //HDP1: component probabilities:
            for( t = 0; t < noppps; t++ )
                fp.write(reinterpret_cast<const char*>(ppprobs+t),sizeof(ppprobs[0]));

            //HDP1: component indices:
            for( t = 0; t < noppps; t++ )
                fp.write(reinterpret_cast<const char*>(pppndxs+t),sizeof(pppndxs[0]));

            if(fp.bad())
                throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");
        }
        //}}

        if( pssm.GetctPsSet()) {
            //HDP ctx
            if(!(Xuninf && res == X && !mildmasking)) {
                bppprob = pssm.GetctBckPPProbAt(m);
                ppprobs = pssm.GetctPPProbsAt(m);
                pppndxs = pssm.GetctPPPIndsAt(m);
                noppps = (int)pssm.GetNoctPPProbsAt(m);
            } else {
                bppprob = 0.0f;
                ppprobs = &prob0;
                pppndxs = &ndxneg;
                noppps = 1;
            }

            fp.write(patstrCT,strlen(patstrCT));
            //HDP ct: number of components:
            fp.write(reinterpret_cast<const char*>(&noppps),sizeof(noppps));
            //HDP ct: background probability:
            fp.write(reinterpret_cast<const char*>(&bppprob),sizeof(bppprob));

            if(fp.bad())
                throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

            if( 0 < noppps ) {
                if( ppprobs == NULL || pppndxs == NULL )
                    throw MYRUNTIME_ERROR("BinaryWriteProfile: Null ct posterior probabilities.");

                //HDP ct: component probabilities:
                for( t = 0; t < noppps; t++ )
                    fp.write(reinterpret_cast<const char*>(ppprobs+t),sizeof(ppprobs[0]));

                //HDP ct: component indices:
                for( t = 0; t < noppps; t++ )
                    fp.write(reinterpret_cast<const char*>(pppndxs+t),sizeof(pppndxs[0]));

                if(fp.bad())
                    throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");
            }
        }

        //{{ context vector
        if( pssm.GetCtxVecSet()) {
            if(!(Xuninf && res == X)) {
                norm2 = pssm.GetCtxVecNorm2At(m);
                lpprb = pssm.GetCtxVecLpprobAt(m);
                ctxvec = pssm.GetCtxVecAt(m);
            } else {
                norm2 = 0.0f;
                lpprb = 0.0f;
                ctxvec = uninfctxvec.GetVector();
            }

            fp.write(patstrCV,strlen(patstrCV));
            //CV: probability:
            fp.write(reinterpret_cast<const char*>(&lpprb),sizeof(lpprb));
            //CV: squared norm:
            fp.write(reinterpret_cast<const char*>(&norm2),sizeof(norm2));
            //CV: size:
            //fp.write(reinterpret_cast<const char*>(&vsize),sizeof(vsize));

            if(fp.bad())
                throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

            if( 0 < vsize ) {
                if( ctxvec == NULL )
                    throw MYRUNTIME_ERROR("BinaryWriteProfile: Null context vector.");

                for( t = 0; t < vsize; t++ )
                    fp.write(reinterpret_cast<const char*>(ctxvec+t),sizeof(ctxvec[0]));

                if(fp.bad())
                    throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");
            }
        }
        //}}

        //{{ SSS
        if( pssm.GetSSSSet()) {
            res = DehashSSCodeLc(pssm.GetSSStateAt(m));
            //use floor in rounding SS state probability
            fp.write(patstrSS,strlen(patstrSS));
            fp.write(&res,1);
            //NOTE: always use SP3 format
            //if( pssm.GetSSSP3()) {
                if(!(Xuninf && res == X)) {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fp.write(reinterpret_cast<const char*>(pssm.sssprob_[m]+t),sizeof(pssm.sssprob_[m][0]));
                } else {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fp.write(reinterpret_cast<const char*>(&prob0),sizeof(prob0));
                }
            //}

            if(fp.bad())
                throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");
        }
        //}}
    }

    //ending:
    fp.write(patstrEND,strlen(patstrEND));

    if(fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Write to file failed.");

    //write the actual size:
    if((fplst = fp.tellp()) == (std::streamoff)-1 || fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Failed to determine file position.");

    totprosize = (int)(fplst - fppos);
    if( totprosize < 1)
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Invalid profile size obtained.");

    fp.seekp(fppos);
    if(fp.fail())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Failed to set file position.");

    fp.write(reinterpret_cast<const char*>(&totprosize),sizeof(totprosize));

    fp.seekp(fplst);
    if(fp.fail() || fp.bad())
        throw MYRUNTIME_ERROR("BinaryWriteProfile: Failed to set file position.");
}

// =========================================================================
// PrintParameters: print profile statistical parameters 
//
void PMProfileModel::PrintParameters( FILE* fp ) const
{
    if( fp == NULL )
        return;

    if( 0.0f <= GetExpectedScore()) {
        fprintf( fp, "%s %.4f!%s", patstrEXPNN, GetExpectedScore(), NL );
        return;
    }

    fprintf( fp, "%-25s  %-6s   %-6s%s", " ", "K", "Lambda", NL );
    fprintf( fp, "%-25s  %6.4f   %6.4f%s", patstrSPCOMP, GetK(), GetLambda(), NL );
    fprintf( fp, "%-25s  %6.4f   %6.4f%s", patstrSPREFR, GetRefK(), GetRefLambda(), NL );
    fprintf( fp, "%s %6.4f; %s %6.4f%s", patstrENTROPY, GetEntropy(), 
             patstrEXPECTED, GetExpectedScore(), NL );
}
// PrintParametersCondensed: print profile statistical parameters in condensed form
//
void PMProfileModel::PrintParametersCondensed( FILE* fp ) const
{
    if( fp == NULL )
        return;

//     if( 0.0f <= GetExpectedScore()) {
//         fprintf( fp, "%s %.4f!%s", patstrEXPNN, GetExpectedScore(), NL );
//         return;
//     }

    fprintf( fp, "%s, %s:%s","K", "Lambda", NL);
    fprintf( fp, "%s %6.4f %6.4f%s", patstrSPCOMP, GetK(), GetLambda(), NL );
    fprintf( fp, "%s %6.4f %6.4f%s", patstrSPREFR, GetRefK(), GetRefLambda(), NL );
    fprintf( fp, "%s %6.4f; %s %6.4f%s", patstrENTROPY, GetEntropy(), 
             patstrEXPECTED, GetExpectedScore(), NL );
}

// =========================================================================
// TextWriteProfileCondensed: write condensed profile data to file
//
void TextWriteProfileCondensed( 
    FILE* fp,
    const PMProfileModel& pssm, const PMTransModel& gaps,
    int scale )
{
    if( !fp )
        return;

    if( pssm.GetSize() != gaps.GetSize())
        throw MYRUNTIME_ERROR( "TextWriteProfileCondensed: Inconsistent Profile data." );

    if( scale < 1 )
        throw MYRUNTIME_ERROR( "TextWriteProfileCondensed: Invalid scale factor." );

    const int noress = NUMAA;//PMProfileModelBase::PVDIM;
    const int showcmdline = MOptions::GetSHOWCMD();

    extspsl::Pslvector uninfctxvec(pmv2DNoCVEls);//uninformative context vector
    const int Xuninf = MOptions::GetX_UNINF();
    constexpr bool mildmasking = true;//use mild masking
    float           bppprob;//background posterior predicitive
    const float*    ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    size_t          noppps;//number of p.p.probability values
    float           norm2;//context vector's squared norm
    float           lpprb;//log of context vector's prior probability
    const float*    ctxvec;//context vector
    int             vsize = pssm.GetCtxVecSize();//size of vector
    const float prob0 = 0.0f;
    const int ndxneg = -1;
    char    res;
    int     m, r, t, n = 0;

    //init each entry to ensure CVS2S score==0 when CVS are not used
    uninfctxvec.AssignTo(0.17302947f);

    fprintf( fp, "%s%s%s", patstrDATVER, PMProfileModel::dataversion, NL );
    fprintf( fp, "%s %s%s", patstrDESC, pssm.GetDescription()? pssm.GetDescription(): "", NL );
    fprintf( fp, "%s %s%s", patstrFILE, pssm.GetName()? pssm.GetName(): "", NL );

    fprintf( fp, "%s ", patstrCMD );
    if( showcmdline )
        print_cmdline( &file_print, fp );
    else
        fprintf( fp, "%s", NL);
    fprintf( fp, "%s %d%s", patstrLEN, pssm.GetSize(), NL );
    fprintf( fp, "%s %d%s", patstrNOS, (int)pssm.GetNoSequences(), NL );
    fprintf( fp, "%s %.1f%s", patstrEFFNOS, pssm.GetEffNoSequences(), NL );
    fprintf( fp, "%s %d%s", patstrSCALE, scale, NL );

    for( r = 0; r < noress; r++ )
        fprintf( fp, " %c", DehashCode(r));

    fprintf( fp, "%s%s", NL,patstrNULL);
    for( r = 0; r < noress; r++ )
        fprintf( fp, " %d", (int)rintf(scale*pssm.GetBackProbsAt(r)));

    fprintf( fp, "%s%s", NL,patstrPOST);
    for( r = 0; r < noress; r++ )
        fprintf( fp, " %d", (int)rintf(scale*pssm.GetPostProbsAt(r)));

    if( 0 < pssm.GetSize()) {
        fprintf( fp, "%s", NL);
        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, " %d", (int)rintf(scale*gaps.GetOrgTrProbsAt(t,-1)));

        fprintf( fp, "%s", NL);
        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, " %d", (int)rintf(INTSCALE*pssm.GetMIDExpNoObservationsBeg(t)));
    }

    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;

        fprintf( fp, "%s%d %c", NL, ++n, DehashCode(pssm.GetResidueAt(m)));

        if(!(Xuninf && res == X && !mildmasking)) {
            for( r = 0; r < noress; r++ )
                fprintf( fp, " %d", (int)rintf(scale*pssm.GetValueAt(m,r)));
        } else {
            for( r = 0; r < noress; r++ )
                fprintf( fp, " %d", (int)rintf(scale*pssm.GetBackProbsAt(r)));
        }

        fprintf( fp, "%s", NL);

//         for( r = 0; r < noress; r++ )//NOTE:!
//             fprintf( fp, "%7d ", ( int )rintf( scale * (*pssm.GetObsFreqsAt(m))[r] ));

//         fprintf( fp, "%s", NL);

        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, " %d", (int)rintf(scale*gaps.GetOrgTrProbsAt(t,m)));

        fprintf( fp, "%s", NL);

        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, " %d", (int)rintf(INTSCALE*pssm.GetMIDExpNoObservationsAt(m,t)));

        //{{HDP1
        if(!(Xuninf && res == X && !mildmasking)) {
            bppprob = pssm.GetBckPPProbAt(m);
            ppprobs = pssm.GetPPProbsAt(m);
            pppndxs = pssm.GetPPPIndsAt(m);
            noppps = (int)pssm.GetNoPPProbsAt(m);
        } else {
            bppprob = 0.0f;
            ppprobs = &prob0;
            pppndxs = &ndxneg;
            noppps = 1;
        }

        fprintf( fp, "%s", NL);
        fprintf( fp, " %d %d", (int)rintf(scale*bppprob), (int)noppps);

        if( 0 < (int)noppps ) {
            if( ppprobs == NULL || pppndxs == NULL )
                throw MYRUNTIME_ERROR("TextWriteProfileCondensed: Null posterior probabilities.");

            fprintf( fp, "%s", NL);

            for( t = 0; t < (int)noppps; t++ )
                fprintf( fp, " %d", (int)rintf(scale*ppprobs[t]));

            fprintf( fp, "%s", NL);

            for( t = 0; t < (int)noppps; t++ )
                fprintf( fp, " %d", pppndxs[t]);
        }
        //}}

        if( pssm.GetctPsSet()) {
            //HDP ctx
            if(!(Xuninf && res == X && !mildmasking)) {
                bppprob = pssm.GetctBckPPProbAt(m);
                ppprobs = pssm.GetctPPProbsAt(m);
                pppndxs = pssm.GetctPPPIndsAt(m);
                noppps = pssm.GetNoctPPProbsAt(m);
            } else {
                bppprob = 0.0f;
                ppprobs = &prob0;
                pppndxs = &ndxneg;
                noppps = 1;
            }

            fprintf( fp, "%s %s%d %d", NL, patstrCT, 
                     (int)noppps, (int)rintf(scale*bppprob));

            if( 0 < (int)noppps ) {
                if( ppprobs == NULL || pppndxs == NULL )
                    throw MYRUNTIME_ERROR("TextWriteProfileCondensed: Null ct posterior probabilities.");

                for( t = 0; t < (int)noppps; t++ )
                    fprintf( fp, " %d", (int)rintf(scale*ppprobs[t]));

                for( t = 0; t < (int)noppps; t++ )
                    fprintf( fp, " %d", pppndxs[t]);
            }
        }

        //{{ context vector
        if( pssm.GetCtxVecSet()) {
            if(!(Xuninf && res == X)) {
                norm2 = pssm.GetCtxVecNorm2At(m);
                lpprb = pssm.GetCtxVecLpprobAt(m);
                ctxvec = pssm.GetCtxVecAt(m);
            } else {
                norm2 = 0.0f;
                lpprb = 0.0f;
                ctxvec = uninfctxvec.GetVector();
            }

            fprintf( fp, "%s %s %d %d %d", NL, patstrCV,
                   (int)rintf(scale*lpprb), (int)rintf(scale*norm2), vsize );

            if( 0 < vsize ) {
                if( ctxvec == NULL )
                    throw MYRUNTIME_ERROR("TextWriteProfileCondensed: Null context vector.");

                for( t = 0; t < vsize; t++ )
                    fprintf( fp, " %d", (int)rintf(scale*ctxvec[t]));
            }
        }
        //}}

        //{{
        if( pssm.GetSSSSet()) {
            //use floor in rounding SS state probability
            fprintf( fp, "%s %s%c", NL, patstrSS, DehashSSCode(pssm.GetSSStateAt(m)));
            if( pssm.GetSSSP3()) {
                if(!(Xuninf && res == X)) {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fprintf( fp, " %d", (int)(scale*pssm.GetSSStateProbAt(m,t)));
                } else {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fprintf( fp, " %d", 0);
                }
            }
            else
                fprintf( fp, " %d", (int)(scale*pssm.GetSSStateProbAt(m)));
        }
        //}}
    }

    fprintf( fp, "%s", NL);

// #ifdef SCALE_PROFILES
    pssm.PrintParametersCondensed(fp);
// #endif
    fprintf( fp, "%s%s", patstrEND, NL );
}

// =========================================================================
// TextWriteProfile: write profile data to file
//
void TextWriteProfile( 
    FILE* fp,
    const PMProfileModel& pssm, const PMTransModel& gaps,
    int scale )
{
    if( !fp )
        return;

    if( pssm.GetSize() != gaps.GetSize())
        throw MYRUNTIME_ERROR( "TextWriteProfile: Inconsistent Profile data." );

    if( scale < 1 )
        throw MYRUNTIME_ERROR( "TextWriteProfile: Invalid scale factor." );

    const int noress = NUMAA;//PMProfileModelBase::PVDIM;
    const int showcmdline = MOptions::GetSHOWCMD();

    extspsl::Pslvector uninfctxvec(pmv2DNoCVEls);//uninformative context vector
    const int Xuninf = MOptions::GetX_UNINF();
    constexpr bool mildmasking = true;//use mild masking
    float           bppprob;//background posterior predicitive
    const float*    ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    size_t          noppps;//number of p.p.probability values
    float           norm2;//context vector's squared norm
    float           lpprb;//log of context vector's prior probability
    const float*    ctxvec;//context vector
    int             vsize = pssm.GetCtxVecSize();//size of vector
    const float prob0 = 0.0f;
    const int ndxneg = -1;
    char    res;
    int     m, r, t, n = 0;

    //init each entry to ensure CVS2S score==0 when CVS are not used
    uninfctxvec.AssignTo(0.17302947f);

    fprintf( fp, "%s%s%s", patstrDATVER, PMProfileModel::dataversion, NL );
    fprintf( fp, "%-9s%s%s", patstrDESC, pssm.GetDescription()? pssm.GetDescription(): "", NL );
    fprintf( fp, "%-9s%s%s", patstrFILE, pssm.GetName()? pssm.GetName(): "", NL );

    fprintf( fp, "%-9s", patstrCMD );
    if( showcmdline )
        print_cmdline( &file_print, fp );
    fprintf( fp, "%-9s%d%s", patstrLEN, pssm.GetSize(), NL );
    fprintf( fp, "%-9s%d%s", patstrNOS, (int)pssm.GetNoSequences(), NL );
    fprintf( fp, "%-9s%.1f%s", patstrEFFNOS, pssm.GetEffNoSequences(), NL );
    fprintf( fp, "%-9s%d%s", patstrSCALE, scale, NL );

    fprintf( fp, "# Target / Observed / ");
    for( t = 0; t < P_NSTATES; t++ ) fprintf( fp, "%s ", gTPTRANS_NAMES[t]);
    fprintf( fp, "/ NexpM NexpI NexpD; Weight; Information /%s", NL );
    fprintf( fp, "# BPProb NoPosts / PProbs / PPNdxs / CT:NoPosts BPProb PProbs PPNdxs ");
    fprintf( fp, "/ CV: Prior Norm2 N Vector / SS:Pred Prob (C,E,H)%s", NL );
// //     fprintf( fp, "/ NexpM NexpI NexpD; Weight; Information; Gap weight, deletion at beg., end, interval\n");
// //     fprintf( fp, "/ Weight, information, exp. no. observations; Gap weight, deletion at beg., end, interval\n");

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < noress; r++ )
        fprintf( fp, " %7c", DehashCode( r ) );


    fprintf( fp, "%s%7s   ", NL, patstrNULL );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", ( int )rintf( scale * pssm.GetBackProbsAt( r )));

    fprintf( fp, "%s%7s   ", NL, patstrPOST );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", ( int )rintf( scale * pssm.GetPostProbsAt( r )));

    if( 0 < pssm.GetSize()) {
        fprintf( fp, "%s%10c", NL, 32 );
        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rintf( scale * gaps.GetOrgTrProbsAt( t, -1 )));

        fprintf( fp, "%s%10c", NL, 32 );
        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rintf( INTSCALE * pssm.GetMIDExpNoObservationsBeg( t )));
    }

    for( m = 0; m < pssm.GetSize(); m++ ) {
        res = pssm.GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;

        fprintf( fp, "%s%5d %c   ", NL, ++n, DehashCode( pssm.GetResidueAt( m )));

        if(!(Xuninf && res == X && !mildmasking)) {
            for( r = 0; r < noress; r++ )
                fprintf( fp, "%7d ", (int)rintf(scale * pssm.GetValueAt(m, r)));
        } else {
            for( r = 0; r < noress; r++ )
                fprintf( fp, "%7d ", (int)rintf(scale * pssm.GetBackProbsAt(r)));
        }

        fprintf( fp, "%s%10c", NL, 32 );

        for( r = 0; r < noress; r++ )//NOTE:!
            fprintf( fp, "%7d ", ( int )rintf( scale * (*pssm.GetObsFreqsAt(m))[r] ));

        fprintf( fp, "%s%10c", NL, 32 );

        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rintf( scale * gaps.GetOrgTrProbsAt( t, m )));

        fprintf( fp, "%s%10c", NL, 32 );

        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rintf( INTSCALE * pssm.GetMIDExpNoObservationsAt( m, t )));

        fprintf( fp, "%7d %7d ",
                ( int )rintf( scale * pssm.GetFrequencyWeightAt( m )),
                ( int )rintf( scale * pssm.GetInformationAt( m )));

// //         fprintf( fp, "%7d %7d %7d %7d ",
// //                 ( int )rint( scale * gaps.GetWeightsAt( m )),
// //                 ( int )rint( scale * gaps.GetDeletesBegAt( m )),
// //                 ( int )rint( scale * gaps.GetDeletesEndAt( m )),
// //                 gaps.GetDeletesIntervalAt( m ));

        //{{HDP1
        if(!(Xuninf && res == X &&! mildmasking)) {
            bppprob = pssm.GetBckPPProbAt(m);
            ppprobs = pssm.GetPPProbsAt(m);
            pppndxs = pssm.GetPPPIndsAt(m);
            noppps = (int)pssm.GetNoPPProbsAt(m);
        } else {
            bppprob = 0.0f;
            ppprobs = &prob0;
            pppndxs = &ndxneg;
            noppps = 1;
        }

        fprintf( fp, "%s%10c", NL, 32 );
        fprintf( fp, "%7d %7d ", (int)rintf( scale*bppprob ), (int)noppps );

        if( 0 < (int)noppps ) {
            if( ppprobs == NULL || pppndxs == NULL )
                throw MYRUNTIME_ERROR("TextWriteProfile: Null posterior probabilities.");

            fprintf( fp, "%s%10c", NL, 32 );

            for( t = 0; t < (int)noppps; t++ )
                fprintf( fp, "%7d ", (int)rintf( scale*ppprobs[t]));

            fprintf( fp, "%s%10c", NL, 32 );

            for( t = 0; t < (int)noppps; t++ )
                fprintf( fp, "%7d ", pppndxs[t]);
        }
        //}}

        if( pssm.GetctPsSet()) {
            //HDP ctx
            if(!(Xuninf && res == X && !mildmasking)) {
                bppprob = pssm.GetctBckPPProbAt(m);
                ppprobs = pssm.GetctPPProbsAt(m);
                pppndxs = pssm.GetctPPPIndsAt(m);
                noppps = (int)pssm.GetNoctPPProbsAt(m);
            } else {
                bppprob = 0.0f;
                ppprobs = &prob0;
                pppndxs = &ndxneg;
                noppps = 1;
            }

            fprintf( fp, "%s%13c%s%d %7d ", NL, 32, patstrCT, 
                     (int)noppps, (int)rintf( scale*bppprob ));

            if( 0 < (int)noppps ) {
                if( ppprobs == NULL || pppndxs == NULL )
                    throw MYRUNTIME_ERROR("TextWriteProfile: Null ct posterior probabilities.");

                for( t = 0; t < (int)noppps; t++ )
                    fprintf( fp, "%7d ", (int)rintf( scale*ppprobs[t]));

                for( t = 0; t < (int)noppps; t++ )
                    fprintf( fp, "%7d ", pppndxs[t]);
            }
        }

        //{{ context vector
        if( pssm.GetCtxVecSet()) {
            if(!(Xuninf && res == X)) {
                norm2 = pssm.GetCtxVecNorm2At(m);
                lpprb = pssm.GetCtxVecLpprobAt(m);
                ctxvec = pssm.GetCtxVecAt(m);
            } else {
                norm2 = 0.0f;
                lpprb = 0.0f;
                ctxvec = uninfctxvec.GetVector();
            }

            fprintf( fp, "%s%14c%s %7d %7d %7d ", NL, 32, patstrCV,
                   (int)rintf( scale*lpprb ), (int)rintf( scale*norm2 ), vsize );

            if( 0 < vsize ) {
                if( ctxvec == NULL )
                    throw MYRUNTIME_ERROR("TextWriteProfile: Null context vector.");

                for( t = 0; t < vsize; t++ )
                    fprintf( fp, "%7d ", (int)rintf( scale*ctxvec[t]));
            }
        }
        //}}

        //{{
        if( pssm.GetSSSSet()) {
            //use floor in rounding SS state probability
            fprintf( fp, "%s%13c%s%c", NL, 32, patstrSS, DehashSSCode( pssm.GetSSStateAt(m)));
            if( pssm.GetSSSP3()) {
                if(!(Xuninf && res == X)) {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fprintf( fp, " %7d", (int)(scale*pssm.GetSSStateProbAt(m,t)));
                } else {
                    for( t = 0; t < SS_NSTATES; t++ )
                        fprintf( fp, " %7d", 0);
                }
            }
            else
                fprintf( fp, " %7d", (int)( scale*pssm.GetSSStateProbAt(m)));
        }
        //}}
    }

    //{{NOTE:THIS BLOCK TO BE REMOVED
    fprintf( fp, "%s%10c%7d %7d", NL, 32,
                ( int )rintf( scale * (-7.f)),
                ( int )rintf( scale * (-1.f)));
    //}}
    fprintf( fp, "%s", NL );

// #ifdef SCALE_PROFILES
    pssm.PrintParameters( fp );
// #endif
    fprintf( fp, "%s%s", patstrEND, NL );
}

// -------------------------------------------------------------------------
// TextReadProfile: read profile data from file
// TODO: reimplement large arrays to use heap (locbuffer)
//
void TextReadProfile( FILE* fp, PMProfileModel& pssm, PMTransModel& gaps )
{
    if( !fp )
        return;

    const mystring preamb = "TextReadProfile: ";
    const float     accuracy = 1.0e-3f;
    const int       maxnocls = 10000;
    const int       noress = NUMAA;
    const int       szdata = PMProfileModelBase::PVDIM;
    size_t          length, rbts;
    bool            lineread = false;
    const size_t    locsize = 10*KBYTE;
    char            locbuffer[locsize+1] = {0};
    char            msgbuf[BUF_MAX];
    char*           p;
    int             emsg;
    mystring    desc, file;
    const char* descp, *filep;
    int         prolen = 0;
    int         nos = 0;
    float       effnos = 0;
    int         scale = 0;
    char        res, sss;
    int         intval;
    float       value, consv;
    float       scores[szdata];
    float       freqns[szdata];
    float       trnpro[P_NSTATES];
    float       midobs[PS_NSTATES];
    float       weight = 0.0f, inform = 0.0f;
// //     float       gapwgt, delbeg, delend, open, extend;
    float       /*expnobs, */sssprob[SS_NSTATES];
// //     int         delint;
    float       bppprob = 0.f, ctbppprob = 0.f;//background posterior predictive
    extspsl::Pslvector ppprobs, ctppprobs;//posterior predictive probabilities
    extspsl::Ivector   pppndxs, ctpppndxs;//indices of posteriors
    size_t      noppps = 0, noctppps = 0;//number of p.p.probability values
    float       norm2 = 0.f, lpprb = 0.f;//squared norm and log of prior probability for context vector
    extspsl::Pslvector ctxvec;//context vector
    int         vsize = 0;//vector size
    float       refunlambda, refunK;
    float       cmpunlambda, cmpunK;
    float       proentropy,  proexpected;
    int         m, r, t, t0 = 0, tn = 0, n = 0;

    memset( scores, 0, szdata * sizeof( float ));
    memset( freqns, 0, szdata * sizeof( float ));

    //read version number
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrDATVER )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrDATVER );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( p, PMProfileModel::dataversion )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Inappropriate profile version number." );


    //read description line
    if(( emsg = skip_comments( fp, desc )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || desc.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( desc.c_str(), patstrDESC )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    descp = p + strlen( patstrDESC );
    for( ; *descp == ' ' || *descp == '\t'; descp++ );
    for( n = (int)desc.length() - 1; 0 <= n && ( desc[n]=='\n' || desc[n]=='\r' ); n-- )
        desc[n] = 0;


    //read filename
    if(( emsg = skip_comments( fp, file )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || file.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( file.c_str(), patstrFILE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    filep = p + strlen( patstrFILE );
    for( ; *filep == ' ' || *filep == '\t' ; filep++ );
    for( n = (int)file.length() - 1; 0 <= n && ( file[n]=='\n' || file[n]=='\r' ); n-- )
        file[n] = 0;


    //read command line
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrCMD )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    //read profile length
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrLEN )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrLEN );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &prolen, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length." );


    //read number of sequences
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrNOS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrNOS );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &nos, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( nos < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of sequences." );


    //read effective number of sequences
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrEFFNOS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrEFFNOS );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_float( p, length - size_t( p - locbuffer ), &effnos, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( effnos <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid effective number of observations." );


    //read scale factor
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrSCALE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrSCALE );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &scale, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( scale < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale factor." );


    //read residue letters; omit checking the order
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));


    //read background probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrNULL )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No background probabilities." );

    p += strlen( patstrNULL );

    consv = 0.0f;
    for( r = 0; r < noress; r++ ) {
        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete background probability data." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || scale < intval )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Invalid background probability value." );

        value = (float)intval / (float)scale;
        consv += value;
        pssm.SetBackProbsAt( r, value );
    }

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid profile data: Invalid background probabilities." );
    //


    //read posterior (generalized target) probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrPOST )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No posterior probabilities." );

    p += strlen( patstrPOST );

    consv = 0.0f;
    for( r = 0; r < noress; r++ ) {
        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete posterior probability data." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || scale < intval )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Invalid posterior probability value." );

        value = (float)intval / (float)scale;
        consv += value;
        pssm.SetPostProbsAt( r, value );
    }

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid profile data: Invalid posterior probabilities." );
    //


    //beginning transition probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    for( t = 0; t < P_NSTATES; t++ )
    {
        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete transition probability data." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        trnpro[t] = (float)intval / (float)scale;
    }
    //


    //beginning MID state expected observations
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    for( t = 0; t < PS_NSTATES; t++ )
    {
        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete MID state obervations data." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        midobs[t] = (float)intval / (float)INTSCALE;
    }
    //

    if( MAXCOLUMNS < prolen )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Too large profile length." );

    pssm.Clear();
    gaps.Clear();

    pssm.Reserve( prolen );
    gaps.Reserve( prolen );

    if( *descp )
        pssm.SetDescription( descp );
    else
        pssm.SetDescription( NULL );

    if( *filep )
        pssm.SetName( filep );
    else
        pssm.SetName( NULL );

    pssm.SetNoSequences( nos );
    pssm.SetEffNoSequences( effnos );
    pssm.SetMIDExpNoObservationsBeg( midobs );

    gaps.SetOrgTrProbsBeg( &trnpro );//set beginning transition probabilities here

    lineread = false;

    for( m = 0; m < prolen; m++ )
    {
        //scores/target probabilities
        if( !lineread ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }
        }

        lineread = false;

        if(( emsg = read_integer( p = locbuffer, length, &n, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( n != m + 1 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: Wrong numbering.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p += rbts;

        if( length <= size_t( p - locbuffer )) {
            sprintf( msgbuf, "Wrong profile format at pos %d: No residue.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_symbol( p, length - size_t( p - locbuffer ), &res, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        p += rbts;

        res = HashAlphSymbol( res );

        for( r = 0; r < noress; r++ )
        {
            if( length <= size_t( p - locbuffer )) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete target probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            scores[r] = (float)intval / (float)scale;
        }


        //observed frequencies
        //NOTE: condensed format does not include this section
//         if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
//         }
// 
//         if( feof( fp ) || !length ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }
// 
//         for( r = 0; r < noress; r++ )
//         {
//             if( length <= size_t( p - locbuffer )) {
//                 sprintf( msgbuf, "Wrong profile format at pos %d: "
//                                  "Incomplete observed frequency data.", m );
//                 throw MYRUNTIME_ERROR( preamb + msgbuf );
//             }
// 
//             if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
//                 sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
//                 throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
//             }
// 
//             p += rbts;
// 
//             freqns[r] = (float)intval / (float)scale;
//         }


        //transition probabilities
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || !length ) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        for( t = 0; t < P_NSTATES; t++ )
        {
            if( length <= size_t( p - locbuffer )) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete transition probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            trnpro[t] = (float)intval / (float)scale;
        }


        //MID state observations
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || !length ) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        for( t = 0; t < PS_NSTATES; t++ )
        {
            if( length <= size_t( p - locbuffer )) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete MID state observations data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( intval < 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Invalid MID state observations value.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p += rbts;

            midobs[t] = (float)intval / (float)INTSCALE;
        }


        //weight, information
        //NOTE: condensed format does not include this section
//         if( length <= size_t( p - locbuffer )) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: "
//                               "No frequency weight.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
//             throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
// 
//         if( intval < 0 ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: "
//                               "Negative frequency weight.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }
// 
//         weight = (float)intval / (float)scale;
// 
//         p += rbts;
// 
//         if( length <= size_t( p - locbuffer )) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: "
//                               "No information content.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
//         }
// 
//         if( intval < 0 ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: "
//                               "Negative information content.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }
// 
//         inform = (float)intval / (float)scale;
// 
//         p += rbts;


        //CLUSTERS: HDP1
        //bck posterior probability, no. posterior probability values
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || !length ) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if( length <= size_t( p - locbuffer )) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        bppprob = (float)intval / (float)scale;

        p += rbts;

        if( length <= size_t( p - locbuffer )) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }
        if( maxnocls < intval ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Too large number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        noppps = (size_t)intval;
        if( noppps ) {
            ppprobs.Allocate((int)noppps );
            pppndxs.Allocate((int)noppps );
            ppprobs.Clear();
            pppndxs.Clear();

            //posterior predictive probabilities
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            for( t = 0; t < (int)noppps; t++ )
            {
                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Incomplete HDP1 posterior probability data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative HDP1 posterior probability value.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                p += rbts;

                ppprobs.Push( (float)intval / (float)scale );
            }

            //indices of posteriors
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                  "No HDP1 cluster indices.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            for( t = 0; t < (int)noppps; t++ )
            {
                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                      "Incomplete HDP1 cluster index data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval != -1 && intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                      "Negative HDP1 cluster index.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                p += rbts;

                pppndxs.Push( intval );
            }
        }//if( noppps )


        //HDP ctx
        if( m < 1 )
            pssm.SetctPsSet( true );//probe for HDP ctx information
        if( pssm.GetctPsSet()) {
            //bck posterior probability, no. posterior probability values
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( locbuffer, patstrCT )) == NULL ) {
                if( m < 1 )
                    pssm.SetctPsSet( false );//no HDP ctx information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No HDP ctx probabilities.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCT );

                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No number of HDP ctx probability values.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative number of HDP ctx probability values.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( maxnocls < intval ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Too large number of HDP ctx probability values.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                noctppps = (size_t)intval;

                p += rbts;

                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No HDP ctx background probability.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative HDP ctx background probability.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                ctbppprob = (float)intval / (float)scale;

                p += rbts;

                if( noctppps ) {
                    ctppprobs.Allocate((int)noctppps );
                    ctpppndxs.Allocate((int)noctppps );
                    ctppprobs.Clear();
                    ctpppndxs.Clear();

                    //posterior predictive probabilities
                    for( t = 0; t < (int)noctppps; t++ )
                    {
                        if( length <= size_t( p - locbuffer )) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                            "Incomplete HDP ctx probability data.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                        }

                        if( intval < 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                            "Negative HDP ctx probability.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        p += rbts;

                        ctppprobs.Push( (float)intval / (float)scale );
                    }

                    //indices of posteriors
                    for( t = 0; t < (int)noctppps; t++ )
                    {
                        if( length <= size_t( p - locbuffer )) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                            "Incomplete HDP ctx cluster index data.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                        }

                        if( intval != -1 && intval < 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                            "Negative HDP ctx cluster index.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        p += rbts;

                        ctpppndxs.Push( intval );
                    }
                }//if( noppps )
            }//else
        }//if( pssm.GetctPsSet())


        //context vector
        if( m < 1 )
            pssm.SetCtxVecSet( true );//probe for context vector data
        if( pssm.GetCtxVecSet()) {
            if( !lineread )
                if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( locbuffer, patstrCV )) == NULL ) {
                if( m < 1 )
                    pssm.SetCtxVecSet( false );//no context vector data
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCV );

                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No prior for context vector.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                lpprb = (float)intval / (float)scale;

                p += rbts;

                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                norm2 = (float)intval / (float)scale;

                p += rbts;

                if( length <= size_t( p - locbuffer )) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( 1000 < intval ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Too large context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( vsize && vsize != intval ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Inconsistent context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                vsize = intval;

                p += rbts;

                if( vsize ) {
                    ctxvec.Allocate( vsize );
                    ctxvec.Clear();

                    //vector elements
                    for( t = 0; t < vsize; t++ )
                    {
                        if( length <= size_t( p - locbuffer )) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "Incomplete context vector data.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                        }

                        p += rbts;

                        ctxvec.Push( (float)intval / (float)scale );
                    }
                }//if( vsize )
            }//else
        }//if( pssm.GetCtxVecSet())


        //SS state
        if( m < 1 ) {
            pssm.SetSSSSet( true );//probe SS information existance
            pssm.SetSSSP3( true );
            t0 = 0; tn = SS_NSTATES-1;
        }
        if( pssm.GetSSSSet()) {
            if( !lineread )
                if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || !length ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( locbuffer, patstrSS )) == NULL ) {
                if( m < 1 )
                    pssm.SetSSSSet( false );//no SS information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: No SS state.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                memset( sssprob, 0, SS_NSTATES*sizeof(float));

                lineread = false;
                p += strlen( patstrSS );

                if(( emsg = read_symbol( p, length - size_t( p - locbuffer ), &sss, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                p += rbts;

                sss = HashSSState( sss );

                if( !pssm.GetSSSP3()) {
                    t0 = sss; tn = sss;
                }

                for( t = t0; t <= tn; t++ )
                {
                    if( length <= size_t( p - locbuffer )) {
                        if( m < 1 && t == t0+1 ) {
                            pssm.SetSSSP3( false );
                            sssprob[(int)sss] = (float)intval / (float)scale;
                            break;
                        }
                        else {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "No SS state probability.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }
                    }

                    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 ) {
                        sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                        throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                    }

                    p += rbts;
                    for( ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n'; p++ );

                    sssprob[t] = (float)intval / (float)scale;
                }
            }
        }//if( pssm.GetSSSSet())


        //SAVE profile position
// //         frqs.PushAt( freqns, res, m );
        pssm.PushAt( scores, freqns, res, weight, inform, midobs, m );
        pssm.SetBckPPProbAt( bppprob, m );
        pssm.SetPPProbsAt( m, ppprobs.GetVector(), pppndxs.GetVector(), (int)noppps );
        if( pssm.GetctPsSet()) {
            pssm.SetctBckPPProbAt( ctbppprob, m );
            pssm.SetctPPProbsAt( m, ctppprobs.GetVector(), ctpppndxs.GetVector(), (int)noctppps );
        }
        if( pssm.GetCtxVecSet()) pssm.SetCtxVecPlusAt( m, norm2, lpprb, ctxvec.GetVector(), vsize );
        if( pssm.GetSSSSet()) pssm.SetSSStateAt( m, sss, sssprob );
// //         gaps.PushAt(&trnpro, gapwgt, delbeg, delend, delint, res, m );
        gaps.PushAt(&trnpro, m );
    }


    gaps.Initialize();


// #ifdef SCALE_PROFILES
    //statistical parameters
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: No statistical parameters." );

    if(( p = strstr( locbuffer, patstrEXPNN )) != NULL ) {
        p += strlen( patstrEXPNN );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No expected value." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( proexpected < 0.0f )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: Wrong expected value." );

        pssm.SetExpectedScore( proexpected );
    }
    else {
        //computed
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameters." );

        if(( p = strstr( locbuffer, patstrSPCOMP )) == NULL )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameters." );

        p += strlen( patstrSPCOMP );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameters." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &cmpunK, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

//         if( cmpunK < 0.0f )
//             throw MYRUNTIME_ERROR( preamb + 
//             "Wrong profile format at the end: Invalid statistical parameters." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameter lambda." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &cmpunlambda, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

//         if( cmpunlambda < 0.0f )
//             throw MYRUNTIME_ERROR( preamb + 
//             "Wrong profile format at the end: Invalid statistical parameter lambda." );


        //reference
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameters." );

        if(( p = strstr( locbuffer, patstrSPREFR )) == NULL )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameters." );

        p += strlen( patstrSPREFR );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameters." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &refunK, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( refunK < 0.0f )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: Invalid reference parameter K." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameter lambda." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &refunlambda, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( refunlambda < 0.0f )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: Invalid reference parameter lambda." );


        //entropy, expected
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No entropy." );

        if(( p = strstr( locbuffer, patstrENTROPY )) == NULL )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No entropy." );

        p += strlen( patstrENTROPY );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No entropy." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &proentropy, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

//         if( proentropy < 0.0f )
//             throw MYRUNTIME_ERROR( preamb + 
//             "Wrong profile format at the end: Invalid entropy." );

        if(( p = strstr( locbuffer, patstrEXPECTED )) == NULL )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No expected value." );

        p += strlen( patstrEXPECTED );

        if( length <= size_t( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No expected value." );

        if(( emsg = read_float( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));


        pssm.SetRefLambda( refunlambda );
        pssm.SetRefK( refunK );
        pssm.SetLambda( cmpunlambda );
        pssm.SetEntropy( proentropy );
        pssm.SetK( cmpunK );
        pssm.SetExpectedScore( proexpected );
    }
// #endif

// //     pssm.Finalize();


    //footer
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( length < lenstrEND )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No ending." );

    for( n = 0; n < (int)lenstrEND; n++ )
        if( locbuffer[n] != patstrEND[n] )
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Invalid ending." );
}

}//namespace pmodel
