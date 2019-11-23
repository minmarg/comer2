/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "liblib/alpha.h"
#include "SUBSTABLE.h"
#include "TRANSPROBS.h"
#include "BMProtoProfile.h"

// -------------------------------------------------------------------------
// CLASS BMProtoProfile
//
float* BMProtoProfile::prexpnores_ = NULL;
const int BMProtoProfile::nSizeOfExpNObservs = 400;

// -------------------------------------------------------------------------
// Constructor
//
BMProtoProfile::BMProtoProfile( size_t reservation )
{
    Init();
    InitPrexpNoDistinctRes();
    Realloc( reservation );
}

// Constructor
//
BMProtoProfile::BMProtoProfile( const BMSequence& sequence )
{
    Init();
    InitPrexpNoDistinctRes();
    if( !sequence.GetCapacity())
        throw MYRUNTIME_ERROR( "BMProtoProfile::BMProtoProfile: Argument's zero capacity." );

    Realloc( sequence.GetCapacity());
    BMSequence::operator=( sequence );
}

// Default constructor
//
BMProtoProfile::BMProtoProfile()
{
    throw MYRUNTIME_ERROR( "BMProtoProfile::BMProtoProfile: Default construction prohibited." );
}

// -------------------------------------------------------------------------
// Destructor
//
BMProtoProfile::~BMProtoProfile()
{
    size_t  n;
    DestroyPrexpNoDistinctRes();

    if( indices_ ) {
        for( n = 0; n < length_; n++ )
            if( indices_[n] )
                free( indices_[n] );
        free( indices_ );
    }
    if( indicesLengths_ )       free( indicesLengths_ );
    if( indicesCapacities_ )    free( indicesCapacities_ );

    if( states_ )       free( states_ );
    if( extents_ )      free( extents_ );
    if( nodifrs_ )      free( nodifrs_ );
    if( MSextents_ )    free( MSextents_ );
    if( counts_ )       free( counts_ );

    FreeSqnWeightsAt();

    if( distribution_)  free( distribution_);
    if( matchWeights_)  free( matchWeights_);
    if( transWeights_)  free( transWeights_);
    if( gapWeights_  )  free( gapWeights_  );
    if( distinctHist_)  free( distinctHist_);
    if( MSdistinctHist_)free( MSdistinctHist_);
    if( targetFreqns_)  free( targetFreqns_);
    if( targetTranst_)  free( targetTranst_);
    if( rawPSSMatrix_)  free( rawPSSMatrix_);
    if( bppprob_ )      free( bppprob_ );
    if( noppps_ )       free( noppps_ );
    if( information_ )  free( information_ );
    if( noseqs_ )       free( noseqs_ );
    if( expnoseqs_ )    free( expnoseqs_ );
    if( expnosmid_ )    free( expnosmid_ );

    if( ppprobs_ ) {
        for( n = 0; n < length_; n++ )
            if( ppprobs_[n] )
                free( ppprobs_[n] );
        free( ppprobs_ );
    }
    if( pppndxs_ ) {
        for( n = 0; n < length_; n++ )
            if( pppndxs_[n] )
                free( pppndxs_[n] );
        free( pppndxs_ );
    }
}

// -------------------------------------------------------------------------
// FreeSqnWeightsAt: free sequence weights for all states
//
void BMProtoProfile::FreeSqnWeightsAt()
{
    int     st;
    size_t  n;
    for( st = 0; st < PS_NSTATES; st++ ) {
        if( sqnweights_[st]) {
            for( n = 0; n < length_; n++ )
                FreeSqnWeightsAt( n, st );
        }
        if( sqnweights_[st]) {
            free( sqnweights_[st]);
            sqnweights_[st] = NULL;
        }
        if( sqnweights_[st+PS_NSTATES]) {
            free( sqnweights_[st+PS_NSTATES]);
            sqnweights_[st+PS_NSTATES] = NULL;
        }
    }
}





// =========================================================================
// InitPrexpNoObservations: initialize expected number of distinct residues:
//     SUM( 1 - ( 1-backp[i]^n )), where n is the number of observations
//
void BMProtoProfile::InitPrexpNoDistinctRes( const float* backprobs )
{
    int noelems = GetSizeOfExpNoDistinctRes();
    int noeffres = NUMAA;

    if( noelems < 1 )
        return;

    DestroyPrexpNoDistinctRes();

    float pterm;//partial sum of probabilities
    int n, r;

    prexpnores_ = ( float* )malloc( noelems * sizeof(float));

    if( prexpnores_ == NULL )
        throw MYRUNTIME_ERROR( "BMProtoProfile::InitPrexpNoDistinctRes: Not enough memory." );

    memset( prexpnores_, 0, noelems * sizeof(float));

    prexpnores_[0] = 0.0f;

    for( n = 1; n < noelems; n++ ) {
        pterm = 0.0f;
        for( r = 0; r < noeffres; r++ )
            pterm += expf( n * logf( 1.0f - ( backprobs? backprobs[r]: STABLE.PROBABility(r)) ));
        prexpnores_[n] = noeffres - pterm;
    }
}

// DestroyPrexpNoObservations: destroy the structure of expected number of 
// distinct residues
//
void BMProtoProfile::DestroyPrexpNoDistinctRes()
{
    if( prexpnores_ ) {
        free( prexpnores_ );
        prexpnores_ = NULL;
    }
}

// GetExpNoObservations: get the expected number of observations
// corresponding to the average number of different observed residues
//
float BMProtoProfile::GetExpNoObservations( float avgnodistres )
{
    int     noelems = GetSizeOfExpNoDistinctRes();
    float   expnobs = 0.0f;
    int     n;

    if( prexpnores_ == NULL )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetExpNoObservations: Memory access error." );

    for( n = 1; n < noelems && prexpnores_[n] <= avgnodistres; n++ );
    expnobs = ( noelems <= n )? n: 
                n - (prexpnores_[n]-avgnodistres) / (prexpnores_[n]-prexpnores_[n-1]);
    return expnobs;
}





// =========================================================================
// Init: reset all values
//
void BMProtoProfile::Init()
{
    int st;
    BMSequence::Init();
    //
    indices_ = NULL;
    indicesLengths_ = NULL;
    indicesCapacities_= NULL;
    //
    memset( backprobs_, 0, sizeof(float) * NUMALPH );
    memset( postprobs_, 0, sizeof(float) * NUMALPH );
    states_ = NULL;
    extents_ = NULL;
    nodifrs_ = NULL;
    MSextents_ = NULL;
    counts_ = NULL;
    distribution_ = NULL;
    for( st = 0; st < PS_NSTATES; st++ ) {
        sqnweights_[st] = NULL;
        sqnweights_[st+PS_NSTATES] = NULL;
    }
    matchWeights_ = NULL;
    transWeights_ = NULL;
    gapWeights_ = NULL;
    distinctHist_ = NULL;
    MSdistinctHist_ = NULL;
    targetFreqns_ = NULL;
    targetTranst_ = NULL;
    rawPSSMatrix_ = NULL;
    information_ = NULL;
    bppprob_ = NULL;
    ppprobs_ = NULL;
    pppndxs_ = NULL;
    noppps_ = NULL;
    noseqs_      = NULL;
    expnoseqs_   = NULL;
    expnosmid_   = NULL;
}

// -------------------------------------------------------------------------
// InitRightExtents: initialize the right boundary of extents
//
void BMProtoProfile::InitRightExtents( size_t from, size_t to )
{
    size_t  n, loclength = GetCapacity();
    int st;

    if( extents_ == NULL || MSextents_ == NULL )
        throw MYRUNTIME_ERROR( "BMProtoProfile::InitRightExtents: Memory access error." );
    if( to )
        loclength = to;

    for( n = from; n < loclength; n++ )
        for( st = 0; st < PS_NSTATES; st++ )
            extents_[n][st][xRight] = (size_t)MYSIZE_MAX;
    for( n = from; n < loclength; n++ )
        MSextents_[n][xRight] = (size_t)MYSIZE_MAX;
}

// -------------------------------------------------------------------------
// assignment operator
//
BMProtoProfile& BMProtoProfile::operator=( const BMSequence& )
{
    throw MYRUNTIME_ERROR("BMProtoProfile::operator=: assignment is disallowed.");
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory
//
void BMProtoProfile::Realloc( size_t newcap )
{
    int st;
    if( newcap < capacity_ )
        return;

    if( capacity_ == 0 ) {
        indices_            = ( size_t** )malloc( sizeof( void* ) * newcap );
        indicesLengths_     = ( size_t* )malloc( sizeof( size_t ) * newcap );
        indicesCapacities_  = ( size_t* )malloc( sizeof( size_t ) * newcap );
        //
        states_ = ( int* )malloc( sizeof( int ) * newcap );
        extents_ = ( size_t(*)[PS_NSTATES][xCount])malloc( sizeof(size_t)*PS_NSTATES*xCount*newcap );
        nodifrs_ = ( float(*)[PS_NSTATES])malloc( sizeof(float)*PS_NSTATES*newcap );
        MSextents_ = ( size_t(*)[xCount] )malloc( sizeof( size_t ) * xCount * newcap );
        counts_  = ( size_t* )malloc( sizeof( size_t ) * newcap );

        for( st = 0; st < PS_NSTATES; st++ ) {
            sqnweights_[st] = ( float** )malloc( sizeof(float*) * newcap );
            sqnweights_[st+PS_NSTATES] = ( float** )malloc( sizeof(float*) * newcap );
        }

        distribution_ = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        matchWeights_ = ( float(*)[NUMALPH] )malloc( sizeof(float)*NUMALPH*newcap );
        transWeights_ = ( float(*)[P_NSTATES] )malloc( sizeof(float)*P_NSTATES*(newcap+1));
        gapWeights_   = ( float* )malloc( sizeof(float) * newcap );
        distinctHist_ = ( size_t(*)[PS_NSTATES][NUMALPH])malloc( sizeof(size_t)*PS_NSTATES*NUMALPH*(newcap+1));
        MSdistinctHist_ = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        targetFreqns_ = ( float(*)[NUMALPH] )malloc( sizeof(float)*NUMALPH*newcap );
        targetTranst_ = ( float(*)[P_NSTATES] )malloc( sizeof(float)*P_NSTATES*(newcap+1));
        rawPSSMatrix_ = ( float(*)[NUMALPH] )malloc( sizeof(float)*NUMALPH*newcap );
        information_  = ( float* )malloc( sizeof(float)*newcap );
        bppprob_     = ( float* )malloc( sizeof(float)*newcap );
        ppprobs_     = ( float** )malloc( sizeof(float*) * newcap );
        pppndxs_     = ( int** )malloc( sizeof( int*) * newcap );
        noppps_      = ( size_t* )malloc( sizeof( size_t ) * newcap );
        noseqs_      = ( size_t* )malloc( sizeof( size_t ) * newcap );
        expnoseqs_   = ( float* )malloc( sizeof(float) * newcap );
        expnosmid_ = ( float(*)[PS_NSTATES])malloc( sizeof(float)*PS_NSTATES*(newcap+1));

    } else {
        indices_ = ( size_t** )realloc( indices_, sizeof( void* ) * newcap );
        indicesLengths_ = ( size_t* )realloc( indicesLengths_, sizeof( size_t ) * newcap );
        indicesCapacities_ = ( size_t* )realloc( indicesCapacities_, sizeof( size_t ) * newcap );
        //
        states_ = ( int* )realloc( states_, sizeof( int ) * newcap );
        extents_ = ( size_t(*)[PS_NSTATES][xCount])realloc( extents_, sizeof(size_t)*PS_NSTATES*xCount*newcap );
        nodifrs_ = ( float(*)[PS_NSTATES])realloc( nodifrs_, sizeof(float)*PS_NSTATES*newcap );
        MSextents_ = ( size_t(*)[xCount] )realloc( MSextents_, sizeof(size_t)*xCount*newcap );
        counts_  = ( size_t* )realloc( counts_, sizeof(size_t)*newcap );

        for( st = 0; st < PS_NSTATES; st++ ) {
            sqnweights_[st] = (float**)realloc( sqnweights_[st], sizeof(float*)* newcap );
            sqnweights_[st+PS_NSTATES] = (float**)realloc( sqnweights_[st+PS_NSTATES], sizeof(float*)*newcap );
        }

        distribution_ = ( size_t(*)[NUMALPH] )realloc( distribution_, sizeof(size_t)*NUMALPH*newcap );
        matchWeights_ = ( float(*)[NUMALPH] )realloc( matchWeights_, sizeof(float)*NUMALPH*newcap );
        transWeights_ = ( float(*)[P_NSTATES] )realloc( transWeights_, sizeof(float)*P_NSTATES*(newcap+1));
        gapWeights_   = ( float* )realloc( gapWeights_, sizeof(float)*newcap );
        distinctHist_ = ( size_t(*)[PS_NSTATES][NUMALPH])realloc( distinctHist_, sizeof(size_t)*PS_NSTATES*NUMALPH*(newcap+1));
        MSdistinctHist_ = ( size_t(*)[NUMALPH] )realloc( MSdistinctHist_, sizeof(size_t)*NUMALPH*newcap );
        targetFreqns_ = ( float(*)[NUMALPH] )realloc( targetFreqns_, sizeof(float)*NUMALPH*newcap );
        targetTranst_ = ( float(*)[P_NSTATES] )realloc( targetTranst_, sizeof(float)*P_NSTATES*(newcap+1));
        rawPSSMatrix_ = ( float(*)[NUMALPH] )realloc( rawPSSMatrix_, sizeof(float)*NUMALPH*newcap );
        information_  = ( float* )realloc( information_, sizeof(float)*newcap );
        bppprob_     = ( float* )realloc( bppprob_, sizeof(float)*newcap );
        ppprobs_     = ( float** )realloc( ppprobs_, sizeof(float*)*newcap );
        pppndxs_     = ( int** )realloc( pppndxs_, sizeof(int*)*newcap );
        noppps_      = ( size_t* )realloc( noppps_, sizeof(size_t)*newcap );
        noseqs_      = ( size_t* )realloc( noseqs_, sizeof(size_t)*newcap );
        expnoseqs_   = ( float* )realloc( expnoseqs_, sizeof(float)*newcap );
        expnosmid_ = ( float(*)[PS_NSTATES])realloc( expnosmid_, sizeof(float)*PS_NSTATES*(newcap+1));
    }

    if( !indices_ || !indicesLengths_ || !indicesCapacities_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    if( !states_ || !extents_ || !nodifrs_ || !MSextents_ || !counts_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    if( !distribution_ || !matchWeights_ || !transWeights_ || !gapWeights_ || 
        !distinctHist_ || !MSdistinctHist_ || !targetFreqns_ || !targetTranst_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    if( !bppprob_ || !ppprobs_ || !pppndxs_ || !noppps_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    if( !information_ || !rawPSSMatrix_ || !noseqs_ || !expnoseqs_ || !expnosmid_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    for( st = 0; st < PS_NSTATES; st++ )
        if( !sqnweights_[st] || !sqnweights_[st+PS_NSTATES])
            throw MYRUNTIME_ERROR( "BMProtoProfile::Realloc: Not enough memory." );

    size_t**            locindices = indices_;
    size_t*             locindicesLengths = indicesLengths_;
    size_t*             locindicesCapacities = indicesCapacities_;
    //
    int*                locstates = states_;
    size_t           ( *locextents )[PS_NSTATES][xCount] = extents_;
    float            ( *locnodifrs )[PS_NSTATES] = nodifrs_;
    size_t           ( *locMSextents )[xCount] = MSextents_;
    size_t*             loccounts = counts_;
    float**             locsqnweights[TIMES2(PS_NSTATES)];
    size_t           ( *locdistribution )[NUMALPH] = distribution_;
    float            ( *locmatchWeights )[NUMALPH] = matchWeights_;
    float            ( *loctransWeights )[P_NSTATES] = transWeights_ + 1;
    float*              locgapWeights = gapWeights_;
    size_t           ( *locdistinctHist )[PS_NSTATES][NUMALPH] = distinctHist_ + 1;
    size_t           ( *locMSdistinctHist )[NUMALPH] = MSdistinctHist_;
    float            ( *loctargetFreqns )[NUMALPH] = targetFreqns_;
    float            ( *loctargetTranst )[P_NSTATES] = targetTranst_ + 1;
    float            ( *locrawPSSMatrix )[NUMALPH] = rawPSSMatrix_;
    float*              locinformation             = information_;
    float*              locbppprob                 = bppprob_;
    float**             locppprobs                 = ppprobs_;
    int**               locpppndxs                 = pppndxs_;
    size_t*             locnoppps                  = noppps_;
    size_t*             locnoseqs                  = noseqs_;
    float*              locexpnoseqs               = expnoseqs_;
    float            ( *locexpnosmid )[PS_NSTATES] = expnosmid_ + 1;

    for( st = 0; st < PS_NSTATES; st++ ) {
        locsqnweights[st] = sqnweights_[st];
        locsqnweights[st+PS_NSTATES] = sqnweights_[st+PS_NSTATES];
    }

    if( capacity_ != 0 ) {
        locindices += capacity_;
        locindicesLengths += capacity_;
        locindicesCapacities += capacity_;
        //
        locstates += capacity_;
        locextents += capacity_;        //!
        locnodifrs += capacity_;
        locMSextents += capacity_;
        loccounts += capacity_;

        for( st = 0; st < PS_NSTATES; st++ ) {
            locsqnweights[st] += capacity_;
            locsqnweights[st+PS_NSTATES] += capacity_;
        }
        locdistribution += capacity_;   //!
        locmatchWeights += capacity_;
        loctransWeights += capacity_;
        locgapWeights   += capacity_;
        locdistinctHist += capacity_;
        locMSdistinctHist += capacity_;
        loctargetFreqns += capacity_;
        loctargetTranst += capacity_;
        locrawPSSMatrix += capacity_;
        locinformation  += capacity_;
        locbppprob      += capacity_;
        locppprobs      += capacity_;
        locpppndxs      += capacity_;
        locnoppps       += capacity_;
        locnoseqs       += capacity_;
        locexpnoseqs    += capacity_;
        locexpnosmid    += capacity_;
    }
    else {
        //beginning transition weights and frequencies
        memset( transWeights_, 0, sizeof( float ) * P_NSTATES );
        memset( targetTranst_, 0, sizeof( float ) * P_NSTATES );
        memset( distinctHist_, 0, sizeof( size_t ) * PS_NSTATES * NUMALPH );
        memset( expnosmid_, 0, sizeof( float ) * PS_NSTATES );
    }


    memset( locindices,             0, sizeof( size_t* ) * ( newcap - capacity_ ));
    memset( locindicesLengths,      0, sizeof( size_t  ) * ( newcap - capacity_ ));
    memset( locindicesCapacities,   0, sizeof( size_t  ) * ( newcap - capacity_ ));
    //
    memset( locstates,  0, sizeof( int ) * ( newcap - capacity_ ));
    memset( locextents, 0, sizeof( size_t ) * ( newcap - capacity_ ) * PS_NSTATES * xCount );
    memset( locnodifrs, 0, sizeof( float ) * ( newcap - capacity_ ) * PS_NSTATES );
    memset( locMSextents, 0, sizeof( size_t ) * ( newcap - capacity_ ) * xCount );
    memset( loccounts,  0, sizeof( size_t ) * ( newcap - capacity_ ));
    InitRightExtents( capacity_, newcap );
    // 
    for( st = 0; st < PS_NSTATES; st++ ) {
        memset( locsqnweights[st], 0, sizeof( float* ) * ( newcap - capacity_ ));
        memset( locsqnweights[st+PS_NSTATES], 0, sizeof( float* ) * ( newcap - capacity_ ));
    }
    //
    memset( locdistribution, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );
    memset( locmatchWeights, 0, sizeof( float ) * ( newcap - capacity_ ) * NUMALPH );
    memset( loctransWeights, 0, sizeof( float ) * ( newcap - capacity_ ) * P_NSTATES );
    memset( locgapWeights,   0, sizeof( float ) * ( newcap - capacity_ ));

    memset( locdistinctHist, 0, sizeof( size_t ) * ( newcap - capacity_ ) * PS_NSTATES * NUMALPH );
    memset( locMSdistinctHist, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );

    memset( loctargetFreqns, 0, sizeof( float ) * ( newcap - capacity_ ) * NUMALPH );
    memset( loctargetTranst, 0, sizeof( float ) * ( newcap - capacity_ ) * P_NSTATES );
    memset( locrawPSSMatrix, 0, sizeof( float ) * ( newcap - capacity_ ) * NUMALPH );
    //
    memset( locinformation,  0, sizeof( float ) * ( newcap - capacity_ ));
    memset( locbppprob,      0, sizeof( float ) * ( newcap - capacity_ ));
    memset( locppprobs,      0, sizeof( float*) * ( newcap - capacity_ ));
    memset( locpppndxs,      0, sizeof( int*) * ( newcap - capacity_ ));
    memset( locnoppps,       0, sizeof( size_t ) * ( newcap - capacity_ ));
    memset( locnoseqs,       0, sizeof( size_t ) * ( newcap - capacity_ ));
    memset( locexpnoseqs,    0, sizeof( float ) * ( newcap - capacity_ ));
    memset( locexpnosmid,    0, sizeof( float ) * ( newcap - capacity_ ) * PS_NSTATES );

    BMSequence::Realloc( newcap );
}

// -------------------------------------------------------------------------
// PushIndexAt: insert the indices of sequences not participating in
//     extent at position n
//
void BMProtoProfile::PushIndexAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::PushIndexAt: Memory access error." );
#endif

    static const size_t scnminsize = 128;

    if( indicesCapacities_[n] <= indicesLengths_[n] ) {
        size_t newindcap = TIMES2( indicesCapacities_[n]);
        if( newindcap <= indicesLengths_[n] )
            newindcap = indicesLengths_[n] + 1;
        if( newindcap <= scnminsize )
            newindcap = scnminsize;
        ReallocIndices( newindcap, n );
    }

    indices_[n][ indicesLengths_[n] ] = value;

    indicesLengths_[n]++;
}

// -------------------------------------------------------------------------
// ReallocIndices: reallocate memory associated with the vector of sequence 
// indices
//
void BMProtoProfile::ReallocIndices( size_t newcap, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || !indicesCapacities_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::ReallocIndices: Memory access error." );
#endif

    if( newcap <= indicesCapacities_[n])
        return;

    if( indicesCapacities_[n] == 0 ) {
        indices_[n] = ( size_t* )malloc( sizeof(size_t) * newcap );
    } else {
        indices_[n] = ( size_t* )realloc( indices_[n], sizeof(size_t) * newcap );
    }

    if( !indices_[n] )
        throw MYRUNTIME_ERROR( "BMProtoProfile::ReallocIndices: Not enough memory." );

    size_t* indarr = indices_[n];

    if( indicesCapacities_[n] != 0 ) {
        indarr += indicesCapacities_[n];
    }

    memset( indarr, 0, sizeof(size_t) * ( newcap - indicesCapacities_[n] ));

    indicesCapacities_[n] = newcap;
}

// -------------------------------------------------------------------------
// clear: clear all buffers associated with this object
//
void BMProtoProfile::clear()
{
    size_t n;
    int st;

    for( n = 0; n < length_; n++ )
        if( indices_[n] )
            free( indices_[n] );

    memset( indices_,           0, sizeof( size_t* ) * capacity_ );
    memset( indicesLengths_,    0, sizeof( size_t  ) * capacity_ );
    memset( indicesCapacities_, 0, sizeof( size_t  ) * capacity_ );
    //
    memset( backprobs_, 0, sizeof( float ) * NUMALPH );
    memset( postprobs_, 0, sizeof( float ) * NUMALPH );
    memset( states_,  0, sizeof( int ) * capacity_ );
    memset( extents_, 0, sizeof( size_t ) * capacity_ * PS_NSTATES * xCount );
    memset( nodifrs_, 0, sizeof( float ) * capacity_ * PS_NSTATES );
    memset( MSextents_, 0, sizeof( size_t ) * capacity_ * xCount );
    memset( counts_,  0, sizeof( size_t ) * capacity_ );
    InitRightExtents();

    if( capacity_ ) {
        for( st = 0; st < PS_NSTATES; st++ ) {
            for( n = 0; n < length_; n++ )
                FreeSqnWeightsAt( n, st );
            memset( sqnweights_[st], 0, sizeof( float* ) * capacity_ );
            memset( sqnweights_[st+PS_NSTATES], 0, sizeof( float* ) * capacity_ );
        }
        for( n = 0; n < length_; n++ ) {
            if( ppprobs_ && ppprobs_[n]) free( ppprobs_[n]);
            if( pppndxs_ && pppndxs_[n]) free( pppndxs_[n]);
        }
        memset( distribution_, 0, sizeof( size_t ) * capacity_ * NUMALPH );
        memset( matchWeights_, 0, sizeof( float ) * capacity_ * NUMALPH );
        memset( transWeights_, 0, sizeof( float ) * ( capacity_+1 ) * P_NSTATES );
        memset( gapWeights_,  0, sizeof( float ) * capacity_ );
        memset( distinctHist_, 0, sizeof( size_t ) * ( capacity_+1 )* PS_NSTATES*NUMALPH );
        memset( MSdistinctHist_, 0, sizeof( size_t ) * capacity_ * NUMALPH );
        memset( targetFreqns_, 0, sizeof( float ) * capacity_ * NUMALPH );
        memset( targetTranst_, 0, sizeof( float ) * ( capacity_+1 ) * P_NSTATES );
        memset( rawPSSMatrix_, 0, sizeof( float ) * capacity_ * NUMALPH );
        memset( information_, 0, sizeof( float ) * capacity_ );
        memset( bppprob_,     0, sizeof( float ) * capacity_ );
        memset( ppprobs_,     0, sizeof( float*) * capacity_ );
        memset( pppndxs_,     0, sizeof( int*) * capacity_ );
        memset( noppps_,      0, sizeof( size_t ) * capacity_ );
        memset( noseqs_,      0, sizeof( size_t ) * capacity_ );
        memset( expnoseqs_,   0, sizeof( float ) * capacity_ );
        memset( expnosmid_,   0, sizeof( float ) * ( capacity_+1 ) * PS_NSTATES );
    }

    BMSequence::clear();
}

// -------------------------------------------------------------------------
// PrintMatchWeights: output match weights, which correspond to observed
//     weighted frequencies
//
void BMProtoProfile::PrintMatchWeights( FILE* fp )
{
    size_t l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c%c", 32, DehashCode( r ) );

    for( size_t p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetMatchWeightsAt(r,p)));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintTransWeights: output transition weights
//
void BMProtoProfile::PrintTransWeights( FILE* fp )
{
    int n;
    size_t l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, " %3s", gTPTRANS_NAMES[n]);

    fprintf( fp, "%s%5zu %c   ", NL, l, 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4d", (int)rintf( 100.0f * GetTransWeightsBeg(n)));

    for( size_t p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetTransWeightsAt(n,p)));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintTransWeights: output target transition weights
//
void BMProtoProfile::PrintTargetTranst( FILE* fp )
{
    int n;
    size_t l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, " %3s", gTPTRANS_NAMES[n]);

    fprintf( fp, "%s%5zu %c   ", NL, l, 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4d", (int)rintf( 100.0f * GetTargetTranstBeg(n)));

    for( size_t p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetTargetTranstAt(n,p)));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintMIDExpNoObservations: output weighted MID state observations
//
void BMProtoProfile::PrintMIDExpNoObservations( FILE* fp )
{
    int n;
    size_t l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < PS_NSTATES; n++ )
        fprintf( fp, " %11s", gTPSTATES_NAMES[n]);

    fprintf( fp, "%s%5zu %c   ", NL, l, 32 );

    for( n = 0; n < PS_NSTATES; n++ )
        fprintf( fp, "%12g", GetMIDExpNoObservationsAt( -1, n ));

    for( size_t p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( n = 0; n < PS_NSTATES; n++ )
            fprintf( fp, "%12g", GetMIDExpNoObservationsAt((int)p, n ));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintPSSMatrix: output PSSM matrix
//
void BMProtoProfile::PrintPSSMatrix( FILE* fp )
{
    size_t l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( size_t p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c  ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%3d", (int)rintf( GetPSSMEntryAt(r,p)));
//             fprintf( fp, "%f ", GetTargetFreqnsAt( r, p ));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintSuppressedPSSMandWeights: output PSSM matrix together with match 
// weights
//
void BMProtoProfile::PrintSuppressedPSSMandWeights( FILE* fp )
{
    int             n;
    size_t          p, l = 0;
    const size_t    effnores = NUMAA;
    unsigned char   r = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%22c Position-specific scoring matrix "
                "%39c Weighted observed frequencies "
                "%18c Target transition frequencies %7c Information%s", 32, 32, 32, 32, NL );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < effnores; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( r = 0; r < effnores; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4s", gTPTRANS_NAMES[n]);

    for( p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( r = 0; r < effnores; r++ )
            fprintf( fp, "%2d ", (int)rintf( GetPSSMEntryAt(r,p)));

        for( r = 0; r < effnores; r++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetMatchWeightsAt(r,p)));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetTargetTranstAt(n,(ssize_t)p)));

        fprintf( fp, " %5.2f", GetInformationAt(p));
    }
    fprintf( fp, "%s", NL );
}

// -------------------------------------------------------------------------
// PrintProfile: output profile information
//
void BMProtoProfile::PrintProfile( FILE* fp )
{
    int             n;
    size_t          p, l = 0;
    const size_t    res_count = NUMALPH - 1;
    unsigned char   r = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%22c Position-specific scoring matrix "
                "%39c Weighted observed frequencies "
                "%18c Target transition frequencies %7c Gap weights %c Information%s", 
                32, 32, 32, 32, 32, NL );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < res_count; r++ )
        fprintf( fp, "%3c", DehashCode(r));

    for( r = 0; r < res_count; r++ )
        fprintf( fp, "%4c", DehashCode(r));

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4s", gTPTRANS_NAMES[n]);

    for( p = 0; p < GetSize(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */GetResidueAt(p) == GAP )
            continue;

        fprintf( fp, "%s%5zu %c   ", NL, ++l, DehashCode( GetResidueAt(p)));

        for( r = 0; r < res_count; r++ )
            fprintf( fp, "%2d ", (int)rintf( GetPSSMEntryAt(r,p)));

        for( r = 0; r < res_count; r++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetMatchWeightsAt(r,p)));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", (int)rintf( 100.0f * GetTargetTranstAt(n,(ssize_t)p)));

        fprintf( fp, " %6d", (int)rintf( 100.0f * GetGapWeightsAt(p)));

        fprintf( fp, " %13.2f", GetInformationAt(p));
    }
    fprintf( fp, "%s", NL );
}
