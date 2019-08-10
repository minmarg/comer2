/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "liblib/alpha.h"
#include "libseg/segdata.h"
#include "libseg/SEGAbstract.h"
#include "SUBSTABLE.h"
#include "PMTransModel.h"
#include "PMProfileModel.h"
#include "SEGProfile.h"

// /////////////////////////////////////////////////////////////////////////
// CLASS SEGProfile
//
//default SEG window length
const size_t SEGProfile::scszDefSegProWinLength = 12;
//default SEG low entropy threshold
const float  SEGProfile::scfpDefSegProLowEntropy = 2.2f;
//default SEG high entropy threshold
const float  SEGProfile::scfpDefSegProHighEntropy = 2.5f;
//default SEG maximum difference between positions
const size_t SEGProfile::scszDefSegProMaxDifference = 100;
//default SEG (Euclidean) distance threshold for profile vectors of probabilities
const float  SEGProfile::scfpDefSegProVecDistance = 12.96f;

// alphabet size for profiles (number of possible columns) --
// Do not specify!
const size_t SEGProfile::sc_sizeproalphabet_ = 0;

float SEGProfile::distthld_ = SEGProfile::scfpDefSegProVecDistance;
float SEGProfile::distsquared_ = SQUARE(SEGProfile::scfpDefSegProVecDistance);

// -------------------------------------------------------------------------
// Constructor:
//
SEGProfile::SEGProfile(
    const pmodel::PMProfileModel& prom,
    size_t  winlen,
    float   lowent,
    float   highent,
    size_t  maxdiff )
:
    SEGAbstract(
        &SEGProfile::ProValidator,
        &SEGProfile::ProVerifier,
        &SEGProfile::ProComparer,
        &SEGProfile::ProEquality,
        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize()/*alphabet size*/
    ),
    addresses_( NULL ),
    length_( 0 )
{
    AllocateAddresses( prom.GetSize());
    Translate( prom );

    if( addresses_ == NULL )
        throw MYRUNTIME_ERROR( "SEGProfile::SEGProfile: Not enough memory." );

    SetRunAddress((void*)GetAddresses(), GetLocalLength());
}

// Alternative constructor: Vectors are constructed as differences between
//     corresponding frequency/probability vectors
//
SEGProfile::SEGProfile(
    const pmodel::PMProfileModel& prom1,
    const pmodel::PMProfileModel& prom2,
    size_t  winlen,
    float   lowent,
    float   highent,
    size_t  maxdiff )
:
    SEGAbstract(
        &SEGProfile::ProValidator,
        &SEGProfile::ProVerifier,
        &SEGProfile::ProComparer,
        &SEGProfile::ProEquality,
        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize()/*alphabet size*/
    ),
    addresses_( NULL ),
    length_( 0 )
{
    if( prom1.GetSize() != prom2.GetSize())
        throw MYRUNTIME_ERROR( "SEGProfile::SEGProfile: Inconsistent profile lengths." );

    AllocateAddresses( prom1.GetSize());
    Translate( prom1, prom2 );

    if( addresses_ == NULL )
        throw MYRUNTIME_ERROR( "SEGProfile::SEGProfile: Not enough memory." );

    SetRunAddress(( void* )GetAddresses(), GetLocalLength());
    SetHighCSearch();
}

// destructor:
//
SEGProfile::~SEGProfile()
{
    Destroy();
}

// -------------------------------------------------------------------------
// AllocateAddresses: allocate memory required by addresses
//
void SEGProfile::AllocateAddresses( size_t newlen )
{
    void* p = NULL;

    Destroy();

    addresses_ = ( provval_t** )malloc( sizeof(provval_t*) * newlen );

    if( addresses_ == NULL )
        throw MYRUNTIME_ERROR( "SEGProfile::AllocateAddresses: Not enough memory." );

    memset( addresses_, 0, sizeof(provval_t*) * newlen );

    length_ = newlen;

    for( size_t n = 0; n < length_; n++ ) {
        p = addresses_[n] = ( provval_t* )malloc( sizeof(provval_t) * GetVectorSize());

        if( p == NULL )
            throw MYRUNTIME_ERROR( "SEGProfile::AllocateAddresses: Not enough memory." );

        memset( p, 0, sizeof(provval_t) * GetVectorSize());
    }
}

// -------------------------------------------------------------------------
// Destroy: destroy allocated addresses
//
void SEGProfile::Destroy()
{
    if( addresses_ ) {
        for( size_t n = 0; n < length_; n++ )
            if( addresses_[n])
                free( addresses_[n] );
        free( addresses_ );
        addresses_ = NULL;
        length_ = 0;
    }
}

// -------------------------------------------------------------------------
// Translate: translate profile columns into a vector of addresses 
// suitable for abstract SEG
//
void SEGProfile::Translate( const pmodel::PMProfileModel& prom )
{
    if( GetLocalLength() < (size_t)prom.GetSize())
        throw MYRUNTIME_ERROR( "SEGProfile::Translate: Memory access error." );

    for( int n = 0; n < prom.GetSize(); n++ )
    {
        SetResidueAt( n, (provval_t)prom.GetResidueAt(n));

        for( size_t r = 0; r < No_frequencies; r++ )
            SetFrequencyAt((size_t)n, r, (provval_t)rintf( FREQUENCY_SUM * prom.GetObsFreqsAt(n,(int)r)) );

        SetFrequencyWeightAt( n, (provval_t)rintf( prom.GetFrequencyWeightAt(n)));
    }
}

// -------------------------------------------------------------------------
// Translate: translate differences of frequency/probability vectors into a
// vector of addresses suitable for abstract SEG
//
void SEGProfile::Translate( const pmodel::PMProfileModel& prom1, const pmodel::PMProfileModel& prom2 )
{
    if( GetLocalLength() < (size_t)prom1.GetSize() || prom1.GetSize() != prom2.GetSize())
        throw MYRUNTIME_ERROR( "SEGProfile::Translate: Inconsistent profile lengths." );

    for( int n = 0; n < prom1.GetSize(); n++ )
    {
        SetResidueAt( n, (provval_t)prom1.GetResidueAt(n));

        for( size_t r = 0; r < No_frequencies; r++ )
            SetFrequencyAt( n, r, 
            (provval_t)rintf( FREQUENCY_SUM * (prom1.GetObsFreqsAt(n,(int)r)-prom2.GetObsFreqsAt(n,(int)r))) );

        SetFrequencyWeightAt( n, (provval_t)0 );
    }
}



// -------------------------------------------------------------------------
// ProEquality: return true if two vectors are supposed to be equal, 
// false otherwise
//
bool SEGProfile::ProEquality( const void* pone, const void* panother )
{
#ifdef __DEBUG__
    if( !pone || !panother )
        throw MYRUNTIME_ERROR( "SEGProfile::ProEquality: Memory access error." );
#endif

    if( pone == panother )
        return true;

    float sqdist = 0.0f;//squared distance >0
    const provval_t* vectone = *(const provval_t**)pone;
    const provval_t* vectanother = *(const provval_t**)panother;
    int lastdiff;//last difference between vector elements

    //compute squared euclidean distance
    for( size_t r = 0; r < No_frequencies; r++ ) {
        lastdiff = (int)GetFrequency(vectone,r) - (int)GetFrequency(vectanother,r);
        sqdist += SEG::PRESQUARES.GetValueOf( abs( lastdiff ));
    }

    //take frequency weights into account as well
//     lastdiff = ( int )GetFrequencyWeight( vectone ) - ( int )GetFrequencyWeight( vectanother );
//     sqdist += SEG::PRESQUARES.GetValueOf( abs( lastdiff ));

    if( sqdist < GetDistance2())
        return true;

    return false;
}

// -------------------------------------------------------------------------
// ProComparer: returns 0 if serial number of two vectors coinsides,
// positive value if the number of the first is greater, negative otherwise
//
int SEGProfile::ProComparer( const void* pone, const void* panother )
{
#ifdef __DEBUG__
    if( !pone || !panother )
        throw MYRUNTIME_ERROR( "SEGProfile::ProComparer: Memory access error." );
#endif

    return ( int )(( size_t )pone - ( size_t )panother );
//     return memcmp(  vectone + start_of_frequencies,
//                     vectanother + start_of_frequencies,
//                 ( No_frequencies + 1 ) * sizeof( provval_t ));
}

// -------------------------------------------------------------------------
// MaskSeggedPositions: mask profile positions identified by the SEG 
// algorithm
//
void SEGProfile::MaskSeggedPositions( pmodel::PMProfileModel& prom, pmodel::PMTransModel& gaps ) const
{
    size_t  n = 0;
    size_t  p = 0;
    size_t  r = 0;
    size_t  left, right;
    const SEG::Segments& loc_segments = GetSegments();

    char    maskX = X;
    float   freqweight =  0.0f;
    float   information = 0.0f;
    float   expMIDs[PS_NSTATES] = {1.0f,0.0f,0.0f};
    float   obsvalues[pmodel::PMProfileModelBase::PVDIM];
    float   trgvalues[pmodel::PMProfileModelBase::PVDIM];
    float   trnprobs[P_NSTATES];
    const int effnos = NUMAA;// effective number of residues

    for( r = 0; r < effnos; r++ ) {
        obsvalues[r] = STABLE.PROBABility((int)r);
        trgvalues[r] = STABLE.PROBABility((int)r);
    }

    for( ; r < pmodel::PMProfileModelBase::PVDIM; r++ ) {
        obsvalues[r] = 0.0f;
        obsvalues[r] = 0.0f;
    }

    memset( trnprobs, 0, sizeof(float) * P_NSTATES );
    if( TRANSPROBS.GetPriors())
        memcpy( trnprobs, *TRANSPROBS.GetPriors(), sizeof(float) * P_NSTATES );

    for( n = 0; n < loc_segments.GetSize(); n++ ) {
        left = loc_segments.GetLeftAt( n );
        right = loc_segments.GetRightAt( n );

        for( p = left; p <= right; p++ ) {
            prom.PushAt( trgvalues, obsvalues, maskX, freqweight, information, expMIDs, (int)p );
            gaps.PushAt( &trnprobs, (int)p );
        }
    }
}

// -------------------------------------------------------------------------
// PrintSequence: print formatted profile sequence 
//
void SEGProfile::PrintSequence( FILE* fp, size_t width )
{
    SEG::SEGAbstract::PrintSequence( fp, &SEGProfile::GetProResidue, width );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: print formatted profile sequence with segments
// found by the algorithm masked with Xs
//
void SEGProfile::PrintSeggedSequence( FILE* fp, size_t width )
{
    SEG::SEGAbstract::PrintSeggedSequence( fp, &SEGProfile::GetProResidue, width );
}

