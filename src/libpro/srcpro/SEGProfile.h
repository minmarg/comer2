/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SEGProfile_h__
#define __SEGProfile_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include "liblib/alpha.h"
#include "PMTransModel.h"
#include "PMProfileModel.h"
#include "libseg/SEGAbstract.h"

typedef ssize_t provval_t;

// _________________________________________________________________________
// Class SEGProfile
//
// implementation of the SEG algorithm applied to profiles
//
class SEGProfile: public SEG::SEGAbstract
{
    enum {
        No_frequencies          = pmodel::PMProfileModelBase::PVDIM
    };
    enum {
        start_of_residue        = 0,
        start_of_frequencies    = start_of_residue + 1,
        start_of_weight         = start_of_frequencies + No_frequencies,
        Overall_size            = start_of_weight + 1
    };

public:
    SEGProfile(
        const pmodel::PMProfileModel& prom,
        size_t  winlen = scszDefSegProWinLength,
        float   lowent = scfpDefSegProLowEntropy,
        float   highent = scfpDefSegProHighEntropy,
        size_t  maxdiff = scszDefSegProMaxDifference
    );

    SEGProfile(
        const pmodel::PMProfileModel& prom1,
        const pmodel::PMProfileModel& prom2,
        size_t  winlen = scszDefSegProWinLength,
        float   lowent = scfpDefSegProLowEntropy,
        float   highent = scfpDefSegProHighEntropy,
        size_t  maxdiff = scszDefSegProMaxDifference
    );

    virtual ~SEGProfile();

    static float    GetDistance() { return distthld_; }
    static float    GetDistance2() { return distsquared_; }

    static void     SetDistance( float value ) { distthld_ = value; distsquared_ = SQUARE(distthld_); }

    void            MaskSeggedPositions( pmodel::PMProfileModel&, pmodel::PMTransModel& ) const;

    void            PrintSequence( FILE*, size_t width = SEGAbstract::scszDefSegSeqPrintWidth );
    void            PrintSeggedSequence( FILE*, size_t width = SEGAbstract::scszDefSegSeqPrintWidth );

protected:
    static bool     ProValidator( const void* );
    static bool     ProVerifier( const void* );
    static int      ProComparer( const void* , const void* );
    static bool     ProEquality( const void* , const void* );

protected:
    size_t          GetLocalLength() const  { return length_; }
    provval_t**     GetAddresses() const    { return addresses_; }

    provval_t       GetResidueAt( size_t pos ) const;
    void            SetResidueAt( size_t pos, provval_t value );

    provval_t       GetFrequencyAt( size_t pos, size_t r ) const;
    void            SetFrequencyAt( size_t pos, size_t r, provval_t value );

    provval_t       GetFrequencyWeightAt( size_t pos ) const;
    void            SetFrequencyWeightAt( size_t pos, provval_t value );

    static provval_t    GetResidue( const provval_t* );
    static provval_t    GetFrequency( const provval_t*, size_t r );
    static provval_t    GetFrequencyWeight( const provval_t* );

    static size_t   GetProResidue( const void* pvectaddr );

    void            AllocateAddresses( size_t newlen );
    void            Destroy();


    static size_t   GetSeqAlphabetSize()    { return sc_sizeproalphabet_; }
    static size_t   GetVectorSize()         { return Overall_size; }

    void            Translate( const pmodel::PMProfileModel& prom );
    void            Translate( const pmodel::PMProfileModel& prom1, const pmodel::PMProfileModel& prom2 );

public:
    static const size_t scszDefSegProWinLength;//default SEG window length
    static const float  scfpDefSegProLowEntropy;//default SEG low entropy threshold
    static const float  scfpDefSegProHighEntropy;//default SEG high entropy threshold
    static const size_t scszDefSegProMaxDifference;//default SEG maximum difference between positions
    static const float  scfpDefSegProVecDistance;//default SEG distance threshold for profile vectors of probabilities

private:
    provval_t**     addresses_;
    size_t          length_;

    static float    distthld_/*equalitydist_*/;//threshold of Euclidean distance between two profile vectors
    static float    distsquared_;//distthld_ squared

    static const size_t sc_sizeproalphabet_;
};

////////////////////////////////////////////////////////////////////////////
// Class SEGProfile INLINES
//
// ProValidator: get the indication of whether the vector is valid
//
inline
bool SEGProfile::ProValidator( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw MYRUNTIME_ERROR( "SEGProfile::ProValidator: Memory access error." );
#endif
    provval_t res = GetResidue(*( const provval_t**)pvectaddr);
    return res != GAP;
}

// -------------------------------------------------------------------------
// ProVerifier: get the indication of whether the vector is to be considered
//
inline
bool SEGProfile::ProVerifier( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw MYRUNTIME_ERROR( "SEGProfile::ProVerifier: Memory access error." );
#endif
    provval_t res = GetResidue(*( const provval_t**)pvectaddr);
    return res != GAP && res != X;
}

// =========================================================================
// GetProResidue: extract residue at the given profile vector address
//
inline
size_t SEGProfile::GetProResidue( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw MYRUNTIME_ERROR( "SEGProfile::GetProResidue: Memory access error." );
#endif
    return (size_t)GetResidue(*(const provval_t**)pvectaddr);
}

// GetResidue: get the residue present in the given profile vector 
//
inline
provval_t SEGProfile::GetResidue( const provval_t* vectaddr )
{
#ifdef __DEBUG__
    if( !vectaddr )
        throw MYRUNTIME_ERROR( "SEGProfile::GetResidue: Memory access error." );
#endif
    return vectaddr[start_of_residue];
}

// GetFrequency: get frequency/probability r present in the given profile
// vector
//
inline
provval_t SEGProfile::GetFrequency( const provval_t* vectaddr, size_t r )
{
#ifdef __DEBUG__
    if( !vectaddr || GetVectorSize() <= r + start_of_frequencies )
        throw MYRUNTIME_ERROR( "SEGProfile::GetFrequency: Memory access error." );
#endif
    return vectaddr[ r + start_of_frequencies ];
}

// GetFrequencyWeight: get the frequency weight present in the given
// profile vector
//
inline
provval_t SEGProfile::GetFrequencyWeight( const provval_t* vectaddr )
{
#ifdef __DEBUG__
    if( !vectaddr )
        throw MYRUNTIME_ERROR( "SEGProfile::GetFrequencyWeight: Memory access error." );
#endif
    return vectaddr[start_of_weight];
}

// =========================================================================
// GetResidueAt: get the residue at the specified position
//
inline
provval_t SEGProfile::GetResidueAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::GetResidueAt: Memory access error." );
#endif
    return addresses_[pos][start_of_residue];
}

// -------------------------------------------------------------------------
// SetResidueAt: set a residue at the specified profile position
//
inline
void SEGProfile::SetResidueAt( size_t pos, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::SetResidueAt: Memory access error." );
#endif
    addresses_[pos][start_of_residue] = value;
}


// -------------------------------------------------------------------------
// GetFrequencyAt: get frequency r at the specified profile position
//
inline
provval_t SEGProfile::GetFrequencyAt( size_t pos, size_t r ) const
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::GetFrequencyAt: Memory access error." );

    if( !addresses_[pos] || GetVectorSize() <= r + start_of_frequencies )
        throw MYRUNTIME_ERROR( "SEGProfile::GetFrequencyAt: Memory access error." );
#endif
    return addresses_[pos][ r + start_of_frequencies ];
}

// -------------------------------------------------------------------------
// SetFrequencyAt: set a frequency for residue r at the specified profile
// position
//
inline
void SEGProfile::SetFrequencyAt( size_t pos, size_t r, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::SetFrequencyAt: Memory access error." );

    if( !addresses_[pos] || GetVectorSize() <= r + start_of_frequencies )
        throw MYRUNTIME_ERROR( "SEGProfile::SetFrequencyAt: Memory access error." );
#endif
    addresses_[pos][ r + start_of_frequencies ] = value;
}


// -------------------------------------------------------------------------
// GetFrequencyWeightAt: get the frequency weight at the given position
//
inline
provval_t SEGProfile::GetFrequencyWeightAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::GetFrequencyWeightAt: Memory access error." );
#endif
    return addresses_[pos][start_of_weight];
}

// -------------------------------------------------------------------------
// SetFrequencyWeightAt: set a frequency weight for position pos
//
inline
void SEGProfile::SetFrequencyWeightAt( size_t pos, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGProfile::SetFrequencyWeightAt: Memory access error." );
#endif
    addresses_[pos][start_of_weight] = value;
}

#endif//__SEGProfile_h__
