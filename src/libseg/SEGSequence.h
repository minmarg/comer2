/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SEGSequence_h__
#define __SEGSequence_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include "liblib/alpha.h"
#include "SEGAbstract.h"

namespace SEG {

// _________________________________________________________________________
// Class SEGSequence
// for application of the SEG algorithm to sequences
//
class SEGSequence: public SEGAbstract
{
public:
    SEGSequence(
        const unsigned char* residues,
        size_t          seqlength,
        bool            hashed = false,
        size_t          winlen = scszDefSegSeqWinLength,
        float           lowent = scfpDefSegSeqLowEntropy,
        float           highent = scfpDefSegSeqHighEntropy,
        size_t          maxdiff = scszDefSegSeqMaxDifference
    );

    virtual ~SEGSequence();

    void            PrintSequence( FILE*, size_t width = SEGAbstract::scszDefSegSeqPrintWidth );
    void            PrintSeggedSequence( FILE*, size_t width = SEGAbstract::scszDefSegSeqPrintWidth );

    static bool     SeqValidator( const void* );
    static bool     SeqVerifier( const void* );
    static int      SeqComparer( const void* , const void* );

protected:
    size_t          GetLocalLength() const { return length_; }

    void            Translate( const unsigned char*, size_t len );

    size_t          GetAddressAt( size_t pos ) const;
    void            SetAddressAt( size_t pos, size_t value );

    static size_t*  AllocateAddresses( size_t newlen );
    static void     Destroy( size_t* );

    static size_t   GetSeqAlphabetSize() { return sc_sizealphabet_; }
    bool            GetHashed() const { return hashed_; }

    static size_t   GetSeqResidue( const void* );

public:
    static const size_t scszDefSegSeqWinLength;//default SEG window length
    static const float  scfpDefSegSeqLowEntropy;//default SEG low entropy threshold
    static const float  scfpDefSegSeqHighEntropy;//default SEG high entropy threshold
    static const size_t scszDefSegSeqMaxDifference;//default SEG maximum difference between positions

private:
    size_t*     addresses_;
    size_t      length_;
    bool        hashed_;

    static const size_t sc_sizealphabet_;
};

////////////////////////////////////////////////////////////////////////////
// Class SEGSequence INLINES
//
// SeqValidator: indicate whether the given letter is valid
//
inline
bool SEGSequence::SeqValidator( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw MYRUNTIME_ERROR( "SEGSequence::SeqValidator: Memory access error." );
#endif
    return *(size_t*)letter != GAP;//  &&  *(size_t*)letter != ASTERISK;
}

// -------------------------------------------------------------------------
// SeqVerifier: indicate whether the letter is to be considered
//
inline
bool SEGSequence::SeqVerifier( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw MYRUNTIME_ERROR( "SEGSequence::SeqVerifier: Memory access error." );
#endif
    return *(size_t*)letter < GetSeqAlphabetSize();
}

// -------------------------------------------------------------------------
// SeqComparer: return 0 if two letters are equal, 1 if the first is
// greater than the second, -1 otherwise
//
inline
int SEGSequence::SeqComparer( const void* one, const void* another )
{
#ifdef __DEBUG__
    if( !one || !another )
        throw MYRUNTIME_ERROR( "SEGSequence::SeqComparer: Memory access error." );
#endif
    return (int)(*(size_t*)one) - (int)(*(size_t*)another);
}

// -------------------------------------------------------------------------
// GetAddressAt: return the address (translated letter) at the given 
// position
//
inline
size_t SEGSequence::GetAddressAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGSequence::GetAddressAt: Memory access error." );
#endif
    return addresses_[pos];
}

// -------------------------------------------------------------------------
// SetAddressAt: set the address (translated letter) at the given position
//
inline
void SEGSequence::SetAddressAt( size_t pos, size_t value )
{
#ifdef __DEBUG__
    if( !addresses_ || length_ <= pos )
        throw MYRUNTIME_ERROR( "SEGSequence::SetAddressAt: Memory access error." );
#endif
    addresses_[pos] = value;
}

// -------------------------------------------------------------------------
// GetResidue: get residue given its address
//
inline
size_t SEGSequence::GetSeqResidue( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw MYRUNTIME_ERROR( "SEGSequence::GetSeqResidue: Memory access error." );
#endif
    return *( size_t*)letter;
}

}//namespace SEG

#endif//__SEGSequence_h__
