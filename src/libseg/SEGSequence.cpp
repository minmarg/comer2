/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdlib.h>
#include <string.h>

#include "liblib/alpha.h"
#include "SEGSequence.h"

namespace SEG {

//alphabet size for sequences
const size_t SEGSequence::sc_sizealphabet_ = NUMAA;

const size_t SEGSequence::scszDefSegSeqWinLength = 12;//default SEG window length
const float  SEGSequence::scfpDefSegSeqLowEntropy = 2.2f;//default SEG low entropy threshold
const float  SEGSequence::scfpDefSegSeqHighEntropy = 2.5f;//default SEG high entropy threshold
const size_t SEGSequence::scszDefSegSeqMaxDifference = 100;//default SEG maximum difference between positions

// /////////////////////////////////////////////////////////////////////////
// CLASS SEGSequence
//
// constructor:
//
// NOTE: address is supposed to point to a vector of pointers void*
//
//
SEGSequence::SEGSequence(
    const unsigned char* residues,
    size_t          seqlength,
    bool            hashed,
    size_t          winlen,
    float           lowent,
    float           highent,
    size_t          maxdiff )
:
    //ATTENTION!
    SEGAbstract(
        &SEGSequence::SeqValidator,
        &SEGSequence::SeqVerifier,
        &SEGSequence::SeqComparer,
        NULL,

        addresses_ = AllocateAddresses( seqlength ),
        seqlength,

        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize()/*alphabet size*/
    ),
    length_( seqlength ),
    hashed_( hashed )
{
    if( addresses_ == NULL )
        throw MYRUNTIME_ERROR( "SEGSequence::SEGSequence: Not enough memory." );

    Translate( residues, seqlength );
}

// destructor:
//
SEGSequence::~SEGSequence()
{
    Destroy( addresses_ );
}

// -------------------------------------------------------------------------
// AllocateAddresses: allocate memory for addresses
//
size_t* SEGSequence::AllocateAddresses( size_t newlen )
{
    size_t* localvector = (size_t*)malloc( sizeof(size_t) * newlen );

    if( localvector == NULL )
        return localvector;

    memset( localvector, 0, sizeof(size_t) * newlen );

    return localvector;
}

// -------------------------------------------------------------------------
// Destroy: deallocate memory occupied by addresses
//
void SEGSequence::Destroy( size_t* addrvector )
{
    if( addrvector )
        free( addrvector );
}

// -------------------------------------------------------------------------
// Translate: translate residue codes into the vector of addresses 
// suitable for abstract SEG
//
void SEGSequence::Translate( const unsigned char* residues, size_t len )
{
    if( GetLocalLength() < len )
        throw MYRUNTIME_ERROR( "SEGSequence::Translate: Memory access error." );

    for( size_t n = 0; n < len; n++ )
        if( GetHashed())
            SetAddressAt( n, residues[n]);
        else
            SetAddressAt( n, HashAlphSymbol(residues[n]));
}

// -------------------------------------------------------------------------
// PrintSequence: print formatted sequence
//
void SEGSequence::PrintSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSequence( fp, &SEGSequence::GetSeqResidue, width );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: print formatted sequence with segments found by
//     running the algorithm and masked with Xs
//
void SEGSequence::PrintSeggedSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSeggedSequence( fp, &SEGSequence::GetSeqResidue, width );
}

}//namespace SEG
