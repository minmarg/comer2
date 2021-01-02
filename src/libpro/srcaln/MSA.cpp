/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
 
// #include <time.h>
// #include <math.h>
#include <cmath>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <memory>
#include <vector>

#include "extsp/psl.h"
#include "liblib/msg.h"
#include "liblib/mysort.h"
#include "liblib/logitnormal.h"
#include "liblib/BinarySearchStructure.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/SUBSTABLE.h"
#include "libpro/srcpro/BMProtoProfile.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libpro/srcsco/ProfileMatrix.h"
#include "libseg/SEGSequence.h"
#include "MSA.h"

//global file variables
//
static const bool   gbTOLXs = false;//include X into calculation of sequence weights
static const float  gdWCPC = 0.5f;//constant power coefficient in weight calculation
static const size_t gszMINEXTINTRV = 10;//minimum extent interval in calculation of seq. weights
static const float  gszMINPRCSEQNS = 0.5f;//minimum percentage of sequences in extent at a position
static const bool   gbUSEXTENTS = true;//use extents in calculation of seq. weights

////////////////////////////////////////////////////////////////////////////
// CLASS MSA
//
const size_t MSA::scnDefNoSeqs = 128;//default number of sequences
const size_t MSA::scnDefNoPoss = 1024;//default number of positions

// Constructor
//

MSA::MSA( size_t noseqs )
:   sequences_( NULL ),
    protoprofile_( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    identity_level_( MOptions::GetIDENTITY()),
    effnoseqs_( 0.0f ),
    name_( NULL ),
    description_( NULL ),
    keepseqdescs_( false ),
    ignoregapsinquery_( false/*true*/ ),
    deletestateson_( true/*false*/ ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
//     HDPctbase_( NULL ),

    usingsegfilt_( false ),
    segfiltwinlenval_( MOptions::GetHCWINDOW()),
    segfiltlowentval_( MOptions::GetHCLOWENT()),
    segfilthighentval_( MOptions::GetHCHIGHENT()),

    usingseqseg_( false ),
    seqsegwinlenval_( MOptions::GetLCWINDOW()),
    seqseglowentval_( MOptions::GetLCLOWENT()),
    seqseghighentval_( MOptions::GetLCHIGHENT()),

    extminwindow_( MOptions::GetMINALNPOS()),
    extminseqperc_( MOptions::GetMINALNFRN()),
    pseudocntweight_( MOptions::GetPCFWEIGHT())
{
    Realloc( noseqs );
}

// -------------------------------------------------------------------------
// Destructor
//
MSA::~MSA()
{
    if( protoprofile_ )
        delete protoprofile_;
    protoprofile_ = NULL;

    clear();

    if( sequences_ )
        free( sequences_ );
    sequences_ = NULL;
}

// -------------------------------------------------------------------------
// PlainPreprocess: preprocess MSA
//
void MSA::PlainPreprocess()
{
    InitProtoProfile();
    SetStates();
}

// -------------------------------------------------------------------------
// SelectSequences: select non-redundant set of sequences from the MSA
//
void MSA::SelectSequences()
{
    PreprocessAlignment();

    if( !GetSize())
        throw MYRUNTIME_ERROR( "MSA::SelectSequences: Memory unallocated." );

    PurgeAtSequenceIdentity();

    InitProtoProfile();
    SetStates();

    if( GetUsingSeqSEGFilter())
        FilterSequencesWithSEG();

    if( GetUsingSEGFilter())
        RefineWithSEG();

// //{{TEST
// if( GetUsingSEGFilter() || GetUsingSeqSEGFilter())
//     PurgeAtSequenceIdentity();
// //}}

    SetBackgroundProbabilities();
    SetPosNoSequences();
}

// -------------------------------------------------------------------------
// CalculateXCovMatrices: calculate cross-covariance matrices betwee the 
// positions of the profile
//
void MSA::CalculateXCovMatrices(const char* filename)
{
    float avgpeseq;

    SelectSequences();

    ComputeGlobSequenceWeights( &avgpeseq );

#if 1
    if( gbUSEXTENTS )
        ComputeExtents();

    if( gbUSEXTENTS )
          ComputeMIDstateSequenceWeights();
    else  ComputeMIDstateSequenceWeightsNoExtents();
#else
    ComputeGWMIDstateSequenceWeights( avgpeseq );
    ComputeTransitionFrequencies( true/*gwghts*/, false );
#endif

    AdjustWeights();

    ComputeTargetFrequenciesMDLVar();

    CalculateAndPrintXCovMatrices(filename);
}

// -------------------------------------------------------------------------
// ConstructProfile: the main procedure for 
// profile construction/model inference
//
void MSA::ConstructProfile()
{
    float avgpeseq;

    SelectSequences();

    ComputeGlobSequenceWeights( &avgpeseq );

#if 1
    if( gbUSEXTENTS )
        ComputeExtents();

//     ComputeMIDstateSequenceWeights();
//     ComputeTransitionFrequencies( true );
//     DeriveExpectedNoObservations();

    if( gbUSEXTENTS )
          ComputeMIDstateSequenceWeights();
    else  ComputeMIDstateSequenceWeightsNoExtents();
    ComputeTransitionFrequencies( true/*gwghts*/, false );
#else
    ComputeGWMIDstateSequenceWeights( avgpeseq );
    ComputeTransitionFrequencies( true/*gwghts*/, false );
#endif

    CalculateEffNoSequences();

    AdjustWeights();
    ComputeTargetTransFrequencies();

//     ComputeTargetFrequencies();
    ComputeTargetFrequenciesMDLVar();
    //do not mix target frequencies here; 
    //they will be mixed, if applied, for query just before searching
    if( GetScoAdjmentHDPCtx())
        //NOTE:calculate posterior predictives once profile has been scaled
        ;//CalcTFPosteriorPredictives();
    else if( GetTarFrMixHDPCtx())
        MixTargetFrequenciesHDPCtx();
    RecalcBackgroundProbabilities();
    CalcPosteriorProbabilities();
    ComputePSSM();
}

// -------------------------------------------------------------------------
// InitProtoProfile: Initialize prototype profile model
//
void MSA::InitProtoProfile()
{
    BMSequence* sequence = GetSequenceAt(0);
    if( !sequence )
        throw MYRUNTIME_ERROR( "MSA::InitProtoProfile: Null sequence 0." );

    if( protoprofile_ )
        delete protoprofile_;
    protoprofile_ = new BMProtoProfile(*sequence);
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory
//
void MSA::Realloc( size_t newcap )
{
    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        sequences_ = (BMSequence**)malloc( sizeof(void*) * newcap );
    } else {
        sequences_ = (BMSequence**)realloc( sequences_, sizeof(void*) * newcap );
    }

    if( !sequences_ )
        throw MYRUNTIME_ERROR( "MSA::Realloc: Not enough memory." );

    BMSequence** seqs = sequences_;

    if( capacity_ != 0 ) {
        seqs = sequences_ + capacity_;
    }

    memset( seqs, 0, sizeof(void*) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: push sequence into the MSA
//
void MSA::push( BMSequence* seq )
{
    if( capacity_ <= length_ ) {
        size_t newcap = TIMES2( capacity_ );
        if( newcap <= length_ )
            newcap = length_ + 1;
        Realloc( newcap );
    }

    sequences_[length_] = seq;

    length_++;
}

// -------------------------------------------------------------------------
// push: clear all the sequences in the MSA
//
void MSA::clear()
{
    if( name_ ) { free( name_ ); name_ = NULL; }
    if( description_ ) { free( description_ ); description_ = NULL; }

    if( sequences_ == NULL )
        return;

    for( size_t n = 0; n < length_; n++ )
        if( sequences_[n])
            delete sequences_[n];

    memset( sequences_, 0, sizeof(void*) * capacity_ );

    length_ = 0;
    effnoseqs_ = 0.0f;
}

// -------------------------------------------------------------------------
// PreprocessAlignment: make terminal sequence gaps unused
//
void MSA::PreprocessAlignment()
{
    BMSequence* one = NULL;
    size_t i, p;

    for( i = 0; i < GetSize(); i++ )
    {
        one = GetSequenceAt(i);
        if( !one || !one->GetUsed())
            continue;

        for( p = 0; p < one->GetSize() && !IsValidResSym( one->GetResidueAt(p)); p++ )
            if(i)
                one->SetUnusedAt(p);

        if( p < one->GetSize() && IsValidResSym( one->GetResidueAt(p)))
            one->SetFirstUsed(p);

        for( p = one->GetSize(); 1 <= p && !IsValidResSym( one->GetResidueAt(p-1)); p-- )
            if(i)
                one->SetUnusedAt(p-1);

        if( 1 <= p && IsValidResSym( one->GetResidueAt(p-1)))
            one->SetLastUsed(p-1);
    }
}

// -------------------------------------------------------------------------
// LenCompare: compare sequence lenths
static inline
int LenCompare( const void* lenvec, size_t n1, size_t n2 )
{
    const size_t* vec = ( const size_t* )lenvec;
    const size_t  nn1 = TIMES2(n1);
    const size_t  nn2 = TIMES2(n2);
    if( vec == NULL )
        return 0;
    //[n2]-[n1] -- to sort in descending order
//     return vec[nn2] - vec[nn1];
    return (int)(( vec[nn2] == vec[nn1])? vec[nn1+1] - vec[nn2+1]: vec[nn2] - vec[nn1] );
}

// -------------------------------------------------------------------------
// IdlCompare: compare sequence identity levels
static inline
int IdlCompare( const void* key1, const void* key2 )
{   //key2-key1 -- to sort in descending order
    return (int)((ssize_t)key1 - (ssize_t)key2 );
}

// -------------------------------------------------------------------------
// PurgeAtSequenceIdentityObs: remove sequences that share similarity above 
// certain degree of sequence identity
//
void MSA::PurgeAtSequenceIdentityObs()
{
    mystring preamb = "MSA::PurgeAtSequenceIdentityObs";
    const bool  bPASS1 = true; 
    const bool  bPASS2 = false; 
    const float dIDTHR2 = 0.9f;//identity threshold in 2nd pass
    const size_t sMINSEG = 5;//minimum segment length in 2nd pass
    bool   simthr;//similar under threshold
    size_t segment;//length of aligned segment
    size_t matches;//number of residue matches between sequences
    size_t n, q;
    size_t i, j, e, p, ith, jth;
    bool            u1, u2;
    unsigned char   r1, r2;
    BMSequence* fst, *sec;

    if( !GetSize())
        throw MYRUNTIME_ERROR( preamb + "No sequences." );

    BMSequence* query = GetSequenceAt(0);

    if( query == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null sequence 0." );

    const size_t effsize = query->GetEffectiveSize();

    if( effsize == 0 )
        throw MYRUNTIME_ERROR( preamb + "Sequence 0 contains no residues." );

    //size_t  locindices[effsize];
	//size_t  lens[GetSize()];
	//size_t  lindsrt[GetSize()];
	std::unique_ptr<size_t[]> locindices(new size_t[effsize]);
	std::unique_ptr<size_t[]> lens(new size_t[GetSize()]);
	std::unique_ptr<size_t[]> lindsrt(new size_t[GetSize()]);

    //{{sort lengths of sequences
    for( i = 0; i < GetSize(); i++ ) {
        BMSequence* sqn = GetSequenceAt(i);
        if( sqn )
              lens[i] = sqn->GetEffectiveSize();
        else  lens[i] = 0;
    }
    lindsrt[0] = 0;
    //sort; the first position is reserved for query
    HeapSortInd( lindsrt.get()+1, lens.get()+1, GetSize()-1, LenCompare );
    for( i = 1; i < GetSize(); i++ )
        lindsrt[i]++;
    //}}

    // collect indices of match positions in query
    for( q = 0, n = 0; q < query->GetSize(); q++ ) {
        if( query->GetResidueAt(q) == GAP )
            continue;
        if( effsize <= n ) {
            n = 0;
            break;
        }
        locindices[n++] = q;
    }

    if( effsize != n )
        throw MYRUNTIME_ERROR( preamb + "Wrong effective size of query." );

    if( bPASS1 ) {
        // 1st pass...
        for( i = 0; i + 1 < GetSize(); i++ )
        {
            fst = GetSequenceAt( ith=lindsrt[i] );
            if( fst == NULL || !fst->GetUsed())
                continue;

            for( j = i + 1; j < GetSize(); j++ )
            {
                sec = GetSequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;

                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < effsize; e++ ) {
                    // iterate over match positions
                    p = locindices[e];

                    u1 = fst->IsUsedAt( p );
                    u2 = sec->IsUsedAt( p );
                    r1 = fst->GetResidueAt( p );
                    r2 = sec->GetResidueAt( p );

                    if( !u1 && !u2 )
                        continue;
                    if( r1 == GAP || r2 == GAP )
                        continue;
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( segment && GetIdentityLevel() <= (float)matches / (float)segment )
                    sec->SetUsed( false );//set this sequence not to be used
            }
        }
    }

    if( bPASS2 ) {
        // 2nd pass...
        for( i = 0; i + 1 < GetSize(); i++ )
        {
            fst = GetSequenceAt( ith=lindsrt[i] );
            if( fst == NULL || !fst->GetUsed())
                continue;

            for( j = i + 1; j < GetSize(); j++ )
            {
                sec = GetSequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;

                simthr = true;
                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < effsize; e++ ) {
                    // iterate over match positions
                    p = locindices[e];

                    u1 = fst->IsUsedAt( p );
                    u2 = sec->IsUsedAt( p );
                    r1 = fst->GetResidueAt( p );
                    r2 = sec->GetResidueAt( p );

                    if( !u1 && !u2 )
                        continue;
                    if( r1 == GAP || r2 == GAP ) {
                        if( segment && sMINSEG <= segment && 
                          (float)matches / (float)segment < dIDTHR2 ) {
                            simthr = false;
                            break;
                        }
                        segment = 0;
                        matches = 0;
                    }
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( simthr )
                    sec->SetUsed( false );
            }
        }
    }

// for( i = 0; i < size(); i++ )
//     if( !GetSequenceAt(i)->GetUsed())
//         fprintf( stderr, "s %d\n", i );
}

// -------------------------------------------------------------------------
// PurgeAtSequenceIdentity: retain diverse subset of sequences
// TODO: reimplement large arrays to use heap
//
void MSA::PurgeAtSequenceIdentity()
{
    mystring preamb = "MSA::PurgeAtSequenceIdentity";
    const bool      bRESORT = false;
    const ssize_t   nDETNOSEQS = 100;
    const size_t    szMINNOSEQS = 50;
    float           dLOSID = 0.2f;
    const float     dHISID = GetIdentityLevel();
    float           dININC = 0.01f;//initial increment
    const float     dUPINC = 0.05f;//upper limit of initial increment
    const float     dCPC = 20.0f;//power coefficient in map function
    const size_t    szFRAGLEN = 55;//length for fragment-based checking
    const size_t    szFRAGLENh = szFRAGLEN >> 1;
//     bool   simthr;//similar under threshold
    size_t segment;//length of aligned segment
    size_t matches;//number of residue matches in the sequences
    size_t q;
    size_t efst, eset, elst;
    size_t i, ii, j, c, e, p, ith, jth;
    bool            u1, u2;
    unsigned char   r1, r2;
    BMSequence* fst, *sec;

    if( !GetSize())
        throw MYRUNTIME_ERROR( preamb + "No sequences." );

    BMSequence* query = GetSequenceAt(0);

    if( query == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null sequence 0." );

    const size_t noseqs = GetSize();
    const size_t effsize = query->GetEffectiveSize();

    if( effsize == 0 )
        throw MYRUNTIME_ERROR( preamb + "Sequence 0 contains no residues." );

    size_t  sum, nomps = 0;
    //size_t  locindices[effsize];
    //size_t  lens[TIMES2(noseqs)];
    //size_t  lindsrt[noseqs];
	std::unique_ptr<size_t[]> locindices(new size_t[effsize]);
	std::unique_ptr<size_t[]> lens(new size_t[TIMES2(noseqs)]);
	std::unique_ptr<size_t[]> lindsrt(new size_t[noseqs]);
	size_t  rdist[NUMALPH];
    size_t  prefnoseqs = noseqs;
    BinarySearchStructure idlevs( IdlCompare, noseqs, true/*keep*/);
    const size_t idlsc = 1000;
    size_t  noids;
    ssize_t idl = 0;
    float   qdist[NUMALPH];
    float   re, val, sfrq;//fraction of set of sequences
    float   medidl = 0.0f;
    float   dsid, dthr, dinc;
    //float   idthrs[effsize];
    //size_t  nosqns[effsize];
	std::unique_ptr<float[]> idthrs(new float[effsize]);
	std::unique_ptr<size_t[]> nosqns(new size_t[effsize]);
	bool	passed, seld;

    //{{collect indices of match positions in query
    for( q = 0, nomps = 0; q < query->GetSize(); q++ ) {
        r1 = query->GetResidueAt(q);
        if( !IsValidResSym(r1) || r1 == X )
            continue;
        if( effsize <= nomps )
            throw MYRUNTIME_ERROR( preamb + "Invalid query size." );
        locindices[nomps++] = q;
    }
//     lens[0] = nomps;
//     lens[1] = idlsc;
    lens[0] = idlsc;
    lens[1] = 0;
    //}}

    //{{calculate sequence identity levels w/r to query
    for( i = 1, ii = 2; i < noseqs; i++, ii+=2 )
    {
        lens[ii] = lens[ii+1] = 0;
        if(( sec = GetSequenceAt(i)) == NULL )
            continue;
        segment = 0;
        matches = 0;
        for( e = 0, p = 0; e < nomps; e++ ) {
            // iterate over match positions
            p = locindices[e];

            u1 = query->IsUsedAt( p );
            u2 = sec->IsUsedAt( p );
            r1 = query->GetResidueAt( p );
            r2 = sec->GetResidueAt( p );

            if( !u1 || !u2 )
                continue;
            if( r1 == GAP || r2 == GAP )
                continue;
            if( r1 == X || r2 == X )
                continue;

            segment++;
            if( u1 && u2 && r1 == r2 )
                matches++;
        }
        if( segment ) {
            idl = (ssize_t)rintf( (float)(matches*idlsc) / (float)segment );
            idlevs.Push((const void*)idl );
//             lens[ii] = segment;
//             lens[ii+1] = ( size_t )idl;
            lens[ii] = (size_t)idl;
            lens[ii+1] = i;
        }
    }
    //}}

    if( bRESORT ) {
        //{{calculate rel. entropy for each sequence
        for( i = 0, ii = 0; i < noseqs; i++, ii+=2 )
        {
            lens[ii+1] = 0;
            if(( sec = GetSequenceAt(i)) == NULL )
                continue;
            sum = 0;
            memset( rdist, 0, NUMALPH * sizeof(size_t));
            for( e = 0, p = 0; e < nomps; e++ ) {
                // iterate over match positions
                p = locindices[e];

                u2 = sec->IsUsedAt( p );
                r2 = sec->GetResidueAt( p );

                if(!( u2 && IsValidResSym(r2) && r2 != X ))
                    continue;
                rdist[r2]++;
                sum++;
            }
            re = 0.0f;
            for( e = 0; e < NUMAA; e++ ) {
                val = ( (float)rdist[e]+STABLE.PROBABility((int)e)) / (float)(sum+1);
                if( !i )
                    qdist[e] = val;
                else if( val && qdist[e])
                    re += val * logf(val/qdist[e]);
            }
            lens[ii+1] = (size_t)rintf( (float)idlsc*re );
        }
        //}}
    }

    //{{sort lengths of sequences
    lindsrt[0] = 0;
    //sort; the first position is reserved for query
    HeapSortInd( lindsrt.get()+1, lens.get()+2, noseqs-1, LenCompare );
    for( i = 1; i < noseqs; i++ )
        lindsrt[i]++;
    //}}


    if( 0 < nDETNOSEQS ) {
        prefnoseqs = (size_t)nDETNOSEQS;
    }
    else if( !nDETNOSEQS ) {
        //get median value of identity levels
        noids = idlevs.GetSize();
        if( noids ) {
            if( noids & 1 )
                medidl = (float)(ssize_t)idlevs.GetValueAt(noids>>1) / (float)idlsc;
            else {
                medidl = (float)((ssize_t)idlevs.GetValueAt((noids>>1)-1) + 
                                 (ssize_t)idlevs.GetValueAt(noids>>1));
                medidl = medidl / (float)(TIMES2(idlsc));
            }
            medidl = SLC_MIN( 1.0f, medidl );
            //map sid level
    //         sfrq = expf( dCPC * medidl * logf( 1.001 - medidl ));
            sfrq = expf( ( dCPC*(1.0f-SQUARE(medidl)) + medidl ) * logf(1.001f-medidl) );
    //         sfrq = 1.001f - medidl;
            prefnoseqs = SLC_MAX( szMINNOSEQS, (size_t)rintf( sfrq*(float)(noids+1) ));
            if( noids + 1 < prefnoseqs )
                prefnoseqs = noids + 1;
        }
    }

    //{{incrementaly process mutual sequence identity levels
    if( noseqs <= prefnoseqs )
        dLOSID = dHISID;
    dININC = SLC_MAX( dININC, SLC_MIN( dUPINC, (float)prefnoseqs*dININC*0.01f ));
    for( i = 1; i < noseqs; i++ ) {
        sec = GetSequenceAt( ith=lindsrt[i] );
        if( sec )
            sec->SetUsed( false );//NOTE:initially unused
    }
    for( e = 0; e < nomps; e++ ) {
        idthrs[e] = 0.0f;
        nosqns[e] = 1;
    }
    //
    for( dsid = dLOSID, dinc = dININC, seld = false; dsid <= dHISID && !seld; dsid += dinc )
    {
        seld = true;
        for( e = 0, p = 0; e < nomps; e++ ) {
            passed = false;
            //borders
            i = ( e <= szFRAGLENh )? 0: e - szFRAGLENh;
            j = i + szFRAGLEN;
            if( nomps <= j ) {
                j = nomps - 1;
                i = 0;
                if( szFRAGLEN < j )
                    i = j - szFRAGLEN;
            }
            for( c = i; c <= j; c++ ) {
                if( prefnoseqs <= nosqns[c]) {
                    passed = true;
                    break;
                }
            }
            if( !passed ) {
                seld = false;
                idthrs[e] = dsid;
            }
        }
//         for( e = 0, p = 0; e < nomps; e++ ) {
//             p = locindices[e];
//             if( nosqns[e] < prefnoseqs ) {
//                 seld = false;
//                 idthrs[e] = dsid;
//             }
//         }
        if( seld )
            break;
        for( i = 0; i < noseqs; i++ )
        {
            fst = GetSequenceAt( ith=lindsrt[i] );
            if( fst == NULL || fst->GetUsed())
                continue;//cont. if processed

            passed = true;
            dthr = 0.0f;
            for( e = 0, p = 0; e < nomps; e++ ) {
                p = locindices[e];
                if( p < fst->GetFirstUsed()) continue;
                if( fst->GetLastUsed() < p ) break;
                if( dthr < idthrs[e])
                    dthr = idthrs[e];
                if( nosqns[e] < prefnoseqs )
                    passed = false;
            }
            if( dthr <= 0.0f )
                continue;
//             if( passed && dLOSID < dthr )
//                 continue;

            seld = false;
            passed = true;
            for( j = 0; j < i; j++ )
            {
                sec = GetSequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;//cont. if unused

                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < nomps; e++ ) {
                    //iterate over match positions
                    p = locindices[e];

                    u1 = fst->GetFirstUsed() <= p && p <= fst->GetLastUsed();//fst->IsUsedAt( p );
                    u2 = sec->GetFirstUsed() <= p && p <= sec->GetLastUsed();//sec->IsUsedAt( p );
                    r1 = fst->GetResidueAt( p );
                    r2 = sec->GetResidueAt( p );

                    if( !u1 || !u2 )
                       continue;
                    if( r1 == GAP || r2 == GAP )
                        continue;
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( segment && dthr <= (float)matches / (float)segment ) {
                    passed = false;
                    break;
                }
            }
            if( !passed )
                continue;
            fst->SetUsed( true );//set used
            //use fragment-based instead of pos.-specific mode
            for( e = efst = eset = elst = 0, p = 0; e < nomps; e++ ) {
                p = locindices[e];
                if( p < fst->GetFirstUsed()) continue;
                if( fst->GetLastUsed() < p ) break;
                if( !eset ) { efst = e; eset = 1; }
                elst = e;
                nosqns[e]++;
            }
//             for( e = efst-1, c = 0; e + 1 && c < szFRAGLENh; e--, c++ )
//                 nosqns[e]++;
//             if( elst )
//                 for( e = elst+1, c = 0; e < nomps && c < szFRAGLENh; e++, c++ )
//                     nosqns[e]++;
        }
    }
    //}}
}

// -------------------------------------------------------------------------
// CalculateNeff: calculate the effective number of sequences using the 
// formula: SUM_i w_i, where w_i = 1/k_i and k_i = SUM_j H(S(i,j)-S_thr), 
// with S(i,j) sequence identity between sequences i and j, S_thr a sequence 
// identity threshold, and H(.) the unit step function;
// idnt_thlds, sequence identity thresholds to calculate Neff at;
// neffs, output Neff for each identity threshold idnt_thlds;
//
void MSA::CalculateNeff(const std::vector<float>& idnt_thlds, std::vector<float>& neffs, 
    size_t& nmsaseqs, size_t& length1)
{
    mystring preamb = "MSA::CalculateNeff";
    if(idnt_thlds.size() < 1)
        return;

    //PreprocessAlignment();

    const size_t noseqs = GetSize();

    nmsaseqs = noseqs;
    length1 = (noseqs && GetSequenceAt(0))? GetSequenceAt(0)->GetEffectiveSize(): 0;

    if(noseqs < 1)
        return;

    size_t segment;//length of aligned segment
    size_t matches;//number of residue matches in the sequences
    size_t i, j, t, e1, e2, p1, p2;
    unsigned char r1;
    BMSequence* fst, *sec;
    std::vector<std::vector<size_t>> rawseqs(noseqs);
    std::vector<std::vector<int>> kseqs(idnt_thlds.size(), std::vector<int>(noseqs, 1));

    for(std::vector<size_t>& v: rawseqs)
        v.reserve(1024);

    //make condense representation of sequences
    for( i = 0; i < noseqs; i++ )
    {
        if((fst = GetSequenceAt(i)) == NULL)
            continue;
        for(e1 = 0; e1 < fst->GetSize(); e1++) {
            r1 = fst->GetResidueAt(e1);
            if(r1 == GAP || r1 == X)
                continue;
            rawseqs[i].push_back(e1);
        }
    }

    neffs.reserve(idnt_thlds.size());

    for( i = 0; i < noseqs; i++ )
    {
        if((fst = GetSequenceAt(i)) == NULL)
            continue;

        for( j = 0; j < i; j++ )
        {
            if((sec = GetSequenceAt(j)) == NULL)
                continue;

            matches = 0;
            segment = SLC_MIN(rawseqs[i].size(), rawseqs[j].size());

            for(e1 = e2 = 0; e1 < rawseqs[i].size() && e2 < rawseqs[j].size();)
            {
                //{{unnecessary but avoids compiler "over-"optimization
                p1 = rawseqs[i][e1];
                p2 = rawseqs[j][e2];
                //}}
                for(; e1 < rawseqs[i].size() && (p1 = rawseqs[i][e1]) < p2; e1++);
                for(; e2 < rawseqs[j].size() && (p2 = rawseqs[j][e2]) < p1; e2++);

                if(p1 == p2 && e1 < rawseqs[i].size() && e2 < rawseqs[j].size()) {
                    if( fst->GetResidueAt(p1) == sec->GetResidueAt(p2))
                        matches++;
                    e1++;
                    e2++;
                }
            }

            float fidn = segment? (float)matches/(float)segment: 0.f;

            for(t = 0; t < idnt_thlds.size(); t++) {
                if(idnt_thlds[t] <= fidn) {
                    kseqs[t][i]++;
                    kseqs[t][j]++;
                }
            }
        }
    }

    //calculate neffs for each given threshold
    for(t = 0; t < kseqs.size(); t++) {
        float nf = 0.f;
        for(i = 0; i < kseqs[t].size(); i++)
            nf += 1.f/(float)kseqs[t][i];
        neffs[t] = nf;
    }
}

// -------------------------------------------------------------------------
// RefineWithSEG: apply the SEG algorithm to MSA columns corresponding to 
// match positions
//
void MSA::RefineWithSEG()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::RefineWithSEG: Null prototype profile." );

    const size_t    szFRAGLEN = 40;//length subtracted from both ends
    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const size_t    noseqs = GetSize();
    size_t  segwinlenval    = GetSEGWindow();
    float   seglowentval    = GetSEGLowEntropy();
    float   seghighentval   = GetSEGHighEntropy();

    if( noseqs <= 1 || noseqs < segwinlenval )
        return;

    myruntime_error mre;
    unsigned char*  validseq = ( unsigned char* )malloc( noseqs * sizeof( unsigned char ));
    unsigned char*  wholeseq = ( unsigned char* )malloc( noseqs * sizeof( unsigned char ));
    char*           omitmask = ( char* )malloc( noseqs * sizeof( char ));

    size_t          validlen = 0;
    size_t          wholelen = 0;
    size_t          p, pm, i;

    BMSequence*     seqn = NULL;
    unsigned char   res;
    int             pstate;

    try {
        if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
            throw MYRUNTIME_ERROR( "MSA::RefineWithSEG: Not enough memory." );

//         memset( validseq, 0, sizeof( unsigned char ) * size());
//         memset( wholeseq, 0, sizeof( unsigned char ) * size());
//         memset( omitmask, 0, sizeof( char ) * size());

        for( p = 0, pm =(size_t)-1; p < protoprofile_->GetSize(); p++ )
        {
            pstate = protoprofile_->GetStateAt( p );
            if( !pstate )
                continue;

            pm++;
            if( szFRAGLEN <= pm && pm + szFRAGLEN < nomstates )
                continue;

            validlen = 0;
            wholelen = 0;

            //iterate over all sequences excluding the first (query)
            for( i = 1; i < noseqs; i++ )
            {
                seqn = GetSequenceAt( i );
                if( seqn == NULL )
                    throw MYRUNTIME_ERROR( "MSA::RefineWithSEG: Memory access error." );
                res = seqn->GetResidueAt( p );

                wholelen++;
                omitmask[wholelen-1] = 0;
                wholeseq[wholelen-1] = res;

                if( !seqn->GetUsed() || !seqn->IsUsedAt( p ) || !IsValidResSym( res ) || res == X ) {
                    omitmask[wholelen-1] = 1;
                    continue;
                }

                validlen++;
                validseq[validlen-1] = res;
            }

            if( validlen < segwinlenval )
                continue;

            //{{SEG logic
            SEG::SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlenval,
                seglowentval,
                seghighentval
            );
            segseq.SetHighCSearch();//Important!
            segseq.Run();
            segseq.MaskSequence((char*)wholeseq, wholelen, X, omitmask, (char)255/*symmask1*/, (char)255/*symmask2*/);
            //}}

            wholelen = 0;
            //set segged symbols as unused
            for( i = 1; i < noseqs; i++ ) {
                seqn = GetSequenceAt( i );

                wholelen++;
                if( omitmask[wholelen-1] )
                    continue;

                if( wholeseq[wholelen-1] == X )
//                     seqn->SetUnusedAt( p );
                    seqn->SetResidueAt( X, p );//sequence modification!
            }
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( validseq ) free( validseq );
    if( wholeseq ) free( wholeseq );
    if( omitmask ) free( omitmask );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// FilterSequencesWithSEG: filter all sequences in the MSA with SEG
//
void MSA::FilterSequencesWithSEG()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::FilterSequencesWithSEG: Null prototype profile." );

    size_t  seqsize = protoprofile_->GetSize();
    size_t  segwinlen   = GetSeqSEGWindow();
    float   seglowent   = GetSeqSEGLowEntropy();
    float   seghighent  = GetSeqSEGHighEntropy();

    if( seqsize <= 1 || seqsize < segwinlen )
        return;

    myruntime_error mre;
    unsigned char*  validseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    unsigned char*  wholeseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    char*           omitmask = ( char* )malloc( sizeof( char ) * seqsize );

    size_t          validlen = 0;
    size_t          wholelen = 0;

    BMSequence*     sequence = NULL;
    unsigned char   residue;

    try {
        if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
            throw MYRUNTIME_ERROR( "MSA::FilterSequencesWithSEG: Not enough memory." );

//         memset( validseq, 0, sizeof( unsigned char ) * seqsize );
//         memset( wholeseq, 0, sizeof( unsigned char ) * seqsize );
//         memset( omitmask, 0, sizeof( char ) * seqsize );

        //iterate over all sequences
        for( size_t i = 0; i < GetSize(); i++ ) {
            sequence = GetSequenceAt( i );
            if( sequence == NULL )
                throw MYRUNTIME_ERROR( "MSA::FilterSequencesWithSEG: Memory access error." );
            if( !sequence->GetUsed())
                continue;

            validlen = 0;
            wholelen = 0;

            //iterate over all positions
            for( size_t p = 0; p < seqsize; p++ )
            {
                residue = sequence->GetResidueAt( p );

                wholelen++;
                omitmask[wholelen-1] = 0;
                wholeseq[wholelen-1] = residue;

                if( !sequence->IsUsedAt( p ) || residue == GAP ) {
                    omitmask[wholelen-1] = 1;
                    continue;
                }

                validlen++;
                validseq[validlen-1] = residue;
            }

            if( validlen < segwinlen )
                continue;

            //{{SEG logic
            SEG::SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlen,
                seglowent,
                seghighent
            );
            segseq.Run();
            segseq.MaskSequence((char*)wholeseq, wholelen, X, omitmask, (char)255/*symmask1*/, (char)255/*symmask2*/ );
            //}}

            wholelen = 0;
            //set segged symbols as unused
            for( size_t p = 0; p < seqsize; p++ )
            {
                wholelen++;
                if( omitmask[wholelen-1] )
                    continue;

                if( wholeseq[wholelen-1] == X )
    //                 sequence->SetUnusedAt( p );
                    sequence->SetResidueAt( X, p );//sequence modification!
            }
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( validseq ) free( validseq );
    if( wholeseq ) free( wholeseq );
    if( omitmask ) free( omitmask );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// SetPosNoSequences: set the number of used sequences per column; 
// the overall number of sequences is also set
//
void MSA::SetPosNoSequences()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::SetPosNoSequences: Null prototype profile." );

    size_t  posnos = 0;//number of sequences per column
    size_t  nos = 0;//number of sequences
    bool    allgaps;

    for( size_t p = 0; p < protoprofile_->GetSize(); p++ ) {
        posnos = 0;
        allgaps = true;

        //iterate over all sequences
        for( size_t i = 0; i < GetSize(); i++ )
        {
            BMSequence* sequence = GetSequenceAt( i );
            if( !sequence )
                continue;
            if( !sequence->GetUsed())
                continue;

            unsigned char residue = sequence->GetResidueAt( p );

            if( allgaps && IsValidResSym( residue ))
                allgaps = false;
            if(!( sequence->IsUsedAt( p ) && IsValidResSym( residue )))
                continue;

            posnos++;
        }
        protoprofile_->SetNoSequencesAt( posnos, p );
        if( allgaps )
            protoprofile_->SetUnusedAt(p);
    }

    //set the number of sequences...
    for( size_t i = 0; i < GetSize(); i++ )
        if( GetSequenceAt(i) && GetSequenceAt(i)->GetUsed())
            nos++;

    SetNoSequences( nos );
}

// -------------------------------------------------------------------------
// SetStates: set match states
//
void MSA::SetStates()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::SetStates: Null prototype profile." );

    unsigned char   residue;
    size_t          p, nn;

    for( p = 0, nn = 0; p < protoprofile_->GetSize(); p++ ) {
        residue = protoprofile_->GetResidueAt( p );
        if( !protoprofile_->IsUsedAt( p ) || !IsValidResSym( residue ))
            continue;
        //match states are determined by query  
        protoprofile_->SetStateAt( p );
        nn++;
    }

    protoprofile_->SetEffectiveSize( nn );
}

// -------------------------------------------------------------------------
// SetBackgroundProbabilities: compute background probabilities for the 
// query; apply pseudo counts
//
void MSA::SetBackgroundProbabilities()
{
    mystring preamb = "MSA::SetBackgroundProbabilities: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    const float     backpseudocounts = 120.0f;

    const float     accuracy = 1.0e-4f;
    float           prob, consv;
    float           weight;
    unsigned char   residue;
    const size_t    effobslen = NUMAA;
    const size_t    obslen = NUMALPH;
    size_t          observs[obslen];
    size_t          noobs;
    size_t          p;

    static int  symB = HashAlphSymbol('B');//N,D
    static int  symZ = HashAlphSymbol('Z');//Q,E
    static int  symJ = HashAlphSymbol('J');//I,L
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');
    static int  resI = HashAlphSymbol('I');
    static int  resL = HashAlphSymbol('L');

    noobs = 0;
    memset( observs, 0, obslen * sizeof(size_t));

    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        residue = protoprofile_->GetResidueAt( p );
        if( !protoprofile_->IsUsedAt( p ) ||
            residue == GAP || residue == X )//|| residue == ASTERISK )
            continue;

        noobs += 2;

        if( residue == symB ) {
            observs[resN]++;
            observs[resD]++;
            continue;
        }
        if( residue == symZ ) {
            observs[resQ]++;
            observs[resE]++;
            continue;
        }
        if( residue == symJ ) {
            observs[resI]++;
            observs[resL]++;
            continue;
        }
        if( NUMAA <= residue )
            throw MYRUNTIME_ERROR( preamb + "Unrecognized residue." );

        observs[residue] += 2;
    }

    if( backpseudocounts <= 0.0f && noobs < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid counts." );

    consv = 0.0f;
    weight = backpseudocounts / (noobs+backpseudocounts);

    if( noobs < 1 )
        noobs = 1;

    for( p = 0; p < effobslen/*obslen*/; p++ ) {
        prob = (1.0f-weight) * (float)observs[p] / (float)noobs + 
                weight * STABLE.PROBABility((int)p);
        consv += prob;
        protoprofile_->SetBackProbsAt((unsigned char)p, prob );
    }
	for (; p < obslen; p++) {
		protoprofile_->SetBackProbsAt((unsigned char)p, 0.0f);
	}

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + "Probabilities are not conserved." );

#ifdef USEPROFBACKPROBS
    STABLE.StoreProbabilities( protoprofile_->GetBackProbs());
#endif
}

// -------------------------------------------------------------------------
// ComputeExtents: compute extents (left-right boundaries) for each match
// position;
// extent is a reduced multiple alignment constructed for each position
// separately so that sequences in the constructed alignment contribute a 
// residue or internal gap symbol
//
void MSA::ComputeExtents()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::ComputeExtents: Null prototype profile." );

    myruntime_error mre;
    unsigned char residue;
    ssize_t matchbeg = -1;//first match state position
    ssize_t matchend = -1;//last match state position

    size_t  extwinlen = GetExtentMinWindow();
    size_t  extcentre =   extwinlen >> 1;
    size_t  extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust
    size_t  i, p;

//     //do not use percenatages as this anyway does not relate to domain boundaries
//     extwinlen = (size_t)((float)protoprofile_->GetEffectiveSize() * GetExtentMinSeqPercentage());
    if( extwinlen < GetExtentMinWindow())
        extwinlen = GetExtentMinWindow();
    extcentre =   extwinlen >> 1;
    extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust

// fprintf(stderr,"W=%d P=%f w=%d c=%d a=%d\n",GetExtentMinWindow(),GetExtentMinSeqPercentage(),extwinlen,extcentre,extadjust);
    size_t*     matchpos = NULL;    //accumulated number of match positions
    size_t*     msvector = NULL;    //accumulated match state vector

    matchpos = (size_t*)malloc( sizeof(size_t) * protoprofile_->GetSize());
    msvector = (size_t*)malloc( sizeof(size_t) * protoprofile_->GetSize());

    try {
        if( !matchpos || !msvector )
            throw MYRUNTIME_ERROR( "MSA::ComputeExtents: Not enough memory." );

        memset( matchpos, 0, sizeof(size_t) * protoprofile_->GetSize());
        memset( msvector, 0, sizeof(size_t) * protoprofile_->GetSize());

        for( p = 0; p < protoprofile_->GetSize(); p++ ) {
            if( IsValidResSym( protoprofile_->GetResidueAt( p ))) {
                if( matchbeg < 0 ) matchbeg = p;
                matchend = p;
            }
        }
        if( matchbeg < 0 || matchend < 0 )
            throw MYRUNTIME_ERROR( "MSA::ComputeExtents: Invalid match state sequence boundaries." );

        //iterate over all sequences in the MSA
        for( i = 0; i < GetSize(); i++ )
        {
            BMSequence* sequence = GetSequenceAt(i);
            if( sequence == NULL || !sequence->GetUsed())
                continue;

            size_t  seqsize = protoprofile_->GetSize();
            size_t  left = MYSIZE_MAX,
                    right = MYSIZE_MAX;

            //set match pos.: cumulative values
            memset( matchpos, 0, seqsize * sizeof(size_t));
            for( p = 0; p < seqsize; p++ ) {
                if( p )
                    matchpos[p] = matchpos[p-1];
                if( sequence->IsUsedAt( p ) &&
                    IsValidResSym( sequence->GetResidueAt( p )) && 
                    IsValidResSym( protoprofile_->GetResidueAt( p ))) {
                        if( p ) matchpos[p] = matchpos[p-1] + 1;
                        else    matchpos[p] = 1; 
                }
            }
            //set support match pos.
            if( i == 0 )
                memcpy( msvector, matchpos, seqsize * sizeof(size_t));

            for( p = 0; p < seqsize; p++ ) {
                if( !sequence->IsUsedAt( p )) {
                    continue;
                }
                residue = sequence->GetResidueAt( p );

                protoprofile_->IncCountAt( p );
                protoprofile_->IncDistributionAt( residue, p );

                if( residue != GAP ) {
                    if( left == MYSIZE_MAX ) left = p;
                    right = p;
                }
            }

            if( left == MYSIZE_MAX || right == MYSIZE_MAX )
                continue;

            size_t  sres, mres, pos, lres;//, rres;
            ssize_t lmsbder = 0, rmsbder = 0;
            ssize_t lborder = 0, rborder = 0;
            ssize_t leftadj = 0, rightadj = 0;

            if((ssize_t)left < matchbeg ) left = matchbeg;
            if((size_t)matchend < right ) right = matchend;

            lres = 0; if(( left? matchpos[left-1]: 0 ) < matchpos[left]) lres = 1;
//             rres = 0; if(( right? matchpos[right-1]: 0 ) < matchpos[right]) rres = 1;

            for( p = 0; p < seqsize; p++ ) {
                if( !sequence->IsUsedAt( p )) {
                    continue;
                }

                leftadj = extcentre - extadjust;
                rightadj = extcentre;

                lborder = matchpos[p] - 1 - leftadj;
                rborder = matchpos[p] - 1 + rightadj;

                lmsbder = msvector[p] - 1 - leftadj;
                rmsbder = msvector[p] - 1 + rightadj;

                //process cases when window is out of boundaries of match-state sequences
                if( lmsbder < 0 ) { 
                    pos = p;
                    if((ssize_t)p < matchbeg ) pos = matchbeg;
                    sres = 0;
                    mres = 0;
                    if( 0 < matchbeg ) {
                        if( matchpos[matchbeg-1] < matchpos[matchbeg]) sres = 1;
                        if( msvector[matchbeg-1] < msvector[matchbeg] ) mres = 1;
                    } else {
                        if( 0 < matchpos[matchbeg]) sres = 1;
                        if( 0 < msvector[matchbeg] ) mres = 1;
                    }
//                     lborder = matchpos[pos] - msvector[pos]; //exactly the same as below
                    lborder = matchpos[matchbeg] - sres + 
                            //difference of beginnings of this and match sequences
                            ( matchpos[pos] - matchpos[matchbeg] + sres ) - 
                            ( msvector[pos] - msvector[matchbeg] + mres );
                    rborder = lborder + extwinlen - 1; 
                }
                if(( ssize_t )( msvector[matchend] - 1 ) < rmsbder ) {
                    pos = p;
                    if((size_t)matchend < p ) pos = matchend;
                    sres = 0;
                    mres = 0;
                    if( 0 < pos ) {
                        if( matchpos[pos-1] < matchpos[pos]) sres = 1;
                        if( msvector[pos-1] < msvector[pos] ) mres = 1;
                    } else {
                        if( 0 < matchpos[pos]) sres = 1;
                        if( 0 < msvector[pos] ) mres = 1;
                    }

                    rborder = matchpos[matchend] - 1 + //for -1 assume there is always at least 1 res. in sequence
                            //difference of tails of this and match sequences
                            ( msvector[matchend] - msvector[pos] + mres ) - 
                            ( matchpos[matchend] - matchpos[pos] + sres ); 
                    lborder = rborder - extwinlen + 1;
                }

                //include sequence in the extent
                if( left <= p && p <= right ) {
//                     protoprofile_->PushIndexAt( i, p );//avoid because of the need for additional memory
                    sequence->SetUsedInExtentAt( p ); //this makes computations faster
                    protoprofile_->IncNoSequencesInExtentAt( p );//one more sequence in the extent at the position
                    protoprofile_->IncNoSequencesInMSExtentAt( p );
                }

//     fprintf(stderr,"%d,%d:: %d - %d <= %d  &&  %d <= %d - %d::   %d <= %d && %d <= %d:: %d\n",
//     i,p,matchpos[left],lres,lborder,rborder,matchpos[right],1,
//     (ssize_t)(matchpos[left]-lres),lborder,rborder,(ssize_t)(matchpos[right]-1),
//     (ssize_t)(matchpos[left]-lres)<=lborder && rborder<=(ssize_t)(matchpos[right]-1));

                //set extent boundaries
                if((( ssize_t ) extwinlen <= ( ssize_t )( msvector[matchend] - msvector[matchbeg] + 1 )) ?
                  (( ssize_t )( matchpos[left] - lres )<= lborder && 
                    rborder <= ( ssize_t )( matchpos[right] - 1 ))
                    :
                  ( left <= p && p <= right ))
                {
                    // omit positions that are unsued or gaps in query
                    if( !protoprofile_->IsUsedAt( p ))
                        continue;

                    if( protoprofile_->GetCountAt( p ) == 1 ||
///                         left < protoprofile_->GetLeftExtentAt( p ))         //was worth to verify but worked a bit worse
                        protoprofile_->GetLeftExtentAt( p ) < left ) {
                        protoprofile_->SetLeftExtentAt( left, p );
                        protoprofile_->SetLeftMSExtentAt( left, p );
                    }

                    if( protoprofile_->GetCountAt( p ) == 1 ||
///                         protoprofile_->GetRightExtentAt( p ) < right )      //was worth to verify but worked a bit worse
                        right < protoprofile_->GetRightExtentAt( p )) {
                        protoprofile_->SetRightExtentAt( right, p );
                        protoprofile_->SetRightMSExtentAt( right, p );
                    }
                }
            }
        }
        //save the intervals of extents for each position
        for( p = 0; p < protoprofile_->GetSize(); p++ ) {
            // omit positions that are unsued or gaps in query
            if( !protoprofile_->IsUsedAt( p ))
                continue;

            size_t  left = protoprofile_->GetLeftExtentAt( p );
            size_t  right = protoprofile_->GetRightExtentAt( p );
            size_t  interval;

            if( right < left || right == MYSIZE_MAX || left == MYSIZE_MAX )
                continue;

//             // adjust left boundary values to properly compute intervals
//             size_t  leftaccpos = left? notused[0][left-1]: 0;
//             size_t  rightaccpos = notused[0][right];
//             interval =  right - left + 1 - ( rightaccpos - leftaccpos );

            left = protoprofile_->GetLeftMSExtentAt( p );
            right = protoprofile_->GetRightMSExtentAt( p );
            interval = msvector[right] - msvector[left] + 1;
            protoprofile_->SetExtentIntervalAt( interval, p );
            protoprofile_->SetMSExtentIntervalAt( interval, p );
//         fprintf( stderr, "%lu %c %lu,%lu  %lu\n", p, DehashCode(protoprofile_->GetResidueAt(p)), left, right, interval );
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( matchpos ) free( matchpos );
    if( msvector ) free( msvector );

    if( mre.isset())
        throw mre;
}





// =========================================================================
// ComputeGlobSequenceWeights: calculate globally sequence weights;
// avgpeseq, average fraction of different residues per one matched position
// TODO: reimplement large arrays to use heap
//
void MSA::ComputeGlobSequenceWeights( float* pavgsq )
{
    const mystring preamb = "MSA::ComputeGlobSequenceWeights: ";
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");
    if( !GetSize())
        throw MYRUNTIME_ERROR( preamb + "No sequences.");

    const size_t    MINUSDLEN = 30;
//     const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;

    const size_t    noseqs = GetSize();
    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const size_t    left = 0;
    const size_t    right = nomstates;
    const size_t    interval = 1;//no avgs: histogram obtained for one column

    size_t  mss;//number of sequences in M state
    size_t  column[NUMALPH];
    size_t  diffsyms;
    size_t  nadiffsyms;
    //size_t  usdlen[noseqs];//length used in weight calculation
	std::unique_ptr<size_t[]> usdlen(new size_t[noseqs]);//length used in weight calculation
	float   avgpeseq;//average fraction of different residues per one match position
    float   w;
    float   gwtsum;//sum of global weights
    //float   gwghts[noseqs];
	std::unique_ptr<float[]> gwghts(new float[noseqs]);

    unsigned char   residue;
    int             pstate;
    size_t  p, i;//position indices


    gwtsum = 0.0f;
    avgpeseq = 0.0f;
    memset( usdlen.get(), 0, sizeof(size_t) * noseqs );
    memset( gwghts.get(), 0, sizeof(float) * noseqs );

    //calculate global weights
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        //reset boundaries
        protoprofile_->SetLeftMSExtentAt( left, p );
        protoprofile_->SetRightMSExtentAt( right, p );
        protoprofile_->SetMSExtentIntervalAt( interval, p );

        pstate = protoprofile_->GetStateAt( p );
        if( !pstate )
            continue;

        mss = 0;
        diffsyms = 0;
        nadiffsyms = 0;
        memset( column, 0, sizeof(size_t) * NUMALPH );

        for( i = 0; i < noseqs; i++ ) {
            // check for valid sequence position
            if( !GetSequenceAt(i)->GetUsed())
                continue;
            if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            residue = GetSequenceAt(i)->GetResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( IsValidResSym( residue )) {
                usdlen[i]++;
                mss++;
            }
            if( column[residue]++ == 0 ) {
                diffsyms++;
                if( IsValidResSym( residue ))
                    nadiffsyms++;
            }
        }

        if( noeffress < nadiffsyms )
            nadiffsyms = noeffress;

        if( pstate )
            if( mss && nadiffsyms )
                //fraction of different residues per one matched position
                avgpeseq += (float)nadiffsyms / (float)mss;

        if( nadiffsyms )
            for( i = 0; i < noseqs; i++ ) {
                // check for valid sequence position
                if( !GetSequenceAt(i)->GetUsed())
                    continue;
                if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                    continue;
                residue = GetSequenceAt(i)->GetResidueAt( p );
                if( !gbTOLXs && residue == X )
                    continue;
                if( IsValidResSym( residue )) {
                    w = 1.0f / (float)(column[residue]*nadiffsyms);
                    gwghts[i] += w;
//                     gwtsum += w;
                }
            }
    }//1:for(p...)

    //average no. different residues per one matched position
    if( avgpeseq && nomstates )
        avgpeseq /= (float)nomstates;
    if( pavgsq )
        *pavgsq = avgpeseq;

    //adjust weights by sequence lengths;
    // each column has total weight =1
    for( i = 0; i < noseqs; i++ )
        if( gwghts[i]) {
            gwghts[i] /= (float)SLC_MAX( usdlen[i] + 1, MINUSDLEN );
            gwtsum += gwghts[i];
        }
    if( !gwtsum )
        throw MYRUNTIME_ERROR( preamb + "Null weights." );
    //normalize weights
    for( i = 0; i < noseqs; i++ ) {
        if( gwghts[i])
            gwghts[i] /= gwtsum;
        if( GetSequenceAt(i))
            GetSequenceAt(i)->SetGlbWeight( gwghts[i]);
    }
}

// -------------------------------------------------------------------------
// CalculateExpNoResiduesAt: calculate expected number of different 
// residues over the extent of the given position 
// TODO: reimplement large arrays to use heap (expnoress)
//
void MSA::CalculateExpNoResiduesAt( 
    size_t p, bool usdgwght, size_t left, size_t right, 
    const float* wghts, float* express, size_t scale )
{
    const mystring preamb = "MSA::CalculateExpNoResiduesAt: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");
    if( !wghts || !express )
        throw MYRUNTIME_ERROR( preamb + "Null arguments.");

    const bool      bMED = true;
    const size_t    noseqs = GetSize();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const size_t    rscale = SLC_MIN( 1000, SLC_MAX( 10, scale ));
    const size_t    maxres = NUMAA;//20;
    const size_t    scsz = maxres * rscale;
    //size_t          expnoress[scsz+1];
	std::unique_ptr<size_t[]> expnoress(new size_t[scsz + 1]);
	size_t          nosinext;//no. sequences in extent
    float           wghtsum[PS_NSTATES] = {0.0f,0.0f,0.0f};
    float           mtwghts[NUMALPH];//match weights
    float           w, ww = 0.0f, wwsum = 0.0f;
    size_t          mss;//number of sequences in M state
    unsigned char   residue;
    int             kstate;
    size_t          k, i, r, n;

    if( usdgwght ) {
        left = 0;
        right = protoprofile_->GetSize();
        if( right )
            right--;
    }

    memset( expnoress.get(), 0, (scsz+1) * sizeof(size_t));
    nosinext = usdgwght?protoprofile_->GetNoSequencesAt( p ):
                        protoprofile_->GetNoSequencesInMSExtentAt( p );

    //iterate over the extent
    for( n = 0, k = left; 0 <= (ssize_t)right && k <= right; k++ ) {
        if( !protoprofile_->IsUsedAt(k))
            continue;
        kstate = protoprofile_->GetStateAt(k);
        if( !kstate )
            continue;
        //use sequence weights calculated for this position
        mss = 0;
        wghtsum[PS_M] = 0.0f;
        memset( mtwghts, 0, sizeof(float) * NUMALPH );
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences absent from the extent
            if( usdgwght? !GetSequenceAt(i)->IsUsedAndInfAt( p ):
                          !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            residue = GetSequenceAt(i)->GetResidueAt( k );
            if( !gbTOLXs && residue == X )
                continue;
            if( IsValidResSym( residue )) {
                mtwghts[residue] += wghts[i];
                wghtsum[PS_M] += wghts[i];
                mss++;
            }
        }
        if( mss < (size_t)( (float)nosinext * gszMINPRCSEQNS ))
            continue;
        if( wghtsum[PS_M]) {
            for( r = 0; r < noress; r++ )
                if( mtwghts[r])
                    mtwghts[r] /= wghtsum[PS_M];
            AdjustWeightsAt((size_t)-1, &mtwghts );//perform adjustment in `mtwghts'
            ww = 0.0f;
            for( r = 0; r < noeffress; r++ ) {
                w = mtwghts[r];
                if( 0.0f < w )
                    ww -= w * logf(w);
            }
            ww = expf( ww );//value in [1,20]
            if( bMED )
                expnoress[ (size_t)rintf(SLC_MIN(maxres,ww)) * rscale ]++;
            else {
                wwsum += ww;
                n++; 
            }
        }
    }//for(k...)
    if( bMED )
        MedianValue( expnoress.get(), scsz, &ww, rscale );
    else if( n )
        ww = wwsum / (float)n;
    *express = ww;
}

// -------------------------------------------------------------------------
// ComputeMIDstateSequenceWeights: calculate MID state sequence weights in 
// one run
// TODO: reimplement large arrays to use heap (gwgflgs, expnors, gwghts, 
// press)
//
void MSA::ComputeMIDstateSequenceWeights()
{
    const mystring preamb = "MSA::ComputeMIDstateSequenceWeights: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const bool      bLGWEIGHTS = true;//use locally global weights
    const size_t    noseqs = GetSize();
    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const float     err = 1.e-6f;

//     size_t  histval;//histogram data value
    size_t  nosinext;//no. sequences in extent
    bool    newset = true;
    size_t  Ihist[NUMALPH];
    size_t  extentIdN[PS_NSTATES] = {0,0,0};
    size_t  column[PS_NSTATES][NUMALPH];
    size_t  diffsyms[PS_NSTATES] = {0,0,0};
    size_t  nadiffsyms[PS_NSTATES] = {0,0,0};
	//     float   gexprss = 0.0f;//expected number of different residues computed globally
	//bool    gwgflgs[nomstates];//flags of using global weights
    //float   expnors[nomstates];//recorded expected number of different residues
    //float   gwghts[noseqs];
	std::unique_ptr<bool[]> gwgflgs(new bool[nomstates]);//flags of using global weights
	std::unique_ptr<float[]> expnors(new float[nomstates]);//recorded expected number of different residues
	std::unique_ptr<float[]> gwghts(new float[noseqs]);
	float   wghtsum[PS_NSTATES] = {0.0f,0.0f,0.0f};
    float*  weights[PS_NSTATES] = {NULL,NULL,NULL};
    bool    usdgwght = true;

    //unsigned char   press[noseqs];
	std::unique_ptr<unsigned char[]> press(new unsigned char[noseqs]);

    float           w, ww;
    float emss, eiss, edss;//expected numbers of observations
    unsigned char   residue;
    int             pstate, kstate;
    char            errbuf[KBYTE];
    size_t  mss = 0;//number of sequences in M state
    size_t  iss, ilen;//the number of sequences in and the number of recurrence/length of I state
    size_t  dss = 0;//number of sequences in D state
    size_t  p, pM, k, i, r;//, rdif;//position indices
    ssize_t ppM, pp, ppp;

    memset( Ihist, 0, sizeof(size_t) * NUMALPH );
    memset( expnors.get(), 0, sizeof(float) * nomstates );
    memset( gwgflgs.get(), 0, sizeof(bool) * nomstates );

    for( i = 0; i < noseqs; i++ ) {
        if( GetSequenceAt(i))
            gwghts[i] = GetSequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0f;
    }

//     //{{calculate exp. no. residues along the entire length
//     for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
//         if( !protoprofile_->GetStateAt(p))
//             continue;
//         if( nomstates>>1 <= ++pp )
//             break;
//     }
//     CalculateExpNoResiduesAt( p, true, 0, protoprofile_->GetSize()-1, gwghts, &gexprss, 100 );
//     //}}

    //1:calculate sequence weights, match weights, expected no. residues
    for( p = 0, pp = ppp = -1, pM = ppM = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );
        if( !pstate )
            continue;

        //{{record residues
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences absent from the extent
            if( !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                    press[i] = GAP;
            else    press[i] = GetSequenceAt(i)->GetResidueAt( p );
        }
//             if( pM )
//                 //save accumulated histogram data of I state
//                 for( rdif = 0; rdif < noress; rdif++ ) {
//                     protoprofile_->SetDistinctHistAt( pM, PS_I, rdif, Ihist[rdif]);
//                     Ihist[rdif] = 0;
//                 }
        pM = p;
        pp++;
        if((ssize_t)nomstates <= pp ) {
            throw MYRUNTIME_ERROR( preamb + "Inconsistent number of match positions.");
        }
        //}}
        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            if( p <= (size_t)ppM ) {
                sprintf( errbuf, "Invalid position %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }

            for( i = 0; i < noseqs; i++ ) {
                if( GetSequenceAt(i)->IsUsedInExtentAt( p ) !=
                    GetSequenceAt(i)->IsUsedInExtentAt( ppM ))
                        break;
                residue = GetSequenceAt(i)->GetResidueAt( ppM );
                if( IsValidResSym( press[i]) != IsValidResSym( residue ))
                    break;
                if( press[i] != residue &&( press[i] == X || residue == X ))
                    break;
            }
            if( i == noseqs ) {
                newset = false;
                protoprofile_->SetSqnWeightsAt( p, PS_M, ppM );
//                 protoprofile_->SetSqnWeightsAt( p, PS_D, ppM );
//                 protoprofile_->SetSqnWeightsAt( p, PS_I, ppM );
                weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );
//                 weights[PS_D] = protoprofile_->GetSqnWeightsAt( p, PS_D );
//                 weights[PS_I] = protoprofile_->GetSqnWeightsAt( p, PS_I );
                if( protoprofile_->GetLeftMSExtentAt(p) != protoprofile_->GetLeftMSExtentAt(ppM) ||
                    protoprofile_->GetRightMSExtentAt(p) != protoprofile_->GetRightMSExtentAt(ppM))
                {
                    newset = true;
                    weights[PS_M] = weights[PS_I] = weights[PS_D] = NULL;
                }
            }

            if( !newset && !( weights[PS_M] /*&& weights[PS_D] && weights[PS_I]*/)) {
//                 newset = true;
                sprintf( errbuf, "Null weights obtained for pos. %d.", (int)(ssize_t)ppM);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }
        }

        if( newset ) {
            ppM = p;
            ppp = pp;

            extentIdN[PS_M] = extentIdN[PS_I] = extentIdN[PS_D] = 0;
            wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0f;

            protoprofile_->NewSqnWeightsAt( p, PS_M, GetSize());
//             protoprofile_->NewSqnWeightsAt( p, PS_D, GetSize());
//             protoprofile_->NewSqnWeightsAt( p, PS_I, GetSize());
            weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );
//             weights[PS_D] = protoprofile_->GetSqnWeightsAt( p, PS_D );
//             weights[PS_I] = protoprofile_->GetSqnWeightsAt( p, PS_I );
            if( !weights[PS_M] /*|| !weights[PS_D] || !weights[PS_I]*/) {
                throw MYRUNTIME_ERROR( preamb + "Null weights." );
            }

            size_t  left = protoprofile_->GetLeftMSExtentAt( p );
            size_t  right = protoprofile_->GetRightMSExtentAt( p );
            size_t  interval = protoprofile_->GetMSExtentIntervalAt( p );
            size_t  nusdposs = 0;

            usdgwght = interval < gszMINEXTINTRV;

            //{{calculate sequence weights
            if( usdgwght ) {
                gwgflgs[pp] = true;
                for( i = 0; i < noseqs; i++ ) {
                    weights[PS_M][i] = gwghts[i];
//                     weights[PS_I][i] = gwghts[i];
//                     weights[PS_D][i] = gwghts[i];
                }
            }
            else {
                nusdposs = 0;
                nosinext = protoprofile_->GetNoSequencesInMSExtentAt( p );
                for( k = left; 0 <= (ssize_t)right && k <= right; k++ ) {
                    //omit unsued and unmatched positions in the extent
                    if( !protoprofile_->IsUsedAt( k ))
                        continue;

                    kstate = protoprofile_->GetStateAt( k );
                    if( !kstate )
                        continue;

                    mss = 0;
                    diffsyms[PS_M] = diffsyms[PS_I] = diffsyms[PS_D] = 0;
                    nadiffsyms[PS_M] = nadiffsyms[PS_I] = nadiffsyms[PS_D] = 0;

                    memset( column[PS_M], 0, sizeof(size_t) * NUMALPH );
                    memset( column[PS_I], 0, sizeof(size_t) * NUMALPH );
                    memset( column[PS_D], 0, sizeof(size_t) * NUMALPH );

                    for( i = 0; i < noseqs; i++ ) {
                        //omit sequences absent from the extent
                        if( !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                            continue;
                        residue = GetSequenceAt(i)->GetResidueAt( k );
                        if( !gbTOLXs && residue == X )
                            continue;
                        if( IsValidResSym( residue ))
                            mss++;
                        if( column[PS_M][residue]++ == 0 ) {
                            diffsyms[PS_M]++;
                            if( IsValidResSym( residue ))
                                nadiffsyms[PS_M]++;
                        }
                    }

                    if( mss < (size_t)( (float)nosinext * gszMINPRCSEQNS ))
                        continue;

                    nusdposs++;
                    if( noeffress < nadiffsyms[PS_M] ) nadiffsyms[PS_M] = noeffress;
//                     if( noeffress < nadiffsyms[PS_D] ) nadiffsyms[PS_D] = noeffress;
//                     if( noeffress < nadiffsyms[PS_I] ) nadiffsyms[PS_I] = noeffress;

                    //update histogram
                    extentIdN[PS_M] += nadiffsyms[PS_M];
//                     extentIdN[PS_D] += nadiffsyms[PS_D];
//                     extentIdN[PS_I] += nadiffsyms[PS_I];
    //                 protoprofile_->IncDistinctHistAt( p, PS_M, nadiffsyms[PS_M]);
    //                 protoprofile_->IncDistinctHistAt( p, PS_D, nadiffsyms[PS_D]);
    //                 protoprofile_->IncDistinctHistAt( p, PS_I, nadiffsyms[PS_I]);

                    if( nadiffsyms[PS_M]) {
                        for( i = 0; i < noseqs; i++ ) {
                            // omit sequences absent from the extent
                            if( !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                                continue;
                            residue = GetSequenceAt(i)->GetResidueAt( k );
                            if( !gbTOLXs && residue == X )
                                continue;
                            if( IsValidResSym( residue )) {
                                w = 1.0f / (float)(column[PS_M][residue] * nadiffsyms[PS_M]);
                                weights[PS_M][i] += w;
                                wghtsum[PS_M] += w;
                            }
                        }
                    }
                }//extent: for(k...)

                if( !wghtsum[PS_M] && gbTOLXs ) {
                    sprintf( errbuf, "Sum of sequence weights is 0 at pos. %d.", (int)(ssize_t)p);
                    throw MYRUNTIME_ERROR( preamb + errbuf );
                }
                //check again
                usdgwght = nusdposs < gszMINEXTINTRV;
                //copy sequence weights
                if( usdgwght ) {
                    gwgflgs[pp] = true;
                    for( i = 0; i < noseqs; i++ ) {
                        weights[PS_M][i] = gwghts[i];
//                         weights[PS_I][i] = gwghts[i];
//                         weights[PS_D][i] = gwghts[i];
                    }
                } else if( wghtsum[PS_M])
                    for( i = 0; i < noseqs; i++ ) {
                        if( weights[PS_M][i])
                            weights[PS_M][i] /= wghtsum[PS_M];
                    }
            }
            //}}
            //{{calculate expected number of different residues
            if( bLGWEIGHTS )
//                 ww = gexprss;
                CalculateExpNoResiduesAt( p, true, left, right, gwghts.get(), &ww, 10 );
            else
                CalculateExpNoResiduesAt( p, usdgwght, left, right, weights[PS_M], &ww, 10 );
            expnors[pp] = SLC_MAX( 1.0f, ww );
            //}}
        }//newset

//         if( newset && !pstate )
//             //accumulate histogram data for I state
//             for( rdif = 0; rdif < noress; rdif++ ) {
//                 histval = protoprofile_->GetDistinctHistAt( p, PS_I, rdif );
//                 Ihist[rdif] += histval;
//             }

        if( !newset ) {
            if( ppp < 0 )
                throw MYRUNTIME_ERROR( preamb + "Invalid source position." );
            //copy expected number of residues
            expnors[pp] = expnors[ppp];
            gwgflgs[pp] = gwgflgs[ppp];
//             //copy histogram data
//             for( rdif = 0; rdif < noress; rdif++ ) {
//                 histval = protoprofile_->GetDistinctHistAt( ppM, PS_M, rdif );
//                 protoprofile_->SetDistinctHistAt( p, PS_M, rdif, histval );
//                 histval = protoprofile_->GetDistinctHistAt( ppM, PS_D, rdif );
//                 protoprofile_->SetDistinctHistAt( p, PS_D, rdif, histval );
//                 histval = protoprofile_->GetDistinctHistAt( ppM, PS_I, rdif );
//                 protoprofile_->SetDistinctHistAt( p, PS_I, rdif, histval );
//             }
        }

        //{{set match weights
        wghtsum[PS_M] = 0.0f;
        for( i = 0; i < noseqs; i++ ) {
            if( usdgwght? !GetSequenceAt(i)->IsUsedAndInfAt( p ):
                          !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            if( !gbTOLXs && press[i] == X )
                continue;
            if( IsValidResSym( press[i] )) {
                wghtsum[PS_M] += weights[PS_M][i];
                protoprofile_->IncMatchWeightsAt( weights[PS_M][i], press[i], p );
            }
        }
        if( !wghtsum[PS_M] && gbTOLXs ) {
            sprintf( errbuf, "Sum of sequence weights is 0 at pos. %d.", (int)(ssize_t)p);
            throw MYRUNTIME_ERROR( preamb + errbuf );
        }
        if( wghtsum[PS_M]) {
            for( r = 0; r < noress; r++ ) {
                w = protoprofile_->GetMatchWeightsAt((unsigned char)r, p );
                if( w )
                    protoprofile_->SetMatchWeightsAt( w/wghtsum[PS_M], (unsigned char)r, p );
            }
        }
        AdjustWeightsAt( p );
        //}}
    }//1:for(p...)

    iss = ilen = 0;
    eiss = 0.0f;
    pM = -1;

    //set expectations at the beginning position
    protoprofile_->SetMIDExpNoObservationsAt( 1.0f, -1, PS_M );
    protoprofile_->SetMIDExpNoObservationsAt( 0.0f, -1, PS_D );
    weights[PS_M] = weights[PS_D] = weights[PS_I] = NULL;

    //2:calculate accumulated weights for each state, expected no. observations
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );

        if( pstate )
            weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );

        mss = dss = 0;
        wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0f;

        //calculate accumulated weights for each state
        for( i = 0; i < noseqs; i++ ) {
            if( !GetSequenceAt(i)->GetUsed())
                continue;
            if( !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            //
            residue = GetSequenceAt(i)->GetResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    wghtsum[PS_M] += gwghts[i];//weights[PS_M][i];
                    mss++;
                }
                else if( residue == GAP ) {
                    wghtsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    wghtsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //calculate expected no. observations (set histogram data)
        if( pstate ) {
            if( !wghtsum[PS_M]) {
                if( !gbTOLXs )
                    wghtsum[PS_M] = 1.0f;
                else {
                    sprintf( errbuf, "Sum of match weights is 0 at pos. %d.", (int)(ssize_t)p);
                    throw MYRUNTIME_ERROR( preamb + errbuf );
                }
            }

            if( iss && ilen ) {
                //average no. observations
                eiss /= (float)iss;
                iss = (size_t)rintf( (float)iss / (float)ilen );
            }

            if( wghtsum[PS_M] < 0.0f || 1.0f+err < wghtsum[PS_M] || eiss < 0.0f || 1.0f+err < eiss ||
                wghtsum[PS_D] < 0.0f || 1.0f+err < wghtsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights at pos. %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }

//             emss = SLC_MIN((float)mss * avgpeseq, (float)noeffress );
//             eiss = SLC_MIN((float)iss * avgpeseq, (float)noeffress );
//             edss = SLC_MIN((float)dss * avgpeseq, (float)noeffress );

//             emss = SLC_MIN( expnors[pp], (float)noeffress );
//             eiss = SLC_MIN( expnors[pp] * eiss / wghtsum[PS_M], (float)noeffress );
//             edss = SLC_MIN( expnors[pp] * wghtsum[PS_D] / wghtsum[PS_M], (float)noeffress );

            //convex function to avoid of rapid expectation drop as weight decreases
            if( bLGWEIGHTS || gwgflgs[pp])
                emss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-wghtsum[PS_M])*logf(wghtsum[PS_M]) ),
                                (float)noeffress);
            else
                emss = SLC_MIN( expnors[pp], (float)noeffress );
            if( eiss )
                eiss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-eiss)*logf(eiss) ),
                                (float)noeffress);
            if( wghtsum[PS_D])
                edss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-wghtsum[PS_D])*logf(wghtsum[PS_D]) ),
                                (float)noeffress);
            else
                edss = 0.0f;

            //for M state, no. different residues is required
            protoprofile_->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0f, emss ));

///             emss = BMProtoProfile::GetExpNoObservations( emss );
///             eiss = BMProtoProfile::GetExpNoObservations( eiss );
///             edss = BMProtoProfile::GetExpNoObservations( edss );

            if((float)( mss + dss ) < emss ) emss = (float)( mss + dss );
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            if((float)( mss + dss ) < edss ) edss = (float)( mss + dss );

            if( emss < 1.0f ) emss = 1.0f;
            if( eiss < 1.0f && eiss ) eiss = 1.0f;
            if( edss < 1.0f && edss ) edss = 1.0f;

            //set expectations directly, avoid using histogram data
            protoprofile_->SetMIDExpNoObservationsAt( emss, (int)p, PS_M );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
            protoprofile_->SetMIDExpNoObservationsAt( edss, (int)p, PS_D );
//             protoprofile_->IncDistinctHistAt( p, PS_M, mss );
//             protoprofile_->IncDistinctHistAt( pM, PS_I, iss );
//             protoprofile_->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0f;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;

    }//2:for(p...)

    if( iss && ilen ) {
        //average no. observations
        eiss /= (float)iss;
        iss = (size_t)rintf( (float)iss / (float)ilen );
        if( /*wghtsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( expnors[nomstates-1] * eiss / wghtsum[PS_M], (float)noeffress );
            eiss = SLC_MIN( expnors[nomstates-1]* expf( gdWCPC*(1.0f-eiss)*logf(eiss) ),
                            (float)noeffress);
///             eiss = BMProtoProfile::GetExpNoObservations( eiss );
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
        }
//         protoprofile_->IncDistinctHistAt( pM, PS_I, (int)rintf( (float)iss * avgpeseq ));
    }
}





// -------------------------------------------------------------------------
// ComputeMIDstateSequenceWeightsNoExtents: calculate MID state sequence 
// weights without using of extents
// TODO: reimplement large arrays to use heap (expnors, gwghts, press)
//
void MSA::ComputeMIDstateSequenceWeightsNoExtents()
{
    const mystring preamb = "MSA::ComputeMIDstateSequenceWeightsNoExtents: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const size_t    noseqs = GetSize();
    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const float     err = 1.e-6f;

    size_t  nosinext;//no. sequences in extent
    bool    newset = true;
    size_t  extentIdN[PS_NSTATES] = {0,0,0};
    size_t  column[PS_NSTATES][NUMALPH];
    size_t  diffsyms[PS_NSTATES] = {0,0,0};
    size_t  nadiffsyms[PS_NSTATES] = {0,0,0};
    //float   expnors[nomstates];//expected number of residues
    //float   gwghts[noseqs];
	std::unique_ptr<float[]> expnors(new float[nomstates]);//expected number of residues
	std::unique_ptr<float[]> gwghts(new float[noseqs]);
	float   wghtsum[PS_NSTATES] = {0.0f,0.0f,0.0f};
    float*  weights[PS_NSTATES] = {NULL,NULL,NULL};
//     float   mtwghts[NUMALPH];//match weights

    //unsigned char   press[noseqs];
	std::unique_ptr<unsigned char[]> press(new unsigned char[noseqs]);

    float           w, ww;
    float emss, eiss, edss;//expected numbers of observations
    unsigned char   residue;
    int             pstate, kstate;
    char			errbuf[KBYTE];
    size_t  mss = 0;//number of sequences in M state
    size_t  iss, ilen;//the number of sequences in and the number of recurrence/length of I state
    size_t  dss = 0;//number of sequences in D state
    size_t  p, pM, k, i, r;//, rdif;//position indices
    ssize_t ppM, pp, ppp;

    memset( expnors.get(), 0, sizeof(float) * nomstates );

    for( i = 0; i < noseqs; i++ ) {
        if( GetSequenceAt(i))
            gwghts[i] = GetSequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0f;
    }

    //1:calculate sequence weights, match weights, expected no. residues
    for( p = 0, pp = ppp = -1, pM = ppM = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );
        if( !pstate )
            continue;

        //{{record residues
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences absent from the extent
            if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                    press[i] = GAP;
            else    press[i] = GetSequenceAt(i)->GetResidueAt( p );
        }
        pM = p;
        pp++;
        if((ssize_t)nomstates <= pp ) {
            throw MYRUNTIME_ERROR( preamb + "Inconsistent number of match positions.");
        }
        //}}
        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            if( p <= (size_t)ppM ) {
                sprintf( errbuf, "Invalid position %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }

            for( i = 0; i < noseqs; i++ ) {
                if( GetSequenceAt(i)->IsUsedAndInfAt(p) != GetSequenceAt(i)->IsUsedAndInfAt(ppM))
                    break;
                residue = GetSequenceAt(i)->GetResidueAt( ppM );
                if( IsValidResSym( press[i]) != IsValidResSym( residue ))
                    break;
            }
            if( i == noseqs ) {
                newset = false;
                protoprofile_->SetSqnWeightsAt( p, PS_M, ppM );
                weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );
            }

            if( !newset && !weights[PS_M]) {
                sprintf( errbuf, "Null weights obtained for pos. %d.", (int)(ssize_t)ppM);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }
        }

        //NEWSET...
        if( newset ) {
            ppM = p;
            ppp = pp;

            extentIdN[PS_M] = extentIdN[PS_I] = extentIdN[PS_D] = 0;
            wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0f;

            protoprofile_->NewSqnWeightsAt( p, PS_M, GetSize());
            weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );
            if( !weights[PS_M]) {
                throw MYRUNTIME_ERROR( preamb + "Null weights." );
            }

            size_t  nusdposs = 0;
            bool    usdgwght = false;

            //{{calculate sequence weights
            nosinext = protoprofile_->GetNoSequencesAt( p );
            for( k = 0; k < protoprofile_->GetSize(); k++ ) {
                //omit unsued and unmatched positions in the extent
                if( !protoprofile_->IsUsedAt( k ))
                    continue;

                kstate = protoprofile_->GetStateAt( k );
                if( !kstate )
                    continue;

                mss = 0;
                diffsyms[PS_M] = diffsyms[PS_I] = diffsyms[PS_D] = 0;
                nadiffsyms[PS_M] = nadiffsyms[PS_I] = nadiffsyms[PS_D] = 0;
                memset( column[PS_M], 0, sizeof(size_t) * NUMALPH );

                for( i = 0; i < noseqs; i++ ) {
                    //omit sequences absent from the extent
                    if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                        continue;
                    residue = GetSequenceAt(i)->GetResidueAt( k );
                    if( !gbTOLXs && residue == X )
                        continue;
                    if( IsValidResSym( residue ))
                        mss++;
                    if( column[PS_M][residue]++ == 0 ) {
                        diffsyms[PS_M]++;
                        if( IsValidResSym( residue ))
                            nadiffsyms[PS_M]++;
                    }
                }

                if( mss < (size_t)( (float)nosinext * gszMINPRCSEQNS ))
                    continue;

                nusdposs++;
                if( noeffress < nadiffsyms[PS_M] ) nadiffsyms[PS_M] = noeffress;

                extentIdN[PS_M] += nadiffsyms[PS_M];

                if( nadiffsyms[PS_M]) {
                    for( i = 0; i < noseqs; i++ ) {
                        // omit sequences absent from the extent
                        if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                            continue;
                        residue = GetSequenceAt(i)->GetResidueAt( k );
                        if( !gbTOLXs && residue == X )
                            continue;
                        if( IsValidResSym( residue )) {
                            w = 1.0f / (float)(column[PS_M][residue] * nadiffsyms[PS_M]);
                            weights[PS_M][i] += w;
                            wghtsum[PS_M] += w;
                        }
                    }
                }
            }//extent: for(k...)

            if( !wghtsum[PS_M] && gbTOLXs ) {
                sprintf( errbuf, "Sum of sequence weights is 0 at pos. %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }
            //check again
            usdgwght = nusdposs < gszMINEXTINTRV;
            //copy sequence weights
            if( usdgwght ) {
                for( i = 0; i < noseqs; i++ ) {
                      weights[PS_M][i] = gwghts[i];
                }
            } else if( wghtsum[PS_M])
                for( i = 0; i < noseqs; i++ ) {
                    if( weights[PS_M][i])
                        weights[PS_M][i] /= wghtsum[PS_M];
                }
            //}}
            //{{calculate expected number of different residues
            CalculateExpNoResiduesAt( p, true, 0, protoprofile_->GetSize()-1, weights[PS_M], &ww, 10 );
            expnors[pp] = SLC_MAX( 1.0f, ww );
            //}}
        }//newset

        if( !newset ) {
            if( ppp < 0 )
                throw MYRUNTIME_ERROR( preamb + "Invalid source position." );
            //copy expected number of residues
            expnors[pp] = expnors[ppp];
        }

        //{{set match weights
        wghtsum[PS_M] = 0.0f;
        for( i = 0; i < noseqs; i++ ) {
            if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            if( !gbTOLXs && press[i] == X )
                continue;
            if( IsValidResSym( press[i] )) {
                wghtsum[PS_M] += weights[PS_M][i];
                protoprofile_->IncMatchWeightsAt( weights[PS_M][i], press[i], p );
            }
        }
        if( !wghtsum[PS_M] && gbTOLXs ) {
            sprintf( errbuf, "Sum of sequence weights is 0 at pos. %d.", (int)(ssize_t)p);
            throw MYRUNTIME_ERROR( preamb + errbuf );
        }
        if( wghtsum[PS_M]) {
            for( r = 0; r < noress; r++ ) {
                w = protoprofile_->GetMatchWeightsAt((unsigned char)r, p );
                if( w )
                    protoprofile_->SetMatchWeightsAt( w/wghtsum[PS_M], (unsigned char)r, p );
            }
        }
        AdjustWeightsAt( p );
        //}}
    }//1:for(p...)

    iss = ilen = 0;
    eiss = 0.0f;
    pM = -1;

    //set expectations at the beginning position
    protoprofile_->SetMIDExpNoObservationsAt( 1.0f, -1, PS_M );
    protoprofile_->SetMIDExpNoObservationsAt( 0.0f, -1, PS_D );
    weights[PS_M] = weights[PS_D] = weights[PS_I] = NULL;

    //2:calculate accumulated weights for each state, expected no. observations
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );

        if( pstate )
            weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );

        mss = dss = 0;
        wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0f;

        //calculate accumulated weights for each state
        for( i = 0; i < noseqs; i++ ) {
            if( !GetSequenceAt(i)->GetUsed())
                continue;
            if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            //
            residue = GetSequenceAt(i)->GetResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    wghtsum[PS_M] += gwghts[i];//weights[PS_M][i];
                    mss++;
                }
                else if( residue == GAP ) {
                    wghtsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    wghtsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //calculate expected no. observations (set histogram data)
        if( pstate ) {
            if( !wghtsum[PS_M]) {
                if( !gbTOLXs )
                    wghtsum[PS_M] = 1.0f;
                else {
                    sprintf( errbuf, "Sum of match weights is 0 at pos. %d.", (int)(ssize_t)p);
                    throw MYRUNTIME_ERROR( preamb + errbuf );
                }
            }

            if( iss && ilen ) {
                //average no. observations
                eiss /= (float)iss;
                iss = (size_t)rintf( (float)iss / (float)ilen );
            }

            if( wghtsum[PS_M] < 0.0f || 1.0f+err < wghtsum[PS_M] || eiss < 0.0f || 1.0f+err < eiss ||
                wghtsum[PS_D] < 0.0f || 1.0f+err < wghtsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights at pos. %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }

//             emss = SLC_MIN((float)mss * avgpeseq, (float)noeffress );
//             eiss = SLC_MIN((float)iss * avgpeseq, (float)noeffress );
//             edss = SLC_MIN((float)dss * avgpeseq, (float)noeffress );

//             emss = SLC_MIN( expnors[pp], (float)noeffress );
//             eiss = SLC_MIN( expnors[pp] * eiss / wghtsum[PS_M], (float)noeffress );
//             edss = SLC_MIN( expnors[pp] * wghtsum[PS_D] / wghtsum[PS_M], (float)noeffress );

            //convex function to avoid of rapid expectation drop as weight decreases
            emss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-wghtsum[PS_M])*logf(wghtsum[PS_M]) ),
                            (float)noeffress);
            if( eiss )
                eiss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-eiss)*logf(eiss) ),
                                (float)noeffress);
            if( wghtsum[PS_D])
                edss = SLC_MIN( expnors[pp] * expf( gdWCPC*(1.0f-wghtsum[PS_D])*logf(wghtsum[PS_D]) ),
                                (float)noeffress);
            else
                edss = 0.0f;

            //for M state, no. different residues is required
            protoprofile_->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0f, emss ));

///             emss = BMProtoProfile::GetExpNoObservations( emss );
///             eiss = BMProtoProfile::GetExpNoObservations( eiss );
///             edss = BMProtoProfile::GetExpNoObservations( edss );

            if((float)( mss + dss ) < emss ) emss = (float)( mss + dss );
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            if((float)( mss + dss ) < edss ) edss = (float)( mss + dss );

            if( emss < 1.0f ) emss = 1.0f;
            if( eiss < 1.0f && eiss ) eiss = 1.0f;
            if( edss < 1.0f && edss ) edss = 1.0f;

            //set expectations directly, avoid using histogram data
            protoprofile_->SetMIDExpNoObservationsAt( emss, (int)p, PS_M );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
            protoprofile_->SetMIDExpNoObservationsAt( edss, (int)p, PS_D );
//             protoprofile_->IncDistinctHistAt( p, PS_M, mss );
//             protoprofile_->IncDistinctHistAt( pM, PS_I, iss );
//             protoprofile_->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0f;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;

    }//2:for(p...)

    if( iss && ilen ) {
        //average no. observations
        eiss /= (float)iss;
        iss = (size_t)rintf( (float)iss / (float)ilen );
        if( /*wghtsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( expnors[nomstates-1] * eiss / wghtsum[PS_M], (float)noeffress );
            eiss = SLC_MIN( expnors[nomstates-1] * expf( gdWCPC*(1.0f-eiss)*logf(eiss) ),
                            (float)noeffress);
///             eiss = BMProtoProfile::GetExpNoObservations( eiss );
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
        }
//         protoprofile_->IncDistinctHistAt( pM, PS_I, (int)rintf( (float)iss * avgpeseq ));
    }
}





// =========================================================================
// ComputeGWMIDstateSequenceWeights: calculate globally sequence weights for
// all states;
// avgpeseq, average fraction of different residues per one matched position
// TODO: reimplement large arrays to use heap (gwghts, ndr, avgndr)
//
void MSA::ComputeGWMIDstateSequenceWeights( float /*avgpeseq*/ )
{
    const mystring preamb = "MSA::ComputeGWMIDstateSequenceWeights: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");
    if( !GetSize())
        throw MYRUNTIME_ERROR( preamb + "No sequences.");

    //for the distribution of positionally normalized weights
    const bool      bPOSNORMAL = false;//positional normalization of seq. weights
    const size_t    NODRLEN = 40;//length for observed no. different residues
    const size_t    NODRLENh = NODRLEN >> 1;
//     const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;

    const size_t    noseqs = GetSize();
    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const float     err = 1.e-6f;

    size_t  mss = 0;//number of sequences in M state
    size_t  iss, ilen;//the number of sequences in and the number of recurrence/length of I state
    size_t  dss = 0;//number of sequences in D state
    float   emss, eiss, edss;//expected numbers of observations
    float   w, ww/*, ndr[nomstates]*/, avglim/*, avgndr[nomstates]*/;
    //float   gwghts[noseqs];
	std::unique_ptr<float[]> ndr(new float[nomstates]);
	std::unique_ptr<float[]> avgndr(new float[nomstates]);
	std::unique_ptr<float[]> gwghts(new float[noseqs]);
	float   poswsum[PS_NSTATES];//positional weight sums
    float*  weights[PS_NSTATES] = {NULL,NULL,NULL};

    unsigned char   residue;
    bool            sused, extused;
    int             pstate;
    char            errbuf[KBYTE];
    size_t  p, pp, i, r, n, nn;//position indices
    ssize_t pM;

    for( i = 0; i < noseqs; i++ ) {
        if( GetSequenceAt(i))
            gwghts[i] = GetSequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0f;
    }

    //distribute global weights at all positions
    if( !bPOSNORMAL ) {
        protoprofile_->NewSqnWeightsAt( 0, PS_M, GetSize());
        protoprofile_->NewSqnWeightsAt( 0, PS_I, GetSize());
        protoprofile_->NewSqnWeightsAt( 0, PS_D, GetSize());
        weights[PS_M] = protoprofile_->GetSqnWeightsAt( 0, PS_M );
        weights[PS_I] = protoprofile_->GetSqnWeightsAt( 0, PS_I );
        weights[PS_D] = protoprofile_->GetSqnWeightsAt( 0, PS_D );
        if( !weights[PS_M] || !weights[PS_I] || !weights[PS_D]) {
            throw MYRUNTIME_ERROR( preamb + "Null weights.");
        }

        for( i = 0; i < noseqs; i++ ) {
            weights[PS_M][i] = gwghts[i];
            weights[PS_I][i] = gwghts[i];
            weights[PS_D][i] = gwghts[i];
        }
    }

    memset( avgndr.get(), 0, sizeof(float) * nomstates );
    avglim = 0.0f;
    nn = NODRLENh;

    //2.adjust and recalculate weights, set MID weights
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );

        if( bPOSNORMAL ) {
            if( pstate ) {
                protoprofile_->NewSqnWeightsAt( p, PS_M, GetSize());
                protoprofile_->NewSqnWeightsAt( p, PS_D, GetSize());

                weights[PS_M] = protoprofile_->GetSqnWeightsAt( p, PS_M );
                weights[PS_D] = protoprofile_->GetSqnWeightsAt( p, PS_D );
                if( !weights[PS_M] || !weights[PS_D]) {
                    throw MYRUNTIME_ERROR( preamb + "Null weights.");
                }
            } else {
                protoprofile_->NewSqnWeightsAt( p, PS_I, GetSize());
                weights[PS_I] = protoprofile_->GetSqnWeightsAt( p, PS_I );
                if( !weights[PS_I]) {
                    throw MYRUNTIME_ERROR( preamb + "Null weights.");
                }
            }
        }//not bPOSNORMAL
        else {
            if( p ) {
                protoprofile_->SetSqnWeightsAt( p, PS_M, 0 );
                protoprofile_->SetSqnWeightsAt( p, PS_I, 0 );
                protoprofile_->SetSqnWeightsAt( p, PS_D, 0 );
            }
        }

        poswsum[PS_M] = poswsum[PS_I] = poswsum[PS_D] = 0;

        //sequence weights at each position are normalized by definition
        for( i = 0; i < noseqs; i++ ) {
            if( !GetSequenceAt(i)->GetUsed())
                continue;
            sused = GetSequenceAt(i)->IsUsedAndInfAt( p );
            extused = GetSequenceAt(i)->IsUsedInExtentAt( p );
            // check for valid sequence position and 
            // modify the flag of belonging to the extent
            if( !sused ) {
                if( extused )
                    GetSequenceAt(i)->UnsetUsedInExtentAt( p );
                continue;
            }
            if( !extused )
                GetSequenceAt(i)->SetUsedInExtentAt( p );
            //
            protoprofile_->IncNoSequencesInMSExtentAt( p );
            //
            residue = GetSequenceAt(i)->GetResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    poswsum[PS_M] += gwghts[i];
                }
                else if( residue == GAP )
                    if( bPOSNORMAL )
                        poswsum[PS_D] += gwghts[i];
            } else
                if( IsValidResSym( residue ))
                    if( bPOSNORMAL )
                        poswsum[PS_I] += gwghts[i];
        }

        if( pstate ) {
            if( !poswsum[PS_M]) {
                if( !gbTOLXs )
                    poswsum[PS_M] = 1.0f;
                else {
                    sprintf( errbuf, "Sum of match weights is 0 at pos. %d.", (int)(ssize_t)p);
                    throw MYRUNTIME_ERROR( preamb + errbuf );
                }
            }
            //set match weights
            for( i = 0; i < noseqs; i++ ) {
                if( !GetSequenceAt(i)->GetUsed())
                    continue;
                if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                    continue;
                residue = GetSequenceAt(i)->GetResidueAt( p );
                if( !gbTOLXs && residue == X )
                    continue;
                if( bPOSNORMAL ) {
                    if( IsValidResSym( residue ))
                        weights[PS_M][i] = gwghts[i] / poswsum[PS_M];
                    else if( residue == GAP && poswsum[PS_D])
                        weights[PS_D][i] = gwghts[i] / poswsum[PS_D];
                    if( IsValidResSym( residue ))
                        protoprofile_->IncMatchWeightsAt( weights[PS_M][i], residue, p );
                } else // !bPOSNORMAL
                    if( IsValidResSym( residue ))
                        protoprofile_->IncMatchWeightsAt( weights[PS_M][i]/poswsum[PS_M], residue, p );
            }
        } else if( bPOSNORMAL ) {
            if( poswsum[PS_I]) {
                for( i = 0; i < noseqs; i++ ) {
                    if( !GetSequenceAt(i)->GetUsed())
                        continue;
                    if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                        continue;
                    residue = GetSequenceAt(i)->GetResidueAt( p );
                    if( !gbTOLXs && residue == X )
                        continue;
                    if( IsValidResSym( residue ))
                        weights[PS_I][i] = gwghts[i] / poswsum[PS_I];
                }
            }
        }

        //calculate observed number of different residues
        if( pstate ) {
            ww = 0.0f;
            for( r = 0; r < noeffress; r++ ) {
                w = protoprofile_->GetMatchWeightsAt((unsigned char)r, p );
                if( 0.0f < w )
                    ww -= w * logf(w);
            }
            //exp(ww) varies between 1 and 20
            ww = expf( ww );
            ndr[pp] = ww;
            avglim += ww;
            pp++;
            if( NODRLEN <= pp ) {
                avgndr[nn] = avglim / (float)NODRLEN;
                if( NODRLEN == pp )
                    for( n = 0; n < nn; n++ )
                        avgndr[n] = avgndr[nn];
                if( nomstates == pp )
                    for( n = nn+1; n < nomstates; n++ )
                        avgndr[n] = avgndr[nn];
                avglim -= ndr[pp-NODRLEN];
                nn++;
            }
        }
    }//2.for(p...)

    if( pp < NODRLEN && nomstates ) {
        avgndr[0] = avglim / (float)nomstates;
        for( n = 1; n < nomstates; n++ )
            avgndr[n] = avgndr[0];
    }
    iss = ilen = 0;
    eiss = 0.0f;
    pM = -1;

    //set expectations at the beginning state
    protoprofile_->SetMIDExpNoObservationsAt( 1.0f, -1, PS_M );
    protoprofile_->SetMIDExpNoObservationsAt( 0.0f, -1, PS_D );

    //3.set histogram data for MID states
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unsued positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;

        pstate = protoprofile_->GetStateAt( p );

        mss = dss = 0;
        poswsum[PS_M] = poswsum[PS_I] = poswsum[PS_D] = 0;

        //sequence weights at each position are normalized by definition
        for( i = 0; i < noseqs; i++ ) {
            if( !GetSequenceAt(i)->GetUsed())
                continue;
            if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            residue = GetSequenceAt(i)->GetResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    poswsum[PS_M] += gwghts[i];
                    mss++;
                }
                else if( residue == GAP ) {
                    poswsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    poswsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //set histogram data
        if( pstate ) {
            if( !poswsum[PS_M]) {
                if( !gbTOLXs )
                    poswsum[PS_M] = 1.0f;
                else {
                    sprintf( errbuf, "Sum of match weights is 0 at pos. %d.", (int)(ssize_t)p);
                    throw MYRUNTIME_ERROR( preamb + errbuf );
                }
            }

            if( iss && ilen ) {
                //average no. observations
                eiss /= (float)iss;
                iss = (size_t)rintf( (float)iss / (float)ilen );
            }

            if( poswsum[PS_M] < 0.0f || 1.0f+err < poswsum[PS_M] || eiss < 0.0f || 1.0f+err < eiss ||
                poswsum[PS_D] < 0.0f || 1.0f+err < poswsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights at pos. %d.", (int)(ssize_t)p);
                throw MYRUNTIME_ERROR( preamb + errbuf );
            }

//             emss = SLC_MIN((float)mss * avgpeseq, (float)noeffress );
//             eiss = SLC_MIN((float)iss * avgpeseq, (float)noeffress );
//             edss = SLC_MIN((float)dss * avgpeseq, (float)noeffress );

//             emss = SLC_MIN( avgndr[pp], (float)noeffress );
//             eiss = SLC_MIN( avgndr[pp] * eiss / poswsum[PS_M], (float)noeffress );
//             edss = SLC_MIN( avgndr[pp] * poswsum[PS_D] / poswsum[PS_M], (float)noeffress );

            //convex function to avoid of rapid expectation drop as weight decreases
            emss = SLC_MIN( avgndr[pp] * expf( gdWCPC*(1.0f-poswsum[PS_M])*logf(poswsum[PS_M]) ),
                            (float)noeffress);
            if( eiss )
                eiss = SLC_MIN( avgndr[pp] * expf( gdWCPC*(1.0f-eiss)*logf(eiss)),
                                (float)noeffress);
            if( poswsum[PS_D])
                edss = SLC_MIN( avgndr[pp] * expf( gdWCPC*(1.0f-poswsum[PS_D])*logf(poswsum[PS_D]) ),
                                (float)noeffress);
            else
                edss = 0.0f;

            //for M state, no. different residues is required
            protoprofile_->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0f, emss ));

            emss = BMProtoProfile::GetExpNoObservations( emss );///NOTE:check commenting this again
            eiss = BMProtoProfile::GetExpNoObservations( eiss );///NOTE:check commenting this again
            edss = BMProtoProfile::GetExpNoObservations( edss );///NOTE:check commenting this again

            if((float)( mss + dss ) < emss ) emss = (float)( mss + dss );
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            if((float)( mss + dss ) < edss ) edss = (float)( mss + dss );

            if( emss < 1.0f ) emss = 1.0f;
            if( eiss < 1.0f && eiss ) eiss = 1.0f;
            if( edss < 1.0f && edss ) edss = 1.0f;

            //set expectations directly, avoid using histogram data
            protoprofile_->SetMIDExpNoObservationsAt( emss, (int)p, PS_M );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
            protoprofile_->SetMIDExpNoObservationsAt( edss, (int)p, PS_D );
//             protoprofile_->IncDistinctHistAt( p, PS_M, mss );
//             protoprofile_->IncDistinctHistAt( pM, PS_I, iss );
//             protoprofile_->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0f;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;
    }//3:for(p...)

    if( iss && ilen ) {
        //average no. observations
        eiss /= (float)iss;
        iss = (size_t)rintf( (float)iss / (float)ilen );
        if( /*poswsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( avgndr[nomstates-1] * eiss / poswsum[PS_M], (float)noeffress );
            eiss = SLC_MIN( avgndr[nomstates-1] * expf( gdWCPC*(1.0f-eiss)*logf(eiss) ),
                            (float)noeffress);
            eiss = BMProtoProfile::GetExpNoObservations( eiss );///NOTE:check commenting this again
            if((float)( mss + dss ) < eiss ) eiss = (float)( mss + dss );
            protoprofile_->SetMIDExpNoObservationsAt( eiss, (int)pM, PS_I );
        }
//         protoprofile_->IncDistinctHistAt( pM, PS_I, (int)rintf( (float)iss * avgpeseq ));
    }
}





// =========================================================================
// CalcTW: calculate combined transition weight given two weights of 
// adjacent positions
//
float CalcTW( float wp, float w )
{
    const bool clbPROD = false;
    float result = wp;
    if( clbPROD && wp && w )
        result = wp * w;
    return result;
}

// -------------------------------------------------------------------------
// ComputeTransitionFrequencies: compute observed weighted transition 
// frequencies for each position
// TODO: reimplement large arrays to use heap (gwghts)
//
void MSA::ComputeTransitionFrequencies( bool usegwgts, bool expmids )
{
    const mystring preamb = "MSA::ComputeTransitionFrequencies: ";
    myruntime_error mre;
    const size_t    noseqs = GetSize();
    const int       maxtr = gTPTRANS_NTPS;// number of transitions per state
    int             states, st;         // state values
    //float           gwghts[noseqs];     // sequence weights computed globally
	std::unique_ptr<float[]> gwghts(new float[noseqs]);// sequence weights computed globally
	float*          Mweights = NULL;    // M state sequence weights at a position
    float*          Iweights = NULL;    // I state sequence weights
    float*          Dweights = NULL;    // D state sequence weights
    float*          ppmMwghs = NULL;    // M state sequence weights at previous position
    float*          ppmIwghs = NULL;    // I state sequence weights at previous position
    float*          ppmDwghs = NULL;    // D state sequence weights at previous position
    float           wM, wI, wD;         // weights of single sequence at a position
    float           wppmM, wppmI, wppmD;
    size_t*         inserts = NULL;     // lengths of insertions between adjacent match states
    unsigned char   residue, ppmres;
    size_t          p, i;//, nosx;
    ssize_t         ppp, ppm;           //indices of previous positions
    size_t          ppm_t;              //indices of previous positions
    int             pppstate;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    inserts = (size_t*)malloc( sizeof(size_t) * (noseqs+1));

     if( inserts == NULL )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory." );

    try {
        memset( inserts, 0, sizeof(size_t) * (noseqs+1));

        ppp = ppm = -1;
        pppstate = 1;//match

        if( usegwgts ) {
            for( i = 0; i < noseqs; i++ ) {
                if( GetSequenceAt(i))
                    gwghts[i] = GetSequenceAt(i)->GetGlbWeight();
                else
                    gwghts[i] = 0.0f;
            }
            Mweights = Iweights = Dweights = gwghts.get();
            ppmMwghs = ppmIwghs = ppmDwghs = gwghts.get();
        }

        // iterate over all query positions
        for( p = 0; p < protoprofile_->GetSize(); p++ ) {
            // omit unsued positions
            if( !protoprofile_->IsUsedAt( p ))
                continue;

            if( !usegwgts ) {
                Mweights = Iweights = Dweights = protoprofile_->GetSqnWeightsAt( p, PS_M );
//                 Iweights = protoprofile_->GetSqnWeightsAt( p, PS_I );
//                 Dweights = protoprofile_->GetSqnWeightsAt( p, PS_D );
            }

            // weights may be not computed if there were no extents for the position
//             if( weights == NULL )
//                 continue;

            ppm_t = p;
            if( 0 <= ppm )
                ppm_t = ppm;

            if( !usegwgts ) {
                ppmMwghs = ppmIwghs = ppmDwghs = protoprofile_->GetSqnWeightsAt( ppm_t, PS_M );
//                 ppmIwghs = protoprofile_->GetSqnWeightsAt( ppm_t, PS_I );
//                 ppmDwghs = protoprofile_->GetSqnWeightsAt( ppm_t, PS_D );
            }

//             if( ppmwghs == NULL )
//                 continue;
//             nosx = protoprofile_->GetNoSequencesInMSExtentAt( ppm_t );
//             if( !nosx )
//                 continue;

            for( i = 0; i < noseqs; i++ ) {
                if( !GetSequenceAt(i)->GetUsed())
                    continue;
//                 if( !GetSequenceAt(i)->IsUsedInExtentAt( p ))
                if( !GetSequenceAt(i)->IsUsedAndInfAt( p ))
                    continue;
                residue = GetSequenceAt(i)->GetResidueAt( p );

                wM = wI = wD = 0.0f;
                if( Mweights ) wM = Mweights[i];
                if( Iweights ) wI = Iweights[i];
                if( Dweights ) wD = Dweights[i];

                //{{PROCESS TRANSITION WEIGHTS
                ppmres = X;//match

                wppmM = wppmI = wppmD = 0.0;
                if( ppmMwghs ) wppmM = ppmMwghs[i];
                if( ppmIwghs ) wppmI = ppmIwghs[i];
                if( ppmDwghs ) wppmD = ppmDwghs[i];

                if( 0 <= ppm ) {
//                     if( GetSequenceAt(i)->IsUsedInExtentAt( ppm ))
                    if( GetSequenceAt(i)->IsUsedAndInfAt( ppm ))
                        ppmres = GetSequenceAt(i)->GetResidueAt( ppm );
                    else
                        continue;//require both positions to be in use
                }

                if( protoprofile_->GetStateAt( p )) {
                    if( pppstate ) {
                        if( IsValidResSym( residue )) {
                            if( IsValidResSym( ppmres )) {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmM, wM ), P_MM, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmM, (int)ppm, PS_M );
                            } else {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmD, wM ), P_DM, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmD, (int)ppm, PS_D );
                            }
                        }
                        else {
                            if( IsValidResSym( ppmres )) {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmM, wD ), P_MD, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmM, (int)ppm, PS_M );
                            } else {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmD, wD ), P_DD, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmD, (int)ppm, PS_D );
                            }
                        }
                    }
                    else {
                        if( IsValidResSym( residue )) {
                            if( 0 < inserts[i]) {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmI, wM ), P_IM, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmI, (int)ppm, PS_I );
                            } else {
                                if( IsValidResSym( ppmres )) {
                                    protoprofile_->IncTransWeightsAt( CalcTW( wppmM, wM ), P_MM, ppm );
                                    if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmM, (int)ppm, PS_M );
                                } else {
                                    protoprofile_->IncTransWeightsAt( CalcTW( wppmD, wM ), P_DM, ppm );
                                    if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmD, (int)ppm, PS_D );
                                }
                            }
                        }
                        else {
                            if( 0 < inserts[i]) {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmI, wD ), P_ID, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmI, (int)ppm, PS_I );
                            } else {
                                if( IsValidResSym( ppmres )) {
                                    protoprofile_->IncTransWeightsAt( CalcTW( wppmM, wD ), P_MD, ppm );
                                    if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmM, (int)ppm, PS_M );
                                } else {
                                    protoprofile_->IncTransWeightsAt( CalcTW( wppmD, wD ), P_DD, ppm );
                                    if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmD, (int)ppm, PS_D );
                                }
                            }
                        }
                    }
                }
                else {
                    if( IsValidResSym( residue ))
                        inserts[i]++;

                    if( IsValidResSym( residue )) {
                        if( 1 < inserts[i]) {
                            protoprofile_->IncTransWeightsAt( CalcTW( wppmI, wI ), P_II, ppm );
                            if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmI, (int)ppm, PS_I );
                        } else {//the first time when an insertion residue is met
                            if( IsValidResSym( ppmres )) {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmM, wI ), P_MI, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmM, (int)ppm, PS_M );
                            } else {
                                protoprofile_->IncTransWeightsAt( CalcTW( wppmD, wI ), P_DI, ppm );
                                if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wppmD, (int)ppm, PS_D );
                            }
                        }
                    }
                }
                //}}PROCESS TRANSITION WEIGHTS
            }

            ppp = p;//last position
            pppstate = protoprofile_->GetStateAt( p );
            if( pppstate ) {
                //reset inserts: this is match state
                memset( inserts, 0, sizeof(size_t) * (noseqs+1));
                ppm = p;//last match state position
            }
        }

        //{{PROCESS END-STATE
        //to assign transition probabilities to the end state,
        //back-iterate over query positions down to a match state
        for( ppp = protoprofile_->GetSize() - 1; 0 <= ppp ; ppp-- ) {
            // omit unsued positions
            if( !protoprofile_->IsUsedAt( ppp ))
                continue;

            bool state = protoprofile_->GetStateAt( ppp );
            if( !state ) {
                if( 0 <= ppm ) {
                    if( ppm <= ppp ) {
                        if( !protoprofile_->IsUsedAt( ppm ))
                            warning(( preamb + 
                            "End-state transitions: Unused last match state position." ).c_str());
                        else
                            ppp = ppm;
                    } else {
                        warning(( preamb + 
                        "End-state transitions: Invalid match state position." ).c_str());
                        break;
                    }
                }
            }

            if( !usegwgts ) {
                Mweights = Iweights = Dweights = protoprofile_->GetSqnWeightsAt( ppp, PS_M );
//                 Iweights = protoprofile_->GetSqnWeightsAt( ppp, PS_I );
//                 Dweights = protoprofile_->GetSqnWeightsAt( ppp, PS_D );
            }

            if( Mweights == NULL && Iweights == NULL && Dweights == NULL ) {
                if( state ) {
                    warning(( preamb + 
                    "End-state transitions: Null sequence weights." ).c_str());
                    break;
                }
                continue;
            }

            for( i = 0; i < noseqs; i++ ) {
                if( !GetSequenceAt(i)->GetUsed())
                    continue;
//                 if( !GetSequenceAt(i)->IsUsedInExtentAt( ppp ))
                if( !GetSequenceAt(i)->IsUsedAndInfAt( ppp ))
                    continue;
                residue = GetSequenceAt(i)->GetResidueAt( ppp );

                wM = wI = wD = 0.0f;
                if( Mweights ) wM = Mweights[i];
                if( Iweights ) wI = Iweights[i];
                if( Dweights ) wD = Dweights[i];

                //{{PROCESS TRANSITION WEIGHTS
                if( protoprofile_->GetStateAt( ppp )) {
                    if( IsValidResSym( residue )) {
                        protoprofile_->IncTransWeightsAt( CalcTW( wM, wM ), P_MM, ppp );
                        if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wM, (int)ppp, PS_M );
                    } else {
                        protoprofile_->IncTransWeightsAt( CalcTW( wD, wM ), P_DM, ppp );
                        if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wD, (int)ppp, PS_D );
                    }
                }
                else {
                    if( 0 < inserts[i]) {
                        protoprofile_->IncTransWeightsAt( CalcTW( wI, wM ), P_IM, ppp );
                        if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wI, (int)ppp, PS_I );
                    } else {
                        if( IsValidResSym( residue ) || ppm < 0 ) {
                            protoprofile_->IncTransWeightsAt( CalcTW( wM, wM ), P_MM, ppm < 0? ppm: ppp );
                            if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wM, (int)(ppm < 0? ppm: ppp), PS_M );
                        } else {
                            protoprofile_->IncTransWeightsAt( CalcTW( wD, wM ), P_DM, ppp );
                            if( expmids ) protoprofile_->IncMIDExpNoObservationsAt( wD, (int)ppp, PS_D );
                        }
                    }
                }
                //}}PROCESS TRANSITION WEIGHTS
            }
            break;//leave the loop
        }
        //}}PROCESS END

        //NORMALIZE Transition weights and observed state frequencies
        for( ppp = -1; ppp < (ssize_t)protoprofile_->GetSize(); ppp++ ) {
            // omit unsued positions
            if( 0 <= ppp && !protoprofile_->IsUsedAt( ppp ))
                continue;
            if( 0 <= ppp && !protoprofile_->GetStateAt( ppp ))
                continue;

            //process MID states
            for( states = 0; states < P_NSTATES; states += maxtr ) {
                wM = 0.0f;//for norm.
                for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
//!                   if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
                    wM += protoprofile_->GetTransWeightsAt( st, ppp );
                }
                if( wM )
                    for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
                        //if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
                        wppmM = protoprofile_->GetTransWeightsAt(st,ppp) / wM;
                        protoprofile_->SetTransWeightsAt( wppmM, st, ppp );
                    }
            }
            if( expmids ) {
                //normalize observed state frequencies
                wM = 0.0f;//for norm.
                for( st = 0; st < PS_NSTATES; st++ )
                    wM += protoprofile_->GetMIDExpNoObservationsAt((int)ppp, st );
                if( wM )
                    for( st = 0; st < PS_NSTATES; st++ ) {
                        wppmM = protoprofile_->GetMIDExpNoObservationsAt((int)ppp,st) / wM;
                        protoprofile_->SetMIDExpNoObservationsAt( wppmM, (int)ppp, st );
                    }
            }
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    free( inserts );

    if( mre.isset())
        throw mre;

// // protoprofile_->PrintTransWeights( stderr );
// // protoprofile_->PrintMIDExpNoObservations( stderr );
}





// =========================================================================
// ComputeTargetTransFrequencies: compute target transition frequencies for
// each position
//
void MSA::ComputeTargetTransFrequencies()
{
    const mystring preamb = "MSA::ComputeTargetTransFrequencies: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const float  ( *ptargs )[P_NSTATES];
    float           trns[P_NSTATES];
//     float           expn;
    float           nm, ni, nd;
    ssize_t         ppp;

    if((ssize_t)protoprofile_->GetSize() < 1 )
        return;

    trns[P_MM] = 1.0f; trns[P_MI] = 0.0f; trns[P_MD] = 0.0f;
    trns[P_IM] = 1.0f; trns[P_II] = 0.0f; trns[P_ID] = 0.0f;
    trns[P_DM] = 1.0f; trns[P_DI] = 0.0f; trns[P_DD] = 0.0f;

    protoprofile_->SetTargetTranstAt( &trns, ppp = -1 );

    // iterate over all query positions; the beginning state is -1
    for( ppp = 0; ppp < (ssize_t)protoprofile_->GetSize(); ppp++ ) {
        if( 0 <= ppp && !protoprofile_->IsUsedAt( ppp ))
            continue;
        if( 0 <= ppp && !protoprofile_->GetStateAt( ppp ))
            continue;

//         if( ppp < 0 )
//               expn = protoprofile_->GetExpNoObservationsAt( 0 );
//         else  expn = protoprofile_->GetExpNoObservationsAt( ppp );

        nm = protoprofile_->GetMIDExpNoObservationsAt((int)ppp, PS_M );
        ni = protoprofile_->GetMIDExpNoObservationsAt((int)ppp, PS_I );
        nd = protoprofile_->GetMIDExpNoObservationsAt((int)ppp, PS_D );

        //one exp. sequence (observation) bears no inform.
        nm = SLC_MAX( 0.0f, nm - 1.0f );

        TRANSPROBS.SetEffNoSequences( nm, ni, nd );
        TRANSPROBS.PME( protoprofile_->GetTransWeightsAt( ppp ));//conserv. there
        ptargs = TRANSPROBS.GetPMEstimators();
        if( ptargs == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null posterior mean estimates.");
        protoprofile_->SetTargetTranstAt( ptargs, ppp );
    }

    if( nomstates )
        protoprofile_->SetTargetTranstAt( &trns, ppp = (ssize_t)nomstates-1 );

// // protoprofile_->PrintTargetTranst( stderr );
}





// =========================================================================
// ExpectedNoObservationsAt: calculate expected number of independent
//  observations; 
//
float MSA::ExpectedNoObservationsAt( size_t pos, int st ) const
{
    const mystring preamb = "MSA::ExpectedNoObservationsAt: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    if( protoprofile_->GetSize() <= pos || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");

    if( !protoprofile_->IsUsedAt( pos ))
        return 0.0f;

    ssize_t interval = (ssize_t)protoprofile_->GetMSExtentIntervalAt( pos );
    ssize_t noeffres = NUMAA;
    ssize_t halfint = 0;
    ssize_t estnocols;  //number of columns having a certain number of different residues
    ssize_t sumnocols;  //accumulated number of columns
    ssize_t distinct;   //accumulated number of different residues over half of the extent
    float   avgdistinct;//average number of different residues
    float   expnobs;    //expected number of observations
    ssize_t d;

    if( interval < 0 || (ssize_t)protoprofile_->GetEffectiveSize() < interval ) {
        warning(( preamb + "Extent interval is out of range." ).c_str());
        return 0.0f;
    }
    if( st == PS_I && interval < 1 )
        return 0.0f;
    if( interval < 1 ) {
        warning(( preamb + "Extent interval is 0." ).c_str());
        return 0.0f;
    }

    halfint = ( interval + 1 ) >> 1;

    sumnocols = 0;
    distinct = 0;

    //calculate average number of different residues over half of the extent,
    //starting with the most diverse columns
    for( d = noeffres; d && sumnocols < halfint; d-- ) {
        estnocols = protoprofile_->GetDistinctHistAt( pos, st, (unsigned char)d );
        sumnocols += estnocols;
        distinct += estnocols * d;

        if( halfint < sumnocols ) {
            distinct -= ( sumnocols - halfint ) * d;
            sumnocols = halfint;
        }
    }

    if( sumnocols < 1 )
//         throw MYRUNTIME_ERROR( preamb + "No observed counts." );
        return 0.0f;

    avgdistinct = (float)distinct / (float)sumnocols;

    //to calculate expected number of independent observations
    //corresponding to the average number of different residues,
    //use the Altschul et al. method as in NAR, 2009

    expnobs = BMProtoProfile::GetExpNoObservations( avgdistinct );

    if( protoprofile_->GetNoSequencesInMSExtentAt( pos ) < expnobs )
        expnobs = (float)protoprofile_->GetNoSequencesInMSExtentAt( pos );

    //gaps are not included in calculating the histogram of different residues
    if( expnobs < 1.0f )
        expnobs = 1.0f;

    return expnobs;
}

// -------------------------------------------------------------------------
// DeriveExpectedNoObservations: calculate expected numbers of 
// observations for each position
//
void MSA::DeriveExpectedNoObservations()
{
    const mystring preamb = "MSA::DeriveExpectedNoObservations: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    float   wm;//, wi, wd;//MID state weights 
    float   expn[PS_NSTATES];
    size_t  p;

    //iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if(!protoprofile_->IsUsedAt( p ))
            continue;
        if(!protoprofile_->GetStateAt( p ))
            continue;

        //MID state observed frequencies are precalculated
        wm = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_M );
//         wi = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_I );
//         wd = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_D );
        if( wm <= 0.0f ) {
            if( protoprofile_->GetStateAt( p ))
                warning(( preamb + "Expected number of observations is 0." ).c_str());
            wm = 1.0f;
        }
        expn[PS_M] = expn[PS_I] = expn[PS_D] = 0.0f;
        expn[PS_M] = ExpectedNoObservationsAt( p, PS_M );
        expn[PS_I] = ExpectedNoObservationsAt( p, PS_I );
        expn[PS_D] = ExpectedNoObservationsAt( p, PS_D );
        protoprofile_->SetExpNoObservationsAt( expn[PS_M], p );
        //overwrite MID observed frequencies with MID state expected observations
        protoprofile_->SetMIDExpNoObservationsAt( expn[PS_M], (int)p, PS_M );
        protoprofile_->SetMIDExpNoObservationsAt( expn[PS_I], (int)p, PS_I );
        protoprofile_->SetMIDExpNoObservationsAt( expn[PS_D], (int)p, PS_D );
//         protoprofile_->SetMIDExpNoObservationsAt( expn, (int)p, PS_M );
//         protoprofile_->SetMIDExpNoObservationsAt( expn * wi / wm, (int)p, PS_I );
//         protoprofile_->SetMIDExpNoObservationsAt( expn * wd / wm, (int)p, PS_D );
    }
}

// -------------------------------------------------------------------------
// MedianValue: calculate median value of the given distribution;
// scale, histogram values scaled by this value
//
void MedianValue( size_t* histdata, size_t szdat, float* median, size_t scale )
{
    if( histdata == NULL )
        return;

    size_t  sumoccs = 0;//sum of occurences
    size_t  novals = 0;//number of values
    size_t  hnovals = 0;//half of that number
    size_t  mass = 0;
    float   medn = 0.0f;
    size_t  e;

    for( e = 0; e <= szdat; e++ ) {
        novals += histdata[e];
        mass += histdata[e] * e;
    }

    hnovals = novals >> 1;
    if( hnovals && ( novals & 1 ))
        hnovals++;//1-based further sum

    for( e = 0; e <= szdat; e++ ) {
        sumoccs += histdata[e];
        if( hnovals <= sumoccs ) {
            medn = (float)e;//median value
            //if the size is even and it is a boundary of bins or bin contains <2 elems.
            if( e < szdat  && !( novals & 1 ) && ( hnovals == sumoccs || histdata[e] < 2 )) { 
                for( e++; e <= szdat && !histdata[e]; e++ );
                if( e <= szdat )
                    //take average value of adjacent expected values
                    medn = ( medn + (float)e ) * 0.5f;
            }
            break;
        }
    }
    if( medn <= 0.0f && novals )
        medn = (float)mass / (float)novals;//average then
    medn /= (float)scale;
    if( median )
        *median = medn;
}

// -------------------------------------------------------------------------
// CalculateEffNoSequences: assign effective number of sequences to
//  expected number of observations calculated over all positions
// TODO: reimplement (potentially) large arrays to use heap (expnhist)
//
void MSA::CalculateEffNoSequences()
{
    const mystring preamb = "MSA::CalculateEffNoSequences: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    float   expn[PS_NSTATES];
    float   effnos = 0.0f;
    bool    pstate;
    size_t  p;

    const int   maxexpval = BMProtoProfile::GetSizeOfExpNoDistinctRes();
    //size_t      expnhist[maxexpval+1];
	std::unique_ptr<size_t[]> expnhist(new size_t[maxexpval+1]);

    memset( expnhist.get(), 0, sizeof(size_t)*(maxexpval+1));

    // iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->IsUsedAt( p ))
            continue;
        pstate = protoprofile_->GetStateAt( p );
        if( !pstate )
            continue;

        expn[PS_M] = expn[PS_I] = expn[PS_D] = 0.0f;
        expn[PS_M] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_M );
        expn[PS_I] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_I );
        expn[PS_D] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_D );
        if( expn[PS_M] <= 0.0f ) {
            if( pstate )
                warning(( preamb + "Expected number of observations is 0." ).c_str());
            expn[PS_M] = 1.0f;
        }
        expnhist[ (int)rintf(SLC_MIN(maxexpval,expn[PS_M])) ]++;
    }

    //calculate the median value of expected nos. observations
    MedianValue( expnhist.get(), maxexpval, &effnos );
    SetEffNoSequences( effnos );
}


    


// =========================================================================
// AdjustWeights: adjusts weights
//
void MSA::AdjustWeights()
{
    size_t p;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::AdjustWeights: Null prototype profile." );

    //iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        // omit unused or non-match positions
        if( !protoprofile_->IsUsedAt( p ))
            continue;
        if( !protoprofile_->GetStateAt( p ))
            continue;
        AdjustWeightsAt( p );
    }
}

// -------------------------------------------------------------------------
// AdjustWeights: adjusts weights at position p
//
void MSA::AdjustWeightsAt( size_t p, float (*mtwghts)[NUMALPH] )
{
    const mystring preamb = "MSA::AdjustWeightsAt: ";
    const float accuracy = 1.0e-4f;
//     const int   efective_number = NUMAA;// effective number of residues
    float   posum = 0.0f,   // sum of weights in a column
            gapw = 0.0f,    // gap weight at a position
            Xw = 0.0f,      // weight for symbol X
            Bw = 0.0f,      // weight for symbol B
            Zw = 0.0f,      // weight for symbol Z
            Jw = 0.0f;      // weight for symbol J
    float   mw;
    char    strbuf[KBYTE];
    unsigned char r;

    static int  symB = HashAlphSymbol('B');//N,D
    static int  symZ = HashAlphSymbol('Z');//Q,E
    static int  symJ = HashAlphSymbol('J');//I,L
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');
    static int  resI = HashAlphSymbol('I');
    static int  resL = HashAlphSymbol('L');

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile." );

    float ( *pmw )[NUMALPH] = mtwghts;
    if( pmw == NULL )
        pmw = const_cast<float(*)[NUMALPH]>( protoprofile_->GetMatchWeightsAt( p ));

    if( pmw == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null match weight vector." );

    posum = 0.0f;
//     for( r = 0; r < NUMALPH; r++ )
    for( r = 0; IsValidResSym( r ); r++ )
        posum += (*pmw)[r];

    if( !posum )
        for( r = 0; r < NUMAA; r++ ) {
            mw = STABLE.PROBABility( r );
            (*pmw)[r] = mw;
            posum += mw;
        }

    if( posum < 1.0f - accuracy || posum > 1.0f + accuracy ) {
        sprintf( strbuf, "Weights are not normalized, %-12g (pos. %d)", posum, (int)(ssize_t)p);
        throw MYRUNTIME_ERROR( preamb + strbuf );
    }

    gapw = (*pmw)[GAP];
    Xw = (*pmw)[X];
    Bw = (*pmw)[symB];
    Zw = (*pmw)[symZ];
    Jw = (*pmw)[symJ];

    // spread out X weight!
    for( r = 0; r < NUMALPH; r++ ) {
        if( 0.0f < STABLE.PROBABility( r ) && r != X )
            (*pmw)[r] += STABLE.PROBABility( r ) * Xw;
    }

    (*pmw)[X] = 0.0f;
    (*pmw)[GAP] = 0.0f;
    if( mtwghts == NULL )
        protoprofile_->SetGapWeightsAt( gapw, p );

    static float Bprob = STABLE.PROBABility( resN ) + STABLE.PROBABility( resD );
    static float Zprob = STABLE.PROBABility( resQ ) + STABLE.PROBABility( resE );
    static float Jprob = STABLE.PROBABility( resI ) + STABLE.PROBABility( resL );

    // spread out B weight: N or D
    if( 0.0f < Bw  ) {
        (*pmw)[resN] += STABLE.PROBABility( resN ) * Bw / Bprob;
        (*pmw)[resD] += STABLE.PROBABility( resD ) * Bw / Bprob;
        (*pmw)[symB] = 0.0f;
    }
    // spread out Z weight: Q or E
    if( 0.0f < Zw  ) {
        (*pmw)[resQ] += STABLE.PROBABility( resQ ) * Zw / Zprob;
        (*pmw)[resE] += STABLE.PROBABility( resE ) * Zw / Zprob;
        (*pmw)[symZ] = 0.0f;
    }
    // spread out J weight: I or L
    if( 0.0f < Jw  ) {
        (*pmw)[resI] += STABLE.PROBABility( resI ) * Jw / Jprob;
        (*pmw)[resL] += STABLE.PROBABility( resL ) * Jw / Jprob;
        (*pmw)[symJ] = 0.0f;
    }
}

// -------------------------------------------------------------------------
// ComputeEstimatedFrequencies: compute estimated probabilities for each
//  residue at each query position;
//  estimated probabilities are computed using pseudo frequencies g and 
//  weighted frequecies f;
//
void MSA::ComputeTargetFrequencies()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::ComputeTargetFrequencies: Null prototype profile.");

    float   information = 0.0f; // information content per position
    float   pseudoFreqn = 0.0f; // pseudocount frequency
    float   obsFrequencyWeight = 0.0f;// weight for observed frequencies (weighted), alpha
    float   denominator = 0.0f; // sum of weights of pseudocounts and observed frequencies
    float   estimated = 0.0f;   // estimated probability
    unsigned char a, b;
    size_t  p;

    // iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        // weight computed as the expected value of different residues per extent column
        obsFrequencyWeight = protoprofile_->ComputeObsFrequencyWeightAt( p );
        information = 0.0f;
        denominator = obsFrequencyWeight + GetPseudoCountWeight();  // alpha + beta

        for( a = 0; a < NUMAA; a++ ) {
            if( STABLE.PROBABility( a ) <= 0.0f )
                continue;

            pseudoFreqn = 0.0f;

            for( b = 0; b < NUMAA; b++ )
                pseudoFreqn += protoprofile_->GetMatchWeightsAt(b,p) * STABLE.FreqRatio(a,b);

            pseudoFreqn *= STABLE.PROBABility( a );
            estimated = ( obsFrequencyWeight * protoprofile_->GetMatchWeightsAt(a,p) + 
                          GetPseudoCountWeight() * pseudoFreqn ) / denominator;
            protoprofile_->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0f < estimated )
                // Q[a] log2( Q[a]/p[a] )
                information += estimated * logf(estimated/STABLE.PROBABility(a)) / SLC_LN2;
        }
        //save information content
        protoprofile_->SetInformationAt(( 0.0f < information )? information: 0.0f, p );
    }
}

// =========================================================================
// ComputeTargetFrequenciesMDL: pseudocounts are computed by the MDL theory
// (see Altschul et al. theory in NAR, 2009)
//
void MSA::ComputeTargetFrequenciesMDL()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::ComputeTargetFrequenciesMDL: Null prototype profile.");

    const float     backcount = 2.0f;       //pseudocount for initial mixing of frequencies
    float           expobscount;            //expected number of observations at a position
    float           pseudomweights[NUMALPH];//obs. frequencies mixed with background probabilities
    float           pseudomsum;             //sum of mixed counts
    float           pseudominfo;            //relative entropy of frequencies wrt background probabilities
    float           pseudofreqweight;       //weight of pseudo frequencies (alpha)
    const float     psnumerator = 0.0806f;  //numerator a in expression a/relent^b
    const float     psexponent = 0.53f;     //exponent b of denominator in a/relent^b
    float           psdenominator;          //denominator in expression a/relent^b
    const float     lowvalue = 0.0001f;

    float   information = 0.0f; //information content per position
    float   pseudoFreqn = 0.0f; //pseudocount frequency
    float   estimated = 0.0f;   //estimated probability
    unsigned char a, b;
    size_t p;

    //iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        //{{
        memset( pseudomweights, 0, NUMALPH * sizeof(float));
//         expobscount = protoprofile_->GetMIDExpNoObservationsAt( p, PS_M );
        expobscount = protoprofile_->GetNoSymbolsInExtentAt( p, PS_M );
        pseudomsum = 0.0f;
        pseudominfo = 0.0f;
        pseudofreqweight = 1.0f;

        for( a = 0; a < NUMAA; a++ ) {
            pseudomweights[a] = ( protoprofile_->GetMatchWeightsAt(a,p ) * expobscount ) +
                                ( STABLE.PROBABility(a) * backcount );
            pseudomsum += pseudomweights[a];
        }
        if( pseudomsum )
            for( a = 0; a < NUMAA; a++ ) {
                pseudomweights[a] /= pseudomsum;
                if( 0.0f < pseudomweights[a] && 0.0f < STABLE.PROBABility(a))
                    pseudominfo += pseudomweights[a] * logf(pseudomweights[a]/STABLE.PROBABility(a));
            }

        pseudominfo /= SLC_LN2;
        if( pseudominfo < 0.0f )
            pseudominfo = 0.0f;
        if( lowvalue < pseudominfo ) {
            psdenominator = powf( pseudominfo, psexponent );
            if( psdenominator )
                pseudofreqweight = psnumerator / psdenominator;
        }

        if( pseudofreqweight < 0.0f ) pseudofreqweight = 0.0f;
        if( 1.0f < pseudofreqweight ) pseudofreqweight = 1.0f;
        //}}

        //pseudocounts do not depend on observations n (m=n*alpha/(1-alpha))
        information = 0.0f;
        for( a = 0; a < NUMAA; a++ ) {
            if( STABLE.PROBABility(a) <= 0.0f )
                continue;

            pseudoFreqn = 0.0f;
            if( pseudofreqweight )
                for( b = 0; b < NUMAA; b++ )
                    pseudoFreqn += protoprofile_->GetMatchWeightsAt(b,p) * STABLE.FreqRatio(a,b);

            pseudoFreqn *= STABLE.PROBABility( a );
            estimated = pseudofreqweight * pseudoFreqn + 
                        (1.0f-pseudofreqweight) * protoprofile_->GetMatchWeightsAt(a,p);
            protoprofile_->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0f < estimated )
                information += estimated * logf(estimated/STABLE.PROBABility(a));
        }

        information /= SLC_LN2;
        if( information < 0.0f )
            information = 0.0f;
        // save information content
        protoprofile_->SetInformationAt( information, p );
    }
}

// =========================================================================
// ComputeTargetFrequenciesMDLVar: pseudocounts are computed by the MDL 
// theory
//
void MSA::ComputeTargetFrequenciesMDLVar()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::ComputeTargetFrequenciesMDLVar: Null prototype profile.");

    float       expobscount;        //expected number of observations at a position
    float       pseudofreqweight;   //weight of pseudo frequencies (alpha)
    const float psnumerator = 0.65f;//numerator a in expression a/n^b
    const float psexponent = 0.90f; //exponent b of denominator in a/n^b
    float       psdenominator;      //denominator in expression a/n^b
//     const float lowvalue = 0.0001f;

    float   information = 0.0f; //information content per position
    float   pseudoFreqn = 0.0f; //pseudocount frequency
    float   estimated = 0.0f;   //estimated probability
    unsigned char a, b;
    size_t p;

    //iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        //{{
//         expobscount = protoprofile_->GetMIDExpNoObservationsAt( p, PS_M );
        expobscount = protoprofile_->GetNoSymbolsInExtentAt( p, PS_M );
        pseudofreqweight = 1.0f;

        psdenominator = powf( expobscount, psexponent );
        if( psdenominator )
            pseudofreqweight = psnumerator / psdenominator;

        if( pseudofreqweight < 0.0f ) pseudofreqweight = 0.0f;
        if( 1.0f < pseudofreqweight ) pseudofreqweight = 1.0f;
        //}}

        //pseudocounts do not depend on observations n (m=n*alpha/(1-alpha))
        information = 0.0f;
        for( a = 0; a < NUMAA; a++ ) {
            if( STABLE.PROBABility(a) <= 0.0f )
                continue;

            pseudoFreqn = 0.0f;
            if( pseudofreqweight )
                for( b = 0; b < NUMAA; b++ )
                    pseudoFreqn += protoprofile_->GetMatchWeightsAt(b,p) * STABLE.FreqRatio(a,b);

            pseudoFreqn *= STABLE.PROBABility( a );
            estimated = pseudofreqweight * pseudoFreqn + 
                        (1.0f-pseudofreqweight) * protoprofile_->GetMatchWeightsAt(a,p);
            protoprofile_->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0f < estimated )
                information += estimated * logf(estimated/STABLE.PROBABility(a));
        }

        information /= SLC_LN2;
        if( information < 0.0f )
            information = 0.0f;
        // save information content
        protoprofile_->SetInformationAt( information, p );
    }
}





// =========================================================================
// CalculateAndPrintXCovMatrices: calculate and print cross-covariance 
// matrices between the positions of the profile
//
void MSA::CalculateAndPrintXCovMatrices(const char* filename)
{
    const mystring preamb = "MSA::CalculateAndPrintXCovMatrices: ";
    myruntime_error mre;
    FILE* fp = NULL;

    //{{local run-time control variables:
    //flag of using extents when determining whether sequence positions 
    //contribute to calculations:
    static const bool cbExtentsUsed = CLOptions::GetXCOV_USEEXTENTS();
    //mix observed frequencies with target probabilities when calculating
    //cross-covariance matrix:
    static const bool cbXCovMixTrgFrqs = CLOptions::GetXCOV_MIXTRGFRQS();
    //compute cross-correlation matrix instead of cross-covariance matrix:
    static const bool cbXCorrelation = CLOptions::GetXCOV_CORR();
    //scale sequence weights by the number of sequences in the respective extent:
    static const bool cbScaleSeqWeights = CLOptions::GetXCOV_SCALEWGTS();
    //additionally calculate mutual information and add to calculated xcov values:
    static const bool cbMI = CLOptions::GetXCOV_MI();
    static constexpr size_t noeffress = NUMAA;
    //}}

    //number of elements in cross-covariance matrix
    static constexpr size_t nxcovels = noeffress * noeffress;
    float xcov[nxcovels], mutinf;
    const float (*trgfrqs)[NUMALPH];//target frequencies at position i
    const float (*trgfrqs2)[NUMALPH];//target frequencies at position j
    float frqs[noeffress], frqs2[noeffress], fct1;
    float avgfrqs[noeffress], avgfrqs2[noeffress];//average frequencies

    const size_t    noseqs = GetSize();
    const size_t    nomstates = protoprofile_->GetEffectiveSize();

    float* weights[PS_NSTATES] = {NULL,NULL,NULL};
    float* weights2[PS_NSTATES] = {NULL,NULL,NULL};

    unsigned char res1, res2, a1, a2;
    std::unique_ptr<char,MyMSAStrDestroyer> errbuf((char*)std::malloc(20*KBYTE));
    size_t  p, p2, pM, pM2, i, n, nseqs;

    fp = fopen( filename, "w" );
    if( fp == NULL )
        throw MYRUNTIME_ERROR(preamb + "Failed to open file for writing.");

    try {
        if( fprintf(fp, "## COMER xcov format v2.2\nLength= %zu xcov_size= %zu\n",
                nomstates, cbMI? nxcovels+1: nxcovels) < 0 )
            throw MYRUNTIME_ERROR(preamb + "Write to file failed.");

        for( p = pM = 0; p < protoprofile_->GetSize(); p++ ) {
            // omit unsued positions
            if( !protoprofile_->IsUsedAt(p))
                continue;
            if( !protoprofile_->GetStateAt(p))
                continue;

            pM++;

            weights[PS_M] = protoprofile_->GetSqnWeightsAt(p, PS_M);
            if( !weights[PS_M]) {
                sprintf(errbuf.get(), "Null weights2 at %zu.", p);
                throw MYRUNTIME_ERROR(preamb + errbuf.get());
            }

            trgfrqs = protoprofile_->GetTargetFreqnsAt(p);
            if(trgfrqs == NULL) {
                sprintf(errbuf.get(), "Null target probabilities at %zu.", p);
                throw MYRUNTIME_ERROR(preamb + errbuf.get());
            }

            for( p2 = p, pM2 = pM-1; p2 < protoprofile_->GetSize(); p2++ ) {
                // omit unsued positions
                if( !protoprofile_->IsUsedAt(p2))
                    continue;
                if( !protoprofile_->GetStateAt(p2))
                    continue;

                pM2++;

                weights2[PS_M] = protoprofile_->GetSqnWeightsAt(p2, PS_M);
                if( !weights2[PS_M]) {
                    sprintf(errbuf.get(), "Null weights2 at %zu.", p2);
                    throw MYRUNTIME_ERROR(preamb + errbuf.get());
                }

                trgfrqs2 = protoprofile_->GetTargetFreqnsAt(p2);
                if(trgfrqs2 == NULL) {
                    sprintf(errbuf.get(), "Null target2 probabilities at %zu.", p2);
                    throw MYRUNTIME_ERROR(preamb + errbuf.get());
                }

                mutinf = 0.0f;
                memset(xcov, 0, nxcovels * sizeof(xcov[0]));
                memset(avgfrqs, 0, noeffress * sizeof(avgfrqs[0]));
                memset(avgfrqs2, 0, noeffress * sizeof(avgfrqs2[0]));

                //{{cross-covariance between positions p and p2
                for( i = 0, nseqs = 0; i < noseqs; i++ ) {
                    if(!GetSequenceAt(i)->GetUsed())
                        continue;
                    if(cbExtentsUsed) {
                        // omit sequences absent from the extent
                        if( !GetSequenceAt(i)->IsUsedInExtentAt(p) ||
                            !GetSequenceAt(i)->IsUsedInExtentAt(p2))
                            continue;
                    }
                    res1 = GetSequenceAt(i)->GetResidueAt(p);
                    res2 = GetSequenceAt(i)->GetResidueAt(p2);
                    if( res1 == GAP || res2 == GAP )
                        continue;

                    if(cbXCovMixTrgFrqs) {
                        memcpy(frqs, *trgfrqs, noeffress * sizeof(float));
                        memcpy(frqs2, *trgfrqs2, noeffress * sizeof(float));
                        if(res1 < noeffress)
                            frqs[res1] = 1.0f;
                        if(res2 < noeffress)
                            frqs2[res2] = 1.0f;
                    }

                    nseqs++;

                    //averages for position p2
                    if(cbXCovMixTrgFrqs) {
                        for(a2 = 0; a2 < noeffress; a2++ )
                            avgfrqs2[a2] += frqs2[a2] * weights2[PS_M][i];
                    } else if(res1 < noeffress && res2 < noeffress)
                        avgfrqs2[res2] += weights2[PS_M][i];

                    //cross-covariance; TODO: unroll
                    if(cbXCovMixTrgFrqs) {
                        for(a1 = 0, n = 0; a1 < noeffress; a1++ ) {
                            fct1 = frqs[a1] * weights[PS_M][i];
                            avgfrqs[a1] += fct1;
                            for(a2 = 0; a2 < noeffress; a2++, n++ )
                                xcov[n] += fct1 * frqs2[a2] * weights2[PS_M][i];
                        }
                    } else if(res1 < noeffress && res2 < noeffress) {
                        avgfrqs[res1] += weights[PS_M][i];
                        n = res1 * noeffress + res2;
                        xcov[n] += weights[PS_M][i] * weights2[PS_M][i];
                    }
                }
                //}}//xcov over sequences
                if(nomstates < pM2)
                    throw MYRUNTIME_ERROR(preamb + "Invalid position 2.");
                //{{adjust cross-covariance for this pair of positions
                if(nseqs >= 1) {
                    if( cbMI) {
                        float nrm1 = 0.0f, nrm2 = 0.0f, nrmp = 0.0f;
                        for(a1 = 0, n = 0; a1 < noeffress; a1++ ) {
                            nrm1 += avgfrqs[a1];
                            nrm2 += avgfrqs2[a1];
                            for(a2 = 0; a2 < noeffress; a2++, n++ )
                                nrmp += xcov[n];
                        }
                        if( nrmp) {
                            for(a1 = 0, n = 0; a1 < noeffress; a1++ )
                                for(a2 = 0; a2 < noeffress; a2++, n++ ) {
                                    //if xcov[n]>0, then both avgfrqs>0 too
                                    if(!xcov[n])
                                        continue;
                                    mutinf += 
                                        xcov[n] * logf(xcov[n]*nrm1*nrm2/(avgfrqs[a1]*avgfrqs2[a2]*nrmp));
                                }
                            mutinf /= nrmp;
                        }
                        //remove error effects
                        if( mutinf < 0.0f)
                            mutinf = 0.0f;
                    }
                    //NOTE: when scaling is used, every element has to be 
                    //scaled by the number of sequences; the scaling then 
                    //reduces for the mean vectors (because of normalization) 
                    //and it becomes nseqs for the cross-covariance matrix 
                    //(scaling=nseqs*nseqs, norm_constant=nseqs)
                    fct1 = cbScaleSeqWeights? (float)nseqs: 1.0f/(float)nseqs;
                    for(a1 = 0, n = 0; a1 < noeffress; a1++ ) {
                        if( !cbScaleSeqWeights) {
                            avgfrqs[a1] *= fct1;
                            avgfrqs2[a1] *= fct1;
                        }
                        if(cbXCorrelation) {
                            for(a2 = 0; a2 < noeffress; a2++, n++ )
                                xcov[n] = fct1 * xcov[n];
                        }
                        else {
                            for(a2 = 0; a2 < noeffress; a2++, n++ )
                                xcov[n] = fct1 * xcov[n] - avgfrqs[a1] * avgfrqs2[a2];
                        }
                    }
                }
                //}}
                //{{write to file
                char *pbuf = errbuf.get();
                for(n = 0; n < nxcovels; n++) {
                    sprintf(pbuf, " %.3g", xcov[n]);
                    pbuf += strlen(pbuf);
                }
                if( cbMI)
                    sprintf(pbuf, " %.3g", mutinf);
                if( fprintf(fp, "%zu %zu  %s\n", pM, pM2, errbuf.get()) < 0 )
                    throw MYRUNTIME_ERROR(preamb + "Write to file failed.");
                //}}
            }//for(p2)
        }//for(p)
    } catch( myexception const& ex ) {
        mre = ex;
    }

    if(fp) 
        fclose(fp);
    if( mre.isset())
        throw mre;
}





// =========================================================================
// CalcTFPosteriorPredictives: using the HDP mixture model, calculate 
// cluster membership posterior predictive probabilities for further 
// processing of target frequencies
//
void MSA::CalcTFPosteriorPredictives()
{
    const mystring preamb = "MSA::CalcTFPosteriorPredictives: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const int       noeffress = NUMAA;
    float           ppr = 0.0f;//posterior probability
    float           tfr = 0.0f;//target frequency
    const HDPbase*  hdpbase = GetHDPbase();

    if( hdpbase == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDPbase object.");

    int     nosupcls = hdpbase->GetNoSupClusters();
    int     ctxtsz = hdpbase->GetCtxtSize();
    size_t  psize = protoprofile_->GetEffectiveSize();
    size_t  p, pp;
    int     left, right, r;
    int     hlf = ctxtsz >> 1;
//     int     parity = ( ctxtsz & 1 ) ^ 1;
//     int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of HDP support clusters.");
    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP context size.");
    if((int)psize < ctxtsz ) {
        warning(( preamb + "Profile length < context size." ).c_str());
        return;
    }

    extspsl::Pslmatrix  promtx( (int)psize, noeffress );//profile matrix
    extspsl::Pslmatrix  ctxmtx;//context matrix
    extspsl::Pslvector  ctxsck;//context stack
    extspsl::Pslvector  ctxnrm( noeffress*ctxtsz );//normal transform
    extspsl::Pslvector  ppprobs;//vector of posterior probabilities
    extspsl::Ivector    cindcs;//indices of clusters

    //make matrix representation of the profile
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tfr = protoprofile_->GetTargetFreqnsAt( r, p );
            promtx.SetValueAt((int)pp, (int)r, tfr );
        }
        pp++;
    }

    //iterate over all query positions
    for( p = 0, pp = (size_t)-1; p < protoprofile_->GetSize(); p++ )
    {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        pp++;

        right = (int)SLC_MIN( (int)psize-1, (int)pp+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //calculate posterior predictive probabilities;
        //input vector ctxnrm will be transformed!
        hdpbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs );

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw MYRUNTIME_ERROR( preamb + "Invalid number of posterior predictive probabilities.");

        //save posteriors
        protoprofile_->SetBckPPProbAt( ppr, p );
        protoprofile_->SetPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}

// =========================================================================
// MixTargetFrequenciesHDPCtx: mix target frequencies using the HDP mixture 
// model
//
void MSA::MixTargetFrequenciesHDPCtx()
{
    const mystring preamb = "MSA::MixTargetFrequenciesHDPCtx: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const int       noeffress = NUMAA;
    float           pme = 0.0f;//estimated probability
    float           tfr = 0.0f;//target frequency
    const HDPbase*  hdpbase = GetHDPbase();

    if( hdpbase == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null HDPbase object.");

    int     ctxtsz = hdpbase->GetCtxtSize();
    size_t  psize = protoprofile_->GetEffectiveSize();
    size_t  p, pp;
    int     left, right, r;
    int     hlf = ctxtsz >> 1;
//     int     parity = ( ctxtsz & 1 ) ^ 1;
//     int     mid = hlf - parity;

    if( ctxtsz < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid HDP context size.");
    if((int)psize < ctxtsz ) {
        warning(( preamb + "Profile length < context size." ).c_str());
        return;
    }

    extspsl::Pslmatrix  promtx(( int)psize, noeffress );//profile matrix
    extspsl::Pslmatrix  ctxmtx;//context matrix
    extspsl::Pslvector  ctxsck;//context stack
    extspsl::Pslvector  ctxnrm( noeffress*ctxtsz );//normal transform
    extspsl::Pslvector  mixed( noeffress );//mixed vector
    float infrm;

    //make matrix representation of the profile
    for( p = 0, pp = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tfr = protoprofile_->GetTargetFreqnsAt( r, p );
            promtx.SetValueAt((int)pp, (int)r, tfr );
        }
        pp++;
    }

    //iterate over all query positions
    for( p = 0, pp = (size_t)-1; p < protoprofile_->GetSize(); p++ )
    {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        pp++;

        right = (int)SLC_MIN( (int)psize-1, (int)pp+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //mix central vector of context;
        //input vector ctxnrm will be transformed!;
        //mixed will be returned in logit-normal space
        hdpbase->MixCLNVar( ctxnrm, &mixed );

        //write PME back to target frequencies
        infrm = 0.0f;
        for( r = 0; r < noeffress; r++ ) {
            pme = mixed.GetValueAt(r);
            protoprofile_->SetTargetFreqnsAt( pme, r, p );
            //calculate relative entropy
            if( STABLE.PROBABility(r) <= 0.0f )
                continue;
            if( 0.0f < pme )
                infrm += pme * logf( pme/ STABLE.PROBABility(r) );
        }
        //save relative entropy
        infrm /= SLC_LN2;
        infrm = SLC_MAX( 0.0f, infrm );
        protoprofile_->SetInformationAt( infrm, p );
    }
}





// -------------------------------------------------------------------------
// RecalcBackgroundProbabilities: recalculate background probabilities by 
// mixing background and target probabilities;
// NOTE: PARAMETERS: tadj
// NOTE: bg probabilities are also recalculated in 
// ProGenerator::RecalcBgProbs
//
void MSA::RecalcBackgroundProbabilities()
{
    const mystring preamb = "MSA::RecalcBackgroundProbabilities: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
//     const size_t    nomstates = protoprofile_->GetEffectiveSize();
    const float     effnos = GetEffNoSequences();
    const float     maxtgtobss = 20.0f;
    //const float     tadj = 1.+4.*exp(-0.5*1.e-6*pow(effnos-18.,6.));//5.0
    const float     tadj = 1.f;//5.f
    const float     bckpseudocounts = SLC_MAX( 1.f, tadj * effnos );//10.0;
    const float     errtol = 1.0e-3f;
    char            errbuf[KBYTE];

    float           bpp[noeffress], mixp[noeffress];
    float           denm, consv;
    float           bckp, tgtp;
    float           expM, tgtpc;
    size_t          p, r;//, nn;

//     consv = 0.0;
//     for( r = 0; r < noeffress; r++ ) {
//         bckp = STABLE.PROBABility( r );
//         bpp[r] = 100.0f / effnos * bckp;
//         mixp[r] = bpp[r];
//         consv += bpp[r];
//     }
//     for( p = 0; p < protoprofile_->GetSize(); p++ ) {
//         if( !protoprofile_->GetStateAt( p ))
//             continue;
//         for( r = 0; r < noeffress; r++ ) {
//             tgtp = protoprofile_->GetTargetFreqnsAt( r, p );
//             mixp[r] += tgtp;
//             consv += tgtp;
//         }
//     }

    for( r = 0; r < noeffress; r++ ) {
        bckp = STABLE.PROBABility((int)r);
        bpp[r] = bckpseudocounts * bckp;
        mixp[r] = 0.0f;
    }

    consv = 0.0f;
    for( p = 0/*, nn = 0*/; p < protoprofile_->GetSize(); p++ ) 
    {
        if( !protoprofile_->GetStateAt( p ))
            continue;

//         nn++;
//         expM = protoprofile_->GetMIDExpNoObservationsAt( p, PS_M );
        expM = protoprofile_->GetNoSymbolsInExtentAt( p, PS_M );
        tgtpc = SLC_MIN( maxtgtobss, expM );
        denm = 1.0f /( bckpseudocounts + tgtpc );

        for( r = 0; r < noeffress; r++ ) {
            tgtp = protoprofile_->GetTargetFreqnsAt((unsigned char)r, p );
            tgtp = ( bpp[r] + tgtpc * tgtp ) * denm;
            mixp[r] += tgtp;
            consv += tgtp;
        }
    }

    denm = 1.0f;
    if( consv )
        denm /= consv;

    consv = 0.0f;
    for( r = 0; r < noeffress; r++ ) {
        mixp[r] *= denm;
        consv += mixp[r];
        protoprofile_->SetBackProbsAt((unsigned char)r, mixp[r]);
    }
    for( ; r < noress; r++ )
        protoprofile_->SetBackProbsAt((unsigned char)r, 0.0f );
    if( consv < 1.0f - errtol || consv > 1.0f + errtol ) {
        sprintf( errbuf, "Probabilities are not conserved: %g", consv );
        throw MYRUNTIME_ERROR( preamb + errbuf );
    }
    return;

//     denm = 1.0f;
//     if( nomstates )
//         denm /= (float)nomstates;
//     for( r = 0; r < noeffress; r++ )
//         mixp[r] *= denm;
// 
//     LogitNormalErrCorr( mixp, noeffress, errtol );//will throw on error
// 
//     for( r = 0; r < noeffress; r++ )
//         protoprofile_->SetBackProbsAt( r, mixp[r]);
//     for( ; r < noress; r++ )
//         protoprofile_->SetBackProbsAt( r, 0.0f );
//     return;
}

// -------------------------------------------------------------------------
// CalcPosteriorProbabilities: calculate posterior (generalized target) 
// probabilities
//
void MSA::CalcPosteriorProbabilities()
{
    const mystring preamb = "MSA::CalcPosteriorProbabilities: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
//     const size_t    nomstates = protoprofile_->GetEffectiveSize();
//     const float     effnos = GetEffNoSequences();
    const float     errtol = 1.0e-3f;
    char            errbuf[KBYTE];

    float   ppp[noeffress];
    float   denm, consv;
    float   tgtp;
    size_t  p, r;

    for( r = 0; r < noeffress; r++ )
        ppp[r] = 0.0f;

    consv = 0.0f;
    for( p = 0; p < protoprofile_->GetSize(); p++ )
    {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tgtp = protoprofile_->GetTargetFreqnsAt((unsigned char)r, p );
            ppp[r] += tgtp;
            consv += tgtp;
        }
    }

    denm = 1.0f;
    if( consv )
        denm /= consv;

    consv = 0.0f;
    for( r = 0; r < noeffress; r++ ) {
        ppp[r] *= denm;
        consv += ppp[r];
        protoprofile_->SetPostProbsAt((unsigned char)r, ppp[r]);
    }
    for( ; r < noress; r++ )
        protoprofile_->SetPostProbsAt((unsigned char)r, 0.0f );
    if( consv < 1.0f - errtol || consv > 1.0f + errtol ) {
        sprintf( errbuf, "Probabilities are not conserved: %g", consv );
        throw MYRUNTIME_ERROR( preamb + errbuf );
    }
}

// -------------------------------------------------------------------------
// ComputePSSM: compute PSSM matrix given estimated target frequencies
//
void MSA::ComputePSSM()
{
    if( !protoprofile_ )
        throw MYRUNTIME_ERROR("MSA::ComputePSSM: Null prototype profile.");

    float   estimated = 0.0f;   // estimated probability
    float   pssmvalue = 0.0f;   // PSSM value to be calculated
    float   ratiosval = 0.0f;   // STABLE frequency ratio value
//     int     scoresval = 0;      // STABLE matrix value
    float   lambda = STABLE.StatisParam( Ungapped, Lambda );//reference lambda
    unsigned char r, residue;
    size_t  p;

    // iterate over all query positions
    for( p = 0; p < protoprofile_->GetSize(); p++ )
    {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        residue = protoprofile_->GetResidueAt( p );

        for( r = 0; r < NUMAA/*NUMALPH*/; r++ ) {
            estimated = protoprofile_->GetTargetFreqnsAt( r, p );
            if( STABLE.PROBABility(r) <= 0.0f || estimated <= 0.0f ) {
                //unable to estimate probabilities if background probabilities are zero
//                 scoresval = STABLE.Entry( residue, r );
                ratiosval = STABLE.FreqRatio( residue, r );

                if( /*scoresval == SCORE_MIN || */ratiosval <= 0.0f )
                    //unable to compute PSSM value
                    pssmvalue = SCORE_MIN;
                else
                    // this is how STABLE matrix values are obtained;
                    // with frequency ratios, STABLE matrix values are 1/s * log2(ratio)
                    // (s is the scaling constant for STABLE values);
                    // pssmvalue = log( ratiosval ) / LN2 * ScaleFactor;
                    // use precomputed values instead...
                    pssmvalue = STABLE.PrecomputedEntry( residue, r );

            } else {
                pssmvalue = logf(estimated/STABLE.PROBABility(r)) / lambda;
            }
            protoprofile_->SetPSSMEntryAt( pssmvalue, r, p );
        }
		for (; r < NUMALPH; r++)
			protoprofile_->SetPSSMEntryAt(0.0f, r, p);
    }
}



// -------------------------------------------------------------------------
// ExportProtoProfile: export profile model
//
void MSA::ExportProtoProfile( pmodel::PMProfileModel& pmpm, bool appfilename ) const
{
    const mystring preamb = "MSA::ExportProtoProfile: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    if( !protoprofile_->GetSize() || !protoprofile_->GetEffectiveSize())
        return;

    size_t  noress = NUMALPH;
	//size_t  noeffress = NUMAA;
	size_t  queffsize = protoprofile_->GetEffectiveSize();
    float   expMIDs[PS_NSTATES];
    const float ( *scores )[NUMALPH];//pssm scores
    const float ( *mwghts )[NUMALPH];//observed weighted match frequencies
    const float ( *tgtfqs )[NUMALPH];//target frequencies
    float           freqweight;
    float           information;
    //float           expnobs;
    //float           bppprob;//background posterior predicitive
    //const float*    ppprobs;//posterior predictive probabilities
    //const int*      pppndxs;//indices of posteriors
    //size_t          noppps;//number of p.p.probability values
    //size_t          noseqs;
    unsigned char   residue;
    size_t  p;

    pmpm.Reserve((int)queffsize);

    if( queffsize ) {
        p = -1;
        expMIDs[PS_M] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_M );
        expMIDs[PS_I] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_I );
        expMIDs[PS_D] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_D );
        pmpm.SetMIDExpNoObservationsBeg( expMIDs );
    }

    for( p = 0; p <noress; p++ )
        pmpm.SetBackProbsAt((char)p, protoprofile_->GetBackProbsAt((unsigned char)p));

    for( p = 0; p <noress; p++ )
        pmpm.SetPostProbsAt((char)p, protoprofile_->GetPostProbsAt((unsigned char)p));

    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        residue = protoprofile_->GetResidueAt( p );

        //save PSSM scores
        scores = protoprofile_->GetPSSMVectorAt( p );
        mwghts = protoprofile_->GetMatchWeightsAt( p );
        tgtfqs = protoprofile_->GetTargetFreqnsAt( p );
        freqweight = protoprofile_->ComputeObsFrequencyWeightAt( p );
        information = protoprofile_->GetInformationAt( p );
        //expnobs = protoprofile_->GetExpNoObservationsAt( p );
        //noseqs = protoprofile_->GetNoSequencesInExtentAt( p );//queryDescription->GetNoSequencesAt( p );

        if( !scores )
            throw MYRUNTIME_ERROR( preamb + "Null scores.");
        if( !mwghts )
            throw MYRUNTIME_ERROR( preamb + "Null match weights.");
        if( !tgtfqs )
            throw MYRUNTIME_ERROR( preamb + "Null target frequencies.");

        expMIDs[PS_M] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_M );
        expMIDs[PS_I] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_I );
        expMIDs[PS_D] = protoprofile_->GetMIDExpNoObservationsAt((int)p, PS_D );

        //bppprob = protoprofile_->GetBckPPProbAt(p);
        //ppprobs = protoprofile_->GetPPProbsAt(p);
        //pppndxs = protoprofile_->GetPPPIndsAt(p);
        //noppps = protoprofile_->GetNoPPProbsAt(p);

        //first, push scores for scaling:
        pmpm.Push( *scores/*vals*/, *mwghts/*frqs*/, (char)residue, freqweight, information, expMIDs );
        //pmpm.Push( *tgtfqs/*vals*/, *mwghts/*frqs*/, (char)residue, freqweight, information, expMIDs );
    }


    Configuration config[NoCTypes];//no filename; read or write attempt issues an error
    SetUngappedParams( config[CTUngapped]);//ungapped configuration
    SetUngappedParams( config[CTGapped]);//not using gapped configuration; make it identical to ungapped

    //table to scale scores and calculate composition-based statistics;
    //NOTE: the first argument must be scores!
    ProfileMatrix<int> matrix(
            (( const pmodel::PMProfileModel& )pmpm).GetVector(),
            pmpm.GetResidues(),
            pmpm.GetSize(),
            config
    );

    if( matrix.GetSupportOptimFreq())
        matrix.OptimizeTargetFrequencies();

#ifdef SCALE_PROFILES
    matrix.ScaleScoreMatrix();
#else
    matrix.ComputeStatisticalParameters();
#endif


    for( int c = 0; c < pmpm.GetSize(); c++ ) {
        residue = pmpm.GetResidueAt( c );
        scores = matrix.GetVectorAt( c );//scaled scores
        if( !scores )
            throw MYRUNTIME_ERROR( preamb + "Null scaled scores.");
        pmpm.PushAt( *scores/*vals*/, (char)residue, c );
        //target frequencies are calculated from scores below...
        //pmpm.PushAt( *tgtfqs/*vals*/, (char)residue, c );
    }

    pmpm.SetNoSequences( GetNoSequences());
    pmpm.SetEffNoSequences( GetEffNoSequences());

    pmpm.SetRefLambda(      matrix.GetRefLambda());
    pmpm.SetRefK(           matrix.GetRefK());
    pmpm.SetLambda(         matrix.GetLambda());
    pmpm.SetEntropy(        matrix.GetEntropy());
    pmpm.SetK(              matrix.GetK());
    pmpm.SetExpectedScore(  matrix.GetExpectedScore());

    if( appfilename )
        pmpm.SetName( GetName());
    pmpm.SetDescription( GetDescription());

    //copnvert scores back to target frequencies to correctly 
    // calculate posterior predictive probabilities and context vector
    pmpm.ConvertToTrgFreqs();//pmpm.Finalize();

    if( GetScoAdjmentHDPCtx()) {
        //calculate posteriors after pmpm.Finalize(); 
        // using target frequencies
        pmpm.CalcTFPostPredictives( GetHDPbase());
//         if( GetHDPctBase())
//             pmpm.Calc_ctTFPostPredictives( GetHDPctBase());
    }
    pmpm.CalcCtxVector();

    //{{NOTE:for TEST only
    //pmpm.ConvertToScores();
    //}}
}

// -------------------------------------------------------------------------
// ExportProtoProfile: export transition probabilities (gap model) 
//
void MSA::ExportProtoProfile( pmodel::PMTransModel& pmtm ) const
{
    const mystring preamb = "MSA::ExportProtoProfile: ";

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( preamb + "Null prototype profile.");

    if( !protoprofile_->GetSize() || !protoprofile_->GetEffectiveSize())
        return;

    size_t  queffsize = protoprofile_->GetEffectiveSize();
    const float ( *ttps )[P_NSTATES] = NULL;//target transition probabilities
    //unsigned char   residue;
    size_t          p;

    pmtm.Reserve((int)queffsize);

    pmtm.SetOrgTrProbsBeg( protoprofile_->GetTargetTranstBeg());

    for( p = 0; p < protoprofile_->GetSize(); p++ ) {
        if( !protoprofile_->GetStateAt( p ))
            continue;

        //residue = protoprofile_->GetResidueAt( p );
        ttps = protoprofile_->GetTargetTranstAt( p );

        pmtm.Push( ttps );
    }
}
