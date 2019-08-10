/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __BMProtoProfile_h__
#define __BMProtoProfile_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>

#include "liblib/alpha.h"
#include "TRANSPROBS.h"
#include "BMSequence.h"

// _________________________________________________________________________
// Class BMProtoProfile
// representing MSA by prototype profile model
//
class BMProtoProfile: public BMSequence {
public:
    // typedefs ...
    //
    enum Direction {    //extent direction type
        xLeft,          //left boundary
        xRight,         //right boundary
        xInterval,      //interval of valid positions in the extent
        xNoSequences,   //number of sequences participating in the extent
        xNoSymbols,     //number of different symbols (sum) occuring in the extent;NOTE:previously
        xCount
    };
    BMProtoProfile( size_t reservation );
    BMProtoProfile( const BMSequence& sequence );
    virtual ~BMProtoProfile();

    virtual BMProtoProfile& operator=( const BMSequence& sequence );

    //output methods
    void            PrintMatchWeights( FILE* );
    void            PrintTransWeights( FILE* );
    void            PrintTargetTranst( FILE* );
    void            PrintMIDExpNoObservations( FILE* );
    void            PrintPSSMatrix( FILE* );
    void            PrintSuppressedPSSMandWeights( FILE* );
    void            PrintProfile( FILE* );

    //get/set methods...
    size_t*         GetIndicesAt( size_t n ) const;
    size_t          GetIndicesSizeAt( size_t n ) const;
    void            PushIndexAt( size_t value, size_t n );

    const float*    GetBackProbs() const { return backprobs_; }
    float           GetBackProbsAt( unsigned char res ) const;
    void            SetBackProbsAt( unsigned char res, float value );

    const float*    GetPostProbs() const { return postprobs_; }
    float           GetPostProbsAt( unsigned char res ) const;
    void            SetPostProbsAt( unsigned char res, float value );

    int             GetStateAt( size_t n ) const;
    void            SetStateAt( size_t n );



    size_t  GetLeftExtentAt( size_t n ) const {return GetLeftExtentAt( n, PS_M );}
    void    SetLeftExtentAt( size_t value, size_t n ) {SetLeftExtentAt( n, PS_M, value );}

    size_t  GetRightExtentAt( size_t n ) const {return GetRightExtentAt( n, PS_M );}
    void    SetRightExtentAt( size_t value, size_t n ) {SetRightExtentAt( n, PS_M, value );}

    size_t  GetExtentIntervalAt( size_t n ) const {return GetExtentIntervalAt( n, PS_M );}
    void    SetExtentIntervalAt( size_t value, size_t n ) {SetExtentIntervalAt( n, PS_M, value );}

    size_t  GetNoSequencesInExtentAt( size_t n ) const {return GetNoSequencesInExtentAt( n, PS_M );}
    void    IncNoSequencesInExtentAt( size_t n ) {IncNoSequencesInExtentAt( n, PS_M );}

    float   GetNoSymbolsInExtentAt( size_t n ) const {return GetNoSymbolsInExtentAt( n, PS_M );}
    void    SetNoSymbolsInExtentAt( float value, size_t n ) {SetNoSymbolsInExtentAt( n, PS_M, value );}


    size_t          GetLeftExtentAt( size_t n, int st ) const;
    void            SetLeftExtentAt( size_t n, int st, size_t value );

    size_t          GetRightExtentAt( size_t n, int st ) const;
    void            SetRightExtentAt( size_t n, int st, size_t value );

    size_t          GetExtentIntervalAt( size_t n, int st ) const;
    void            SetExtentIntervalAt( size_t n, int st, size_t value );

    size_t          GetNoSequencesInExtentAt( size_t n, int st ) const;
    void            IncNoSequencesInExtentAt( size_t n, int st );

    float           GetNoSymbolsInExtentAt( size_t n, int st ) const;
    void            SetNoSymbolsInExtentAt( size_t n, int st, float value );



    size_t          GetLeftMSExtentAt( size_t n ) const;
    void            SetLeftMSExtentAt( size_t value, size_t n );

    size_t          GetRightMSExtentAt( size_t n ) const;
    void            SetRightMSExtentAt( size_t value, size_t n );

    size_t          GetMSExtentIntervalAt( size_t n ) const;
    void            SetMSExtentIntervalAt( size_t value, size_t n );

    size_t          GetNoSequencesInMSExtentAt( size_t n ) const;
    void            IncNoSequencesInMSExtentAt( size_t n );

    size_t          GetNoSymbolsInMSExtentAt( size_t n ) const;
    void            SetNoSymbolsInMSExtentAt( size_t value, size_t n );


    size_t          GetCountAt( size_t n ) const;
    void            IncCountAt( size_t n );

    float           ComputeObsFrequencyWeightAt( size_t n ) const;

    size_t          GetDistributionAt( unsigned char res, size_t n ) const;
    void            IncDistributionAt( unsigned char res, size_t n );

    const float*    GetSqnWeightsAt( size_t p, int st ) const;
    float*          GetSqnWeightsAt( size_t p, int st );
    void            SetSqnWeightsAt( size_t p, int st, size_t pp );
    void            NewSqnWeightsAt( size_t n, int st, size_t size );
    void            FreeSqnWeightsAt( size_t n, int st );
    void            FreeSqnWeightsAt();

    float           GetMatchWeightsAt( unsigned char res, size_t n ) const;
    const float  ( *GetMatchWeightsAt( size_t n ) const )[NUMALPH];
    void            SetMatchWeightsAt( float value, unsigned char res, size_t n );
    void            IncMatchWeightsAt( float value, unsigned char res, size_t n );

    float           GetTransWeightsBeg( int trans ) const;
    const float  ( *GetTransWeightsBeg() const )[P_NSTATES];
    void            SetTransWeightsBeg( float value, int trans );
    void            IncTransWeightsBeg( float value, int trans );

    float           GetTransWeightsAt( int trans, ssize_t n ) const;
    const float  ( *GetTransWeightsAt( ssize_t n ) const )[P_NSTATES];
    void            SetTransWeightsAt( float value, int trans, ssize_t n );
    void            IncTransWeightsAt( float value, int trans, ssize_t n );

    float           GetGapWeightsAt( size_t n ) const;
    void            SetGapWeightsAt( float value, size_t n );


    size_t  GetDistinctHistAt( unsigned char res, ssize_t n ) const {return GetDistinctHistAt(n,PS_M,res);}
    void    SetDistinctHistAt( size_t value, unsigned char res, ssize_t n ) {SetDistinctHistAt(n,PS_M,res,value);}
    void    IncDistinctHistAt( unsigned char res, ssize_t n ) {IncDistinctHistAt(n,PS_M,res);}

    size_t          GetDistinctHistAt( ssize_t n, int st, unsigned char res ) const;
    void            SetDistinctHistAt( ssize_t n, int st, unsigned char res, size_t value );
    void            IncDistinctHistAt( ssize_t n, int st, unsigned char res );

    size_t          GetMSDistinctHistAt( unsigned char res, size_t n ) const;
    void            SetMSDistinctHistAt( size_t value, unsigned char res, size_t n );
    void            IncMSDistinctHistAt( unsigned char res, size_t n );


    float           GetTargetFreqnsAt( unsigned char res, size_t n ) const;
    const float  ( *GetTargetFreqnsAt( size_t n ) const )[NUMALPH];
    void            SetTargetFreqnsAt( float value, unsigned char res, size_t n );
    void            NullTargetFreqnsAt(size_t n );

    float           GetTargetTranstBeg( int trans ) const;
    const float  ( *GetTargetTranstBeg() const )[P_NSTATES];
    void            SetTargetTranstBeg( float value, int trans );
    void            SetTargetTranstBeg( const float(*)[P_NSTATES]);
    void            NullTargetTranstBeg();

    float           GetTargetTranstAt( int trans, ssize_t n ) const;
    const float  ( *GetTargetTranstAt( ssize_t n ) const )[P_NSTATES];
    void            SetTargetTranstAt( float value, int trans, ssize_t n );
    void            SetTargetTranstAt( const float(*)[P_NSTATES], ssize_t n );
    void            NullTargetTranstAt(ssize_t n );

    const float  ( *GetPSSMVector() const )[NUMALPH] { return rawPSSMatrix_; }
    float           GetPSSMEntryAt( unsigned char res, size_t n ) const;
    void            SetPSSMEntryAt( float value, unsigned char res, size_t n );
    const float  ( *GetPSSMVectorAt( size_t n ) const )[NUMALPH];

    float           GetInformationAt( size_t n ) const;
    void            SetInformationAt( float value, size_t n );

    
    float           GetBckPPProbAt( size_t n ) const;
    void            SetBckPPProbAt( float value, size_t n );

    const float*    GetPPProbsAt( size_t n ) const;
    const int*      GetPPPIndsAt( size_t n ) const;
    void            SetPPProbsAt( size_t n, const float* probs, const int* ndxs, int size );

    size_t          GetNoPPProbsAt( size_t n ) const;


    size_t          GetNoSequencesAt( size_t n ) const;
    void            SetNoSequencesAt( size_t value, size_t n );

    float           GetExpNoObservationsAt( size_t n ) const;
    void            SetExpNoObservationsAt( float value, size_t n );


    float           GetMIDExpNoObservationsAt( int n, int st ) const;
    void            SetMIDExpNoObservationsAt( float value, int n, int st );
    void            IncMIDExpNoObservationsAt( float value, int n, int st );
    void            MulMIDExpNoObservationsAt( float value, int n, int st );


    virtual void    clear();//clear all data

    static int      GetSizeOfExpNoDistinctRes() { return nSizeOfExpNObservs; };
    static void     InitPrexpNoDistinctRes( const float* = NULL );
    static void     DestroyPrexpNoDistinctRes();
    static float    GetExpNoObservations( float avgnodistres );

protected:
    explicit        BMProtoProfile();
    virtual void    Realloc( size_t newcap );
    void            ReallocIndices( size_t newcap, size_t n );
    virtual void    Init();

    void            InitRightExtents( size_t from = 0, size_t to = 0 );

private:
    size_t**            indices_;                   //indices of sequences outside the extent
    size_t*             indicesLengths_;            //lengths of indices for each position
    size_t*             indicesCapacities_;         //capacities of indices for each position
    //
    int                *states_;                    //encodes states for each position
    float               backprobs_[NUMALPH];        //background probabilities of this profile model
    float               postprobs_[NUMALPH];        //posterior probabilities of this profile model
    size_t           ( *extents_ )[PS_NSTATES][xCount];//left and right boundaries of extents
    float            ( *nodifrs_ )[PS_NSTATES];     //number of different residues in extents
    size_t           ( *MSextents_ )[xCount];       //match-state version of extents
    size_t*             counts_;                    //number of matched residues (gap inc.) at each position
    size_t           ( *distribution_ )[NUMALPH];   //distribution of residues for each position
    float**             sqnweights_[TIMES2(PS_NSTATES)];//MID sequence weights with allocation mark vectors
    float            ( *matchWeights_ )[NUMALPH];   //match weights of residues
    float            ( *transWeights_ )[P_NSTATES]; //transition weights
    float              *gapWeights_;                //gap weights
    size_t           ( *distinctHist_ )[PS_NSTATES][NUMALPH];//histogram of distinct residues for each position
    size_t           ( *MSdistinctHist_ )[NUMALPH]; //histogram of distinct residues for each match-state position
    float            ( *targetFreqns_ )[NUMALPH];   //estimated frequencies
    float            ( *targetTranst_ )[P_NSTATES]; //estimated transition frequencies
    float            ( *rawPSSMatrix_ )[NUMALPH];   //unscaled PSSM matrix
    float*              bppprob_;                   //background posterior predictive
    float**             ppprobs_;                   //posterior predictive probabilities for each cluster
    int**               pppndxs_;                   //indices of clusters p.p.probabilities were calculated for
    size_t*             noppps_;                    //number of p.p.probability values
    float*              information_;               //information content for each position
    size_t*             noseqs_;                    //number of sequences for each position
    float*              expnoseqs_;                 //expected number of observations (sequences) for each position
    float            ( *expnosmid_ )[PS_NSTATES];   //expected number of observations for each state

    static float*       prexpnores_;                //precomputed expected numbers of (distinct) residues
    static const int nSizeOfExpNObservs;            //size of epxected number of observations
};

////////////////////////////////////////////////////////////////////////////
// Class BMProtoProfile inlines
//
// GetBackProbsAt: get background probability for residue res
//
inline
float BMProtoProfile::GetBackProbsAt( unsigned char res ) const
{
#ifdef __DEBUG__
    if( backprobs_ == NULL || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetBackProbsAt: Memory access error." );
#endif
    return backprobs_[res];
}
// SetBackProbsAt: set background probability for residue res
//
inline
void BMProtoProfile::SetBackProbsAt( unsigned char res, float value )
{
#ifdef __DEBUG__
    if( backprobs_ == NULL || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetBackProbsAt: Memory access error." );
#endif
    backprobs_[res] = value;
}


// -------------------------------------------------------------------------
// GetPostProbsAt: get posterior probability for residue res
//
inline
float BMProtoProfile::GetPostProbsAt( unsigned char res ) const
{
#ifdef __DEBUG__
    if( postprobs_ == NULL || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetPostProbsAt: Memory access error." );
#endif
    return postprobs_[res];
}
// SetPostProbsAt: set posterior probability for residue res
//
inline
void BMProtoProfile::SetPostProbsAt( unsigned char res, float value )
{
#ifdef __DEBUG__
    if( postprobs_ == NULL || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetPostProbsAt: Memory access error." );
#endif
    postprobs_[res] = value;
}


// -------------------------------------------------------------------------
// GetIndicesAt: return array of indices at the given position
//
inline
size_t* BMProtoProfile::GetIndicesAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !indices_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetIndicesAt: Memory access error." );
#endif
    return indices_[n];
}
// GetIndicesSizeAt: return number of entries in the array of indices at 
// position n
//
inline
size_t BMProtoProfile::GetIndicesSizeAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !indicesLengths_ || length_ <= n )
        throw myruntime_error( "BMProtoProfile::GetIndicesSizeAt: Memory access error." );
#endif
    return indicesLengths_[n];
}


// -------------------------------------------------------------------------
// GetStateAt: get state at position n
inline
int BMProtoProfile::GetStateAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !states_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetStateAt: Memory access error." );
#endif
    return states_[n];
}
// SetStateAt: set match state at position n
//
inline
void BMProtoProfile::SetStateAt( size_t n )
{
#ifdef __DEBUG__
    if( states_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetStateAt: Memory access error." );
#endif
    states_[n] = 1;
}


// -------------------------------------------------------------------------
// GetLeftExtentAt: get the left boundary of extent computed at position n
//
inline
size_t BMProtoProfile::GetLeftExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetLeftExtentAt: Memory access error." );
#endif
    return extents_[n][st][xLeft];
}
// SetLeftExtentAt: set left extent boundary at position n
//
inline
void BMProtoProfile::SetLeftExtentAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetLeftExtentAt: Memory access error." );
#endif
    extents_[n][st][xLeft] = value;
}


// -------------------------------------------------------------------------
// GetRightExtentAt: get the right boundary of extent computed at position n
//
inline
size_t BMProtoProfile::GetRightExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetRightExtentAt: Memory access error." );
#endif
        return extents_[n][st][xRight];
}
// SetRightExtentAt: set right extent boundary at position n
//
inline
void BMProtoProfile::SetRightExtentAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetRightExtentAt: Memory access error." );
#endif
    extents_[n][st][xRight] = value;
}


// -------------------------------------------------------------------------
// GetExtentIntervalAt: get the interval size of extent computed at the 
// given position
//
inline
size_t BMProtoProfile::GetExtentIntervalAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetExtentIntervalAt: Memory access error." );
#endif
    return extents_[n][st][xInterval];
}
// SetExtentIntervalAt: set the interval size of extent computed at 
// position n
//
inline
void BMProtoProfile::SetExtentIntervalAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetExtentIntervalAt: Memory access error." );
#endif
    extents_[n][st][xInterval] = value;
}


// -------------------------------------------------------------------------
// GetNoSequencesInExtentAt: get the number of sequences in extent at 
// position n
//
inline
size_t BMProtoProfile::GetNoSequencesInExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetNoSequencesInExtentAt: Memory access error." );
#endif
    return extents_[n][st][xNoSequences];
}
// IncNoSequencesInExtentAt: increment the number of sequences in extent at 
// position n
//
inline
void BMProtoProfile::IncNoSequencesInExtentAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( extents_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncNoSequencesInExtentAt: Memory access error." );
#endif
    extents_[n][st][xNoSequences]++;
}


// -------------------------------------------------------------------------
// GetNoSymbolsInExtentAt: the get number of symbols present in the extent 
// computed at position n
//
inline
float BMProtoProfile::GetNoSymbolsInExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !nodifrs_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetNoSymbolsInExtentAt: Memory access error." );
#endif
    return nodifrs_[n][st];
}
// SetNoSymbolsInExtentAt: set the number of symbols in extent at position n
//
inline
void BMProtoProfile::SetNoSymbolsInExtentAt( size_t n, int st, float value )
{
#ifdef __DEBUG__
    if( nodifrs_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetNoSymbolsInExtentAt: Memory access error." );
#endif
    nodifrs_[n][st] = value;
}



// -------------------------------------------------------------------------
// GetLeftMSExtentAt: get the left boundary of match-state extent 
// computed at position n
//
inline
size_t BMProtoProfile::GetLeftMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetLeftMSExtentAt: Memory access error." );
#endif
    return MSextents_[n][xLeft];
}
// SetLeftMSExtentAt: set the left boundary of match-state extent at 
// position n
//
inline
void BMProtoProfile::SetLeftMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetLeftMSExtentAt: Memory access error." );
#endif
    MSextents_[n][xLeft] = value;
}
// -------------------------------------------------------------------------
// GetRightMSExtentAt: get the right boundary of match-state extent at the 
// given position
//
inline
size_t BMProtoProfile::GetRightMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetRightMSExtentAt: Memory access error." );
#endif
    return MSextents_[n][xRight];
}
// SetRightMSExtentAt: set the right boundary of match-state extent at 
// position n
//
inline
void BMProtoProfile::SetRightMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetRightMSExtentAt: Memory access error." );
#endif
    MSextents_[n][xRight] = value;
}
// -------------------------------------------------------------------------
// GetMSExtentIntervalAt: get the interval size of match-state extent 
// computed at the given position
//
inline
size_t BMProtoProfile::GetMSExtentIntervalAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetMSExtentIntervalAt: Memory access error." );
#endif
    return MSextents_[n][xInterval];
}
// SetMSExtentIntervalAt: set the interval size of match-state extent at 
// position n
//
inline
void BMProtoProfile::SetMSExtentIntervalAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetMSExtentIntervalAt: Memory access error." );
#endif
    MSextents_[n][xInterval] = value;
}
// -------------------------------------------------------------------------
// GetNoSequencesInMSExtentAt: get the number of sequences in match-state 
// extent computed at position n
//
inline
size_t BMProtoProfile::GetNoSequencesInMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetNoSequencesInMSExtentAt: Memory access error." );
#endif
    return MSextents_[n][xNoSequences];
}
// IncNoSequencesInMSExtentAt: increment the number of sequences in 
// match-state extent at position n
//
inline
void BMProtoProfile::IncNoSequencesInMSExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( MSextents_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncNoSequencesInMSExtentAt: Memory access error." );
#endif
    MSextents_[n][xNoSequences]++;
}
// -------------------------------------------------------------------------
// GetNoSymbolsInMSExtentAt: get the number of symbols in match-state 
// extent computed at position n
//
inline
size_t BMProtoProfile::GetNoSymbolsInMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetNoSymbolsInMSExtentAt: Memory access error." );
#endif
    return MSextents_[n][xNoSymbols];
}
// SetNoSymbolsInMSExtentAt: set teh number of symbols in match-state 
// extent at position n
//
inline
void BMProtoProfile::SetNoSymbolsInMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetNoSymbolsInMSExtentAt: Memory access error." );
#endif
    MSextents_[n][xNoSymbols] = value;
}




// -------------------------------------------------------------------------
// ComputeObsFrequencyWeightAt: Compute observed frequency weight known as 
// alpha coefficient;
// it represents the mean value of different residues per column (of extent)
//
inline
float BMProtoProfile::ComputeObsFrequencyWeightAt( size_t n ) const
{
//     size_t  interval = GetExtentIntervalAt( n );
//     size_t  diffsyms = GetNoSymbolsInExtentAt( n );
//     size_t  interval = GetMSExtentIntervalAt( n );
    float   diffsyms = GetNoSymbolsInExtentAt( n, PS_M );
    float   weight = diffsyms;

//     if( !interval || !diffsyms )
//         return 0.0f;
//     //NOTE:already processed over all interval
//     if( 1 < interval )
//         weight /= (float)interval;
    if( weight < 0.0f )
        weight = 0.0f;
    return weight;
}




// -------------------------------------------------------------------------
// GetCountAt: get the number of residues evaluated at position n
//
inline
size_t BMProtoProfile::GetCountAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !counts_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetCountAt: Memory access error." );
#endif
    return counts_[n];
}
// IncCountAt: increase the number of residues observed at the position
//
inline
void BMProtoProfile::IncCountAt( size_t n )
{
#ifdef __DEBUG__
    if( counts_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncCountAt: Memory access error." );
#endif
    counts_[n]++;
}


// -------------------------------------------------------------------------
// GetDistributionAt: get the observed frequency of residue `res' at the 
// given position
//
inline
size_t BMProtoProfile::GetDistributionAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !distribution_ || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetDistributionAt: Memory access error." );
#endif
    return distribution_[n][res];
}
// IncDistributionAt: increase the observed frequency of residue `res' at 
// position n
//
inline
void BMProtoProfile::IncDistributionAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( distribution_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::: Memory access error." );
#endif
    distribution_[n][res]++;
}





// -------------------------------------------------------------------------
// GetSqnWeightsAt: get weights of sequences for position `n' and state `st'
//
inline
const float* BMProtoProfile::GetSqnWeightsAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetSqnWeightsAt: Memory access error." );
#endif
    return sqnweights_[st][n];
}
// GetSqnWeightsAt: get weights of sequences for position `n' and state `st'
//
inline
float* BMProtoProfile::GetSqnWeightsAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetSqnWeightsAt: Memory access error." );
#endif
    return sqnweights_[st][n];
}
// SetSqnWeightsAt: set weights of sequences for position `n' and state `st'
//  equal to those for position `pn'
//
inline
void BMProtoProfile::SetSqnWeightsAt( size_t n, int st, size_t pn )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || n <= pn || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetSqnWeightsAt: Memory access error." );
#endif
    sqnweights_[st][n] = sqnweights_[st][pn];
}
// SetSqnWeightsAt: allocate space for weights of sequences for 
//  position `n' and state `st'
//
inline
void BMProtoProfile::NewSqnWeightsAt( size_t n, int st, size_t size )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::NewSqnWeightsAt: Memory access error." );
#endif
    if( size < 1 )
        return;
    float* wghts = NULL;
    //check if not occupied
    FreeSqnWeightsAt( n, st );
    wghts = sqnweights_[st][n] = ( float* )malloc( sizeof(float) * size );
    if( wghts == NULL )
        throw MYRUNTIME_ERROR( "BMProtoProfile::NewSqnWeightsAt: Not enough memory." );
    //initialize memory and mark the allocated block as belonging to this position;
    //the latter is for deallocation
    memset( wghts, 0, sizeof(float) * size );
    sqnweights_[st+PS_NSTATES][n] = wghts;
}
// FreeSqnWeightsAt: deallocate memory for weights of sequences for 
//  position `n' and state `st'
//
inline
void BMProtoProfile::FreeSqnWeightsAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::FreeSqnWeightsAt: Memory access error." );
#endif
    float* marks = sqnweights_[st+PS_NSTATES][n];
    if( marks == NULL )
        return;
    free( marks );
    sqnweights_[st+PS_NSTATES][n] = NULL;
    sqnweights_[st][n] = NULL;
}





// -------------------------------------------------------------------------
// GetMatchWeightsAt: get match weight computed for residue of type 
// `res' at the given position
//
inline
float BMProtoProfile::GetMatchWeightsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !matchWeights_ || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetMatchWeightsAt: Memory access error." );
#endif
    return matchWeights_[n][res];
}
// GetMatchWeightsAt: get match weight vector at the given position
//
inline
const float ( *BMProtoProfile::GetMatchWeightsAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !matchWeights_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetMatchWeightsAt: Memory access error." );
#endif
    return matchWeights_ + n;
}
// SetMatchWeightsAt: set weight for residue of type `res' at position n
//
inline
void BMProtoProfile::SetMatchWeightsAt( float value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( matchWeights_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetMatchWeightsAt: Memory access error." );
#endif
    matchWeights_[n][res] = value;
}
// IncMatchWeightsAt: increment weight for residue of type `res' at 
// position n
//
inline
void BMProtoProfile::IncMatchWeightsAt( float value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( matchWeights_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncMatchWeightsAt: Memory access error." );
#endif
    matchWeights_[n][res] += value;
}





// -------------------------------------------------------------------------
// GetTransWeightsBeg: get the beginning transition weight: 0th pos.
//
inline
float BMProtoProfile::GetTransWeightsBeg( int trans ) const
{
#ifdef __DEBUG__
    if( !transWeights_ || !length_ || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTransWeightsBeg: Memory access error." );
#endif
    return transWeights_[0][trans];
}
// GetTransWeightsBeg: get the beginning transition weight vector: 0th pos.
//
inline
const float ( *BMProtoProfile::GetTransWeightsBeg() const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( !transWeights_ || !length_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTransWeightsBeg: Memory access error." );
#endif
    return transWeights_;
}
// SetTransWeightsBeg: set the beginning transition weight: 0th pos.
//
inline
void BMProtoProfile::SetTransWeightsBeg( float value, int trans )
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || !length_ || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTransWeightsBeg: Memory access error." );
#endif
    transWeights_[0][trans] = value;
}
// IncTransWeightsBeg: increment the beginning transition weight: 0th pos.
//
inline
void BMProtoProfile::IncTransWeightsBeg( float value, int trans )
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || !length_ || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncTransWeightsBeg: Memory access error." );
#endif
    transWeights_[0][trans] += value;
}


// -------------------------------------------------------------------------
// GetTransWeightsAt: get the transition weight computed at position n; 
// 0th pos. is reserved
//
inline
float BMProtoProfile::GetTransWeightsAt( int trans, ssize_t n ) const
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || (ssize_t)length_+1 <= n+1 || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTransWeightsAt: Memory access error." );
#endif
    if( n < 0 )
        return GetTransWeightsBeg( trans );
    return transWeights_[n+1][trans];
}
// GetTransWeightsAt: get transition weight vector at position n
//
inline
const float ( *BMProtoProfile::GetTransWeightsAt( ssize_t n ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || (ssize_t)length_+1 <= n+1 )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTransWeightsAt: Memory access error." );
#endif
    if( n < 0 )
        return GetTransWeightsBeg();
    return transWeights_ +( n+1 );
}
// SetTransWeightsAt: set transition weight at position n; 0th pos. is 
// reserved
//
inline
void BMProtoProfile::SetTransWeightsAt( float value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || (ssize_t)length_+1 <= n+1 || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTransWeightsAt: Memory access error." );
#endif
    if( n < 0 )
        SetTransWeightsBeg( value, trans );
    else
        transWeights_[n+1][trans] = value;
}
// IncTransWeightsAt: increment transition weight at position n; 
// 0th pos. is reserved
//
inline
void BMProtoProfile::IncTransWeightsAt( float value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( transWeights_ == NULL || (ssize_t)length_+1 <= n+1 || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncTransWeightsAt: Memory access error." );
#endif
    if( n < 0 )
        IncTransWeightsBeg( value, trans );
    else
        transWeights_[n+1][trans] += value;
}





// -------------------------------------------------------------------------
// GetGapWeightsAt: get the gap weight computed at position n
//
inline
float BMProtoProfile::GetGapWeightsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !gapWeights_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetGapWeightsAt: Memory access error." );
#endif
    return gapWeights_[n];
}
// SetGapWeightsAt: set gap weight at position n
//
inline
void BMProtoProfile::SetGapWeightsAt( float value, size_t n )
{
#ifdef __DEBUG__
    if( gapWeights_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetGapWeightsAt: Memory access error." );
#endif
    gapWeights_[n] = value;
}





// -------------------------------------------------------------------------
// GetDistinctHistAt: get the number of occurences of residue `res'
//
inline
size_t BMProtoProfile::GetDistinctHistAt( ssize_t n, int st, unsigned char res ) const
{
#ifdef __DEBUG__
    if( !distinctHist_ || n < -1 || (ssize_t)length_+1 <= n+1 || 
        st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetDistinctHistAt: Memory access error." );
#endif
    return distinctHist_[n+1][st][res];
}
// SetDistinctHistAt: set the number of occurrences of residue `res'
//
inline
void BMProtoProfile::SetDistinctHistAt( ssize_t n, int st, unsigned char res, size_t value )
{
#ifdef __DEBUG__
    if( distinctHist_ == NULL || n < -1 || (ssize_t)length_+1 <= n+1 || 
        st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetDistinctHistAt: Memory access error." );
#endif
    distinctHist_[n+1][st][res] = value;
}
// IncDistinctHistAt: increase the number of occurrences of residue `res'
//
inline
void BMProtoProfile::IncDistinctHistAt( ssize_t n, int st, unsigned char res )
{
#ifdef __DEBUG__
    if( distinctHist_ == NULL || n < -1 || (ssize_t)length_+1 <= n+1 || 
        st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncDistinctHistAt: Memory access error." );
#endif
    distinctHist_[n+1][st][res]++;
}


// -------------------------------------------------------------------------
// GetMSDistinctHistAt: get the number of occurences of residue `res' at a
//  match-state position
//
inline
size_t BMProtoProfile::GetMSDistinctHistAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !MSdistinctHist_ || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetMSDistinctHistAt: Memory access error." );
#endif
    return MSdistinctHist_[n][res];
}
// SetMSDistinctHistAt: set the number of occurrences of residue `res' at a
//  match-state position
//
inline
void BMProtoProfile::SetMSDistinctHistAt( size_t value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( MSdistinctHist_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetMSDistinctHistAt: Memory access error." );
#endif
    MSdistinctHist_[n][res] = value;
}
// IncMSDistinctHistAt: increase the number of occurrences of residue 
// `res' at a match-state position
//
inline
void BMProtoProfile::IncMSDistinctHistAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( MSdistinctHist_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncMSDistinctHistAt: Memory access error." );
#endif
    MSdistinctHist_[n][res]++;
}





// -------------------------------------------------------------------------
// GetTargetFreqnsAt: get the estimated probability of residue of type 
// `res' at position n
//
inline
float BMProtoProfile::GetTargetFreqnsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !targetFreqns_ || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetFreqnsAt: Memory access error." );
#endif
    return targetFreqns_[n][res];
}
// GetTargetFreqnsAt: get target probabilites vector at the given position
//
inline
const float ( *BMProtoProfile::GetTargetFreqnsAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !targetFreqns_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetFreqnsAt: Memory access error." );
#endif
    return targetFreqns_ + n;
}
// SetTargetFreqnsAt: set the estimated probability of residue of type 
// `res' at the given position
//
inline
void BMProtoProfile::SetTargetFreqnsAt( float value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( targetFreqns_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTargetFreqnsAt: Memory access error." );
#endif
    targetFreqns_[n][res] = value;
}
// NullTargetFreqnsAt: make all target frequencies at the given position 
// equal zero
//
inline
void BMProtoProfile::NullTargetFreqnsAt( size_t n )
{
#ifdef __DEBUG__
    if( targetFreqns_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::NullTargetFreqnsAt: Memory access error." );
#endif
    memset( targetFreqns_ + n, 0, sizeof(float) * NUMALPH );
}





// -------------------------------------------------------------------------
// GetTargetTranstBeg: get the estimated beginning transition probability 
// (at 0th pos.)
//
inline
float BMProtoProfile::GetTargetTranstBeg( int trans ) const
{
#ifdef __DEBUG__
    if( !targetTranst_ || !length_ || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetTranstBeg: Memory access error." );
#endif
    return targetTranst_[0][trans];
}
// GetTargetTranstBeg: get the estimated beginning transition probabilites: 
// 0th pos.
//
inline
const float ( *BMProtoProfile::GetTargetTranstBeg() const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( !targetTranst_ || !length_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetTranstBeg: Memory access error." );
#endif
    return targetTranst_;
}
// -------------------------------------------------------------------------
// SetTargetTranstBeg: set the estimated beginning transition probability: 
// 0th pos.
//
inline
void BMProtoProfile::SetTargetTranstBeg( float value, int trans )
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || !length_ || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTargetTranstBeg: Memory access error." );
#endif
    targetTranst_[0][trans] = value;
}
// SetTargetTranstBeg: set the estimated beginning transition 
// probabilities: 0th pos.
//
inline
void BMProtoProfile::SetTargetTranstBeg( const float (*values)[P_NSTATES])
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || values == NULL || *values == NULL || !length_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTargetTranstBeg: Memory access error." );
#endif
    memcpy( targetTranst_[0], *values, sizeof(float) * P_NSTATES );
}
// NullTargetTranstBeg: make all beginning target transition frequencies 
// equal zero: 0th pos.
//
inline
void BMProtoProfile::NullTargetTranstBeg()
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || !length_ )
        throw MYRUNTIME_ERROR( "BMProtoProfile::NullTargetTranstBeg: Memory access error." );
#endif
    memset( targetTranst_, 0, sizeof(float) * P_NSTATES );
}


// -------------------------------------------------------------------------
// GetTargetTranstAt: get estimated transition probability at position n;
// 0th pos. is reserved
//
inline
float BMProtoProfile::GetTargetTranstAt( int trans, ssize_t n ) const
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || (ssize_t)length_+1 <= n+1 || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetTranstAt: Memory access error." );
#endif
    if( n < 0 )
        return GetTargetTranstBeg( trans );
    return targetTranst_[n+1][trans];
}
// GetTargetTranstAt: get estimated transition probabilities at position n;
// 0th pos. is reserved
//
inline
const float ( *BMProtoProfile::GetTargetTranstAt( ssize_t n ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || (ssize_t)length_+1 <= n+1 )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetTargetTranstAt: Memory access error." );
#endif
    if( n < 0 )
        return GetTargetTranstBeg();
    return targetTranst_ + ( n+1 );
}
// -------------------------------------------------------------------------
// SetTargetTranstAt: set estimated transition probability at position n; 
// 0th pos. is reserved
//
inline
void BMProtoProfile::SetTargetTranstAt( float value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || (ssize_t)length_+1 <= n+1 || P_NSTATES <= trans )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTargetTranstAt: Memory access error." );
#endif
    if( n < 0 )
        SetTargetTranstBeg( value, trans );
    else
        targetTranst_[n+1][trans] = value;
}
// SetTargetTranstAt: set estimated transition probabilities at position n;
// 0th pos. is reserved
//
inline
void BMProtoProfile::SetTargetTranstAt( const float (*values)[P_NSTATES], ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || values == NULL || *values == NULL || 
        (ssize_t)length_+1 <= n+1 )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetTargetTranstAt: Memory access error." );
#endif
    if( n < 0 )
        SetTargetTranstBeg( values );
    else
        memcpy( targetTranst_[n+1], *values, sizeof(float) * P_NSTATES );
}
// NullTargetTranstAt: make all target transition frequencies at position n
// equal zero; 0th pos. is reserved
//
inline
void BMProtoProfile::NullTargetTranstAt( ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst_ == NULL || (ssize_t)length_+1 <= n+1 )
        throw MYRUNTIME_ERROR( "BMProtoProfile::NullTargetTranstAt: Memory access error." );
#endif
    if( n < 0 )
        NullTargetTranstBeg();
    else
        memset( targetTranst_ +( n+1 ), 0, sizeof(float) * P_NSTATES );
}





// -------------------------------------------------------------------------
// GetPSSMEntryAt: get PSSM matrix entry given residue type and position
//
inline
float BMProtoProfile::GetPSSMEntryAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !rawPSSMatrix_ || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetPSSMEntryAt: Memory access error." );
#endif
    return rawPSSMatrix_[n][res];
}
// GetPSSMVectorAt: get vector of values at position n from the PSSM matrix
//
inline
const float ( *BMProtoProfile::GetPSSMVectorAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !rawPSSMatrix_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetPSSMVectorAt: Memory access error." );
#endif
    return rawPSSMatrix_ + n;
}
// SetPSSMEntryAt: set PSSM matrix entry given residue type and position
//
inline
void BMProtoProfile::SetPSSMEntryAt( float value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( rawPSSMatrix_ == NULL || length_ <= n || NUMALPH <= res )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetPSSMEntryAt: Memory access error." );
#endif
    rawPSSMatrix_[n][res] = value;
}





// -------------------------------------------------------------------------
// GetInformationAt: get information content calculated for position n
//
inline
float BMProtoProfile::GetInformationAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !information_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetInformationAt: Memory access error." );
#endif
    return information_[n];
}
// SetInformationAt: set information content at position n
//
inline
void BMProtoProfile::SetInformationAt( float value, size_t n )
{
#ifdef __DEBUG__
    if( information_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetInformationAt: Memory access error." );
#endif
    information_[n] = value;
}





// -------------------------------------------------------------------------
// GetBckPPProbAt: get background posterior predictive probability at the 
//  given position
//
inline
float BMProtoProfile::GetBckPPProbAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !bppprob_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::GetBckPPProbAt: Memory access error.");
#endif
    return bppprob_[n];
}

// SetBckPPProbAt: set background posterior predictive probability at the 
//  given position
inline
void BMProtoProfile::SetBckPPProbAt( float value, size_t n )
{
#ifdef __DEBUG__
    if( !bppprob_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::SetBckPPProbAt: Memory access error.");
#endif
    bppprob_[n] = value;
}

// GetPPProbsAt: get posterior predictive probabilities at position n
inline
const float* BMProtoProfile::GetPPProbsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !ppprobs_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::GetPPProbsAt: Memory access error.");
#endif
    return ppprobs_[n];
}

// GetPPPIndsAt: get indices of posterior predictive probabilities at 
//  position n
inline
const int* BMProtoProfile::GetPPPIndsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !pppndxs_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::GetPPPIndsAt: Memory access error.");
#endif
    return pppndxs_[n];
}

// SetPPProbsAt: set posterior predictive probabilities and their 
//  indices at position n
inline
void BMProtoProfile::SetPPProbsAt( size_t n, const float* probs, const int* ndxs, int size )
{
    int k;
    if( probs == NULL || ndxs == NULL )
        throw MYRUNTIME_ERROR("BMProtoProfile::SetPPProbsAt: Null parameters.");
    if( !ppprobs_ || !pppndxs_ || !noppps_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::SetPPProbsAt: Memory access error.");
    if( size < 0 )
        throw MYRUNTIME_ERROR("BMProtoProfile::SetPPProbsAt: Invalid size of data.");
    if( ppprobs_[n]) { free( ppprobs_[n]); ppprobs_[n] = NULL; }
    if( pppndxs_[n]) { free( pppndxs_[n]); pppndxs_[n] = NULL; }
    noppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ppprobs_[n] = ( float* )malloc( size * sizeof(float));
    pppndxs_[n] = ( int* )malloc( size * sizeof( int ));
    if( ppprobs_[n] == NULL || pppndxs_[n] == NULL )
        throw MYRUNTIME_ERROR("BMProtoProfile::SetPPProbsAt: Not enough memory.");
    for( k = 0; k < size; k ++ ) {
        ppprobs_[n][k] = probs[k];
        pppndxs_[n][k] = ndxs[k];
    }
}

// GetNoPPProbsAt: get the number of posterior predictive probabilities at 
// position n
inline
size_t BMProtoProfile::GetNoPPProbsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !noppps_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMProtoProfile::GetNoPPProbsAt: Memory access error.");
#endif
    return noppps_[n];
}





// -------------------------------------------------------------------------
// GetNoSequencesAt: get the number of sequences at position n
//
inline
size_t BMProtoProfile::GetNoSequencesAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !noseqs_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetNoSequencesAt: Memory access error." );
#endif
    return noseqs_[n];
}
// SetNoSequencesAt: set the number of sequences at position n
//
inline
void BMProtoProfile::SetNoSequencesAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( noseqs_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetNoSequencesAt: Memory access error." );
#endif
    noseqs_[n] = value;
}





// -------------------------------------------------------------------------
// GetExpNoObservationsAt: get the expected number of observations at 
// position n
//
inline
float BMProtoProfile::GetExpNoObservationsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !expnoseqs_ || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetExpNoObservationsAt: Memory access error." );
#endif
    return expnoseqs_[n];
}
// SetExpNoObservationsAt: set the expected number of observations at 
// position n
//
inline
void BMProtoProfile::SetExpNoObservationsAt( float value, size_t n )
{
#ifdef __DEBUG__
    if( expnoseqs_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetExpNoObservationsAt: Memory access error." );
#endif
    expnoseqs_[n] = value;
}


// -------------------------------------------------------------------------
// GetMIDExpNoObservationsAt: get the expected number of observations for 
// state `st' at position n
//
inline
float BMProtoProfile::GetMIDExpNoObservationsAt( int n, int st ) const
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || (int)length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::GetMIDExpNoObservationsAt: Memory access error." );
#endif
    //beginning position starts at index 0
    return expnosmid_[n+1][st];
}
// SetMIDExpNoObservationsAt: set the expected number of observations for 
// state `st' at position n
//
inline
void BMProtoProfile::SetMIDExpNoObservationsAt( float value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || (int)length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::SetMIDExpNoObservationsAt: Memory access error." );
#endif
    //beginning position is saved at ndex 0
    expnosmid_[n+1][st] = value;
}
// IncMIDExpNoObservationsAt: add the expected number of observations for 
// state `st' at position n
//
inline
void BMProtoProfile::IncMIDExpNoObservationsAt( float value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || (int)length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::IncMIDExpNoObservationsAt: Memory access error." );
#endif
    //beginning position is saved at index 0
    expnosmid_[n+1][st] += value;
}
// MulMIDExpNoObservationsAt: multiply the expected number of 
// observations for state `st' at position `n' by `value'
//
inline
void BMProtoProfile::MulMIDExpNoObservationsAt( float value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || (int)length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "BMProtoProfile::MulMIDExpNoObservationsAt: Memory access error." );
#endif
    //beginning position is saved at index 0
    expnosmid_[n+1][st] *= value;
}

#endif//__BMProtoProfile_h__
