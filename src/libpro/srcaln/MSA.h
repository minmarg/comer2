/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MSA_h__
#define __MSA_h__

#include "liblib/mybase.h"

#include <stdlib.h>

#include <vector>

#include "liblib/alpha.h"
#include "libpro/srcpro/BMProtoProfile.h"
// #include "libpro/srcpro/PMProfileModel.h"
// #include "libpro/srcpro/PMTransModel.h"
#include "libHDP/HDPbase.h"

//{{functional interface
void MedianValue( size_t* histdata, size_t szdat, float* median, size_t scale = 1 );
//}}

// class AlignmentSimulation;
// class RndMSAGenerator;

namespace pmodel {
class PMProfileModel;
class PMTransModel;
}

struct MyMSAStrDestroyer {
    void operator()(char* p) const {
        std::free(p);
    };
};

// _________________________________________________________________________
// CLASS MSA
// Multiple sequence alignment for profile model inference
//
class MSA
{
public:
    MSA( size_t = scnDefNoSeqs );
    ~MSA();

    void    PlainPreprocess();
    void    SelectSequences();
    void    CalculateXCovMatrices(const char* filename);
    void    ConstructProfile();

    void    CalculateNeff(const std::vector<float>& idnt_thlds, std::vector<float>& neffs, size_t&, size_t&);

    size_t  GetSize() const { return length_; };
    size_t  GetCapacity() const { return capacity_; }

    float   GetIdentityLevel() const { return identity_level_; }
    void    SetIdentityLevel( float value ) { identity_level_ = value; };

    const char* GetName() const { return name_; }
    const char* GetDescription() const { return description_; }

    void    SetName( const char* );
    void    SetDescription( const char* );
                                      //append text to the description
    void    AppendDescription( const char*, size_t = 0, size_t = UINT_MAX );

    void    ReadAlignment( const char* );
    void    ReadFASTA( const char* filename );
    void    ReadSTOCKHOLM1( const char* filename, int* = NULL );
    void    PutAlignment( const char* = NULL );

    bool    GetKeepSeqDescriptions() const { return keepseqdescs_; }
    void    SetKeepSeqDescriptions( bool value ) { keepseqdescs_ = value; }

    void    SetIgnoreGapsInQuery( bool value ) { ignoregapsinquery_ = value; }
    bool    GetIgnoreGapsInQuery() const { return ignoregapsinquery_; }

    void    SetComputeDELETEstates( bool value );
    bool    GetComputeDELETEstates() const { return deletestateson_; }

    void    SetTarFrMix( int value ) { tfrmix_ = value; }
    bool    GetTarFrMixHDPCtx() const { return tfrmix_ == tfrmixHDPCtx; }

    void    SetScoAdjment( int value ) { scoadj_ = value; }
    bool    GetScoAdjmentHDPCtx() const {return scoadj_==scoadjHDPCtx || scoadj_==scoadjHDPsco;}
    void    SetHDPbase( const HDPbase* value ) { HDPbase_ = value; }
    const HDPbase* GetHDPbase() const { return HDPbase_; }
//     void    SetHDPctBase( const HDPbase* value ) { HDPctbase_ = value; }
//     const HDPbase* GetHDPctBase() const { return HDPctbase_; }

    bool    GetUsingSEGFilter() const { return usingsegfilt_; }
    void    SetUsingSEGFilter( bool value ) { usingsegfilt_ = value; }

    size_t  GetSEGWindow() const { return segfiltwinlenval_; }
    void    SetSEGWindow( size_t value ) { segfiltwinlenval_ = value; SetUsingSEGFilter( true ); }

    float   GetSEGLowEntropy() const { return segfiltlowentval_; }
    void    SetSEGLowEntropy( float value ) { segfiltlowentval_ = value; SetUsingSEGFilter( true ); }

    float   GetSEGHighEntropy() const { return segfilthighentval_; }
    void    SetSEGHighEntropy( float value ) { segfilthighentval_ = value;SetUsingSEGFilter( true ); }


    bool    GetUsingSeqSEGFilter() const { return usingseqseg_; }
    void    SetUsingSeqSEGFilter( bool value ) { usingseqseg_ = value; }

    size_t  GetSeqSEGWindow() const { return seqsegwinlenval_; }
    void    SetSeqSEGWindow( size_t value ) { seqsegwinlenval_ = value; SetUsingSeqSEGFilter( true ); }

    float   GetSeqSEGLowEntropy() const { return seqseglowentval_; }
    void    SetSeqSEGLowEntropy( float value ) { seqseglowentval_ = value; SetUsingSeqSEGFilter( true ); }

    float   GetSeqSEGHighEntropy() const { return seqseghighentval_; }
    void    SetSeqSEGHighEntropy( float value ){ seqseghighentval_ = value;SetUsingSeqSEGFilter( true ); }


    size_t  GetExtentMinWindow() const { return extminwindow_; }
    void    SetExtentMinWindow( size_t value ) { extminwindow_ = value; }

    float   GetExtentMinSeqPercentage() const { return extminseqperc_; }
    void    SetExtentMinSeqPercentage( float value ) { extminseqperc_ = value; }

    float   GetPseudoCountWeight() const { return pseudocntweight_; }
    void    SetPseudoCountWeight( float value ) { pseudocntweight_ = value; }


    void    ExportProtoProfile( pmodel::PMProfileModel&, bool = true ) const;
    void    ExportProtoProfile( pmodel::PMTransModel& ) const;

    
    void    OutputWeightedFrequencies( const char* = NULL );
    void    OutputPSSM( const char* = NULL );
    void    OutputSuppressedProfile( const char* = NULL );
    void    OutputProfile( const char* = NULL );

    float   GetEffNoSequences() const { return effnoseqs_; }

    BMSequence* GetSequenceAt( size_t n ) const;

protected:
    void    InitProtoProfile();
    void    PreprocessAlignment();
    void    PurgeAtSequenceIdentity();
    void    PurgeAtSequenceIdentityObs();

    void    RefineWithSEG();
    void    FilterSequencesWithSEG();

    void    SetPosNoSequences();

    void    SetStates();
    void    SetBackgroundProbabilities();
    void    RecalcBackgroundProbabilities();
    void    CalcPosteriorProbabilities();

    void    ComputeExtents();

    void    CalculateExpNoResiduesAt( 
                    size_t p, bool usdgwght, size_t left, size_t right, 
                    const float* wghts, float* express, size_t scale = 10 );
    void    ComputeGlobSequenceWeights( float* );
    void    ComputeMIDstateSequenceWeights();
    void    ComputeMIDstateSequenceWeightsNoExtents();
    void    ComputeGWMIDstateSequenceWeights( float );

    void    ComputeTransitionFrequencies( bool usegwgts, bool expmids = false );

    void    AdjustWeights();
    void    AdjustWeightsAt( size_t, float (*)[NUMALPH] = NULL );

    float   ExpectedNoObservationsAt( size_t pos, int st ) const;
    void    DeriveExpectedNoObservations();
    void    CalculateEffNoSequences();

    void    ComputeTargetTransFrequencies();

    void    CalculateAndPrintXCovMatrices(const char* filename);

    void    ComputeTargetFrequenciesMDL();
    void    ComputeTargetFrequenciesMDLVar();
    void    ComputeTargetFrequencies();
    void    MixTargetFrequenciesHDPCtx();
    void    CalcTFPosteriorPredictives();
    void    ComputePSSM();

//     void            Reset() { counter_ = 0; }
//     bool            Eof() { return length_ <= counter_; }
//     void            Inc() { counter_++; }

//     bool                        NotNull() const { return SequenceAt( counter_ ) != NULL; }
//     bool            NotNull( size_t n ) const { return SequenceAt( n ) != NULL; }
//     BMSequence*       Sequence() const { return SequenceAt( counter_ ); }

    void    DestroyAt( size_t n );

    void    SetSequenceAt( size_t, BMSequence* );
    void    push( BMSequence* );
    void    clear();//clear all sequences in the MSA
    void    Realloc( size_t newcap );

    size_t  GetNoSequences() const { return noseqs_; }
    void    SetNoSequences( size_t value ) { noseqs_ = value; }

    void    SetEffNoSequences( float value ) { effnoseqs_ = value; }

//     BMSequence* NewSequence/*NewPositionVector*/( size_t = scnDefNoPoss ) const;
//     void    DeleteSequence/*DeletePositionVector*/( BMSequence* );

//     BMProtoProfile*  NewProtoProfile/*NewExtendedDescriptionVector*/( size_t = scnDefNoPoss ) const;
//     BMProtoProfile*  NewProtoProfile/*NewExtendedDescriptionVector*/( const BMSequence& ) const;
//     void    DeleteProtoProfile/*DeleteExtendedDescriptionVector*/( BMProtoProfile* );

    void TranslateSTOConsSequence( const mystring&, BMSequence* svs );
    void TranslateSTOStates( const mystring&, BMSequence* svs );

private:
    BMSequence**    sequences_/*alignmentMatrix_*/;//aligned sequences
    BMProtoProfile* protoprofile_/*queryDescription*/;//profile model representing the MSA

    size_t          length_;        //number of positions the MSA consists of
    size_t          capacity_;      //current capacity of the MSA
    float           identity_level_;//sequence identity level; sequences with mutual...
                                    //  similarity above this value will be ignored
    size_t          noseqs_;        //number of sequences in the MSA
    float           effnoseqs_;     //effective number of sequences
    //
//     size_t          counter_;       //initial counter to iterate over all sequences
    //
    char*           name_;          //name of the MSA
    char*           description_/*titletext*/;//description of the MSA

    bool            keepseqdescs_;  //whether to keep descriptions of sequences
    bool            ignoregapsinquery_;//whether to ignore gaps in query sequence
    bool            deletestateson_;//whether to compute delete states

    int             tfrmix_;        //the flag of mixing target frequencies
    int             scoadj_;        //the flag of score adjustment
    const HDPbase*  HDPbase_;       //filled HDP base structure
//     const HDPbase*  HDPctbase_;     //filled auciliary HDP base structure

    bool            usingsegfilt_;  //the flag of using filter
    size_t          segfiltwinlenval_;//window length
    float           segfiltlowentval_;//low entropy value
    float           segfilthighentval_;//high entropy value

    bool            usingseqseg_;   //whether using SEG for sequences in the MSA
    size_t          seqsegwinlenval_;//window length of SEG for sequences
    float           seqseglowentval_;//low entropy value of SEG for sequences
    float           seqseghighentval_;//high entropy value of SEG for sequences

    size_t          extminwindow_;  //minimum required window length when computing extents
    float           extminseqperc_; //minimum required sequence length percentage that an extent must cover

    float           pseudocntweight_;//weight of pseudo count frequencies

    static const size_t scnDefNoSeqs;//default number of sequences
    static const size_t scnDefNoPoss;//default number of positions

    //FRIENDS...
//     friend class AlignmentSimulation;
//     friend class RndMSAGenerator;
};

////////////////////////////////////////////////////////////////////////////
// Class MSA inlines
//
// NewSequence: new sequence
//
// inline
// BMSequence* MSA::NewSequence( size_t initcap ) const
// {
//     return new BMSequence( initcap );
// }
// 
// // DeleteSequence: delete sequence
// //
// inline
// void MSA::DeleteSequence( BMSequence* v )
// {
//     if( v )
//         delete v;
// }
// 
// // new default extended description vector
// 
// inline
// ExtendedDescriptionVector* MSA::NewExtendedDescriptionVector( size_t reserve ) const
// {
//     return new ExtendedDescriptionVector( reserve );
// }
// 
// // new extended description vector
// 
// inline
// ExtendedDescriptionVector* MSA::NewExtendedDescriptionVector( const PosDescriptionVector& s ) const
// {
//     return new ExtendedDescriptionVector( s );
// }
// 
// // deallocates extended description vector
// 
// inline
// void MSA::DeleteExtendedDescriptionVector( ExtendedDescriptionVector* v )
// {
//     if( v )
//         delete v;
// }

// -------------------------------------------------------------------------
// SetComputeDELETEstates: set the flag of computing delete states
//
inline
void MSA::SetComputeDELETEstates( bool value )
{
    if( value )
        SetIgnoreGapsInQuery( !value );
    deletestateson_ = value;
}

// -------------------------------------------------------------------------
// GetSequenceAt: get sequence n in the MSA
//
inline
BMSequence* MSA::GetSequenceAt( size_t n ) const
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw MYRUNTIME_ERROR("MSA::GetSequenceAt: Memory access error.");
#endif
    return sequences_[n];
}
// -------------------------------------------------------------------------
// SetSequenceAt: set/replace sequence n in the MSA
//  NOTE: memory management is not performed!
//
inline
void MSA::SetSequenceAt( size_t n, BMSequence* seq )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw MYRUNTIME_ERROR("MSA::SetSequenceAt: Memory access error.");
#endif
    sequences_[n] = seq;
}

#endif//__MSA_h__
