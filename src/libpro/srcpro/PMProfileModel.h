/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMProfileModel_h__
#define __PMProfileModel_h__

#include "liblib/mybase.h"

#include <stdlib.h>
#include <string.h>

#include <fstream>

#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "libHDP/HDPbase.h"
#include "PMTransModel.h"
#include "PMProfileModelBase.h"

#define MAXFILESIZE 1099511627776UL

// // class Serializer;
// // class GapScheme;

// profile model namespace
namespace pmodel {

// _________________________________________________________________________
// Class PMProfileModel
// profile model
// TODO: change the types of the members name_ and description_ to use 
// mystring
//
class PMProfileModel: public PMProfileModelBase
{
public:
    enum {v2_2_NSECTIONS = 1/*pps2DLen*/ + 1/*pps2DENO*/ + NUMAA/*pps2DBkgPrbs*/ +
        P_NSTATES-2/*ptr2DTrnPrbs*/ + NUMAA/*pmv2DTrgFrqs*/ +
        NUMAA-1/*pmv2DCVentrs*/ + 1/*pmv2DCVprior*/+ 1/*pmv2DCVnorm2*/ +
        SS_NSTATES/*pmv2DSSsprbs*/ + 1/*pmv2DHDP1prb*/ + 1/*pmv2DHDP1ind*/ +
        1/*pmv2Daa*/ + 1/*pmv2DSSstate*/ +
        1/*addrdesc_end_addrs*/ + 1/*addrdesc*/
    };

    // COMER profile format version
    static const char* dataversion;
    static const char* bindataversion;
    // 
    enum {
        ReadWriteScale = 10000,
        DescBufferSize = 4096
    };
    // 
    PMProfileModel();
    virtual ~PMProfileModel();

    static const char* GetDataVersion() { return dataversion; }
    static const char* GetBinaryDataVersion() { return bindataversion; }

    void            Push( const float vals[PVDIM], const float frqs[PVDIM], char, float weight, float info, float[PS_NSTATES]);
    virtual void    PushAt( const float vals[PVDIM], const float frqs[PVDIM], char, float weight, float info, float[PS_NSTATES], int pos );
    virtual void    PushAt( const float vals[PVDIM], const float frqs[PVDIM], char, int pos );
    virtual void    PushAt( const float values[PVDIM], char res, int pos ) { 
                        return PMProfileModelBase::PushAt(values, res, pos ); 
    }

    //{{
    float           GetBckPPProbAt( int n ) const;
    void            SetBckPPProbAt( float value, int n );
    const float*    GetPPProbsAt( int n ) const;
    const int*      GetPPPIndsAt( int n ) const;
    void            SetPPProbsAt( int n, const float* probs, const int* ndxs, int size );
    size_t          GetNoPPProbsAt( int n ) const;
    //}}

    //{{
    bool            GetctPsSet() const { return ctpsset_; }
    void            SetctPsSet( bool value ) { ctpsset_ = value; }
    float           GetctBckPPProbAt( int n ) const;
    void            SetctBckPPProbAt( float value, int n );
    const float*    GetctPPProbsAt( int n ) const;
    const int*      GetctPPPIndsAt( int n ) const;
    void            SetctPPProbsAt( int n, const float* probs, const int* ndxs, int size );
    size_t          GetNoctPPProbsAt( int n ) const;
    //}}

    //{{methods for processing context vectors
    bool            GetCtxVecSet() const { return ctxvecset_; }
    void            SetCtxVecSet( bool value ) { ctxvecset_ = value; }
    float           GetLogPProbFactor() const { return lppfact_; }
    void            CalcLogPProbFactor();
    float           GetCtxVecNorm2At( int n ) const;
    float           GetCtxVecLpprobAt( int n ) const;
    const float*    GetCtxVecAt( int n ) const;
    int             GetCtxVecSize() const { return ctxvecsize_; }
    void            SetCtxVecPlusAt( int n, float norm2, float pprob, const float* vec, int vsize );
    //}}

    //{{
    bool            GetSSSSet() const { return sssset_; }
    void            SetSSSSet( bool value ) { sssset_ = value; }
    bool            GetSSSP3() const { return sssp3_; }
    void            SetSSSP3( bool value ) { sssp3_ = value; }
    char            GetSSStateAt( int n ) const;
    float           GetSSStateProbAt( int n ) const;
    float           GetSSStateProbAt( int n, char st ) const;
    void            SetSSStateAt( int n, char st, float prob );
    void            SetSSStateAt( int n, char st, float[SS_NSTATES]);
    //}}

    void            Finalize();
    void            ConvertToTrgFreqs();
    void            ConvertToScores();

    void            MixTrgFrequencies( const HDPbase* );
    void            CalcTFPostPredictives( const HDPbase* );
    void            Calc_ctTFPostPredictives( const HDPbase* );
    void            CalcCtxVector();

// //     virtual void    Serialize( Serializer& ) const;
// //     virtual void    Deserialize( Serializer& );

    virtual void    Print/*OutputMatrix*/( const char* = NULL ) const;

    size_t          GetNoSequences() const { return nosequences_; }
    void            SetNoSequences( size_t value ) { nosequences_ = value; }

    float           GetEffNoSequences() const       { return effnosequences_; }
    void            SetEffNoSequences( float value ) { effnosequences_ = value; }


    float           GetRefLambda() const            { return referenceLambda_; }
    void            SetRefLambda( float value )     { referenceLambda_ = value; }

    float           GetRefK() const                 { return referenceK_; }
    void            SetRefK( float value )          { referenceK_ = value; }

    float           GetLambda() const               { return lambda_; }
    void            SetLambda( float value )        { lambda_ = value; }

    float           GetEntropy() const              { return entropy_; }
    void            SetEntropy( float value )       { entropy_ = value; }

    float           GetK() const                    { return parameterK_; }
    void            SetK( float value )             { parameterK_ = value; }

    float           GetExpectedScore() const        { return expscore_; }
    void            SetExpectedScore( float value ) { expscore_ = value; }


    const float*    GetBackProbs() const { return backprobs_; }
    float           GetBackProbsAt( char res ) const;
    void            SetBackProbsAt( char res, float value );

    const float*    GetPostProbs() const { return postprobs_; }
    float           GetPostProbsAt( char res ) const;
    void            SetPostProbsAt( char res, float value );


    const float   (*GetTrgFreqsAt( int m ) const )[PVDIM];
    const float   (*GetTrgFreqs() const )[PVDIM] { return GetVector(); }

    const float   (*GetObsFreqsAt( int m ) const )[PVDIM];
    const float   (*GetObsFreqs() const )[PVDIM] { return obsfrqs_; }
    float           GetObsFreqsAt( int m, int a ) const;


    float           GetFrequencyWeightAt( int ) const;//alpha


    const float*    GetInformation() const { return information_; }
    float           GetInformationAt( int ) const;


    float           GetMIDExpNoObservationsBeg( int st ) const { return GetMIDExpNoObservationsAt( -1, st ); }
    float           GetMIDExpNoObservationsAt( int, int st ) const;
    void            SetMIDExpNoObservationsBeg( float vv[PS_NSTATES]) { SetMIDExpNoObservationsAt( -1, vv ); }
    void            SetMIDExpNoObservationsAt( int, float[PS_NSTATES]);


    size_t          GetNameSize() const         { return szname_; }
    size_t          GetDescriptionSize() const  { return szdescription_; }

    const char*     GetName() const         { return name_; }
    const char*     GetDescription() const  { return description_; }

// //     void            SetNameSize( size_t sz )        { szname = sz; }
// //     void            SetDescriptionSize( size_t sz ) { szdescription = sz; }

    void            SetName( const char* );                 //sets name
    void            SetDescription( const char* );          //sets description

    void            PrintAnnotation( char* sp ) const;      //print short annotation to string stream
    void            PrintAnnotation( FILE* fp ) const;      //print short annotation to file
    void            PrintAnnotation( TPrintFunction, void* vpn ) const;//print short annotation

    void            PrintDescriptionFirst( FILE* fp ) const;//format and print to file
    void            PrintDescriptionFirst( TPrintFunction, void* ) const;//format and print description

    void            PrintDescription( char* sp ) const;     //format and print to string stream
    void            PrintDescription( FILE* fp ) const;     //format and print to file
    void            PrintDescription( TPrintFunction, void* vpn ) const;//format and print description

    size_t          GetMaxAnnotationWidth() const;          //maximum annotation width
    size_t          GetMinimumRequiredSizeForDescription() const;//minimum size required to contain description

    void            PrintParameters( FILE* ) const;         //print statistical parameters
    void            PrintParametersCondensed( FILE* ) const;//condensed print

    virtual void    Clear();//clear all object's information

protected:
    float         (*GetTrgFreqs())[PVDIM] { return GetVector(); }
    float         (*GetObsFreqs())[PVDIM] { return obsfrqs_; }
    void            SetObsFreqsAt( int n, const float frqs[PVDIM]);
    float*          GetInformation() { return information_; }


    void    PrintDescriptionHelper(//helper method for formating and printing the description
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width,
        size_t max_rows, size_t max_length, bool annotation ) const;


    virtual void    destroy();              //deallocation of memory and reset of values
    virtual void    reallocate( int size ); //memory allocation
    virtual void    init();                 //initialization

    size_t          GetPrivateBufferSize() const    { return DescBufferSize; }
    char*           GetPrivateBuffer() const        { return private_buffer_; }

    void            FormatBuffer( char*&, const char*,//format buffer
                size_t&, size_t&, const size_t, const size_t, const size_t ) const;

private:
// //     float     (*trgfrqs_)[NUMAA];       //target frequencies
    float     (*obsfrqs_)[PVDIM];       //observed weighted frequencies
    float*      freqweights_;           //frequency weights at each position (known as alpha coefficient)
    float*      information_;           //information content vector
    float     (*expMIDs_)[PS_NSTATES];  //expected number of observations in MID states
    //
    float*      bppprob_;               //background posterior predictive
    float**     ppprobs_;               //posterior predictive probabilities for each global cluster
    int**       pppndxs_;               //indices of clusters p.p.probabilities were calculated for
    size_t*     noppps_;                //number of p.p.probability values
    //
    bool        ctpsset_;               //HDP ctx probabilities set
    float*      ctbppprob_;             //background posterior predictive
    float**     ctppprobs_;             //posterior predictive probabilities for each global cluster
    int**       ctpppndxs_;             //indices of clusters p.p.probabilities were calculated for
    size_t*     noctppps_;              //number of p.p.probability values
    //
    bool        ctxvecset_;             //whether ctx vector calculation is on
    float       lppfact_;               //log of prior probability factor used in processing context vectors
    float*      ctxvecnorm2_;           //squared norm of vector at each position
    float*      ctxveclpprb_;           //log of prior probabilities at each position
    float**     ctxvecs_;               //context vectors at each position
    int         ctxvecsize_;            //size of context vectors
    //
    bool        sssset_;                //SS states are set
    bool        sssp3_;                 //probabilities for all SS states are set
    char*       sss_;                   //SS state
    float     (*sssprob_)[SS_NSTATES];  //probabilities for SS states
    //
    char*       name_;                  //profile name
    char*       description_;           //profile description
    size_t      szname_;                //size of name
    size_t      szdescription_;         //size of description

    char*       private_buffer_;        //buffer to contain description

    size_t      nosequences_;           //number of sequences used in model inference
    float       effnosequences_;        //effective number of sequences

    float       backprobs_[PVDIM];      //background probabilities
    float       postprobs_[PVDIM];      //posterior (generalized target) probabilities

    float       referenceLambda_;       //reference lambda parameter
    float       referenceK_;            //reference parameter K

    float       lambda_;                //computed statistical parameter, lambda
    float       entropy_;               //computed entropy given lambda
    float       parameterK_;            //computed Karlin's parameter K
    float       expscore_;              //expected score per column pair

    friend class RndMSAGenerator;
    friend void BinaryWriteProfile_v2_2(
        std::ofstream& fp,
        size_t* addrfields,
        size_t& addrdesc_end_addrs, size_t& addrdesc,
        const PMProfileModel&, const PMTransModel&);
    friend void BinaryWriteProfile( 
        std::ofstream& fp, const PMProfileModel&, const PMTransModel&);
};

// -------------------------------------------------------------------------

void CalculateAndWriteAddressTable_v2_2(
    std::ofstream& fp,
    const size_t szpreamble,
    size_t* addrfields,
    size_t& addrdesc_end_addrs, size_t& addrdesc,
    const size_t nprofiles, const size_t totalposs );

// -------------------------------------------------------------------------
// Functions for saving/reading profile model in text format
//
// // static const int    gcpDMscaling = 65536;

// // void PrintProfile( const char* filename, const PMProfileModel&, const PMTransModel& );
void TextWriteProfileCondensed( FILE*, const PMProfileModel&, const PMTransModel&, int scale = PMProfileModel::ReadWriteScale );
void TextWriteProfile( FILE*, const PMProfileModel&, const PMTransModel&, int scale = PMProfileModel::ReadWriteScale );
void TextReadProfile( FILE*, PMProfileModel&, PMTransModel& );

// // void OutputProfile( const char* filename, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme& );
// // void TextWriteProfile( FILE*, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, int = gcpDMscaling );
// // void TextReadProfile( FILE*, FrequencyMatrix&, LogOddsMatrix&, GapScheme& );

////////////////////////////////////////////////////////////////////////////
// PMProfileModel INLINES
//
// GetTrgFreqsAt: get the vector of target frequencies at the given position
//
inline
const float ( *PMProfileModel::GetTrgFreqsAt( int m ) const )[PVDIM]
{
    return GetVectorAt(m);
}

// -------------------------------------------------------------------------
// GetObsFreqsAt: get the vector of observed weighted frequencies at the 
// given position
//
inline
const float ( *PMProfileModel::GetObsFreqsAt( int m ) const )[PVDIM]
{
#ifdef __DEBUG__
    if( !obsfrqs_ || GetSize() <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetObsFreqsAt: Memory access error." );
#endif
    return obsfrqs_ + m;
}

// GetValue: get observed weighted frequency for the specified 
// residue at the given position 
//
inline
float PMProfileModel::GetObsFreqsAt( int m, int a ) const
{
#ifdef __DEBUG__
    if( GetSize() <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetObsFreqsAt: Memory access error." );

    if( PVDIM <= a || a < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetObsFreqsAt: Memory access error." );
#endif
    return obsfrqs_[m][a];
}

// SetObsFreqsAt: Set a vector of observed weighted frequencies at the 
// given position
//
inline
void PMProfileModel::SetObsFreqsAt( int m, const float frqs[PVDIM])
{
#ifdef __DEBUG__
    if( !obsfrqs_ || GetSize() <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetObsFreqsAt: Memory access error." );
#endif
    for( int k = 0; k < PVDIM; k ++ )
        obsfrqs_[m][k] = frqs[k];
}

// -------------------------------------------------------------------------
// GetFrequencyWeightAt: get frequency weight at position m
//
inline
float PMProfileModel::GetFrequencyWeightAt( int m ) const
{
#ifdef __DEBUG__
    if( GetSize() <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetFrequencyWeightAt: Memory access error." );
#endif
    return freqweights_[m];
}

// -------------------------------------------------------------------------
// GetInformationAt: get information content at position m
//
inline
float PMProfileModel::GetInformationAt( int m ) const
{
#ifdef __DEBUG__
    if( GetSize() <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetInformationAt: Memory access error." );
#endif
    return information_[m];
}

// -------------------------------------------------------------------------
// GetMIDExpNoObservationsAt: get the expected number of observations for a 
// MID state at the given position
//
inline
float PMProfileModel::GetMIDExpNoObservationsAt( int m, int st ) const
{
#ifdef __DEBUG__
    if( m < -1 || GetSize() <= m || st < 0 || PS_NSTATES <= st )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetMIDExpNoObservationsAt: Memory access error." );
#endif
    return expMIDs_[m+1][st];
}

// SetMIDExpNoObservationsAt: set expected numbers of observations for the 
// MID states at the given position
//
inline
void PMProfileModel::SetMIDExpNoObservationsAt( int m, float mids[PS_NSTATES])
{
#ifdef __DEBUG__
    if( m < -1 || GetSize() <= m )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetMIDExpNoObservationsAt: Memory access error." );
#endif
    expMIDs_[m+1][PS_M] = mids[PS_M];
    expMIDs_[m+1][PS_I] = mids[PS_I];
    expMIDs_[m+1][PS_D] = mids[PS_D];
}

// -------------------------------------------------------------------------
// GetBackProbsAt: get the background probability of residue res
//
inline
float PMProfileModel::GetBackProbsAt( char res ) const
{
#ifdef __DEBUG__
    if( res < 0 || PVDIM <= res )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetBackProbsAt: Memory access error." );
#endif
    return backprobs_[(int)res];
}

// SetBackProbsAt: set a background probability for residue res
//
inline
void PMProfileModel::SetBackProbsAt( char res, float value )
{
#ifdef __DEBUG__
    if( res < 0 || PVDIM <= res )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetBackProbsAt: Memory access error." );
#endif
    backprobs_[(int)res] = value;
}

// -------------------------------------------------------------------------
// GetPostProbsAt: get the posterior probability of residue res
//
inline
float PMProfileModel::GetPostProbsAt( char res ) const
{
#ifdef __DEBUG__
    if( res < 0 || PVDIM <= res )
        throw MYRUNTIME_ERROR( "PMProfileModel::GetPostProbsAt: Memory access error." );
#endif
    return postprobs_[(int)res];
}

// SetPostProbsAt: set a posterior probability for residue res
//
inline
void PMProfileModel::SetPostProbsAt( char res, float value )
{
#ifdef __DEBUG__
    if( res < 0 || PVDIM <= res )
        throw MYRUNTIME_ERROR( "PMProfileModel::SetPostProbsAt: Memory access error." );
#endif
    postprobs_[(int)res] = value;
}

// -------------------------------------------------------------------------
// GetBckPPProbAt: get background posterior predictive probability at the 
// given position
//
inline
float PMProfileModel::GetBckPPProbAt( int n ) const
{
#ifdef __DEBUG__
    if( !bppprob_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetBckPPProbAt: Memory access error.");
#endif
    return bppprob_[n];
}

// SetBckPPProbAt: set background posterior predictive probability at the 
// given position
inline
void PMProfileModel::SetBckPPProbAt( float value, int n )
{
#ifdef __DEBUG__
    if( !bppprob_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::SetBckPPProbAt: Memory access error.");
#endif
    bppprob_[n] = value;
}

// GetPPProbsAt: get posterior predictive probabilities at position n
inline
const float* PMProfileModel::GetPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ppprobs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetPPProbsAt: Memory access error.");
#endif
    return ppprobs_[n];
}

// GetPPPIndsAt: get indices of posterior predictive probabilities at the 
// given position
inline
const int* PMProfileModel::GetPPPIndsAt( int n ) const
{
#ifdef __DEBUG__
    if( !pppndxs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetPPPIndsAt: Memory access error.");
#endif
    return pppndxs_[n];
}

// SetPPProbsAt: set posterior predictive probabilities and their 
// indices at position n
inline
void PMProfileModel::SetPPProbsAt( int n, const float* probs, const int* ndxs, int size )
{
    int k;
    if( !ppprobs_ || !pppndxs_ || !noppps_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::SetPPProbsAt: Memory access error.");
    if( size < 0 )
        throw MYRUNTIME_ERROR("PMProfileModel::SetPPProbsAt: Invalid data size.");
    if( ppprobs_[n]) { free( ppprobs_[n]); ppprobs_[n] = NULL; }
    if( pppndxs_[n]) { free( pppndxs_[n]); pppndxs_[n] = NULL; }
    noppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ppprobs_[n] = ( float* )malloc( size * sizeof(float));
    pppndxs_[n] = ( int* )malloc( size * sizeof(int));
    if( ppprobs_[n] == NULL || pppndxs_[n] == NULL ) {
        if( ppprobs_[n]) { free( ppprobs_[n]); ppprobs_[n] = NULL; }
        if( pppndxs_[n]) { free( pppndxs_[n]); pppndxs_[n] = NULL; }
        throw MYRUNTIME_ERROR("PMProfileModel::SetPPProbsAt: Not enough memory.");
    }
    if( probs == NULL || ndxs == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::SetPPProbsAt: Null parameters.");
    for( k = 0; k < size; k ++ ) {
        ppprobs_[n][k] = probs[k];
        pppndxs_[n][k] = ndxs[k];
    }
}

// GetNoPPProbsAt: get the number of posterior predictive probabilities at 
// position n
inline
size_t PMProfileModel::GetNoPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !noppps_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetNoPPProbsAt: Memory access error.");
#endif
    return noppps_[n];
}

// =========================================================================
// GetctBckPPProbAt: get background posterior predictive probability at 
// position n
//
inline
float PMProfileModel::GetctBckPPProbAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctbppprob_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetctBckPPProbAt: Memory access error.");
#endif
    return ctbppprob_[n];
}

// SetctBckPPProbAt: set background posterior predictive probability at 
// position n
inline
void PMProfileModel::SetctBckPPProbAt( float value, int n )
{
#ifdef __DEBUG__
    if( !ctbppprob_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::SetctBckPPProbAt: Memory access error.");
#endif
    ctbppprob_[n] = value;
}

// GetctPPProbsAt: get posterior predictive probabilities at position n
inline
const float* PMProfileModel::GetctPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctppprobs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetctPPProbsAt: Memory access error.");
#endif
    return ctppprobs_[n];
}

// GetctPPPIndsAt: get indices of posterior predictive probabilities at 
// position n
inline
const int* PMProfileModel::GetctPPPIndsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctpppndxs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetctPPPIndsAt: Memory access error.");
#endif
    return ctpppndxs_[n];
}

// SetctPPProbsAt: set posterior predictive probabilities and their 
// indices at position n
inline
void PMProfileModel::SetctPPProbsAt( int n, const float* probs, const int* ndxs, int size )
{
    int k;
    if( !ctppprobs_ || !ctpppndxs_ || !noctppps_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::SetctPPProbsAt: Memory access error.");
    if( size < 0 )
        throw MYRUNTIME_ERROR("PMProfileModel::SetctPPProbsAt: Invalid data size.");
    if( ctppprobs_[n]) { free( ctppprobs_[n]); ctppprobs_[n] = NULL; }
    if( ctpppndxs_[n]) { free( ctpppndxs_[n]); ctpppndxs_[n] = NULL; }
    noctppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ctppprobs_[n] = ( float* )malloc( size * sizeof(float));
    ctpppndxs_[n] = ( int* )malloc( size * sizeof(int));
    if( ctppprobs_[n] == NULL || ctpppndxs_[n] == NULL ) {
        if( ctppprobs_[n]) { free( ctppprobs_[n]); ctppprobs_[n] = NULL; }
        if( ctpppndxs_[n]) { free( ctpppndxs_[n]); ctpppndxs_[n] = NULL; }
        throw MYRUNTIME_ERROR("PMProfileModel::SetctPPProbsAt: Not enough memory.");
    }
    if( probs == NULL || ndxs == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::SetctPPProbsAt: Null parameters.");
    for( k = 0; k < size; k ++ ) {
        ctppprobs_[n][k] = probs[k];
        ctpppndxs_[n][k] = ndxs[k];
    }
}

// GetNoctPPProbsAt: get the number of posterior predictive 
// probabilities at position n
inline
size_t PMProfileModel::GetNoctPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !noctppps_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetNoctPPProbsAt: Memory access error.");
#endif
    return noctppps_[n];
}

// =========================================================================
// GetCtxVecNorm2At: get the squared norm of the vector at position n
//
inline
float PMProfileModel::GetCtxVecNorm2At( int n ) const
{
#ifdef __DEBUG__
    if( !ctxvecnorm2_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetCtxVecNorm2At: Memory access error.");
#endif
    return ctxvecnorm2_[n];
}

// GetCtxVecLpprobAt: get the prior probability of the vector at position n
inline
float PMProfileModel::GetCtxVecLpprobAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctxveclpprb_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel: GetCtxVecLpprobAt: Memory access error.");
#endif
    return ctxveclpprb_[n];
}

// GetCtxVecAt: get the context vector at position n
inline
const float* PMProfileModel::GetCtxVecAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctxvecs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetCtxVecAt: Memory access error.");
#endif
    return ctxvecs_[n];
}

// SetCtxVecPlusAt: set context vector, its norm and log of prior 
// probability at position n
inline
void PMProfileModel::SetCtxVecPlusAt( int n, 
    float norm2, float pprob, const float* vec, int vsize )
{
    int k;
    if( !ctxvecnorm2_ || !ctxveclpprb_ || !ctxvecs_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::SetCtxVecPlusAt: Memory access error.");
    if( vsize < 0 )
        throw MYRUNTIME_ERROR("PMProfileModel::SetCtxVecPlusAt: Invalid vector size.");
    if( ctxvecs_[n]) { free( ctxvecs_[n]); ctxvecs_[n] = NULL; }
    if( ctxvecsize_ && ctxvecsize_ != vsize )
        throw MYRUNTIME_ERROR("PMProfileModel::SetCtxVecPlusAt: Inconsistent vector size.");
    if( !ctxvecsize_ )
        ctxvecsize_ = vsize;
    if( vsize < 1 )
        return;
    if( vec == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::SetCtxVecPlusAt: Null argument vector.");
    ctxvecs_[n] = ( float* )malloc( vsize * sizeof(float));
    if( ctxvecs_[n] == NULL )
        throw MYRUNTIME_ERROR("PMProfileModel::SetCtxVecPlusAt: Not enough memory.");
    for( k = 0; k < vsize; k++ )
        ctxvecs_[n][k] = vec[k];
    ctxvecnorm2_[n] = norm2;
    ctxveclpprb_[n] = pprob;
}

// =========================================================================
// GetSSStateAt: get the SS state at position n
//
inline
char PMProfileModel::GetSSStateAt( int n ) const
{
#ifdef __DEBUG__
    if( !sss_ || n < 0 || GetSize() <= n )
        throw MYRUNTIME_ERROR("PMProfileModel::GetSSStateAt: Memory access error.");
#endif
    return sss_[n];
}

// GetSSStateProbAt: get the SS state probability at position n
//
inline
float PMProfileModel::GetSSStateProbAt( int n ) const
{
#ifdef __DEBUG__
    if( !sss_ || !sssprob_ || n < 0 || GetSize() <= n || sss_[n] < 0 || SS_NSTATES <= sss_[n] )
        throw MYRUNTIME_ERROR("PMProfileModel::GetSSStateProbAt: Memory access error.");
#endif
    return sssprob_[n][(int)sss_[n]];
}
inline
float PMProfileModel::GetSSStateProbAt( int n, char st ) const
{
#ifdef __DEBUG__
    if( !sssprob_ || n < 0 || GetSize() <= n || st < 0 || SS_NSTATES <= st )
        throw MYRUNTIME_ERROR("PMProfileModel::GetSSStateProbAt: Memory access error.");
#endif
    return sssprob_[n][(int)st];
}

// SetSSStateAt: set SS state and probability at position n
inline
void PMProfileModel::SetSSStateAt( int n, char st, float prob )
{
#ifdef __DEBUG__
    if( !sss_ || !sssprob_ || n < 0 || GetSize() <= n || st < 0 || SS_NSTATES <= st )
        throw MYRUNTIME_ERROR("PMProfileModel::SetSSStateAt: Memory access error.");
#endif
    sss_[n] = st;
    sssprob_[n][(int)st] = prob;
}
inline
void PMProfileModel::SetSSStateAt( int n, char st, float probs[SS_NSTATES])
{
    int s;
#ifdef __DEBUG__
    if( !sss_ || !sssprob_ || n < 0 || GetSize() <= n || st < 0 || SS_NSTATES <= st )
        throw MYRUNTIME_ERROR("PMProfileModel::SetSSStateAt: Memory access error.");
#endif
    sss_[n] = st;
    for( s = 0; s < SS_NSTATES; s++ )
        sssprob_[n][s] = probs[s];
}

}//namespace pmodel

#endif//__PMProfileModel_h__
