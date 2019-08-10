/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMTransModel_h__
#define __PMTransModel_h__

#include "liblib/mybase.h"

#include <stdlib.h>
#include <string.h>

#include <fstream>

#include "TRANSPROBS.h"

//NOTE: COMMENTS LIKE THE ONE BELOW TAG LINES TO BE REMOVED!
// // 

// // #include "datapro.h"
// // #include "Configuration.h"
// // #include "DistributionMatrix.h"

// // #define DEFAULT_AC_CORRECTION   ( 1.2 )

// // class Serializer;

// profile model namespace
namespace pmodel {

class PMProfileModel;

// // enumeration for delete state boundaries
// // enum TDeleteBounds {
// //     DBeg,   //beginning
// //     DEnd,   //end
// //     DSlope, //slope of probability
// //     DCnt    //number of terms
// // };

////////////////////////////////////////////////////////////////////////////
// CLASS PMTransModel
// implements transition probabilities (gap model) inherent in profile model
//
class PMTransModel
{
public:
// //     PMTransModel( double gap_open = DEFAULTGAPOPENCOST, double gap_extend = DEFAULTEXTENDCOST );
    PMTransModel();
    ~PMTransModel();

// //     const Configuration*    GetConfiguration() const { return configuration; }
// //     void        SetConfiguration( const Configuration* config );

    int GetSize/*GetColumns*/() const { return length_; }

// //     bool        GetFixed() const            { return fixed; }
// //     void        SetFixed( bool value )      { fixed = value; }

// //     double      GetOpenCost() const         { return openCost; }
// //     double      GetExtendCost() const       { return extendCost; }

// //     bool        GetUsePosACcorrections() const                  { return useposaccorrect; }
// //     void        SetUsePosACcorrections( bool value )            { useposaccorrect = value; }


// //     double          GetGapProbabFactorEvalue() const            { return gapprobfactevalue; }
// //     void            SetGapProbabFactorEvalue( double value )    { gapprobfactevalue = value; }

// //     double          GetGapProbabFactorWeight() const            { return gapprobfactweight; }
// //     void            SetGapProbabFactorWeight( double value )    { gapprobfactweight = value; }

// //     double          GetGapProbabFactorShift() const             { return gapprobfactshift; }
// //     void            SetGapProbabFactorShift( double value )     { gapprobfactshift = value; }


// //     double      GetAutocorrectionNumerator1st() const           { return acorrnumerator1st; }
// //     void        SetAutocorrectionNumerator1st( double value )   { acorrnumerator1st = value; }

// //     double      GetAutocorrectionNumerator2nd() const           { return acorrnumerator2nd; }
// //     void        SetAutocorrectionNumerator2nd( double value )   { acorrnumerator2nd = value; }

// //     double      GetAutocorrectionLogScale() const               { return acorrlogscale; }
// //     void        SetAutocorrectionLogScale( double value )       { acorrlogscale = value; }

// //     double      GetAutocorrectionDenomScale() const             { return acorrdenomscale; }
// //     void        SetAutocorrectionDenomScale( double value )     { acorrdenomscale = value; }

// //     void        AdjustContextByEval( double evalue, double relent, const double* = NULL );
// //     void        ResetContext();

    static float    GetTransPowIndex( int tr );
    static void     SetTransPowIndex( int tr, float value );

    float           GetOrgTrProbsAt( int trans, int pos ) const;
    const float  ( *GetOrgTrProbsAt( int pos ) const )[P_NSTATES];

    float           GetTransProbsAt( int trans, int pos ) const;
    const float  ( *GetTransProbsAt( int pos ) const )[P_NSTATES];

    float           GetLogTransProbsAt( int trans, int pos ) const;
    const float  ( *GetLogTransProbsAt( int pos ) const )[P_NSTATES];

    float           GetNSLogTransProbsAt( int trans, int pos ) const;
    const float  ( *GetNSLogTransProbsAt( int pos ) const )[P_NSTATES];

// //     double      GetDeleteOpenWeightAt( int pos ) const  { return GetDeletesBegAt( pos ); }
// //     double      GetDeleteOpenProbAt( int pos ) const;
// //     double      GetDeleteExtendWeightAt( int pos, int time ) const;
// //     double      GetDeleteExtensionProbAt( int pos, int time ) const;

// //     void        SetOpenCost( double value )                     { openCost = value; }
// //     void        SetExtendCost( double value )                   { extendCost = value; }

// //     double      GetPosOpenAt( int pos ) const;
// //     double      GetPosExtendAt( int pos ) const;

// //     void        SetPosOpenAt( double value, int pos );
// //     void        SetPosExtendAt( double value, int pos );

    void        Initialize();//initialization

// //     double      GetOpenAt( int ) const;             //gap opening cost at the position
// //     double      GetExtendAt( int ) const;           //gap extend cost at the position
// //     double      GetWeightsAt( int ) const;          //gap weight at the position
// //     double      GetProbabilityAt( int ) const;      //gap probability at the position

// //     double      GetDeletesBegAt( int ) const;       //starting deletion weight at the position
// //     double      GetDeletesEndAt( int ) const;       //final deletion weight at the position
// //     double      GetDeletesSlopeAt( int ) const;     //slope of deletion variability at the position

// //     int         GetDeletesIntervalAt( int ) const;  //interval of deletion variability at the position

// //     double      operator[]( int ind ) const;

    void        Clear();
    void        Push( const float (*)[P_NSTATES]);
    void        PushAt( const float (*)[P_NSTATES], int pos );

// //     void        Push( const double (*)[P_NSTATES], double w, double delbeg, double delend, int delta, char );
// //     void        PushAt( const double (*)[P_NSTATES], double w, double delbeg, double delend, int delta, char, int pos );

    void        SetOrgTrProbsBeg( const float (*)[P_NSTATES]);

    float       GetScoresMultiplier() { return scoremult_; }
    void        SetScoresMultiplier( float value ) { scoremult_ = value; }

    void        Prepare( int scale = 1 );//prepare data for using for alignment
// //     void        Prepare( int );                                     //prepare gap costs for using for alignment
// //     void        Prepare( int, const float* posinfcontents );        //prepare augmented with init of positional values of correct.
// //     void        Prepare( float opencost, float extncost, bool );    //prepare: overloaded

// //                                                                 //adjust gap costs with autocorrelation function
// //     void        AdjustCosts( double*, size_t, const LogOddsMatrix&, int window, bool* masks, size_t, double );

// //     void        Serialize( Serializer& ) const;
// //     void        Deserialize( Serializer& );

// //     bool        IsCompatible( const PMTransModel& one ) const;         //is it compatible with other object?
// //     bool        IsComposIdentical( const PMTransModel& one ) const;    //whether it is compositionally identical

    void        Print/*OutputGapScheme*/( const char* = NULL );

// //     char        AAcid( int ind ) const;

    void        Reserve( int amount ) { reallocate( amount ); }

protected:
    void        destroy();//memory deallocation
    void        reallocate( int size );//memory (re)allocation

// //     void        SetColumns( int col )               { length_ = col; }

    void        ReevaluateTransProbs();

// //     void        InitializePosGapCosts();            //initialize positional gap costs (upper bound)
// //     void        InitializePosGapOpenCosts();
// //     void        InitializePosGapExtendCosts();
// //                                                     //initialize positional values of correction
// //     void        InitializePosACcorrection( const double* posinfcontents, double = -1.0 );

    void        Validate/*CheckIntegrity*/() const;//verify probabilities

// //     void        ResetExtendCost();
// //     void        SetExtendCostByEval( double evalue, double relent );

// //     double      GetACcorrection() const             { return autcorrection; }
// //     void        ResetACcorrection();                //reset alternative correction term for autocorrelation function
// //     void        ResetPosACcorrection();             //reset positional values of correction for autocorrelation function
// //     void        SetACcorrection( double value );
// //     void        SetACcorrectionByEval( double evalue, double relent, const double* = NULL );

// //     double      GetScaledACcorrection() const           { return ac_correction; }
// //     void        SetScaledACcorrection( double value )   { ac_correction = value; }

    int         GetScaleFactor() const { return scalefactor_; }
    void        SetScaleFactor( int value ) { scalefactor_ = value; }

// //     double      GetScaledOpenCost() const           { return scaledopenCost; }
// //     void        SetScaledOpenCost( double value )   { scaledopenCost = value; }

// //     double      GetScaledExtendCost() const         { return scaledextendCost; }
// //     void        SetScaledExtendCost( double value ) { scaledextendCost = value; }

    void        SetOrgTrProbsAt( float value, int trans, int pos );
    void        SetOrgTrProbsAt( const float (*)[P_NSTATES], int pos );

    void        SetTransProbsAt( float value, int trans, int pos );
    void        SetTransProbsAt( const float (*)[P_NSTATES], int pos );

    void        SetLogTransProbsAt( float value, int trans, int pos );
    void        SetLogTransProbsAt( const float (*)[P_NSTATES], int pos );

    void        SetNSLogTransProbsAt( float value, int trans, int pos );
    void        SetNSLogTransProbsAt( const float (*)[P_NSTATES], int pos );

// //     void        SetOpenAt( double value, int pos );
// //     void        SetExtendAt( double value, int pos );
// //     void        SetWeightsAt( double value, int pos );

// //     const double*   GetPosACcorrection() const      { return positionalaccorr; }
// //     double          GetPosACcorrectionAt( int pos ) const;          //get correction for autocorrelation function at position pos
// //     void            SetPosACcorrectionAt( double value, int pos );  //set correction for autocorrelation function at position pos

// //     void        ComputeDeleteSlopes();                      //compute delete-state slopes for all positions

// //     void        SetDeletesBegAt( double value, int );       //set starting deletion weight at the position
// //     void        SetDeletesEndAt( double value, int );       //set final deletion weight at the position
// //     void        ComputeDeletesSlopeAt( int );               //compute slope of deletion variability at the position

// //     void        SetDeletesIntervalAt( int value, int );     //set interval of deletion variability at the position
// //                                                             //set delete-state information at the position
// //     void        SetDeleteStateAt( double beg, double end, int delta, int );

// //                                                     //autocorrelation function
// //     double      Autocorrelation( double*, int window, const double* info, const double* posco );
// //     void        ComputeCosts();                     //compute position-specific gap costs

// //     size_t      GetThickness() const                { return thickness; }
// //     void        SetThickness( size_t value )        { thickness = value; }
// //     void        ResetThickness()                    { thickness = 0; }

// //     bool        GetContextOn() const                { return contextadjusted; }
// //     void        SetContextOn()                      { contextadjusted = true; }
// //     void        SetContextOff()                     { contextadjusted = false; }

// //     double      GetContextEvalue() const            { return contextevalue; }
// //     void        SetContextEvalue( double evalue )   { contextevalue = evalue; }
// //     void        ResetContextEvalue()                { contextevalue = -1.0; }

// //     double      GetProbabilityFactor() const;
// //     static double   SigmoidResponse( double );

private:
    int     length_;                //profile length
    int     allocated_;             //positions allocated
    int     scalefactor_;           //the scale of values
    float   scoremult_;             //multiplier of scores (from scaling of a score system)

////////////////////---
// //     bool    fixed;                  //whether gap costs are fixed at the positions
// //     double* posvec;                 //vector of positional gap opening costs (upper bound)
// //     double* posext;                 //vector of positional gap extension costs (upper bound)

    static float s_trnpowind_[P_NSTATES];//transition power indices
    float (*orgtrobs_)[P_NSTATES];  //estimated (target) transition probabilities
    float (*trnprobs_)[P_NSTATES];  //recalculated transition probabilities
    float (*logtrnps_)[P_NSTATES];  //logs of recalculated transition probabilities
    float (*nsltrnps_)[P_NSTATES];  //non-scaled logs of recalculated transition probabilities

// //     double* vector;                 //vector of gap costs for each position of sequence
// //     double* extvec;                 //vector of gap extend costs for each position of sequence
// //     double* weights;                //vector of position-specific gap weights
// //     double* deletes[DCnt];          //position-specific deletion weights
// //     int*    deleteInt;              //intervals of deletions
// //     int     length;                 //length of gapCostVector
// //     double  openCost;               //default gap opening cost
// //     double  extendCost;             //default cost for extending a gap
// //     double  scaledopenCost;         //scaled gap opening cost
// //     double  scaledextendCost;       //scaled cost for extending a gap
// //     char*   aacids;                 //sequence of amino acids
// //     int     allocated;              //how many positions allocated
// //     int     scalefactor;            //factor to scale gap-cost-related values
// //     double  scoremult_;             //multiplier of scores (from scaling of score system)
// //     const Configuration*
// //             configuration;          //system configuration

// //     double  autcorrection;          //logical value of correction used in autocorrelation function
// //     double  ac_correction;          //scaled correction used in autocorrelation function
// //     double* positionalaccorr;       //positional values of correction used in autocorrelation function

// //     double  gapprobfactevalue;      //evalue threshold for computing of gap probability factor
// //     double  gapprobfactweight;      //argument weight in expression of gap probability factor
// //     double  gapprobfactshift;       //argument shift in expression of gap probability factor

// //     double  acorrnumerator1st;      //numerator of expression to compute 1st-pass autocorrection
// //     double  acorrnumerator2nd;      //numerator of expression to compute 2nd-pass upper bound for autocorrection
// //     double  acorrlogscale;          //logarithmic scale to compute 2nd-pass autocorrection
// //     double  acorrdenomscale;        //denominator scale to compute 2nd-pass autocorrection

// //     size_t  thickness;              //effective thickness in sequences used to construct corresponding profile
// //     double  contextevalue;          //context evalue; -1 in case of reset context
// //     bool    contextadjusted;        //whether a context is adjusted
// //     bool    useposaccorrect;        //whether to use positional values of correction to autocorrelation function

    friend void BinaryWriteProfile(
        std::ofstream& fp, const PMProfileModel&, const PMTransModel&);
};

// =========================================================================
// INLINES...
//
// // // SetConfiguration: sets the configuration and resets parameters
// // //     depending on it
// // //
// // inline
// // void PMTransModel::SetConfiguration( const Configuration* config )
// // {
// //     configuration = config;
// //     ResetContext();
// // }

// // // -------------------------------------------------------------------------
// // // SigmoidResponse: gives sigmoid response corresponding to weight which is
// // //     supposed to be in [0,1]. Range of the returned value is (-.5,.5)
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::SigmoidResponse( double weight )
// // {
// //     static double   wshift = 0.5;
// //     static double   rshift = 0.5;
// //     double          expont = 10.0;
// // 
// //     return 1.0 / ( 1.0 + exp( -expont * ( weight - wshift ))) - rshift;
// // }

// // // -------------------------------------------------------------------------
// // // GetProbabilityFactor: returns factor of probability that depends on
// // //     effective thickness of corresponding profile; return value is
// // //     always in [0,1]
// // //
// // inline
// // double PMTransModel::GetProbabilityFactor() const
// // {
// // // return 1.0;
// //     double  evalthresh = GetGapProbabFactorEvalue();
// //     double  weightfact = GetGapProbabFactorWeight();
// //     double  argshift = GetGapProbabFactorShift();
// // 
// //     double  effthick = GetThickness();
// //     double  evalue = GetContextEvalue();
// // #ifdef __DEBUG__
// //     if( effthick < 0 )
// //         throw myruntime_error( mystring( "PMTransModel: GetProbabilityFactor: Negative effective thickness." ));
// // #endif
// //     if( 0.0 <= evalue && evalue < evalthresh )
// //         return 1.0;
// // 
// //     double  factor01 = 1.0 / ( 1.0 + exp( -weightfact * effthick + argshift ));
// //     return  factor01;
// // }

// // // -------------------------------------------------------------------------
// // // operator[]: used to access a gap opening cost at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::operator[]( int n ) const
// // {
// //     return GetOpenAt( n );
// // }

// // // -------------------------------------------------------------------------
// // // AAcid: used to access amino acid at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // char PMTransModel::AAcid( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || aacids == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return aacids[n];
// // }

// -------------------------------------------------------------------------
// GetTransFactor/SetTransPowIndex: get/set transition probability power 
// index
//
inline
float PMTransModel::GetTransPowIndex( int tr )
{
#ifdef __DEBUG__
    if( tr < 0 || P_NSTATES <= tr )
        throw MYRUNTIME_ERROR("PMTransModel::GetTransPowIndex: Memory access error.");
#endif
    return s_trnpowind_[tr];
}

inline
void PMTransModel::SetTransPowIndex( int tr, float value )
{
#ifdef __DEBUG__
    if( tr < 0 || P_NSTATES <= tr )
        throw MYRUNTIME_ERROR("PMTransModel::SetTransPowIndex: Memory access error.");
#endif
    s_trnpowind_[tr] = value;
}

// -------------------------------------------------------------------------
// GetOrgTrProbsAt: get specified target transition probability at the given
// position
//
inline
float PMTransModel::GetOrgTrProbsAt( int trans, int pos ) const
{
#ifdef __DEBUG__
    if( orgtrobs_ == NULL ||
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetOrgTrProbsAt: Memory access error." );
#endif
    return orgtrobs_[pos+1][trans];
}

// SetOrgTrProbsAt: set specified target transition probability at the given
// position
//
inline
void PMTransModel::SetOrgTrProbsAt( float value, int trans, int pos )
{
#ifdef __DEBUG__
    if( orgtrobs_ == NULL || 
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::SetOrgTrProbsAt: Memory access error." );
#endif
    orgtrobs_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// GetOrgTrProbsAt: get target transition probabilities at position pos
//
inline
const float ( *PMTransModel::GetOrgTrProbsAt( int pos ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( orgtrobs_ == NULL || length_ <= pos || pos < -1 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetOrgTrProbsAt: Memory access error." );
#endif
    return orgtrobs_ + ( pos+1 );
}

// SetOrgTrProbsAt: set target transition probabilities at position pos
//
inline
void PMTransModel::SetOrgTrProbsAt( const float (*values)[P_NSTATES], int pos )
{
#ifdef __DEBUG__
    if( orgtrobs_ == NULL || values == NULL || length_ <= pos || pos < -1 ||
        orgtrobs_[pos+1] == NULL || *values == NULL )
        throw MYRUNTIME_ERROR( "PMTransModel::SetOrgTrProbsAt: Memory access error." );
#endif
    memcpy( orgtrobs_[pos+1], *values, sizeof(float) * P_NSTATES );
}

// SetTransProbsBeg: set beginning target transition probabilities
//
inline
void PMTransModel::SetOrgTrProbsBeg( const float (*values)[P_NSTATES])
{
    SetOrgTrProbsAt( values, -1 );
}

// -------------------------------------------------------------------------
// GetTransProbsAt: get specified recalculated target transition 
//     probability at position pos
//
inline
float PMTransModel::GetTransProbsAt( int trans, int pos ) const
{
#ifdef __DEBUG__
    if( trnprobs_ == NULL ||
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetTransProbsAt: Memory access error." );
#endif
    return trnprobs_[pos+1][trans];
}

// SetTransProbsAt: set specified recalculated target transition 
//     probability at position pos
//
inline
void PMTransModel::SetTransProbsAt( float value, int trans, int pos )
{
#ifdef __DEBUG__
    if( trnprobs_ == NULL || 
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::SetTransProbsAt: Memory access error." );
#endif
    trnprobs_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// GetTransProbsAt: get recalculated target transition probabilities at 
// position pos
//
inline
const float ( *PMTransModel::GetTransProbsAt( int pos ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( trnprobs_ == NULL || length_ <= pos || pos < -1 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetTransProbsAt: Memory access error." );
#endif
    return trnprobs_ + ( pos+1 );
}

// SetTransProbsAt: set recalculated target transition probabilities at 
// position pos
//
inline
void PMTransModel::SetTransProbsAt( const float (*values)[P_NSTATES], int pos )
{
#ifdef __DEBUG__
    if( trnprobs_ == NULL || values == NULL || length_ <= pos || pos < -1 ||
        trnprobs_[pos+1] == NULL || *values == NULL )
        throw MYRUNTIME_ERROR( "PMTransModel::SetTransProbsAt: Memory access error." );
#endif
    memcpy( trnprobs_[pos+1], *values, sizeof(float) * P_NSTATES );
}

// -------------------------------------------------------------------------
// GetLogTransProbsAt: get log of specified recalculated target transition 
// probability at position pos
//
inline
float PMTransModel::GetLogTransProbsAt( int trans, int pos ) const
{
#ifdef __DEBUG__
    if( logtrnps_ == NULL ||
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetLogTransProbsAt: Memory access error." );
#endif
    return logtrnps_[pos+1][trans];
}

// SetLogTransProbsAt: set log of specified target transition probability at
// position pos
//
inline
void PMTransModel::SetLogTransProbsAt( float value, int trans, int pos )
{
#ifdef __DEBUG__
    if( logtrnps_ == NULL || 
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::SetLogTransProbsAt: Memory access error." );
#endif
    logtrnps_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// GetLogTransProbsAt: get log of recalculated target transition 
// probabilities at position pos
//
inline
const float ( *PMTransModel::GetLogTransProbsAt( int pos ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( logtrnps_ == NULL || length_ <= pos || pos < -1 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetLogTransProbsAt: Memory access error." );
#endif
    return logtrnps_ + ( pos+1 );
}

// SetLogTransProbsAt: set log of target transition 
//     probabilities at the position
//
inline
void PMTransModel::SetLogTransProbsAt( const float (*values)[P_NSTATES], int pos )
{
#ifdef __DEBUG__
    if( logtrnps_ == NULL || values == NULL || length_ <= pos || pos < -1 ||
        logtrnps_[pos+1] == NULL || *values == NULL )
        throw MYRUNTIME_ERROR( "PMTransModel::SetLogTransProbsAt: Memory access error." );
#endif
    memcpy( logtrnps_[pos+1], *values, sizeof(float) * P_NSTATES );
}

// -------------------------------------------------------------------------
// GetNSLogTransProbsAt: get non-scaled log of specified recalculated target 
// transition probability at position pos
//
inline
float PMTransModel::GetNSLogTransProbsAt( int trans, int pos ) const
{
#ifdef __DEBUG__
    if( nsltrnps_ == NULL ||
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetNSLogTransProbsAt: Memory access error." );
#endif
    return nsltrnps_[pos+1][trans];
}

// SetNSLogTransProbsAt: set non-scaled log of specified target transition 
// probability at position pos
//
inline
void PMTransModel::SetNSLogTransProbsAt( float value, int trans, int pos )
{
#ifdef __DEBUG__
    if( nsltrnps_ == NULL || 
        length_ <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw MYRUNTIME_ERROR( "PMTransModel::SetNSLogTransProbsAt: Memory access error." );
#endif
    nsltrnps_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// GetNSLogTransProbsAt: get non-scaled log of recalculated target 
// transition probabilities at position pos
//
inline
const float ( *PMTransModel::GetNSLogTransProbsAt( int pos ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( nsltrnps_ == NULL || length_ <= pos || pos < -1 )
        throw MYRUNTIME_ERROR( "PMTransModel::GetNSLogTransProbsAt: Memory access error." );
#endif
    return nsltrnps_ + ( pos+1 );
}

// SetNSLogTransProbsAt: set non-scaled log of target transition 
// probabilities at position pos
//
inline
void PMTransModel::SetNSLogTransProbsAt( const float (*values)[P_NSTATES], int pos )
{
#ifdef __DEBUG__
    if( nsltrnps_ == NULL || values == NULL || length_ <= pos || pos < -1 ||
        nsltrnps_[pos+1] == NULL || *values == NULL )
        throw MYRUNTIME_ERROR( "PMTransModel::SetNSLogTransProbsAt: Memory access error." );
#endif
    memcpy( nsltrnps_[pos+1], *values, sizeof(float) * P_NSTATES );
}

// // // -------------------------------------------------------------------------
// // // GetPosOpenAt: returns positional gap opening cost at the position
// // //     (upper bound)
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetPosOpenAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || posvec == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return posvec[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetPosExtendAt: returns positional gap extension cost at the position
// // //     (upper bound)
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetPosExtendAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || posext == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return posext[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetOpenAt: accesses gap opening cost at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetOpenAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || vector == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return vector[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetExtendAt: accesses gap extend cost at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetExtendAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || extvec == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return extvec[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetWeightsAt: returns weight value at the specified position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetWeightsAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || weights == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return weights[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetDeletesBegAt: returns starting deletion weight at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetDeletesBegAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DBeg] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return deletes[DBeg][n];
// // }

// // // -------------------------------------------------------------------------
// // // GetDeletesEndAt: returns final deletion weight at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetDeletesEndAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DEnd] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return deletes[DEnd][n];
// // }

// // // -------------------------------------------------------------------------
// // // GetDeletesSlopeAt: returns slope of deletion variability at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetDeletesSlopeAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DSlope] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return deletes[DSlope][n];
// // }

// // // -------------------------------------------------------------------------
// // // GetDeletesIntervalAt: returns interval of deletion variability at the
// // //     position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // int PMTransModel::GetDeletesIntervalAt( int n ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deleteInt == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return deleteInt[n];
// // }

// // // -------------------------------------------------------------------------
// // // GetDeleteExtendWeightAt: returns extend cost for deletion at the given
// // //     position and time moment of extend
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double PMTransModel::GetDeleteExtendWeightAt( int n, int time ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 ||
// //         deletes[DBeg] == NULL ||
// //         deletes[DEnd] == NULL ||
// //         deletes[DSlope] == NULL ||
// //         deleteInt == NULL
// //     )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     double  beg = deletes[DBeg][n];
// //     double  end = deletes[DEnd][n];
// //     double  slope = deletes[DSlope][n];
// //     int     delta = deleteInt[n];
// // 
// //     if( delta <= 0 )
// //         return 0.0;
// // 
// //     if( time <= 0 )
// //         return beg;
// // 
// //     if( delta <= time )
// //         return end;
// // 
// //     return slope * ( double ) time + beg;
// // }

// // // -------------------------------------------------------------------------
// // // GetProbabilityAt: returns gap probability at the given position
// // // -------------------------------------------------------------------------
// // inline
// // double PMTransModel::GetProbabilityAt( int pos ) const
// // {
// //     double  factor01 = GetProbabilityFactor();
// //     double  gapweight = GetWeightsAt( pos );
// //     double  corweight = gapweight * factor01;   //corrected weight
// //     double  gapprobab = 2.0 * (( corweight <= 0.5 )? corweight: 1.0 - corweight );
// // #ifdef __DEBUG__
// //     if( gapprobab < 0.0 || 1.0 < gapprobab )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Gap probability inconsistent." ));
// // #endif
// //     return  gapprobab;
// // }

// // // -------------------------------------------------------------------------
// // // GetDeleteOpenProbAt: returns deletion opening probability at the given
// // //     position
// // // -------------------------------------------------------------------------
// // inline
// // double PMTransModel::GetDeleteOpenProbAt( int pos ) const
// // {
// //     double  factor01 = GetProbabilityFactor();
// //     double  opnweight = GetDeleteOpenWeightAt( pos ) * factor01;
// //     double  delopenpr = 2.0 * (( opnweight <= 0.5 )? opnweight: 1.0 - opnweight );
// // #ifdef __DEBUG__
// //     if( delopenpr < 0.0 || 1.0 < delopenpr )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Deletion open probability inconsistent." ));
// // #endif
// //     return  delopenpr;
// // }

// // // -------------------------------------------------------------------------
// // // GetDeleteExtensionProbAt: returns extension probability of deletion at
// // //     the given position and time moment
// // // -------------------------------------------------------------------------
// // inline
// // double PMTransModel::GetDeleteExtensionProbAt( int pos, int time ) const
// // {
// //     double  factor01 = GetProbabilityFactor();
// //     double  extweight = GetDeleteExtendWeightAt( pos, time ) * factor01;
// //     double  delprob = 2.0 * (( extweight <= 0.5 )? extweight: 1.0 - extweight );
// // #ifdef __DEBUG__
// //     if( delprob < 0.0 || 1.0 < delprob )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Deletion extension probability inconsistent." ));
// // #endif
// //     return  delprob;
// // }

// // // -------------------------------------------------------------------------
// // // GetPosACcorrectionAt: gets correction value for autocorrelation
// // //     function at position pos
// // // -------------------------------------------------------------------------
// // inline
// // double PMTransModel::GetPosACcorrectionAt( int pos ) const
// // {
// // #ifdef __DEBUG__
// //     if( length <= pos || pos < 0 || positionalaccorr == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     return positionalaccorr[pos];
// // }

// /////////////////////////////////////////////////////////////////////////
// Set methods
//


// // // -------------------------------------------------------------------------
// // // SetPosOpenAt: sets positional gap opening cost at the position
// // //     (upper bound)
// // // -------------------------------------------------------------------------
// // inline
// // void PMTransModel::SetPosOpenAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || posvec == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     posvec[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetPosExtendAt: sets positional gap extension cost at the position
// // //     (upper bound)
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetPosExtendAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || posext == NULL)
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     posext[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetOpenAt: sets gap open cost at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetOpenAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || vector == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     vector[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetExtendAt: sets cost to extend a gap
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetExtendAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || extvec == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     extvec[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetWeightsAt: sets gap weight at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetWeightsAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || weights == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     weights[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetDeletesBegAt: sets starting deletion weight at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetDeletesBegAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DBeg] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     deletes[DBeg][n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetDeletesEndAt: sets final deletion weight at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetDeletesEndAt( double value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DEnd] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     deletes[DEnd][n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // ComputeDeletesSlopeAt: compute slope of deletion variability at the position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::ComputeDeletesSlopeAt( int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deletes[DSlope] == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     double  beg = GetDeletesBegAt( n );
// //     double  end = GetDeletesEndAt( n );
// //     int     delta = GetDeletesIntervalAt( n );
// // 
// //     if( delta <= 0 ) {
// //         deletes[DSlope][n] = 0.0;
// //         return;
// //     }
// // 
// //     deletes[DSlope][n] = ( end - beg ) / ( double ) delta;
// // }

// // // -------------------------------------------------------------------------
// // // SetDeletesIntervalAt: sets interval of deletion variability at the
// // //     position
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetDeletesIntervalAt( int value, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 || deleteInt == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     deleteInt[n] = value;
// // }

// // // -------------------------------------------------------------------------
// // // SetDeleteStateAt: sets delete-state information at the position including
// // //     slope computation
// // // -------------------------------------------------------------------------
// // 
// // inline
// // void PMTransModel::SetDeleteStateAt( double beg, double end, int delta, int n )
// // {
// // #ifdef __DEBUG__
// //     if( length <= n || n < 0 ||
// //         deletes[DBeg] == NULL ||
// //         deletes[DEnd] == NULL ||
// //         deletes[DSlope] == NULL ||
// //         deleteInt == NULL
// //     )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     deletes[DBeg][n] = beg;
// //     deletes[DEnd][n] = end;
// //     deleteInt[n] = delta;
// // 
// //     if( delta <= 0 ) {
// //         deletes[DSlope][n] = 0.0;
// //         return;
// //     }
// // 
// //     deletes[DSlope][n] = ( end - beg ) / ( double ) delta;
// // }

// // // -------------------------------------------------------------------------
// // // SetPosACcorrectionAt: sets correction value for autocorrelation
// // //     function at position pos
// // // -------------------------------------------------------------------------
// // inline
// // void PMTransModel::SetPosACcorrectionAt( double value, int pos )
// // {
// // #ifdef __DEBUG__
// //     if( length <= pos || pos < 0 || positionalaccorr == NULL )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Memory access error." ));
// // #endif
// //     positionalaccorr[pos] = value;
// // }

// // // -------------------------------------------------------------------------
// // // IsCompatible: verifies whether the vector containing gap opening costs is 
// // //     compatible with another one
// // // -------------------------------------------------------------------------
// // 
// // inline
// // bool PMTransModel::IsCompatible( const PMTransModel& one ) const
// // {
// //     return  GetOpenCost()   == one.GetOpenCost() &&
// //             GetExtendCost() == one.GetExtendCost();
// // }

}//namespace pmodel

#endif//__PMTransModel_h__
