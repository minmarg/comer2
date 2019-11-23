/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <math.h>
#include <cmath>

#include "TRANSPROBS.h"
#include "PMTransModel.h"

//NOTE: COMMENTS LIKE THE ONE BELOW TAG LINES TO BE REMOVED!
// // 

namespace pmodel {

//local constants
static const float  lc_tpm = 1.0f;
static const float  lc_tpi = 0.4f;

//static members
//transition power indices
float PMTransModel::s_trnpowind_[P_NSTATES] = 
    {lc_tpm,lc_tpm,lc_tpm, lc_tpi,lc_tpi,lc_tpi, lc_tpi,lc_tpi,lc_tpi};

// -------------------------------------------------------------------------
// Constructor: 
//
PMTransModel::PMTransModel()
:
    length_( 0 ),
    allocated_( 0 ),
    scalefactor_( 1 ),
    scoremult_( 1.0f ),
    orgtrobs_( NULL ),
    trnprobs_( NULL ),
    logtrnps_( NULL ),
    nsltrnps_( NULL )
{
}

// // // -------------------------------------------------------------------------
// // // default constructor: gap opening costs are initialized to zero, no memory
// // //     allocation is performed
// // // -------------------------------------------------------------------------
// // 
// // PMTransModel::PMTransModel( double gap_open, double gap_extend )
// // :
// //     fixed( false ),
// //     posvec( NULL ),
// //     posext( NULL ),
// // 
// //     orgtrobs_( NULL ),
// //     trnprobs_( NULL ),
// //     logtrnps_( NULL ),
// //     nsltrnps_( NULL ),
// //     vector( NULL ),
// //     extvec( NULL ),
// //     weights( NULL ),
// //     deleteInt( NULL ),
// //     length( 0 ),
// //     openCost( gap_open ),
// //     extendCost( gap_extend ),
// //     scaledopenCost( gap_open ),
// //     scaledextendCost( gap_extend ),
// //     aacids( NULL ),
// //     allocated( 0 ),
// //     scalefactor( 1 ),
// //     scoremult_( 1.0 ),
// //     configuration( NULL ),
// // 
// //     positionalaccorr( NULL ),
// // 
// //     gapprobfactevalue( 1.0 ),
// //     gapprobfactweight( 0.0 ),
// //     gapprobfactshift( 0.0 ),
// // 
// //     acorrnumerator1st( -1.0 ),
// //     acorrnumerator2nd( -1.0 ),
// //     acorrlogscale( -1.0 ),
// //     acorrdenomscale( -1.0 ),
// // 
// //     thickness( 0 ),
// //     contextevalue( -1.0 ),
// //     contextadjusted( false ),
// //     useposaccorrect( false )
// // {
// //     ResetContextEvalue();
// //     for( size_t n = 0; n < DCnt; n++ )
// //         deletes[n] = NULL;
// // }

// -------------------------------------------------------------------------
// Destructor: deallocate memory
//
PMTransModel::~PMTransModel()
{
    destroy();
}

// // // -------------------------------------------------------------------------
// // // IsCompatible: verifies whether the vector containing gap opening costs is 
// // //     compositionally identical to another one
// // // -------------------------------------------------------------------------
// // 
// // bool PMTransModel::IsComposIdentical( const PMTransModel& one ) const
// // {
// //     if( GetColumns() != one.GetColumns())
// //         return false;
// // 
// //     for( int n = 0; n < GetColumns(); n++ )
// //         if( aacids[n] != one.aacids[n] )
// //             return false;
// // 
// //     return true;
// // }

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
//
void PMTransModel::destroy()
{
    if( orgtrobs_ ) { free( orgtrobs_ ); orgtrobs_ = NULL; }
    if( trnprobs_ ) { free( trnprobs_ ); trnprobs_ = NULL; }
    if( logtrnps_ ) { free( logtrnps_ ); logtrnps_ = NULL; }
    if( nsltrnps_ ) { free( nsltrnps_ ); nsltrnps_ = NULL; }

// //     if( posext ) { free( posext ); posext = NULL; }
// //     if( posvec ) { free( posvec ); posvec = NULL; }

// //     if( weights) { free( weights); weights= NULL; }
// //     if( extvec ) { free( extvec ); extvec = NULL; }
// //     if( vector ) { free( vector ); vector = NULL; }
// //     if( aacids ) { free( aacids ); aacids = NULL; }

// //     if( positionalaccorr ) { free( positionalaccorr ); positionalaccorr = NULL; }

// //     if( deleteInt ) { free( deleteInt ); deleteInt = NULL; }

// //     for( size_t n = 0; n < DCnt; n++ )
// //         if( deletes[n] ) {
// //             free( deletes[n] ); deletes[n] = NULL;
// //         }

    length_ = 0;
    allocated_ = 0;
    scalefactor_ = 1;
    scoremult_ = 1.0f;
}

// -------------------------------------------------------------------------
// Clear: clear memory and reset values
//
void PMTransModel::Clear()
{
    if( allocated_ <= 0 )
        return;

    if( orgtrobs_ ) memset( orgtrobs_, 0, sizeof(float) * P_NSTATES * ( allocated_ + 1 ));
    if( trnprobs_ ) memset( trnprobs_, 0, sizeof(float) * P_NSTATES * ( allocated_ + 1 ));
    if( logtrnps_ ) memset( logtrnps_, 0, sizeof(float) * P_NSTATES * ( allocated_ + 1 ));
    if( nsltrnps_ ) memset( nsltrnps_, 0, sizeof(float) * P_NSTATES * ( allocated_ + 1 ));

// //     if( posext )    memset( posext, 0, sizeof( double ) * allocated );
// //     if( posvec )    memset( posvec, 0, sizeof( double ) * allocated );

// //     if( weights)    memset( weights, 0, sizeof( double ) * allocated );
// //     if( extvec )    memset( extvec, 0, sizeof( double ) * allocated );
// //     if( vector )    memset( vector, 0, sizeof( double ) * allocated );
// //     if( aacids )    memset( aacids, 0, sizeof( char ) * allocated );

// //     if( positionalaccorr )  memset( positionalaccorr, 0, sizeof( double ) * allocated );

// //     if( deleteInt ) memset( deleteInt, 0, sizeof( int ) * allocated );

// //     for( size_t n = 0; n < DCnt; n++ )
// //         if( deletes[n] )
// //             memset( deletes[n], 0, sizeof( double ) * allocated );

    length_ = 0;
}

// -------------------------------------------------------------------------
// reallocate: reallocate memory
//
void PMTransModel::reallocate( int size )
{
    float (*tmp_orgtrobs)[P_NSTATES] = NULL;
    float (*tmp_trnprobs)[P_NSTATES] = NULL;
    float (*tmp_logtrnps)[P_NSTATES] = NULL;
    float (*tmp_nsltrnps)[P_NSTATES] = NULL;

    if( size <= allocated_ )
        return;

    if( allocated_ <= 0 ) {
        tmp_orgtrobs = ( float(*)[P_NSTATES] )malloc( sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_trnprobs = ( float(*)[P_NSTATES] )malloc( sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_logtrnps = ( float(*)[P_NSTATES] )malloc( sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_nsltrnps = ( float(*)[P_NSTATES] )malloc( sizeof(float) * P_NSTATES * ( size+1 ));

    } else {
        tmp_orgtrobs = ( float(*)[P_NSTATES] )realloc( orgtrobs_, sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_trnprobs = ( float(*)[P_NSTATES] )realloc( trnprobs_, sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_logtrnps = ( float(*)[P_NSTATES] )realloc( logtrnps_, sizeof(float) * P_NSTATES * ( size+1 ));
        tmp_nsltrnps = ( float(*)[P_NSTATES] )realloc( nsltrnps_, sizeof(float) * P_NSTATES * ( size+1 ));
    }

    if( !tmp_orgtrobs ||
        !tmp_trnprobs ||
        !tmp_logtrnps ||
        !tmp_nsltrnps ) 
    {
        if( tmp_orgtrobs ) { free( tmp_orgtrobs ); tmp_orgtrobs = NULL; }
        if( tmp_trnprobs ) { free( tmp_trnprobs ); tmp_trnprobs = NULL; }
        if( tmp_logtrnps ) { free( tmp_logtrnps ); tmp_logtrnps = NULL; }
        if( tmp_nsltrnps ) { free( tmp_nsltrnps ); tmp_nsltrnps = NULL; }
        throw MYRUNTIME_ERROR( "PMTransModel::reallocate: Not enough memory." );
    }

    orgtrobs_ = tmp_orgtrobs;
    trnprobs_ = tmp_trnprobs;
    logtrnps_ = tmp_logtrnps;
    nsltrnps_ = tmp_nsltrnps;

    if( allocated_ <= 0 ) {
        memset( tmp_orgtrobs, 0, sizeof(float) * P_NSTATES );
        memset( tmp_trnprobs, 0, sizeof(float) * P_NSTATES );
        memset( tmp_logtrnps, 0, sizeof(float) * P_NSTATES );
        memset( tmp_nsltrnps, 0, sizeof(float) * P_NSTATES );
    }
    tmp_orgtrobs++;
    tmp_trnprobs++;
    tmp_logtrnps++;
    tmp_nsltrnps++;

    // fill uninitialized memory with zeros
    memset( tmp_orgtrobs + allocated_,  0, sizeof(float) * P_NSTATES *( size - allocated_ ));
    memset( tmp_trnprobs + allocated_,  0, sizeof(float) * P_NSTATES *( size - allocated_ ));
    memset( tmp_logtrnps + allocated_,  0, sizeof(float) * P_NSTATES *( size - allocated_ ));
    memset( tmp_nsltrnps + allocated_,  0, sizeof(float) * P_NSTATES *( size - allocated_ ));

    allocated_ = size;
}

// // // -------------------------------------------------------------------------
// // // reallocate: memory reallocation for the vector of gap costs and vector
// // //     'aacids'
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::reallocate( int howmuch )
// // {
// //     double(*tmp_orgtrobs )[P_NSTATES];
// //     double(*tmp_trnprobs )[P_NSTATES];
// //     double(*tmp_logtrnps )[P_NSTATES];
// //     double(*tmp_nsltrnps )[P_NSTATES];
// //     int*    tmp_deleteInt;
// //     double* tmp_deletes[DCnt];
// //     double* tmp_weights;
// //     double* tmp_posvec;
// //     double* tmp_posext;
// //     double* tmp_vector;
// //     double* tmp_extvec;
// //     char*   tmp_aacids;
// //     double* tmp_positionalaccorr;
// // 
// //     if( howmuch <= allocated )
// //         return;
// // 
// //     if( allocated <= 0 ) {
// //         tmp_orgtrobs = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_trnprobs = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_logtrnps = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_nsltrnps = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_deleteInt = ( int* )malloc( sizeof( int ) * howmuch );
// // 
// //         tmp_deletes[DBeg]   = ( double* )malloc( sizeof( double ) * howmuch );
// //         tmp_deletes[DEnd]   = ( double* )malloc( sizeof( double ) * howmuch );
// //         tmp_deletes[DSlope] = ( double* )malloc( sizeof( double ) * howmuch );
// // 
// //         tmp_weights = ( double* )malloc( sizeof( double ) * howmuch );
// // 
// //         tmp_posext  = ( double* )malloc( sizeof( double ) * howmuch );
// //         tmp_posvec  = ( double* )malloc( sizeof( double ) * howmuch );
// // 
// //         tmp_extvec  = ( double* )malloc( sizeof( double ) * howmuch );
// //         tmp_vector  = ( double* )malloc( sizeof( double ) * howmuch );
// //         tmp_aacids  = ( char*   )malloc( sizeof( char   ) * howmuch );
// // 
// //         tmp_positionalaccorr = ( double* )malloc( sizeof( double ) * howmuch );
// // 
// //     } else {
// //         tmp_orgtrobs = ( double(*)[P_NSTATES] )realloc( orgtrobs_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_trnprobs = ( double(*)[P_NSTATES] )realloc( trnprobs_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_logtrnps = ( double(*)[P_NSTATES] )realloc( logtrnps_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_nsltrnps = ( double(*)[P_NSTATES] )realloc( nsltrnps_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
// //         tmp_deleteInt = ( int* )realloc( deleteInt, sizeof( int ) * howmuch );
// // 
// //         tmp_deletes[DBeg]   = ( double* )realloc( deletes[DBeg],    sizeof( double ) * howmuch );
// //         tmp_deletes[DEnd]   = ( double* )realloc( deletes[DEnd],    sizeof( double ) * howmuch );
// //         tmp_deletes[DSlope] = ( double* )realloc( deletes[DSlope],  sizeof( double ) * howmuch );
// // 
// //         tmp_weights = ( double* )realloc( weights,  sizeof( double  ) * howmuch );
// // 
// //         tmp_posext  = ( double* )realloc( posext,   sizeof( double  ) * howmuch );
// //         tmp_posvec  = ( double* )realloc( posvec,   sizeof( double  ) * howmuch );
// // 
// //         tmp_extvec  = ( double* )realloc( extvec,   sizeof( double  ) * howmuch );
// //         tmp_vector  = ( double* )realloc( vector,   sizeof( double  ) * howmuch );
// //         tmp_aacids  = ( char*   )realloc( aacids,   sizeof( char    ) * howmuch );
// // 
// //         tmp_positionalaccorr = ( double* )realloc( positionalaccorr, sizeof( double ) * howmuch );
// //     }
// // 
// //     if( !tmp_orgtrobs ||
// //         !tmp_trnprobs ||
// //         !tmp_logtrnps ||
// //         !tmp_nsltrnps ||
// //         !tmp_deleteInt ||
// //         !tmp_deletes[DBeg] ||
// //         !tmp_deletes[DEnd] ||
// //         !tmp_deletes[DSlope] ||
// //         !tmp_weights ||
// //         !tmp_posext || !tmp_posvec ||
// //         !tmp_extvec || !tmp_vector || !tmp_aacids ||
// //         !tmp_positionalaccorr
// //     )
// //         throw myruntime_error( mystring( "PMTransModel: Not enough memory." ));
// // 
// //     orgtrobs_ = tmp_orgtrobs;
// //     trnprobs_ = tmp_trnprobs;
// //     logtrnps_ = tmp_logtrnps;
// //     nsltrnps_ = tmp_nsltrnps;
// //     deleteInt = tmp_deleteInt;
// //     deletes[DBeg]   = tmp_deletes[DBeg];
// //     deletes[DEnd]   = tmp_deletes[DEnd];
// //     deletes[DSlope] = tmp_deletes[DSlope];
// //     weights = tmp_weights;
// //     posext  = tmp_posext;
// //     posvec  = tmp_posvec;
// //     extvec  = tmp_extvec;
// //     vector  = tmp_vector;
// //     aacids  = tmp_aacids;
// //     positionalaccorr = tmp_positionalaccorr;
// // 
// //     if( allocated <= 0 ) {
// //         memset( tmp_orgtrobs, 0, sizeof( double ) * P_NSTATES );
// //         memset( tmp_trnprobs, 0, sizeof( double ) * P_NSTATES );
// //         memset( tmp_logtrnps, 0, sizeof( double ) * P_NSTATES );
// //         memset( tmp_nsltrnps, 0, sizeof( double ) * P_NSTATES );
// //     }
// //     tmp_orgtrobs++;
// //     tmp_trnprobs++;
// //     tmp_logtrnps++;
// //     tmp_nsltrnps++;
// // 
// //     // fill uninitialized memory with zeros
// //     memset( tmp_orgtrobs + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
// //     memset( tmp_trnprobs + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
// //     memset( tmp_logtrnps + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
// //     memset( tmp_nsltrnps + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
// //     memset( deleteInt + allocated,  0, sizeof( int ) * ( howmuch - allocated ));
// // 
// //     memset( deletes[DBeg]   + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// //     memset( deletes[DEnd]   + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// //     memset( deletes[DSlope] + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// // 
// //     memset( weights + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// // 
// //     memset( posext  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// //     memset( posvec  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// // 
// //     memset( extvec  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// //     memset( vector  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// //     memset( aacids  + allocated,  0, sizeof( char   ) * ( howmuch - allocated ));
// // 
// //     memset( positionalaccorr + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
// // 
// //     allocated = howmuch;
// // }

// -------------------------------------------------------------------------
// Initialize: initialize
//
void PMTransModel::Initialize()
{
// //     InitializePosGapCosts();
// //     ComputeDeleteSlopes();
    ReevaluateTransProbs();
}

// // // -------------------------------------------------------------------------
// // // Initialize: initializes the gap opening and extension vectors
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::InitializePosGapCosts()
// // {
// //     InitializePosGapOpenCosts();
// //     InitializePosGapExtendCosts();
// // }

// // // InitializePosGapOpenCosts: positional gap open costs
// // //
// // void PMTransModel::InitializePosGapOpenCosts()
// // {
// //     for( int n = 0; n < GetColumns(); n++ )
// //         SetPosOpenAt( GetScaledOpenCost(), n );
// // }

// // // InitializePosGapExtendCosts: positional gap extend costs
// // //
// // void PMTransModel::InitializePosGapExtendCosts()
// // {
// //     for( int n = 0; n < GetColumns(); n++ )
// //         SetPosExtendAt( GetScaledExtendCost(), n );
// // }

// -------------------------------------------------------------------------
// ReevaluateTransProbs: reevaluate transition probabilities
//
void PMTransModel::ReevaluateTransProbs()
{
    int n, s;
    float val, logval, logfct;
    int     scale = GetScaleFactor();
    //multiplier of scores in used the procedure of obtaining required lambda (cbs)
    float   mult = GetScoresMultiplier();

    if( scale <= 0 )
        scale = 1;
    //if( mult <= 0.0f )///cbs does not affect log values when commented
    mult = 1.0f;

    for( n = -1; n < GetSize(); n++ ) {
        for( s = 0; s < P_NSTATES; s++ ) {
            val = GetOrgTrProbsAt( s, n );
            logval = LOG_PROB_MIN;
            if( 0.0f < val ) {
                logval = logf(val) * mult;
                logfct = GetTransPowIndex( s );
                if( P_MM + gTPTRANS_NTPS <= s && logfct != 1.0f )
                    logval *= logfct;
                val = expf( logval );
            }
            else
                val = 0.0f;
            SetTransProbsAt( val, s, n );
            SetLogTransProbsAt( logval * scale, s, n );
            SetNSLogTransProbsAt( logval, s, n );
        }
    }
}

// // // -------------------------------------------------------------------------
// // // InitializePosACcorrection: initialize positional values of correction for
// // //     autocorrelation function given array of information contents
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::InitializePosACcorrection( const double* posinfcontents, double correct )
// // {
// //     double  numerator = GetAutocorrectionNumerator1st();
// //     double  defaultac = DEFAULT_AC_CORRECTION * 10.0;
// //     double  deviation = 0.25;
// //     double  value, entropy, sysrelent;
// // 
// // //     if( posinfcontents == NULL )
// // //         throw myruntime_error( mystring( "PMTransModel: InitializePosACcorrection: Null address of entropies." ));
// // 
// //     defaultac += deviation;
// // 
// //     if( GetConfiguration()) {
// //         sysrelent = GetConfiguration()[ProcomUngapped].GetH();
// //         if( 0.0 < sysrelent )
// //             defaultac = numerator / sqrt( sysrelent );
// //     }
// // 
// //     for( int n = 0; n < GetColumns(); n++ ) {
// //         entropy = -1.0;
// //         if( posinfcontents )
// //             entropy = posinfcontents[n];
// //         if( 0.0 < entropy ) {
// //             if( 0.0 < correct )
// //                 value = correct / sqrt( entropy );
// //             else
// //                 value = numerator / sqrt( entropy );
// //         } else
// //             value = defaultac;
// //         SetPosACcorrectionAt( value * ( double )GetScaleFactor(), n );
// //     }
// // }

// // // -------------------------------------------------------------------------
// // // ResetACcorrection: resets correction for autocorrelation function to its
// // //     default value
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ResetACcorrection()
// // {
// //     double  correct = DEFAULT_AC_CORRECTION;
// //     double  sysrelent = 1.0;
// //     double  deviation = 0.25;
// // 
// //     if( GetAutocorrectionNumerator1st() < 0.0 )
// //         throw myruntime_error( mystring( "PMTransModel: Invalid numerator for autocorrection." ));
// // 
// //     if( GetConfiguration()) {
// //         sysrelent = GetConfiguration()[ProcomUngapped].GetH();
// //         correct = sysrelent + deviation;
// //         if( 0.0 < sysrelent )
// //             correct = GetAutocorrectionNumerator1st() / sqrt( sysrelent );
// //     }
// // 
// //     SetACcorrection( correct );
// // }

// // // -------------------------------------------------------------------------
// // // ResetPosACcorrection: resets positional values of correction for
// // //     autocorrelation function
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ResetPosACcorrection()
// // {
// //     double  numerator = GetAutocorrectionNumerator1st();
// //     double  defaultac = DEFAULT_AC_CORRECTION * 10.0;
// //     double  deviation = 0.25;
// //     double  sysrelent;
// // 
// //     defaultac += deviation;
// // 
// //     if( GetConfiguration()) {
// //         sysrelent = GetConfiguration()[ProcomUngapped].GetH();
// //         if( 0.0 < sysrelent )
// //             defaultac = numerator / sqrt( sysrelent );
// //     }
// // 
// //     defaultac *= ( double )GetScaleFactor();
// // 
// //     for( int n = 0; n < GetColumns(); n++ ) {
// //         SetPosACcorrectionAt( defaultac, n );
// //     }
// // }

// // // -------------------------------------------------------------------------
// // // SetACcorrection: sets correction for autocorrelation function
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::SetACcorrection( double value )
// // {
// //     autcorrection = value;
// //     SetScaledACcorrection( autcorrection * ( double )GetScaleFactor());
// // }

// // // -------------------------------------------------------------------------
// // // SetACcorrectionByEval: sets correction for autocorrelation function
// // //     according to the given evalue
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::SetACcorrectionByEval( double evalue, double relent, const double* posinfcontents )
// // {
// //     double  correct = 0.0;
// //     double  sysrelent = 1.0;
// //     double  evalbound = 1e-5;   //used to be computed
// //     double  logscale = GetAutocorrectionDenomScale();
// //     double  evalscale = GetAutocorrectionLogScale();    //not a typo!
// //     double  upper = GetAutocorrectionNumerator2nd();
// // 
// //     if( GetAutocorrectionNumerator2nd() < 0.0 ||
// //         GetAutocorrectionLogScale() < 0.0 ||
// //         GetAutocorrectionDenomScale() < 0.0 )
// //         throw myruntime_error( mystring( "PMTransModel: Invalid parameter values for autocorrection on 2nd pass." ));
// // 
// //     if( GetConfiguration())
// //         sysrelent = GetConfiguration()[ProcomUngapped].GetH();
// // 
// //     if( 0.0 < sysrelent )
// //         upper /= sqrt( sysrelent );
// // //     if( 0.0 < relent )
// // //         upper /= sqrt( relent );
// // 
// // 
// //     evalbound = exp( -evalscale - 1.0 / ( logscale * upper ));
// //     correct = upper;
// // 
// //     if( evalbound < evalue )
// //         //avoid phase transition of function
// //         correct = upper;
// //     else
// //         if( 0.0 < evalue ) {
// //             correct = -1.0 / (( log( evalue ) + evalscale ) * logscale );
// //             if( upper < correct )
// //                 correct = upper;
// //             else
// //                 if( posinfcontents )
// //                     InitializePosACcorrection( posinfcontents, correct );
// //         }
// // 
// //     SetACcorrection( correct );
// // }

// // // -------------------------------------------------------------------------
// // // ResetExtendCost: resets extend cost to its initial value
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ResetExtendCost()
// // {
// //     SetScaledExtendCost( GetExtendCost() * GetScaleFactor());
// //     InitializePosGapExtendCosts();
// // }

// // // -------------------------------------------------------------------------
// // // SetExtendCostByEval: determines extend cost according to e-value given
// // // TODO: requires more investigation
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::SetExtendCostByEval( double evalue, double relent )
// // {
// //     if( evalue < 0.0 )
// //         return;
// // /// see todo; constants should be reestimated
// //     return;
// // 
// //     double  numerator = -30.0;
// //     double  maxcost = -TIMES2( GetExtendCost());
// //     double  mincost = 0.0;//-GetExtendCost();
// //     double  cost = mincost;
// // 
// //     if( 0.0 < evalue ) {
// //         cost += numerator / log( evalue );
// //         if( maxcost < cost || cost < 0.0 )
// //             cost = maxcost;
// //     }
// // 
// //     SetScaledExtendCost( -cost * GetScaleFactor());
// //     InitializePosGapExtendCosts();
// // }

// // // -------------------------------------------------------------------------
// // // AdjustContextByEval: adjusts context parameters according to e-value
// // //     given
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::AdjustContextByEval( double evalue, double relent, const double* posinfcontents )
// // {
// //     SetACcorrectionByEval( evalue, relent, posinfcontents );
// //     SetExtendCostByEval( evalue, relent );
// //     SetContextOn();
// //     SetContextEvalue( evalue );
// // }

// // // -------------------------------------------------------------------------
// // // ResetContext: resets context parameters
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ResetContext()
// // {
// //     if( GetContextOn())
// //         ResetPosACcorrection();
// //     ResetACcorrection();
// //     ResetExtendCost();
// //     ResetThickness();
// //     SetContextOff();
// // }

// -------------------------------------------------------------------------
// Prepare: prepare data for using it for alignment
//
void PMTransModel::Prepare( int scale )
{
    SetScaleFactor( scale );
    ReevaluateTransProbs();
}

// // // -------------------------------------------------------------------------
// // // Prepare: prepares gap costs for using for alignments
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::Prepare( int scfact )
// // {
// //     SetScaleFactor( scfact );
// //     SetScaledOpenCost( GetOpenCost() * GetScaleFactor());
// //     if( !GetContextOn()) {
// //         SetScaledExtendCost( GetExtendCost() * GetScaleFactor());
// //         SetACcorrection( GetACcorrection());
// //         InitializePosGapExtendCosts();
// //         ResetContextEvalue();
// //     }
// //     InitializePosGapOpenCosts();
// // //     InitializePosGapCosts();
// //     ComputeCosts();
// // }

// // // Prepare: overloaded
// // //
// // void PMTransModel::Prepare( int scfact, const double* posinfcontents )
// // {
// //     Prepare( scfact );
// //     if( !GetContextOn())
// //         InitializePosACcorrection( posinfcontents );
// //     ReevaluateTransProbs();
// // }

// // // Prepare: overloaded
// // //
// // void PMTransModel::Prepare( double opencost, double extncost, bool fixedcosts )
// // {
// //     SetFixed( fixedcosts );
// //     SetOpenCost( opencost );
// //     SetExtendCost( extncost );
// //     Prepare(( int ) 1 );
// //     ReevaluateTransProbs();
// // }

// // // -------------------------------------------------------------------------
// // // ComputeDeleteSlopes: computes delete-state slopes for all positions
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ComputeDeleteSlopes()
// // {
// //     for( int n = 0; n < GetColumns(); n++ )
// //         ComputeDeletesSlopeAt( n );
// // }

// -------------------------------------------------------------------------
// Push: push transition probabilities
//
void PMTransModel::Push( const float (*tps)[P_NSTATES])
{
    if( allocated_ <= length_ ) {
        int newcap = TIMES2( allocated_ );
        if( newcap <= length_ )
            newcap = length_ + 1;
        reallocate( newcap );
    }
    PushAt( tps, GetSize());
}

// -------------------------------------------------------------------------
// PushAt: assign transition probabilities at the given position
//
void PMTransModel::PushAt( const float (*tps)[P_NSTATES], int position )
{
    if( allocated_ <= position )
        throw MYRUNTIME_ERROR( "PMTransModel::PushAt: Memory access error." );

    if( length_ <= position )
        length_ = position + 1;

    SetOrgTrProbsAt( tps, position );
}

// // // -------------------------------------------------------------------------
// // // Push: push values related to the gap processing scheme
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::Push( const double( *tps )[P_NSTATES], double w, double delbeg, double delend, int delta, char aa )
// // {
// //     if( allocated < GetColumns() + 1 ) {
// //         reallocate( allocated + allocated + 1 );
// //     }
// //     PushAt( tps, w, delbeg, delend, delta, aa, GetColumns());
// // }

// // // -------------------------------------------------------------------------
// // // PushAt: pushes values at the given position
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::PushAt( const double( *tps )[P_NSTATES], double w, double delbeg, double delend, int delta, 
// //                         char aa, int position )
// // {
// //     if( allocated < position + 1 )
// //         throw myruntime_error( "PMTransModel: Memory access error." );
// // 
// //     weights[position] = w;
// //     aacids [position] = aa;
// // 
// //     if( GetColumns() < position + 1 )
// //         SetColumns( position + 1 );
// // 
// //     SetOrgTrProbsAt( tps, position );
// //     SetDeleteStateAt( delbeg, delend, delta, position );
// // }

// // // -------------------------------------------------------------------------
// // // Autocorrelation: computes normalized sum of autocorrelations given
// // //     vector of values; i.e. computes all autocorrelation function values
// // //     by varying autocorrelation window size.
// // // -------------------------------------------------------------------------
// // 
// // double PMTransModel::Autocorrelation( double* maxs, int wsize, const double* info, const double* posco )
// // {
// //     if( maxs == NULL || info == NULL || posco == NULL )
// //         throw myruntime_error( mystring( "PMTransModel: Memory access error." ));
// // 
// // #ifdef __DEBUG__
// //     if( wsize <= 0 )
// //         throw myruntime_error( mystring( "PMTransModel: Invalid window paramaters." ));
// // #endif
// // 
// //     double  ro = 0.0;                       //sum of autocorrelations
// //     double  ac = GetScaledACcorrection();   //correction for autocorrelation function
// //     double  aci, acj, dmi, dmj;
// //     double  dmlim = 100.0;                  //1/0.1^2
// //     double  icorrect = 0.00;
// //     double  scale = 0.4;
// //     size_t  no_pairs = 0;   //number of pairs, which is wsize! / ((wsize-2)! 2! )
// //     int     pribeg = 0;
// // 
// //     aci = acj = ac;
// // 
// //     if( wsize == 1 ) {
// //         if( GetUsePosACcorrections()) {
// //             aci = acj = posco[0];
// //         }
// // //         aci = GetScaleFactor() * 1.5 * ( icorrect - info[0] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
// // //         aci = GetScaleFactor() * ( scale - ( scale * exp( scale * ( maxs[0] / GetScaleFactor()))));
// //         ro += ( maxs[0] + aci )*( maxs[0] + aci );
// //         no_pairs++;
// //     }
// //     else
// //         for( int priwin = wsize; 0 < priwin; priwin-- ) {
// //             int winpar = ( priwin  & 1 ) ^ 1;
// //             int midwin = ( priwin >> 1 ) - winpar;
// // 
// //             for( int i = pribeg; i < pribeg + priwin; i++ ) {
// //                 int j = i + midwin;
// //                 if( priwin <= j )
// //                     j -= priwin;
// // 
// //                 if( GetUsePosACcorrections()) {
// //                     aci = posco[i];
// //                     acj = posco[j];
// //                 }
// // 
// // //                 dmi = maxs[i] / GetScaleFactor(); dmj = maxs[j] / GetScaleFactor();
// // //                 if( 0.0 < dmi ) dmi = 1.0 / ( dmi * 0.4 ); else dmi = dmlim;
// // //                 if( 0.0 < dmj ) dmj = 1.0 / ( dmj * 0.4 ); else dmj = dmlim;
// // //                 dmi *= GetScaleFactor(); dmj *= GetScaleFactor();
// // 
// // //                 ro += ( maxs[i] + dmi + aci )*( maxs[j] + dmj + acj );
// // //                 ro += ( TIMES2( maxs[i] )+ aci )*( TIMES2( maxs[j] ) + acj );
// //                 ro += ( maxs[i] + aci )*( maxs[j] + acj );//original
// //                 no_pairs++;
// // 
// // //             for( int j = i + 1; j < wsize; j++ ) {
// // //                 if( !maxs[i] || !maxs[j] )
// // //                     continue;
// // //                 aci = GetScaleFactor() * 1.5 * ( icorrect - info[i] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
// // //                 acj = GetScaleFactor() * 1.5 * ( icorrect - info[j] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
// // //                 aci = GetScaleFactor() * ( scale - ( scale * exp( 0.2 * ( maxs[i] / GetScaleFactor()))));
// // //                 acj = GetScaleFactor() * ( scale - ( scale * exp( 0.2 * ( maxs[j] / GetScaleFactor()))));
// // //                 ro += ( maxs[i] + aci )*( maxs[j] + acj );
// // //                 no_pairs++;
// // //             }
// //             }
// //             if( priwin & 1 )
// //                 pribeg++;
// //         }
// // 
// //     if( no_pairs )
// //         return ro / ( double )no_pairs;
// //     return ro;
// // }

// // // -------------------------------------------------------------------------
// // // AdjustCosts: adjusts gap costs according to the (max) scores obtained for
// // //     query or subject for each position of profile scoring matrix
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::AdjustCosts(
// //     double* maxscores,
// //     size_t  szvalues,
// //     const   LogOddsMatrix& logo,
// //     int     window,
// //     bool*   masks,
// //     size_t  szmasks,
// //     double  relent )
// // {
// //     if( szvalues == 0 )
// //         return;
// // 
// //     const double*   infovect = logo.GetInformation();
// //     const double*   corrects = GetPosACcorrection();
// // 
// //     SetThickness( logo.GetEffNoSequences());
// // 
// //     if( maxscores == NULL || infovect == NULL )
// //         throw myruntime_error( mystring( "PMTransModel: Memory access error." ));
// // 
// //     if( szvalues < window )
// //         window = szvalues;
// // 
// //     if( window <= 0 )
// //         throw myruntime_error( mystring( "PMTransModel: Non-positive autocorrelation window size." ));
// // 
// //     int     winpar = ( window  & 1 ) ^ 1;
// //     int     midwin = ( window >> 1 ) - winpar;
// // 
// //     double  extcost = GetScaledExtendCost();
// //     int     gapsize = GetColumns();
// // 
// //     int     m;
// //     double  value = 1.0, preval;
// //     double  avgmaxsc = 0.0; //average of max scores
// //     double  stddev = 0.0;   //standard deviation
// //     double  ro = 0.0;       //computed sum of autocorrelations
// //     double  sqrtro, roa;
// //     double  rolim = 100.0 * GetScaleFactor();   //1/0.1^2
// //     bool    winmasked;
// // 
// //     if( gapsize != szvalues || szvalues != szmasks ||
// //         gapsize != logo.GetColumns())
// //         throw myruntime_error( mystring( "PMTransModel: Unable to compute gap costs." ));
// // 
// // 
// //     //window must be strictly positive !!
// //     for( m = 0; m < gapsize - window + 1; m++ ) {
// //         winmasked = false;
// //         for( int j = m; j < m + window; j ++ )
// //             if( masks[j] == true ) {
// //                 winmasked = true;
// //                 break;
// //             }
// // 
// //         ro = Autocorrelation( maxscores + m, window, infovect + m, corrects + m );
// //         //set positional gap costs for subject
// //         if( 0.0 < ro ) {
// //             preval = value;
// // //             value = -( int )rint( sqrt( ro ));
// // //             sqrtro = sqrt( ro );
// // //             roa = 0.0;
// // //             roa = sqrtro / GetScaleFactor();
// // //             roa = 1.0 / ( roa * 0.3 );
// // //             roa *= GetScaleFactor();
// // //             value = - TIMES2( sqrtro ) - roa;
// //             value = -sqrt( ro ); //original
// // //             if( extcost < value )
// // //                 value = extcost;
// //             //if simple derivative is 0
// // //             if( abs( value - preval ) <= 1 * GetScaleFactor() && 2 * extcost < value )
// // //                 value = extcost * 2;
// //             SetPosOpenAt( value, m + midwin );
// //         } else
// // //             SetPosOpenAt( preval = rolim, m + midwin );
// //             SetPosOpenAt( preval = 0.0/*extcost*/, m + midwin );//original
// //     }
// // 
// //     if( midwin < gapsize )
// //         for( m = 0; m < midwin && m < gapsize; m++ )
// //             SetPosOpenAt( GetPosOpenAt( midwin ), m );
// // 
// //     if( window - midwin < gapsize )
// //         for( m = gapsize - window + midwin + 1; m < gapsize; m++ )
// //             SetPosOpenAt( GetPosOpenAt( gapsize - window + midwin ), m );
// // 
// //     ComputeCosts();
// //     //reset context for proceeding calls of gap cost adjustment
// //     ResetContext();
// // 
// // //{{TEST***
// // // fprintf( stderr, "\nGap costs:\n" );
// // // for( m = 0; m < gapsize; m++ )
// // //     fprintf( stderr, "%4d  %7.2f  %7.2f  %4.2fi %4.2fd\n", m + 1,
// // //         GetOpenAt( m )/GetScaleFactor(), maxscores[m]/GetScaleFactor(), GetProbabilityAt( m ), GetDeleteOpenProbAt( m ));
// // //     fprintf( stderr, "%7.2f  %7.2f  %7.2f  %7.2f  %5.2f -- %d %5.2f\n",
// // //         maxscores[m], GetPosOpenAt( m ), GetOpenAt( m ), GetExtendAt( m ), GetWeightsAt( m ), masks[m], infovect[m] );
// // // fprintf( stderr, "\n//\n" );
// // //}}
// // 
// // }

// // // -------------------------------------------------------------------------
// // // ComputeCosts: computes position-specific gap costs given gap weights
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::ComputeCosts()
// // {
// //     double  minopen = GetScaledOpenCost();      //this is minimum score value a gap penalty can acquire
// //     double  maxopen = 0; //GetScaledExtendCost();//maximum score value a gap penalty can acquire
// //     double  minextn = GetScaledExtendCost();    //minimum score value a gap extension penalty can acquire
// //     double  maxextn = 0;                //maximum score value a gap extension penalty can acquire
// //     double  midextn = ( double )( maxextn - minextn ) / 2.0;
// //     double  midopen = ( double )( maxopen - minopen ) / 2.0;
// //     double  cost = 0;
// // 
// // #ifdef __DEBUG__
// //     if( posvec == NULL || posext == NULL ||
// //         vector == NULL || extvec == NULL || weights == NULL )
// //         throw myruntime_error( mystring( "PMTransModel: Memory access error." ));
// // #endif
// // 
// //     if( maxopen < minopen )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Invalid gap opening costs." ));
// // 
// //     for( int n = 0; n < GetColumns(); n++ ) {
// //         cost = minopen = GetPosOpenAt( n );
// //         if( !GetFixed()) {
// //                 //do NOT uncomment: gap costs are computed by using insertion and deletion probabilities
// // //             cost =  minopen + ( maxopen - minopen ) * GetWeightsAt( n );
// // //             cost = -midopen + ( maxopen - minopen ) * SigmoidResponse( GetWeightsAt( n ));
// //         }
// //         if( 0.0 < cost )
// //             cost = 0.0;
// //         SetOpenAt( cost, n );
// // 
// //         cost = minextn = GetPosExtendAt( n );
// //         if( !GetFixed()) {
// //                 //do NOT uncomment: gap costs are computed by using insertion and deletion probabilities
// // //             cost =  minextn + ( maxextn - minextn ) * GetWeightsAt( n );
// // //             cost = -midextn + ( maxextn - minextn ) * SigmoidResponse( GetWeightsAt( n ));
// //         }
// //         if( 0.0 < cost )
// //             cost = 0.0;
// // //         if( cost < GetOpenAt( n ))
// // //             cost = GetOpenAt( n );
// //         SetExtendAt( cost, n );
// //     }
// // }

// -------------------------------------------------------------------------
// Print: print transition probabilities to file
//
void PMTransModel::Print( const char* filename )
{
    myruntime_error mre;
    FILE*   fp = stderr;
    int     t, p;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "PMTransModel::Print: Failed to open file for writing." );

    try {
        fprintf( fp, "%9c", 32 );
        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, " %5s", gTPTRANS_NAMES[t]);
        fprintf( fp, "%s", NL );
        for( p = -1; p < GetSize(); p++ ) 
        {
            if( 0 <= p )
                fprintf( fp, "%s%5d     ", NL, p+1 );
            else
                fprintf( fp, "%10c", 32 );
            for( t = 0; t < P_NSTATES; t++ )
                fprintf( fp, "%6.3f", GetTransProbsAt/*GetOrgTrProbsAt*/( t, p ));
        }
        fprintf( fp, "%s", NL );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( fp != stderr )
        fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// Validate: verify transition probabilities
//
void PMTransModel::Validate() const
{
    int st, states;
    int pos;
    float sum;
    const int   maxtr = gTPTRANS_NTPS;
    const float accuracy = 1.0e-6f;

    for( pos = -1; pos < GetSize(); pos++ ) {
        //verify probability conservation
        for( states = 0; states < P_NSTATES; states += maxtr ) {
            sum = 0.0f;
            for( st = states; st < states + maxtr && st < P_NSTATES; st++ )
                sum += GetOrgTrProbsAt( st, pos );
            if( sum < 1.0f - accuracy || sum > 1.0f + accuracy )
                throw MYRUNTIME_ERROR( "PMTransModel::Validate: Invalid probabilities." );
        }
    }
}

// // // -------------------------------------------------------------------------
// // // Serialize: write the class data to binary file for reading them later
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::Serialize( Serializer& serializer ) const
// // {
// //     serializer.Write(( char* )&length, sizeof( length ), 1 );
// // 
// //     if( length <= 0 )
// //         return;
// // 
// //     int tmpopencost = ( int )rint( GetOpenCost());
// //     int tmpextendcost = ( int )rint( GetExtendCost());
// // 
// //     serializer.Write(( char* )orgtrobs_, sizeof( double ) * P_NSTATES, length + 1 );
// // //     serializer.Write(( char* )vector, sizeof( int ), length );
// //     serializer.Write(( char* )weights, sizeof( double ), length );
// //     serializer.Write(( char* )deletes[DBeg], sizeof( double ), length );
// //     serializer.Write(( char* )deletes[DEnd], sizeof( double ), length );
// //     serializer.Write(( char* )deleteInt, sizeof( int ), length );
// //     serializer.Write(( char* )aacids, sizeof( char ), length );
// // 
// // //     serializer.Write(( char* )&openCost, sizeof( openCost ), 1 );
// // //     serializer.Write(( char* )&extendCost, sizeof( extendCost ), 1 );
// // 
// //     serializer.Write(( char* )&tmpopencost, sizeof( tmpopencost ), 1 );
// //     serializer.Write(( char* )&tmpextendcost, sizeof( tmpextendcost ), 1 );
// // }

// // // -------------------------------------------------------------------------
// // // Deserialize: read data from binary file into the class members
// // // -------------------------------------------------------------------------
// // 
// // void PMTransModel::Deserialize( Serializer& serializer )
// // {
// //     serializer.Read(( char* )&length, sizeof( length ), 1 );
// // 
// //     if( MAXCOLUMNS < length )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Number of positions read from file is larger than maximum allowed." ));
// // 
// //     if( length <= 0 )
// //         throw myruntime_error(
// //             mystring( "PMTransModel: Invalid number of positions read from file." ));
// // 
// //     int tmpopencost = SCORE_MIN;
// //     int tmpextendcost = SCORE_MIN;
// // 
// //     Reserve( length );
// // 
// //     serializer.Read(( char* )orgtrobs_, sizeof( double ) * P_NSTATES, length + 1 );
// // //     serializer.Read(( char* )vector, sizeof( int ), length );
// //     serializer.Read(( char* )weights, sizeof( double ), length );
// //     serializer.Read(( char* )deletes[DBeg], sizeof( double ), length );
// //     serializer.Read(( char* )deletes[DEnd], sizeof( double ), length );
// //     serializer.Read(( char* )deleteInt, sizeof( int ), length );
// //     serializer.Read(( char* )aacids, sizeof( char ), length );
// // 
// // //     serializer.Read(( char* )&openCost, sizeof( openCost ), 1 );
// // //     serializer.Read(( char* )&extendCost, sizeof( extendCost ), 1 );
// // 
// //     serializer.Read(( char* )&tmpopencost, sizeof( tmpopencost ), 1 );
// //     serializer.Read(( char* )&tmpextendcost, sizeof( tmpextendcost ), 1 );
// // 
// //     SetOpenCost( tmpopencost );
// //     SetExtendCost( tmpextendcost );
// //     //
// //     CheckIntegrity();
// //     Initialize();
// // }

}//namespace pmodel
