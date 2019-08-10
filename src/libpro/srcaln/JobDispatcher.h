/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __JobDispatcher_h__
#define __JobDispatcher_h__

#include "liblib/mybase.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/Db.h"
#include "libpro/srcpro/VirtualRoDb.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"

// _________________________________________________________________________
// Class JobDispatcher
//
// Implementation of distributing jobs over CPU and GPU threads
//
class JobDispatcher
{
public:
    JobDispatcher(
            const char* configfile,
            const char* input,
            const char* database,
            const char* output,
            float eval_thld,
            int no_hits,
            int no_alns,
            bool showpars
    );
    ~JobDispatcher();

    void            Run();

    const char*     GetConfigFile() const { return configFile_; }
    const char*     GetInput() const            { return input_; }
    const char*     GetDatabase() const         { return database_; }
    const char*     GetOutput() const           { return output_; }

    double          GetEvalueUpper() const  { return eval_upper_; }
    int             GetMaxNoHits() const  { return max_no_hits_; }
    int             GetMaxNoAlns() const  { return max_no_alns_; }
    bool            GetShowPars() const          { return show_pars_; }

    void            SetTarFrMix( int value ) { tfrmix_ = value; }
    bool            GetTarFrMixHDPCtx() const { return tfrmix_ == tfrmixHDPCtx; }

    void            SetScoAdjment( int value ) { scoadj_ = value; }
    bool            GetScoAdjmentHDPCtx() const { return scoadj_ == scoadjHDPCtx || scoadj_ == scoadjHDPsco; }

    void            SetHDPbase( const HDPbase* value ) { HDPbase_ = value; }
    const HDPbase*  GetHDPbase() const { return HDPbase_; }

//     void            SetHDPctbase( const HDPbase* value ) { HDPctbase_ = value; }
//     const HDPbase*  GetHDPctbase() const { return HDPctbase_; }

    void            SetSSEModel( int modndx ) { ssemodel_ = modndx; }
    int             GetSSEModel() const { return ssemodel_; }


    int             GetHSPLength() const                        { return hsplength_; }
    void            SetHSPLength( int value )                   { hsplength_ = value; }

    int             GetHSPScore() const                         { return hspscore_; }
    void            SetHSPScore( int value )                    { hspscore_ = value; }

    int             GetHSPDistance() const                      { return hspdistance_; }
    void            SetHSPDistance( int value )                 { hspdistance_ = value; }

    int             GetHSPNoHSPs() const                        { return hspnohsps_; }
    void            SetHSPNoHSPs( int value )                   { hspnohsps_ = value; }

    void            PrintMethodName( FILE* fp ) const;//scoring method name
    void            PrintParameterTable( FILE* ) const;

protected:
    explicit JobDispatcher();

    void    ReadConfiguration();

//     void                        CreateScoreSystem( const FrequencyMatrix&, const LogOddsMatrix& );
//     void                        CreateScoreSystem();            //create member score system
//     void                        DestroyScoreSystem();           //destroy score system
//     void                        ComputeScoreSystem();           //compute score system if needed
//     void                        ScaleScoreSystem();             //scale score system
//     bool    ScanForHSPs( double minhspscore, int hsplen, int nohsps, int mindist, 
//                         int* possbjct = NULL, int* posquery = NULL );
//     bool                        PreprocessSubject(              //preprocessing of subject profile
//         const FrequencyMatrix&, const LogOddsMatrix&, GapScheme&, bool firstpass );

//     const AbstractScoreMatrix*  GetScoreSystem() const  { return scoreSystem; }
//     AbstractScoreMatrix*        GetScoreSystem()        { return scoreSystem; }

//     void                        ComputationLogicWithProfiles( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
//     void                        PostComputationLogic();


    void    PrintHeader( FILE* );//header of a search
    void    PrintHits( FILE* );//results of a search
    void    PrintFooter( FILE* );//footer of a search


    const Configuration&    GetConfiguration( TConfigType ) const;
    Configuration&          GetConfiguration( TConfigType );
    Configuration*          GetConfiguration() { return config_; }

//     void                        ProcessQuery(); //read query information from the file
//     void                        FillMatrices(   //fill profile matrices by reading and processing information from file
//                                     FrequencyMatrix&, LogOddsMatrix&, GapScheme&,
//                                     const char* filename );


    size_t      GetNoSequences() const { return db_no_seqs_; }
    void        SetNoSequences( size_t value ) { db_no_seqs_ = value; }

    uint64_mt   GetDbSize() const { return db_size_; }
    void        SetDbSize( uint64_mt value ) { db_size_ = value; }

private:
    const char*             configFile_;//configuration file of statistical parameters
    const char*             input_;//input's filename
    const char*             database_;//profile database
    const char*             output_;//pattern for output file (null=standard output)
//     AbstractScoreMatrix<>*    scoreSystem_;//profile scoring system
    Configuration           config_[NoCTypes];//configurations

    float                   eval_upper_;//e-value upper threshold
    int                     max_no_hits_;//maximum number of hits to show in the result list
    int                     max_no_alns_;//maximum number of alignments to show in the output
    bool                    show_pars_;//whether to show statistical parameters below alignments

    VirtualRoDb             inpdb_;//input batch of profiles, which also can be a database
    Db                      prodb_;//profile database
    size_t                  db_no_seqs_;//number of profiles in the database
    uint64_mt               db_size_;//size of the database in positions

    int                     tfrmix_;//mixing of target frequencies
    int                     scoadj_;////score adjustment
    const HDPbase*          HDPbase_;//HDP base data structure
//     const HDPbase*          HDPctbase_;//HDP ctx base structure

    int                     ssemodel_;//index of a model for the estimation of stat. sign.

    int                     hsplength_;//minimum length for HSPs
    int                     hspscore_;//minimum HSP score
    int                     hspdistance_;//maximum distance between the HSPs
    int                     hspnohsps_;//minimum number of HSPs in a diagonal
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// GetConfiguration: returns reference to parameter configuration object
//
inline
const Configuration& JobDispatcher::GetConfiguration( TConfigType ct ) const
{
#ifdef __DEBUG__
    if( NoCTypes <= ct )
        throw MYRUNTIME_ERROR( "JobDispatcher: Wrong configuration index." );
#endif
    return config_[ct];
}
inline
Configuration& JobDispatcher::GetConfiguration( TConfigType ct )
{
#ifdef __DEBUG__
    if( NoCTypes <= ct )
        throw MYRUNTIME_ERROR( "JobDispatcher: Wrong configuration index." );
#endif
    return config_[ct];
}

// // CreateScoreSystem: creates score system of alternative type
// //
// inline
// void JobDispatcher::CreateScoreSystem(
//         const FrequencyMatrix& sbjctfreq,
//         const LogOddsMatrix& sbjctpssm )
// {
//     switch( GetMethod()) {
//         case AbstractScoreMatrix::ProfSpecLSO:
//                 scoreSystem = new LSOSMatrix(
//                     query_freq,
//                     query_pssm,
//                     sbjctfreq,
//                     sbjctpssm,
//                     GetScoAdjmentHDPCtx()? GetHDPbase(): NULL,
//                     GetScoAdjmentHDPCtx()? GetHDPctbase(): NULL,
//                     GetInformationThreshold(),
//                     GetThicknessNumber(),
//                     GetThicknessPercents(),
//                     GetMaskscalePercents(),
//                     GetConfiguration(),
//                     GetBehaviour(), //ComputeStatistics or StatisticsGiven
//                     GetScaling(),   //FPScaling, AutoScalling, or NoScaling
//                     GetMasking()
//                 );
//         break;
//         default: return;
//     };
// 
//     if( scoreSystem == NULL )
//         throw myruntime_error( mystring( "JobDispatcher: Not enough memory." ));
// 
//     scoreSystem->SetName( sbjctpssm.GetName());
// }

// // -------------------------------------------------------------------------
// // DestroyScoreSystem: destroys score system
// //
// inline
// void JobDispatcher::DestroyScoreSystem()
// {
//     if( scoreSystem )
//         delete scoreSystem;
// }

// // -------------------------------------------------------------------------
// // computeProfileScoringMatrix: calls appropriate method of ScoresMatrix object
// //
// inline
// void JobDispatcher::ComputeScoreSystem()
// {
// #ifdef __DEBUG__
//     if( !scoreSystem )
//         throw myruntime_error( mystring( "JobDispatcher: Unable to compute score matrix." ));
// #endif
//     scoreSystem->ComputeProfileScoringMatrix();
// }

// // -------------------------------------------------------------------------
// // ScanForHSPs: run algorithm of hsps
// //
// inline
// bool JobDispatcher::ScanForHSPs(
//     double minhspscore, int hsplen, int nohsps, int mindist,
//     int* possbjct, int* posquery )
// {
// #ifdef __DEBUG__
//     if( !scoreSystem )
//         throw myruntime_error( mystring( "JobDispatcher: Unable to scan score matrix for HSPs." ));
// #endif
//     if( minhspscore <= 0.0 )
//         return true;
// 
//     switch( scoreSystem->GetType()) {
//         case AbstractScoreMatrix::ProfileSpecific:
//         case AbstractScoreMatrix::ProfSpecLSO:
//         case AbstractScoreMatrix::AdjustedProfileSpecific:
//         case AbstractScoreMatrix::HDPProfileSpecific:
//         case AbstractScoreMatrix::HDP0CtxProfileSpecific:
//             return scoreSystem->ScanForHSPs( minhspscore, hsplen, nohsps, mindist, possbjct, posquery );
// //             break;
// 
//         case AbstractScoreMatrix::Universal:
//             break;
// 
//         default:
//             break;
//     };
//     return true;
// }

// // -------------------------------------------------------------------------
// // ScaleScoreSystem: scale score system
// //
// inline
// void JobDispatcher::ScaleScoreSystem()
// {
// #ifdef __DEBUG__
//     if( !scoreSystem )
//         throw myruntime_error( mystring( "JobDispatcher: Unable to scale score matrix." ));
// #endif
//     switch( scoreSystem->GetType()) {
//         case AbstractScoreMatrix::ProfSpecLSO:
//             scoreSystem->ScaleScoringMatrix();
//             break;
// 
//         default:
//             break;
//     };
// }

// // -------------------------------------------------------------------------
// // PreprocessSubject: prepares scoring system for alignment with subject
// //     profile
// //
// inline
// bool JobDispatcher::PreprocessSubject(
//     const FrequencyMatrix&  freq,
//     const LogOddsMatrix&    pssm,
//     GapScheme&              gaps,
//     bool                    firstpass )
// {
//     if( firstpass ) {
//         if( scoreSystem == NULL ||
//             scoreSystem->GetType() == AbstractScoreMatrix::ProfSpecLSO )
//         {
//             DestroyScoreSystem();
//             CreateScoreSystem( freq, pssm );
//             ComputeScoreSystem();
//             if( !ScanForHSPs( GetHSPScore(), GetHSPLength(), GetHSPNoHSPs(), GetHSPDistance()))
//                 return false;
//             if( scoreSystem->GetSupportOptimFreq())///originally commented
//                 scoreSystem->OptimizeTargetFrequencies();
//             ScaleScoreSystem();
//         }
// #ifdef __DEBUG__
//         if( !scoreSystem )
//             throw myruntime_error( mystring( "JobDispatcher: Unable to preprocess subject profile." ));
// #endif
// 
//         if( scoreSystem->GetType() == AbstractScoreMatrix::Universal )
//             dynamic_cast<UniversalScoreMatrix*>( scoreSystem )->PreserveSubject( freq, pssm );
//     }
// 
//     gaps.SetUsePosACcorrections( GetAutocorrectionPositional());
//     query_gaps.SetUsePosACcorrections( GetAutocorrectionPositional());
// 
//     scoreSystem->PostScalingProc(
//         query_pssm, pssm,
//         query_gaps, gaps,
//         GetAutoGapCosts(),
//         GetAutocorrWinsize(),
//         firstpass
//     );
//     return true;
// }

// // -------------------------------------------------------------------------
// // PrintMethodName: print the name of a profile scoring method
// //
// inline
// void JobDispatcher::PrintMethodName( FILE* fp ) const
// {
//     const size_t    padding = OUTPUTINDENT;
//     char            blanks[padding+1];
//     size_t          n = 0;
// 
//     for( n = 0; n < padding; blanks[n++] = 32 );
//     blanks[n] = 0;
// 
//     if( scoreSystem_ ) {
//         fprintf( fp, "Scoring method: %s%s", scoreSystem_->GetMethodName(), NL );
//         fprintf( fp, "%s", NL );
//         return;
//     }
// }

// // -------------------------------------------------------------------------
// // PrintParameterTable: print the parameter table attributable to the 
// // scoring system used
// //
// inline
// void JobDispatcher::PrintParameterTable( FILE* fp ) const
// {
//     if( !scoreSystem_ )
//         return;
// 
//     scoreSystem_->PrintParameterTable( fp );
// }

#endif//__JobDispatcher_h__
