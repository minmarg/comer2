/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "liblib/mybase.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/Db.h"
#include "libpro/srcpro/VirtualRoDb.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cusco/CuBatchScoreMatrix.cuh"
#include "JobDispatcher.h"

// _________________________________________________________________________
// Class JobDispatcher
//
// Constructor
//
JobDispatcher::JobDispatcher(
    const char* configfile,
    const char* input,
    const char* database,
    const char* output,
    float eval_thld,
    int no_hits,
    int no_alns,
    bool showpars )
:
    configFile_( configfile ),
    input_( input ),
    database_( database ),
    output_( output ),
//     scoreSystem_( NULL ),

    eval_upper_( eval_thld ),
    max_no_hits_( no_hits ),
    max_no_alns_( no_alns ),
    show_pars_( showpars ),

    inpdb_( input ),
    prodb_( database ),
    db_no_seqs_( 0 ),
    db_size_( 0ULL ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
//     HDPctbase_( NULL ),

    ssemodel_( 0 ),

    hsplength_( 0 ),
    hspscore_( 0 ),
    hspdistance_( 0 ),
    hspnohsps_( 0 )
{
}

// Default constructor
//
JobDispatcher::JobDispatcher()
:
    configFile_( NULL ),
    input_( NULL ),
    database_( NULL ),
    output_( NULL ),
//     scoreSystem_( NULL ),

    eval_upper_( 0.0f ),
    max_no_hits_( 0 ),
    max_no_alns_( 0 ),
    show_pars_( false ),

    inpdb_( NULL ),
    prodb_( NULL ),
    db_no_seqs_( 0 ),
    db_size_( 0ULL ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
//     HDPctbase_( NULL ),

    ssemodel_( 0 ),

    hsplength_( 0 ),
    hspscore_( 0 ),
    hspdistance_( 0 ),
    hspnohsps_( 0 )
{
    throw MYRUNTIME_ERROR("JobDispatcher::JobDispatcher: Default initialization is prohibited.");
}

// Destructor
//
JobDispatcher::~JobDispatcher()
{
//     DestroyScoreSystem();
}

// -------------------------------------------------------------------------
// ReadConfiguration: read configurations
//
void JobDispatcher::ReadConfiguration()
{
    MYMSG( "JobDispatcher::ReadConfiguration", 3 );
    GetConfiguration(CTUngapped).SetFilename( GetConfigFile());
    GetConfiguration(CTGapped).SetFilename( GetConfigFile());

    //read parameters
    static const bool autogapcosts = true;
    static const int autogapopencostunused = 4;
    static const int autogapextncostunused = 1;
    GetConfiguration(CTUngapped).ReadUngapped();//read ungapped configuration
    GetConfiguration(CTGapped).SetAutoGapOpenCost( autogapcosts );
    GetConfiguration(CTGapped).SetGapOpenCost( autogapopencostunused );
    GetConfiguration(CTGapped).SetGapExtendCost( autogapextncostunused );
    GetConfiguration(CTGapped).Read();//read reference parameters
}

// -------------------------------------------------------------------------
// Run: launch worker threads for searching and aligning profiles
//
void JobDispatcher::Run()
{
    MYMSG( "JobDispatcher::Run", 3 );
    ReadConfiguration();
    CuBatchScoreMatrix* cu = new CuBatchScoreMatrix( GetConfiguration());
    if( cu == NULL )
        throw MYRUNTIME_ERROR( "JobDispatcher::Run: Not enough memory." );
    cu->ComputeScoreMatrix();
    delete cu;
    cu = NULL;

//     if( scoreSystem )
//         throw myruntime_error( mystring( "JobDispatcher: Score system unexpectedly initialized." ));
// 
//     query_gaps.SetConfiguration( GetConfiguration());
//     dbgaps.SetConfiguration( GetConfiguration());
// 
//     ProcessQuery();
// 
//     try{
//         //try read database header next
//         profile_db.Open();
// 
//     } catch( myexception const& ex )
//     {
//         try{
//             //try profile or multiple alignment in fasta
//             FillMatrices( dbfreq, dbpssm, dbgaps, GetDatabase());
//         } catch( myexception const& ex1 ) {
//             throw myruntime_error(( mystring( ex.what()) + "\n" )+ ex1.what(), ex.eclass());
//         }
//         
// 
//         message( "Searching..." );
// 
//         SetNoSequences( 1 );
//         SetDbSize( dbpssm.GetColumns());
// 
//         AbstractScoreMatrix::ComputeLengthAdjustment(
//                         GetConfiguration( ProcomGapped ),
//                         query_pssm.GetColumns(),
//                         GetDbSize(),
//                         GetNoSequences()
//         );
// 
//         ComputationLogicWithProfiles( dbfreq, dbpssm, dbgaps );
// 
//         PostComputationLogic();
//         message( "Finished." );
//         return;
//     }
// 
//     message( "Searching..." );
// 
//     SetNoSequences( profile_db.GetNoSequences());
//     SetDbSize( profile_db.GetDbSize());
// 
//     AbstractScoreMatrix::ComputeLengthAdjustment(
//                     GetConfiguration( ProcomGapped ),
//                     query_pssm.GetColumns(),
//                     GetDbSize(),
//                     GetNoSequences()
//     );
// 
//     //while not having reached the end of database
//     while( profile_db.Next( dbfreq, dbpssm, dbgaps, GetGapOpenCost(), GetGapExtnCost(), GetFixedCosts())) {
//         ComputationLogicWithProfiles( dbfreq, dbpssm, dbgaps );
//     }
// 
//     PostComputationLogic();
// 
//     profile_db.Close();
//     message( "Finished." );
}

// // -------------------------------------------------------------------------
// // PrintHeader: print a header related to a search
// //
// void JobDispatcher::PrintHeader( FILE* fp )
// {
//     if( fp == NULL )
//         return;
// 
//     const size_t    padding = OUTPUTINDENT;
//     char            blanks[padding+1];
//     size_t          n = 0;
// 
//     for( n = 0; n < padding; blanks[n++] = 32 );
//     blanks[n] = 0;
// 
//     progname_and_version( fp );
// 
//     if( GetFixedCosts()) {
//         if( GetAutoGapCosts())
//             fprintf( fp, " -positional gap costs unaffected by probabilities-\n" );
//         else
//             fprintf( fp, " -positionally invariable gap costs-\n" );
//     }
// 
//     fprintf( fp, "\n" );
//     fprintf( fp, " Query (%d positions):\n", query_pssm.GetColumns());
// 
//     query_pssm.PrintDescriptionFirst( fp );
// 
//     fprintf( fp, "\n\n" );
//     fprintf( fp, " Database:\n" );
//     fprintf( fp, "%s\n", my_basename( profile_db.GetDbName()));
// 
//     
//     fprintf( fp, "%s%d profiles\n%s%llu total positions\n\n\n",
//             blanks, GetNoSequences(),
//             blanks, GetDbSize());
// }

// // -------------------------------------------------------------------------
// // PrintHits: print hits found
// //
// void JobDispatcher::PrintHits( FILE* fp )
// {
//     if( fp == NULL )
//         return;
// 
//     size_t      width = OUTPUTINDENT + OUTPUTWIDTH;
//     const char* title = " Profiles found below the e-value threshold:";
//     const char* nohitstitle = " No profiles found below the e-value threshold";
//     size_t      padding = width - strlen( title );
// 
//     if( !GetHitlistSize()) {
//         fprintf( fp, "%s of %g.\n\n\n", nohitstitle, GetEvalueThreshold());
//         return;
//     }
// 
//     fprintf( fp, "%s", title );
//     for( size_t n = 0; n < padding; n++ ) fprintf( fp, " " );
//     fprintf( fp, " %7s %7s\n\n", "Score", "E-value" );
// 
//     for( size_t h = 0; h < GetHitlistSize() && ( int )h < GetNoHitsThreshold(); h++ )
//     {
//         const HitInformation*   hit = GetHitAt( h );
//         //annotation will never be null
//         size_t                  alength = strlen( hit->GetAnnotation());
// 
//         fprintf( fp, "%s", hit->GetAnnotation());
//         for( size_t n = 0; n < width - alength; n++ ) fprintf( fp, " " );
// 
//         if( hit->GetEvalue() < 0.0 )
// //             fprintf( fp, " %7d %7s\n", ( int )hit->GetScore(), "n/a" );
//             fprintf( fp, " %7d %7.1g\n", ( int )hit->GetScore(), hit->GetRefEvalue());
//         else
//             fprintf( fp, " %7d %7.1g\n", ( int )hit->GetScore(), hit->GetEvalue());
//     }
//     fprintf( fp, "\n\n" );
// 
// 
//     for( size_t h = 0; h < GetHitlistSize() && ( int )h < GetNoAlnsThreshold(); h++ )
//     {
//         const HitInformation*   hit = GetHitAt( h );
//         fprintf( fp, "%s\n", hit->GetFullAlignment());
//     }
// }

// // -------------------------------------------------------------------------
// // PrintFooter: print a footer, statistical summary of a search
// //
// void JobDispatcher::PrintFooter( FILE* fp )
// {
//     if( fp == NULL )
//         return;
// 
//     size_t  length = ScoringMatrix::GetDeltaLength();       //length adjustment
//     Uint64  sspace = ScoringMatrix::GetSearchSpace();       //effective search space
// 
//     int     eff_query_length = query_pssm.GetColumns() - length;
//     Int64   eff_db_length = GetDbSize() - GetNoSequences() * length;
// 
//     fprintf( fp, "\n\n" );
//     PrintMethodName( fp );
//     if( UngappedAlignments())
//         fprintf( fp, "No gap costs used: Ungapped alignments\n\n" );
//     else if(0) {
//         if( GetAutoGapCosts()) {
//             fprintf( fp, "Gap open cost, Computed (window size, %d)\n", GetAutocorrWinsize());
//             fprintf( fp, "Initial gap extension cost, %d\n", abs( GetGapExtnCost()));
//         } else {
//             fprintf( fp, "Gap open cost, %d\n", abs( GetGapOpenCost()));
//             fprintf( fp, "Gap extension cost, %d\n", abs( GetGapExtnCost()));
//         }
//         if( !GetFixedCosts())
//             fprintf( fp, "Deletion probability weight, %.2f\n", GetDeletionCoefficient());
//         fprintf( fp, "\n" );
//     }
//     PrintParameterTable( fp );
//     fprintf( fp, "\n" );
//     fprintf( fp, "Length of query, %d\n", query_pssm.GetColumns());
//     fprintf( fp, "Length of database, %llu\n", GetDbSize());
//     fprintf( fp, "Effective length of query, %d\n", eff_query_length );
//     fprintf( fp, "Effective length of database, %lld\n", eff_db_length );
//     fprintf( fp, "Effective search space, %llu\n\n", sspace );
// }
