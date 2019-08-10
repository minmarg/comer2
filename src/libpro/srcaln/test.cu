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
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "JobDispatcher.h"

// -------------------------------------------------------------------------
// Run: launch worker threads for searching and aligning profiles
//
 __global__ void Run()
{
    GetConfiguration(CTUngapped).SetFilename( GetConfigFile());
    GetConfiguration(CTGapped).SetFilename( GetConfigFile());

    int asd;
    __device__(asd);
    Run<<<my,args>>>();
    //read parameters
    static const bool autogapcosts = true;
    static const int autogapopencostunused = 4;
    static const int autogapextncostunused = 1;
    GetConfiguration(CTUngapped).ReadUngapped();//read ungapped configuration
    GetConfiguration(CTGapped).SetAutoGapOpenCost( autogapcosts );
    GetConfiguration(CTGapped).SetGapOpenCost( autogapopencostunused );
    GetConfiguration(CTGapped).SetGapExtendCost( autogapextncostunused );
    GetConfiguration(CTGapped).Read();//read reference parameters

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
// #ifdef  UNIVERSALSCORES
//     message( "Reading frequencies..." );
//     //read the distinct frequency vectors first
//     profile_db.ReadInFrequencies();
// 
//     message( "Computing scores..." );
//     CreateScoreSystem();
// #endif
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
