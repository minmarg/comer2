/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "liblib/mygetopt.h"
#include "liblib/msg.h"
#include "incconf/localconfig.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/JobDispatcher.h"
#include "libHDP/HDPbase.h"
#include "comer.h"

// =========================================================================

int main( int argc, char *argv[] )
{
    int             c;
    char*           p;
//     char            strbuf[BUF_MAX];
    //names of input file and database
    mystring        myoptarg;
    mystring        input;
    mystring        database;
    mystring        output;
    mystring        optfile;
    mystring        verbose;
    mystring        devmaxN;
    mystring        devmaxmem;
    mystring        devcachep;
    mystring        expct_prolen;
    mystring        devpass2memp;
    bool            nofilemap = 0;//do not use file mapping
    bool            cpufreemem = 0;//free unused memory
    bool            listdevices = 0;//whether to list available devices
    int             verblev = 0;//suppress warnings

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],instructs,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"i", my_required_argument, 'i'},
        {"d", my_required_argument, 'd'},
        {"o", my_required_argument, 'o'},
        {"p", my_required_argument, 'p'},
        {"v", my_optional_argument, 'v'},
        {"h", my_no_argument,       'h'},
        {"dev-N", my_required_argument, 'N'},
        {"dev-mem", my_required_argument, 'M'},
        {"dev-cachep", my_required_argument, 'C'},
        {"dev-expected-length", my_required_argument, 'L'},
        {"dev-pass2memp", my_required_argument, '2'},
        {"io-nofilemap", my_no_argument, 'F'},
        {"cpu-freemem",my_no_argument, 'R'},
        {"dev-list", my_no_argument, 'D'},
        { NULL, my_n_targflags, 0 }
    };

    try {
        try {
            MyGetopt mygetopt( long_options, (const char**)argv, argc );
            while(( c = mygetopt.GetNextOption( &myoptarg )) >= 0 ) {
                switch( c ) {
                    case ':':   fprintf( stdout, "Argument missing. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case '?':   fprintf( stdout, "Unrecognized option. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case '!':   fprintf( stdout, "Ill-formed option. Please try option -h for help.%s", NL );
                                return EXIT_FAILURE;
                    case 'h':   fprintf( stdout, "%s", usage(argv[0],instructs,version,verdate).c_str());
                                return EXIT_SUCCESS;
                    case 'i':   input = myoptarg; break;
                    case 'd':   database = myoptarg; break;
                    case 'o':   output = myoptarg; break;
                    case 'p':   optfile = myoptarg; break;
                    case 'v':   verblev = 1; verbose = myoptarg; break;
                    //
                    case 'N':   devmaxN = myoptarg; break;
                    case 'M':   devmaxmem = myoptarg; break;
                    case 'C':   devcachep = myoptarg; break;
                    case 'L':   expct_prolen = myoptarg; break;
                    case '2':   devpass2memp = myoptarg; break;
                    case 'F':   nofilemap = 1; break;
                    case 'R':   cpufreemem = 1; break;
                    case 'D':   listdevices = 1; break;
                    default:    break;
                }
            }
        } catch( myexception const& ex ) {
            error( dynamic_cast<myruntime_error const&>(ex).pretty_format().c_str());
            return EXIT_FAILURE;
        }
    } catch ( ... ) {
        error("Unknown exception caught.");
        return EXIT_FAILURE;
    }


    if( !verbose.empty()) {
        verblev = strtol( verbose.c_str(), &p, 10 );
        if( errno || *p || verblev < 0 ) {
            error( "Invalid verbose mode argument." );
            return EXIT_FAILURE;
        }
    }

    SetVerboseMode( verblev );


    //{{list available devices and exit
    if( listdevices ) {
        TRY
            DEVPROPs.PrintDevices( stdout );
        CATCH_ERROR_RETURN(;);
        return EXIT_SUCCESS;
    }
    //}}



    if( input.empty() && !database.empty()) {
        error( "Profile is not specified." );
        return EXIT_FAILURE;
    }
    if( !input.empty() && database.empty()) {
        error( "Database is not specified." );
        return EXIT_FAILURE;
    }
    if( input.empty() && database.empty()) {
        fprintf( stdout, "%s", usage(argv[0],instructs,version,verdate).c_str());
        return EXIT_FAILURE;
    }

    if( output.empty()) {
        error( "Output directory is not specified." );
        return EXIT_FAILURE;
    }
    if( directory_exists(output.c_str()) == false && 
        mymkdir(output.c_str()) < 0 )
    {
        error(("Failed to create directory: " + output).c_str());
        return EXIT_FAILURE;
    }


    print_dtime(1);


    //{{ COMMAND-LINE OPTIONS
    TRY
        if( !devmaxN.empty()) {
            CLOPTASSIGN( DEV_N, devmaxN );
        }
        if( !devmaxmem.empty()) {
            c = strtol( devmaxmem.c_str(), &p, 10 );
            if( errno || *p ) {
                error( "Invalid argument of option dev-mem." );
                return EXIT_FAILURE;
            }
            CLOPTASSIGN( DEV_MEM, c );
        }
        if( !devcachep.empty()) {
            c = strtol( devcachep.c_str(), &p, 10 );
            if( errno || *p ) {
                error( "Invalid argument of option dev-cachep." );
                return EXIT_FAILURE;
            }
            CLOPTASSIGN( DEV_CACHEP, c );
        }
        if( !expct_prolen.empty()) {
            c = strtol( expct_prolen.c_str(), &p, 10 );
            if( errno || *p ) {
                error( "Invalid argument of option dev-expected-length." );
                return EXIT_FAILURE;
            }
            CLOPTASSIGN( DEV_EXPCT_DBPROLEN, c );
        }
        if( !devpass2memp.empty()) {
            c = strtol( devpass2memp.c_str(), &p, 10 );
            if( errno || *p ) {
                error( "Invalid argument of option dev-pass2memp." );
                return EXIT_FAILURE;
            }
            CLOPTASSIGN( DEV_PASS2MEMP, c );
        }

        CLOPTASSIGN( IO_NOFILEMAP, nofilemap );
        CLOPTASSIGN( CPU_FREEMEM, cpufreemem );

        DEVPROPs.RegisterDevices();

    CATCH_ERROR_RETURN(;);
    //}}



    mystring        insparamfile = GetFullParamFilename();
    mystring        altinsparamfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetParamFilename();

    if( !file_exists( insparamfile.c_str()))
        insparamfile = altinsparamfile;
    //// hdp1
    mystring        inhdp1file = GetFullHDP1Filename();
    mystring        altinhdp1file =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDP1Filename();

    if( !file_exists( inhdp1file.c_str()))
        inhdp1file = altinhdp1file;

    mystring        inhdpscofile = GetFullHDPscoFilename();
    mystring        altinhdpscofile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPscoFilename();

    if( !file_exists( inhdpscofile.c_str()))
        inhdpscofile = altinhdpscofile;
    //// hdpctx
//     mystring        inhdpctxfile = GetFullHDPCtxFilename();
//     mystring        altinhdpctxfile =
//                         mystring( my_dirname( argv[0] )) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetHDPCtxFilename();
// 
//     if( !file_exists( inhdpctxfile.c_str()))
//         inhdpctxfile = altinhdpctxfile;

//     mystring        inhdpctsfile = GetFullHDPctsFilename();
//     mystring        altinhdpctsfile =
//                         mystring( my_dirname( argv[0] )) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetHDPctsFilename();
// 
//     if( !file_exists( inhdpctsfile.c_str()))
//         inhdpctsfile = altinhdpctsfile;
    //// cvs2s
    mystring        incvs2sfile = GetFullCVS2SFilename();
    mystring        altincvs2sfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetCVS2SFilename();

    if( !file_exists( incvs2sfile.c_str()))
        incvs2sfile = altincvs2sfile;


    //// ssss
    mystring        inssssfile = GetFullSSSSFilename();
    mystring        altinssssfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetSSSSFilename();

    if( !file_exists( inssssfile.c_str()))
        inssssfile = altinssssfile;


//     mystring        inhdpssssfile = GetFullHDPSSSSFilename();
//     mystring        altinhdpssssfile =
//                         mystring( my_dirname( argv[0] )) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetHDPSSSSFilename();

//     if( !file_exists( inhdpssssfile.c_str()))
//         inhdpssssfile = altinhdpssssfile;


//     mystring        inihdpssssfile = GetFulliHDPSSSSFilename();
//     mystring        altinihdpssssfile =
//                         mystring( my_dirname( argv[0] )) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetiHDPSSSSFilename();

//     if( !file_exists( inihdpssssfile.c_str()))
//         inihdpssssfile = altinihdpssssfile;


    //// OPTIONS
    mystring        insoptfile = GetFullOptionsFilename();
    mystring        altinsoptfile =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetOptionsFilename();

    if( !file_exists( insoptfile.c_str()))
        insoptfile = altinsoptfile;

    if( optfile.empty())
        optfile = insoptfile;


    TRY
        MOptions::Read( optfile.c_str());
    CATCH_ERROR_RETURN(;);

    float           valEVAL = MOptions::GetEVAL();
    int             valNOHITS = MOptions::GetNOHITS();
    int             valNOALNS = MOptions::GetNOALNS();
    int             valSHOW = MOptions::GetSHOW();

//     float           valIDENTITY = MOptions::GetIDENTITY();
//     int             valDELSTATE = MOptions::GetDELSTATE();

//     float           valPCFWEIGHT = MOptions::GetPCFWEIGHT();
//     float           valMINALNFRN = MOptions::GetMINALNFRN();
//     int             valMINALNPOS = MOptions::GetMINALNPOS();
    mystring        valSUBMAT = MOptions::GetSUBMAT();

    mystring        valTFRMIX = MOptions::GetTFRMIX();
    int             inttfrmix = tfrmixNo;
    float           valMIXWGT = MOptions::GetMIXWGT();
    mystring        valSCOADJ = MOptions::GetSCOADJ();
    int             intscoadj = scoadjNo;
    int             valSUPCLT = MOptions::GetSUPCLT();
    float           valADJWGT = MOptions::GetADJWGT();
//     float           valcADJWGT = MOptions::GetcADJWGT();

    float           valCVSWGT = MOptions::GetCVSWGT();

    float           valSSSWGT = MOptions::GetSSSWGT();
//     float           valSSSHDP = MOptions::GetSSSHDP();

    int             valSSEMODEL = MOptions::GetSSEMODEL();

//     int             valHCFILTER = MOptions::GetHCFILTER();
//     int             valHCWINDOW = MOptions::GetHCWINDOW();
//     float           valHCLOWENT = MOptions::GetHCLOWENT();
//     float           valHCHIGHENT = MOptions::GetHCHIGHENT();

//     int             valINVLCFILTER = MOptions::GetINVLCFILTER();
//     int             valLCFILTEREACH = MOptions::GetLCFILTEREACH();
//     int             valLCWINDOW = MOptions::GetLCWINDOW();
//     float           valLCLOWENT = MOptions::GetLCLOWENT();
//     float           valLCHIGHENT = MOptions::GetLCHIGHENT();
//     float           valDISTANCE = MOptions::GetDISTANCE();

//     mystring        valSCHEME = MOptions::GetSCHEME();
//     float           valMINPP = MOptions::GetMINPP();

    int             valHSPLEN = MOptions::GetHSPLEN();
    int             valHSPMINSCORE = MOptions::GetHSPMINSCORE();
    int             valHSPMAXDIST = MOptions::GetHSPMAXDIST();
    int             valNOHSPS = MOptions::GetNOHSPS();

//     AbstractScoreMatrix::TType      method      = DEFAULT_SCORING_SCHEME;
//     AbstractScoreMatrix::TBehaviour behaviour   = DEFAULT_STATISTICAL_BEHAVIOUR;
//     AbstractScoreMatrix::TScaling   precision   = DEFAULT_PRECISION;
//     ProfileAlignment::TAlgorithm    alnalgo     = DEFAULT_ALNALGORITHM;


    TRY
        //{{ -- sub. matrix --
        SetSTABLE( valSUBMAT );
        //}}
        //{{ SS state scoring parameters
        if( valSSSWGT ) {
            SSSSCORES.ReadScores( inssssfile.c_str());
            SSSSCORES.SetSSSWeight( valSSSWGT );
//             if( valSSSHDP ) {
// #if defined( USEiHDPSSSSCORES )
//                 iHDPSSSSCORES.ReadScores( inihdpssssfile.c_str());
// #else
//                 HDPSSSSCORES.ReadScores( inhdpssssfile.c_str());
// #endif
//                 iHDPSSSSCORES.SetiHDPSSSWeight( valSSSHDP );
//                 HDPSSSSCORES.SetHDPSSSWeight( valSSSHDP );
//             }
        }
        //}}
        //{{ context vector scoring parameters
        if( valCVSWGT ) {
            CVS2SCORES.ReadScores( incvs2sfile.c_str());
            CVS2SCORES.SetCVSWeight( valCVSWGT );
        }
        //}}
        //{{ -- mix. parameters --
        mystring lcval1 = valTFRMIX;
        mystring lcval2 = valSCOADJ;
        lcval1.lower();
        lcval2.lower();
        if( lcval1 == "hdpctx" || 
            lcval2 == "hdpctx" || lcval2 == "hdpsco") {
            SetHDPBASE( inhdp1file.c_str());
//             if( GetHDPctsUsed())
//                 SetHDPctBASE( inhdpctxfile.c_str());
            if( lcval2 == "hdpsco") {
                HDPSCORES.ReadScores( inhdpscofile.c_str());
                HDPBASE.SetScores( &HDPSCORES );
//                 if( GetHDPctsUsed()) {
//                     HDPctSCORES.ReadScores( inhdpctsfile.c_str());
//                     HDPctBASE.SetScores( &HDPctSCORES );
//                 }
                intscoadj = scoadjHDPsco;
            }
            else if( lcval2 == "hdpctx")
                intscoadj = scoadjHDPCtx;
            else if( lcval1 == "hdpctx")
                inttfrmix = tfrmixHDPCtx;
            HDPBASE.SetMixWeight( valMIXWGT );
            HDPBASE.SetNoSupClusters( valSUPCLT );
            HDPBASE.SetAdjWeight( valADJWGT );
            //
//             if( GetHDPctsUsed()) {
//                 HDPctBASE.SetAdjWeight( valcADJWGT );
//                 HDPctBASE.SetNoSupClusters( valSUPCLT );
//             }
        }
        //}}

    CATCH_ERROR_RETURN(;);


    int     ret = EXIT_SUCCESS;
    JobDispatcher* jdisp = NULL;

    TRY
        jdisp = new JobDispatcher(
                    insparamfile.c_str(),
                    input.c_str(),
                    database.c_str(),
                    output.c_str(),
                    valEVAL,
                    valNOHITS,
                    valNOALNS,
                    valSHOW
        );
        if( !jdisp )
            throw MYRUNTIME_ERROR("Not enough memory.");

        jdisp->SetTarFrMix( inttfrmix );
        jdisp->SetScoAdjment( intscoadj );
        jdisp->SetHDPbase( &HDPBASE );
//         if( GetHDPctsUsed())
//             jdisp->SetHDPctbase( &HDPctBASE );

        jdisp->SetSSEModel( valSSEMODEL );
//         jdisp->SetScoringType( valSCHEME );

//         jdisp->SetExtentMinWindow( valMINALNPOS );
//         jdisp->SetExtentMinSeqPercentage( valMINALNFRN );
//         jdisp->SetPseudoCountWeight( valPCFWEIGHT );

        jdisp->SetHSPLength( valHSPLEN );
        jdisp->SetHSPScore( valHSPMINSCORE );
        jdisp->SetHSPDistance( valHSPMAXDIST );
        jdisp->SetHSPNoHSPs( valNOHSPS );

//         if( valHCFILTER )
//             jdisp->SetHCParameters(
//                 valHCWINDOW,
//                 valHCLOWENT,
//                 valHCHIGHENT
//             );
//         if( valINVLCFILTER )
//             jdisp->SetSegParameters(
//                 valLCWINDOW,
//                 valLCLOWENT,
//                 valLCHIGHENT,
//                 valDISTANCE
//             );
//         if( valLCFILTEREACH )
//             jdisp->SetSeqSegParameters(
//                 valLCWINDOW,
//                 valLCLOWENT,
//                 valLCHIGHENT
//             );

        jdisp->Run();

    CATCH_ERROR_RETURN(if(jdisp) delete jdisp);

    if( ret == EXIT_SUCCESS )
        checkforwarnings();

    return ret;
}



// -------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <vector>
#include <chrono>
#include "liblib/mysort.h"
inline int mycom( const void* vec, size_t n1, size_t n2)
{
    std::vector<float>* v = (std::vector<float>*)vec;
//     return (*v)[n1]<(*v)[n2]? -1: 1;
    return (*v)[n1]<(*v)[n2]? -1: ((*v)[n1]==(*v)[n2]? 0: 1);
}
inline int myswp( void* vec, size_t n1, size_t n2)
{
    std::vector<float>* v = (std::vector<float>*)vec;
    float tmp = (*v)[n1]; (*v)[n1]=(*v)[n2]; (*v)[n2]=tmp;
    return 0;
}
inline int mycom3( const void* vec, size_t n1, size_t n2)
{
    float* v = (float*)vec;
//     return v[n1]<v[n2]? -1: 1;
    return v[n1]<v[n2]? -1: (v[n1]==v[n2]? 0: 1);
}
inline int myswp3( void* vec, size_t n1, size_t n2)
{
    float* v = (float*)vec;
    float tmp = v[n1]; v[n1]=v[n2]; v[n2]=tmp;
    return 0;
}
void mytestsort()
{
    std::vector<float> myvec(20000000);
    std::generate(myvec.begin(), myvec.end(), []{return (float)rand()/(float)INT_MAX;});

    std::vector<float> vec1(myvec);

    auto start = std::chrono::high_resolution_clock::now();
    std::sort(vec1.begin(), vec1.end());
    auto finish = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

    std::cout << "std::sort (" << vec1.size() << ") took " << duration << "ms" << std::endl;
    for(int i=0;i<10;i++)
        std::cout << "  " << vec1[i] << std::endl;
    std::cout << std::endl;

    std::vector<float> vec2(myvec);

    start = std::chrono::high_resolution_clock::now();
    HeapSort( &vec2, vec2.size(), mycom, myswp );
    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

    std::cout << "HeapSort took (" << vec2.size() << ") took " << duration << "ms" << std::endl;
    for(int i=0;i<10;i++)
        std::cout << "  " << vec2[i] << std::endl;
    std::cout << std::endl;

    float* vec3 = (float*)malloc(myvec.size()*sizeof(float));
    for(size_t i=0;i<myvec.size();i++)
        vec3[i] = myvec[i];

    start = std::chrono::high_resolution_clock::now();
    HeapSort( vec3, myvec.size(), mycom3, myswp3 );
    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

    std::cout << "HeapSort took (" << myvec.size() << ") took " << duration << "ms" << std::endl;
    for(int i=0;i<10;i++)
        std::cout << "  " << vec3[i] << std::endl;
    std::cout << std::endl;
    free(vec3);
}
