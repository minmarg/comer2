/***************************************************************************
 *   Copyright (C) 2013-2020 Mindaugas Margelevicius                       *
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
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcaln/MSA.h"
#include "makecov.h"



int main( int argc, char *argv[] )
{
    int             c;
    //char            strbuf[BUF_MAX];
    //string values of parameters
    mystring        myoptarg;
    mystring        input;
    mystring        output;
    mystring        optfile;
    bool            suppress = true;//suppress warnings
    bool extentsused = false;//use extents for checking whether positions coontribute to calculation
    bool mixwithtrgfrqs = false;//mix observations with target frequencies
    bool xcorrelation = false;//calculate cross-correlation matrices instead of xcov
    bool scaleseqwgts = false;//scale sequence weights by the number of sequences in the extent
    bool mutinf = false;//calculate mutual information additionally

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"i", my_required_argument, 'i'},
        {"o", my_required_argument, 'o'},
        {"p", my_required_argument, 'p'},
        {"v", my_no_argument,       'v'},
        {"h", my_no_argument,       'h'},
        {"extents", my_no_argument,  'E'},
        {"mix", my_no_argument,      'M'},
        {"xcorr", my_no_argument,    'C'},
        {"scale", my_no_argument,    'S'},
        {"mi", my_no_argument,       'I'},
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
                    case 'h':   fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
                                return EXIT_SUCCESS;
                    case 'i':   input = myoptarg; break;
                    case 'o':   output = myoptarg; break;
                    case 'p':   optfile = myoptarg; break;
                    case 'v':   suppress = false; break;
                    //
                    case 'E':   extentsused = 1; break;
                    case 'M':   mixwithtrgfrqs = 1; break;
                    case 'C':   xcorrelation = 1; break;
                    case 'S':   scaleseqwgts = 1; break;
                    case 'I':   mutinf = 1; break;
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

    SetVerboseMode( !suppress );

    if( input.empty()) {
        error( "Input multiple alignment is not specified." );
        return EXIT_FAILURE;
    }
    if( output.empty()) {
        error( "Output file is not specified." );
        return EXIT_FAILURE;
    }

    print_dtime(1);

    CLOPTASSIGN( XCOV_USEEXTENTS, extentsused );
    CLOPTASSIGN( XCOV_MIXTRGFRQS, mixwithtrgfrqs );
    CLOPTASSIGN( XCOV_CORR, xcorrelation );
    CLOPTASSIGN( XCOV_SCALEWGTS, scaleseqwgts );
    CLOPTASSIGN( XCOV_MI, mutinf );

    mystring        insparamfile = GetFullParamFilename();
    mystring        altinsparamfile =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetParamFilename();

    if( !file_exists( insparamfile.c_str()))
        insparamfile = altinsparamfile;
    //// hdp1
    mystring        inhdp1file = GetFullHDP1Filename();
    mystring        altinhdp1file =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDP1Filename();

    if( !file_exists( inhdp1file.c_str()))
        inhdp1file = altinhdp1file;
    //// hdpctx
//     mystring        inhdpctxfile = GetFullHDPCtxFilename();
//     mystring        altinhdpctxfile =
//                         mystring( my_dirname( argv[0])) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetHDPCtxFilename();

//     if( !file_exists( inhdpctxfile.c_str()))
//         inhdpctxfile = altinhdpctxfile;
    ////
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

    float           valIDENTITY = MOptions::GetIDENTITY();
    int             valDELSTATE = MOptions::GetDELSTATE();

    float           valPCFWEIGHT = MOptions::GetPCFWEIGHT();
    float           valMINALNFRN = MOptions::GetMINALNFRN();
    int             valMINALNPOS = MOptions::GetMINALNPOS();
    mystring        valSUBMAT = MOptions::GetSUBMAT();
    mystring        valTFRMIX = MOptions::GetTFRMIX();
    //int             inttfrmix = tfrmixNo;
    float           valMIXWGT = MOptions::GetMIXWGT();
    mystring        valSCOADJ = MOptions::GetSCOADJ();
    int             intscoadj = scoadjNo;
    int             valSUPCLT = MOptions::GetSUPCLT();
    float           valADJWGT = MOptions::GetADJWGT();
    //float           valcADJWGT = MOptions::GetcADJWGT();

    int             valHCFILTER = MOptions::GetHCFILTER();
    int             valHCWINDOW = MOptions::GetHCWINDOW();
    float           valHCLOWENT = MOptions::GetHCLOWENT();
    float           valHCHIGHENT = MOptions::GetHCHIGHENT();

//     int             valINVLCFILTER = MOptions::GetINVLCFILTER();
    int             valLCFILTEREACH = MOptions::GetLCFILTEREACH();
    int             valLCWINDOW = MOptions::GetLCWINDOW();
    float           valLCLOWENT = MOptions::GetLCLOWENT();
    float           valLCHIGHENT = MOptions::GetLCHIGHENT();
//     float           valDISTANCE = MOptions::GetDISTANCE();


    TRY
        //{{ -- sub. matrix --
        SetSTABLE( valSUBMAT );
        //}}
        //{{ -- mix. parameters --
        mystring lcval1 = valTFRMIX;
        mystring lcval2 = valSCOADJ;
        lcval1.lower();
        lcval2.lower();
        if( lcval1 == "hdpctx" || 
            lcval2 == "hdpctx" || lcval2 == "hdpsco") {
            SetHDPBASE( inhdp1file.c_str());
            if( lcval2 == "hdpsco")
                intscoadj = scoadjHDPsco;
            else if( lcval2 == "hdpctx")
                intscoadj = scoadjHDPCtx;
            //else if( lcval1 == "hdpctx")
            //    inttfrmix = tfrmixHDPCtx;
            HDPBASE.SetMixWeight( valMIXWGT );
            HDPBASE.SetNoSupClusters( valSUPCLT );
            HDPBASE.SetAdjWeight( valADJWGT );
        }
        //}}

#if 0//printing parameters
        //{{ -- Print --
        sprintf( strbuf, "IDENTITY = %.2f", valIDENTITY );
        message( strbuf, false );

        sprintf( strbuf, "PCFWEIGHT = %.1f, MINALNFRN = %.2f, MINALNPOS = %d", valPCFWEIGHT, valMINALNFRN, valMINALNPOS );
        message( strbuf, false );

        sprintf( strbuf, "HCFILTER = %s", valHCFILTER? "yes": "no" );
        message( strbuf, false );
        if( valHCFILTER ){
            sprintf( strbuf, "HCWINDOW = %d, HCLOWENT = %.2f, HCHIGHENT = %.2f",
                        valHCWINDOW, valHCLOWENT, valHCHIGHENT );
            message( strbuf, false );
        }

        sprintf( strbuf, "INVLCFILTER = %s, LCFILTEREACH = %s",
                    valINVLCFILTER? "yes": "no", valLCFILTEREACH? "yes": "no" );
        message( strbuf, false );
        if( valINVLCFILTER || valLCFILTEREACH ){
            sprintf( strbuf, "LCWINDOW = %d, LCLOWENT = %.2f, LCHIGHENT = %.2f",
                        valLCWINDOW, valLCLOWENT, valLCHIGHENT );
            if( valINVLCFILTER )
                sprintf( strbuf + strlen( strbuf ), ", DISTANCE = %.2f", valDISTANCE );
            message( strbuf, false );
        }

        message( NULL );
        //}}
#endif//printing parameters

    CATCH_ERROR_RETURN(;);


    int     ret = EXIT_SUCCESS;
    MSA*    msa = NULL;

    TRY
        msa = new MSA();
        if( !msa )
            throw MYRUNTIME_ERROR("Not enough memory.");

        //target frequencies will be mixed, if required, just before searching
//         msa->SetTarFrMix( inttfrmix );
        msa->SetScoAdjment( intscoadj );
        msa->SetHDPbase( &HDPBASE );
//         if( GetHDPctsUsed())
//             msa->SetHDPctBase( &HDPctBASE );

        if( valHCFILTER ) {
            msa->SetSEGWindow( valHCWINDOW );
            msa->SetSEGLowEntropy( valHCLOWENT );
            msa->SetSEGHighEntropy( valHCHIGHENT );
        }
        else
            msa->SetUsingSEGFilter( false );

        if( valLCFILTEREACH ) {
            msa->SetSeqSEGWindow( valLCWINDOW );
            msa->SetSeqSEGLowEntropy( valLCLOWENT );
            msa->SetSeqSEGHighEntropy( valLCHIGHENT );
        }
        else
            msa->SetUsingSeqSEGFilter( false );


        msa->SetIdentityLevel( valIDENTITY );
        msa->SetComputeDELETEstates( valDELSTATE );
        msa->SetExtentMinWindow( valMINALNPOS );
        msa->SetExtentMinSeqPercentage( valMINALNFRN );
        msa->SetPseudoCountWeight( valPCFWEIGHT );


        msa->ReadAlignment( input.c_str());
        message( "Calculating cross-covariance matrices...");
        msa->CalculateXCovMatrices( output.c_str());
        message( "Done.");

    CATCH_ERROR_RETURN(if(msa) delete msa);

    return ret;
}
