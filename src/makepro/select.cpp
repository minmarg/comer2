/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
// #include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "liblib/mygetopt.h"
#include "liblib/msg.h"
#include "incconf/localconfig.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/SEGProfile.h"
// #include "libpro/srcpro/PMTransModel.h"
// #include "libpro/srcpro/PMProfileModel.h"
// #include "libHDP/HDPbase.h"
#include "libpro/srcaln/MSA.h"
#include "select.h"



int main( int argc, char *argv[] )
{
    int             c;
    //char            strbuf[BUF_MAX];
    //string values of parameters
    mystring        myoptarg;
    mystring        input;
    mystring        output;
    mystring        optfile;
    bool            noselect = false;//do not apply select
    bool            suppress = true;//suppress warnings

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], selinst, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"i", my_required_argument, 'i'},
        {"o", my_required_argument, 'o'},
        {"p", my_required_argument, 'p'},
        {"t", my_no_argument,       't'},
        {"v", my_no_argument,       'v'},
        {"h", my_no_argument,       'h'},
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
                    case 'h':   fprintf( stdout, "%s", usage(argv[0],selinst,version,verdate).c_str());
                                return EXIT_SUCCESS;
                    case 'i':   input = myoptarg; break;
                    case 'o':   output = myoptarg; break;
                    case 'p':   optfile = myoptarg; break;
                    case 't':   noselect = true; break;
                    case 'v':   suppress = false; break;
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


    mystring        insparamfile = GetFullParamFilename();
    mystring        altinsparamfile =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetParamFilename();

    if( !file_exists( insparamfile.c_str()))
        insparamfile = altinsparamfile;

    //// hdpctx
//     mystring        inhdpctxfile = GetFullHDPCtxFilename();
//     mystring        altinhdpctxfile =
//                         mystring( my_dirname( argv[0])) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetHDPCtxFilename();
// 
//     if( !file_exists( inhdpctxfile.c_str()))
//         inhdpctxfile = altinhdpctxfile;
//     mystring        inhdpfile = inhdpctxfile;
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
    //float           valMIXWGT = MOptions::GetMIXWGT();
    mystring        valSCOADJ = MOptions::GetSCOADJ();
    //int             intscoadj = scoadjNo;
    //int             valSUPCLT = MOptions::GetSUPCLT();
    //float           valADJWGT = MOptions::GetADJWGT();
    //float           valcADJWGT = MOptions::GetcADJWGT();

    int             valHCFILTER = MOptions::GetHCFILTER();
    int             valHCWINDOW = MOptions::GetHCWINDOW();
    float           valHCLOWENT = MOptions::GetHCLOWENT();
    float           valHCHIGHENT = MOptions::GetHCHIGHENT();

    //int             valINVLCFILTER = MOptions::GetINVLCFILTER();
    int             valLCFILTEREACH = MOptions::GetLCFILTEREACH();
    int             valLCWINDOW = MOptions::GetLCWINDOW();
    float           valLCLOWENT = MOptions::GetLCLOWENT();
    float           valLCHIGHENT = MOptions::GetLCHIGHENT();
    //float           valDISTANCE = MOptions::GetDISTANCE();


    TRY
        //{{ -- sub. matrix --
        SetSTABLE( valSUBMAT );
        //}}

#if 0//printing parameters
        //{{ -- Print --
        if( noselect ) {
            sprintf( strbuf, "IDENTITY = %.2f, DELSTATE = %s", valIDENTITY, valDELSTATE? "yes": "no" );
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
        }
        //}}
#endif//printing parameters

    CATCH_ERROR_RETURN(;);


    int     ret = EXIT_SUCCESS;
    MSA*    msa = NULL;
    FILE*   fp = NULL;

    TRY
        msa = new MSA();
        if( !msa )
            throw MYRUNTIME_ERROR("Not enough memory.");

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


        msa->SetKeepSeqDescriptions( true );
        msa->SetIdentityLevel( valIDENTITY );
        msa->SetComputeDELETEstates( valDELSTATE );
        msa->SetExtentMinWindow( valMINALNPOS );
        msa->SetExtentMinSeqPercentage( valMINALNFRN );
        msa->SetPseudoCountWeight( valPCFWEIGHT );


        msa->ReadAlignment( input.c_str());
        if( noselect ) {
            msa->PlainPreprocess();
        } else {
            message( "Selecting sequences...");
            msa->SelectSequences();
        }
        msa->PutAlignment( output.c_str());

        message( "Done.");

    CATCH_ERROR_RETURN(if(fp) fclose(fp);if(msa) delete msa);

    return ret;
}
