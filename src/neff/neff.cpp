/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <numeric>
#include <algorithm>

#include "liblib/mygetopt.h"
#include "liblib/msg.h"
#include "incconf/localconfig.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcaln/MSA.h"
#include "neff.h"


// =========================================================================
// declarations of functions:
void GetIdntThresholds(mystring, std::vector<float>&);
// =========================================================================


int main( int argc, char *argv[] )
{
    int c;
    //string values of parameters
    mystring myoptarg;
    mystring input;
    mystring optfile;
    mystring idnlist;//list of sequence identity thresholds
    std::vector<float> vidnthlds(1, 0.62f);//vector of unique values of identity thresholds (def=62)
    bool suppress = true;//suppress warnings

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"i", my_required_argument, 'i'},
        {"v", my_no_argument,       'v'},
        {"h", my_no_argument,       'h'},
        {"idn", my_required_argument, 'L'},
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
                    case 'v':   suppress = false; break;
                    //
                    case 'L':   idnlist = myoptarg; break;
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


    TRY
        if( !idnlist.empty())
            GetIdntThresholds(idnlist, vidnthlds);
        if(vidnthlds.size() < 1) {
            error( "No sequence identity thresholds given." );
            return EXIT_FAILURE;
        }
        //fprintf(stderr,"Identities:\n");
        //std::for_each(vidnthlds.begin(), vidnthlds.end(), [](const float& n){fprintf(stderr,"%f\n",n);});
    CATCH_ERROR_RETURN(;);


//     mystring        insoptfile = GetFullOptionsFilename();
//     mystring        altinsoptfile =
//                         mystring( my_dirname( argv[0])) + DIRSEP +
//                         UPDIR + DIRSEP +
//                         GetParamDirectory() + DIRSEP +
//                         GetOptionsFilename();
// 
//     if( !file_exists( insoptfile.c_str()))
//         insoptfile = altinsoptfile;
// 
//     if( optfile.empty())
//         optfile = insoptfile;
// 
// 
//     TRY
//         MOptions::Read( optfile.c_str());
//     CATCH_ERROR_RETURN(;);


    int     ret = EXIT_SUCCESS;
    MSA*    msa = NULL;

    TRY
        std::vector<float> neffs(vidnthlds.size(), 0.f);
        size_t nmsaseqs = 0, length1 = 0;

        msa = new MSA();
        if( !msa )
            throw MYRUNTIME_ERROR("Not enough memory.");

        msa->ReadAlignment( input.c_str());

        message("Calculating...");
        msa->CalculateNeff(vidnthlds, neffs, nmsaseqs, length1);

        //print result
        fprintf(stdout, "\nTotal_number_of_sequences= %zu  Length_of_1st_sequence= %zu\n\n",
                nmsaseqs, length1);
        for(size_t t = 0; t < vidnthlds.size(); t++)
            fprintf(stdout, "Neff= %9.2f @ %.3g%s\n", neffs[t], vidnthlds[t], t? "": " sequence identity");
        fprintf(stdout, "\n");

    CATCH_ERROR_RETURN(;);

    return ret;
}

// -------------------------------------------------------------------------
// GetIdntThresholds: get sequence identity thresholds written in one line;
// idlist, string of identity thresholds separated by commas;
// thldvec, sorted vector of unique sequence identity thresholds
//
void GetIdntThresholds(mystring idlist, std::vector<float>& thldvec)
{
    MYMSG( "Main::GetIdntThresholds", 4 );
    char* p;
    size_t pos;
    thldvec.clear();
    for( pos = idlist.find(','); ; pos = idlist.find(',')) {
        mystring strval = idlist.substr( 0, pos );
        errno = 0;
        int c = (int)strtol( strval.c_str(), &p, 10 );
        if( errno || *p || strval.empty())
            throw MYRUNTIME_ERROR("Invalid value specified by option idn.");
        if( c < 1 || 100 < c )
            throw MYRUNTIME_ERROR(
            "A threshold specified by command-line option idn is outside the interval [1,100].");
        thldvec.push_back((float)c/100.f);
        if( pos == mystring::npos )
            break;
        idlist = idlist.substr(pos+1);
    }
    std::sort(thldvec.begin(), thldvec.end());
    auto last = std::unique(thldvec.begin(), thldvec.end());
    thldvec.erase(last, thldvec.end());
}

