/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
// #include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/mygetopt.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcpro/Db.h"
#include "makedb.h"



int main( int argc, char *argv[] )
{
    int             c;
    //string values of the parameters
    mystring        myoptarg;
    mystring        output;
    mystring        directory;
    mystring        strmaxopenfiles;
    mystring        segwindow;
    mystring        seglowent;
    mystring        seghighent;
    mystring        segdistance;
    bool            suppress = true;    //suppress output
    bool            usingseg = false;   //whether using seg

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"d", my_required_argument, 'd'},
        {"o", my_required_argument, 'o'},
        {"n", my_required_argument, 'n'},
        //
        {"U", my_no_argument,       'U'},
        {"w", my_required_argument, 'w'},
        {"f", my_required_argument, 'f'},
        {"F", my_required_argument, 'F'},
        {"D", my_required_argument, 'D'},
        //
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
                    case 'h':   fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
                                return EXIT_SUCCESS;
                    case 'v':   suppress = false; break;
                    case 'd':   directory = myoptarg; break;
                    case 'o':   output = myoptarg; break;
                    case 'n':   strmaxopenfiles = myoptarg; break;
                    case 'U':   usingseg = true; break;
                    case 'w':   segwindow = myoptarg; usingseg = true; break;
                    case 'f':   seglowent = myoptarg; usingseg = true; break;
                    case 'F':   seghighent = myoptarg; usingseg = true; break;
                    case 'D':   segdistance = myoptarg; usingseg = true; break;
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

    char*   p;
    float   tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    mystring valSUBMAT = MOptions::GetSUBMAT();

    int nmaxopenfiles = Db::cMaxNOpenedFiles;

    size_t  segwinlenval = SEGProfile::scszDefSegProWinLength;
    float   seglowentval = SEGProfile::scfpDefSegProLowEntropy;
    float   seghighentval = SEGProfile::scfpDefSegProHighEntropy;
    float   segdistanceval = SEGProfile::scfpDefSegProVecDistance;

    if( output.empty()) {
        error( "Output database name is not specified." );
        return EXIT_FAILURE;
    }
    if( !strmaxopenfiles.empty()) {
        intvalue = strtol( strmaxopenfiles.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Maximum number of open files is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 0 ) {
            error( "Maximum number of open files is negative." );
            return EXIT_FAILURE;
        }
        nmaxopenfiles = intvalue;
    }
    //{{SEG options
    if( !segwindow.empty()) {
        intvalue = strtol( segwindow.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Window length is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 2 ) {
            error( "Window length must be >1." );
            return EXIT_FAILURE;
        }
        segwinlenval = intvalue;
    }
    if( !seglowent.empty()) {
        tmpval = strtof( seglowent.c_str(), &p );
        if( errno || *p ) {
            error( "Low entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0f ) {
            error( "Low entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seglowentval = tmpval;
    }
    if( !seghighent.empty()) {
        tmpval = strtof( seghighent.c_str(), &p );
        if( errno || *p ) {
            error( "High entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0f ) {
            error( "High entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seghighentval = tmpval;
    }
    if( !segdistance.empty()) {
        tmpval = strtof( segdistance.c_str(), &p );
        if( errno || *p ) {
            error( "Distance of equivalence is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0f ) {
            error( "Distance of equivalence is negative." );
            return EXIT_FAILURE;
        }
        segdistanceval = tmpval;
    }
    //}}SEG

    int     ret = EXIT_SUCCESS;
    Db*     db = NULL;

    TRY
        //{{ -- sub. matrix --
        SetSTABLE( valSUBMAT );
        //}}

        if( !directory.empty())
            db = new Db( output.c_str(), directory.c_str());
        else
            db = new Db( output.c_str(), argv + 1, argc - 1 );

        if( !db )
            throw MYRUNTIME_ERROR("Not enough memory.");

        db->SetMaxNumberOfOpenFiles(nmaxopenfiles);

        if( usingseg ) {
            db->SetSegParameters(
                segwinlenval,
                seglowentval,
                seghighentval,
                segdistanceval
            );
            sprintf( strbuf, "SEG:  %4zu, %.2f, %.2f; Distance %.2f",
                    segwinlenval, seglowentval, seghighentval, segdistanceval );
            message( strbuf );
        }

        db->Make();

    CATCH_ERROR_RETURN(if(db) delete db);

    return ret;
}
