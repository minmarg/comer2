/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "liblib/msg.h"
#include "liblib/mygetopt.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "libpro/srcpro/Db.h"
#include "db2bin.h"



int main( int argc, char *argv[] )
{
    int c;
	char* p;
	//string values of the parameters
    mystring myoptarg;
    mystring input;
	mystring verbose;
	int verblev = 0;//suppress warnings

    SetArguments( &argc, &argv );
    SetProgramName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage(argv[0],makeinst,version,verdate).c_str());
        return EXIT_SUCCESS;
    }

    static struct myoption long_options[] = {
        {"d", my_required_argument, 'd'},
        //
        {"v", my_optional_argument, 'v'},
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
					case 'v':   verblev = 1; verbose = myoptarg; break;
					case 'd':   input = myoptarg; break;
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


	if (!verbose.empty()) {
		verblev = strtol(verbose.c_str(), &p, 10);
		if (errno || *p || verblev < 0) {
			error("Invalid verbose mode argument.");
			return EXIT_FAILURE;
		}
	}

	SetVerboseMode(verblev);


    mystring valSUBMAT = MOptions::GetSUBMAT();

    if( input.empty()) {
        error( "Input database name is not specified." );
        return EXIT_FAILURE;
    }

    int     ret = EXIT_SUCCESS;
    Db*     db = NULL;

    TRY
        //{{ -- sub. matrix --
        SetSTABLE( valSUBMAT );
        //}}

        db = new Db( input.c_str());

        if( !db )
            throw MYRUNTIME_ERROR("Not enough memory.");

        db->MakeBin();

    CATCH_ERROR_RETURN(if(db) delete db);

    return ret;
}
