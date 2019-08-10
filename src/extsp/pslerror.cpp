/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "pslerror.h"

using namespace extspsl;

// -------------------------------------------------------------------------
// priverror: trivial error handler
//
void extspsl::priverror( const char* errstr, 
    const char* file, unsigned int line, const char* func )
{
    fprintf( stderr, "ERROR: %s.\n", errstr );
    fprintf( stderr, "    (File: %s; Line: %u; Function: %s)\n", file, line, func );
    abort();
}

