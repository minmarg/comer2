/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "liblib/mybase.h"
#include "CLOptions.h"

using namespace CLOptions;

// =========================================================================
// MACROs
//
#define CLDEFINEASSIGNMENT( NAME, TYPE, DEFVALUE, COND ) \
  TYPE CLOptions::val##NAME##_ = DEFVALUE; \
  void CLOptions::AssignCLOpt##NAME( TYPE newvalue ) { \
    TYPE value = newvalue; \
    if(!( COND )) \
        throw MYRUNTIME_ERROR(mystring("Invalid command-line option " TOSTR(NAME))); \
    val##NAME##_ = value; \
  }

// -------------------------------------------------------------------------
// command-line options for comer
//
CLDEFINEASSIGNMENT( DEV_N, mystring, "1", 1 )
CLDEFINEASSIGNMENT( DEV_MEM, int, -1, value>=100 && value<=1000000 )
CLDEFINEASSIGNMENT( DEV_NGRIDS, int, 1, value>=1 && value<=20 );
CLDEFINEASSIGNMENT( DEV_EXPCT_DBPROLEN, int, 50, value>=20 && value<=200 );
CLDEFINEASSIGNMENT( DEV_PASS2MEMP, int, 10, value>=1 && value<=100 )
CLDEFINEASSIGNMENT( IO_NBUFFERS, int, 4, value>=1 && value<=10 );
CLDEFINEASSIGNMENT( IO_FILEMAP, int, 0, value==0 || value==1 );
CLDEFINEASSIGNMENT( IO_UNPINNED, int, 0, value==0 || value==1 );

// command-line options for makecov
//
CLDEFINEASSIGNMENT( XCOV_USEEXTENTS, int, 0, value==0 || value==1 );
CLDEFINEASSIGNMENT( XCOV_MIXTRGFRQS, int, 0, value==0 || value==1 );
CLDEFINEASSIGNMENT( XCOV_CORR, int, 0, value==0 || value==1 );
CLDEFINEASSIGNMENT( XCOV_SCALEWGTS, int, 0, value==0 || value==1 );
CLDEFINEASSIGNMENT( XCOV_MI, int, 0, value==0 || value==1 );
