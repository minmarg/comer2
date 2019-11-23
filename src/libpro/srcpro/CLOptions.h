/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CLOptions__
#define __CLOptions__

#include "liblib/preproc.h"
#include "MOptions.h"

class mystring;

#define CLDECLAREOPTION( NAME, TYPE, TYPECAST, OPER ) \
    DECLAREOPTION( NAME, TYPE, TYPECAST, OPER ); \
    void AssignCLOpt##NAME( TYPE );

#define CLOPTASSIGN( NAME, VALUE ) \
    CLOptions::AssignCLOpt##NAME( VALUE );

//command-line options
namespace CLOptions {

CLDECLAREOPTION( DEV_N, mystring, mystring, );
CLDECLAREOPTION( DEV_MEM, int, int, );
CLDECLAREOPTION( DEV_NGRIDS, int, int, );
CLDECLAREOPTION( DEV_EXPCT_DBPROLEN, int, int, );
CLDECLAREOPTION( DEV_PASS2MEMP, int, float, *0.01f );
CLDECLAREOPTION( IO_NBUFFERS, int, int, );
CLDECLAREOPTION( IO_FILEMAP, int, int, );
CLDECLAREOPTION( IO_UNPINNED, int, int, );

}//namespace CLOptions

#endif//__CLOptions__
