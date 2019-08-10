/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuRoDb_h__
#define __CuRoDb_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include "libpro/srcpro/Db.h"

// _________________________________________________________________________
// Class CuRoDb
//
// Read-only profile database for batch processing using CUDA or pthreads
//
class CuRoDb: virtual public Db
{
public:
    CuRoDb( const char* name );
    virtual ~CuRoDb();

    bool NextProfileHeader( mystring* desc, mystring* file,
            int* prolen, int* scale, char** pmdata );
    bool Next( unsigned int distn,
            int profnr, int prolen, int scale, 
            char** pmorgaddr, char** pmdata );

protected:
    explicit CuRoDb();

private:
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//

#endif//__CuRoDb_h__
