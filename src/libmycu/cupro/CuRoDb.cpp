/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libpro/srcpro/Db.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "CuRoDb.h"

// -------------------------------------------------------------------------
// Constructor
//
CuRoDb::CuRoDb( const char* name )
:   Db( name )
{
    MYMSG("CuRoDb::CuRoDb",4);
}

// Default constructor
//
CuRoDb::CuRoDb()
:   Db()
{
    throw MYRUNTIME_ERROR( "CuRoDb::CuRoDb: Default initialization is not allowed." );
}

// -------------------------------------------------------------------------
// Destructor
//
CuRoDb::~CuRoDb()
{
    MYMSG("CuRoDb::~CuRoDb",4);
}

// -------------------------------------------------------------------------
// NextProfileHeader: read the header of the next profile from the database; 
// return a flag indicating whether the database end has been reached; 
// eof would mean an error
//
bool CuRoDb::NextProfileHeader( 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata )
{
    MYMSG("CuRoDb::NextProfileHeader",5);

    int nsym = 0;

    if( GetNextSym() == EOF )
        return false;

    if( GetDbDesc() == NULL )
        throw MYRUNTIME_ERROR( "CuRoDb::NextProfileHeader: Db is not opened." );

    try {
        TextReadProfileHeader( GetDbDesc(), desc, file, prolen, scale, pmdata );

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error(ex);
    }

    if( GetNextSym() != EOF ) {
        nsym = fgetc( GetDbDesc());
        ungetc( nsym, GetDbDesc());
        SetNextSym( nsym );
    }

    return true;
}

// -------------------------------------------------------------------------
// Next: read the data of the next profile from the database; 
// return a flag indicating whether the database end has been reached
//
bool CuRoDb::Next( 
    unsigned int distn,
    int profnr, int prolen, int scale, 
    char** pmorgaddr, char** pmdata )
{
    MYMSG("CuRoDb::Next",5);

    int nsym = 0;

    if( GetNextSym() == EOF )
        return false;

    if( GetDbDesc() == NULL )
        throw MYRUNTIME_ERROR( "Db::Next: Db is not opened." );

    try {
        TextReadProfileData( 
            GetDbDesc(), distn, profnr, prolen, scale, 
            pmorgaddr, pmdata );

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error(ex);
    }

    //NOTE: Prepare() for initializating transitions

    if( GetNextSym() != EOF ) {
        nsym = fgetc( GetDbDesc());
        ungetc( nsym, GetDbDesc());
        SetNextSym( nsym );
    }

    return true;
}
