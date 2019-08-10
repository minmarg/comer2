/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "libpro/srcpro/PMTransModel.h"
// #include "libpro/srcpro/PMProfileModel.h"
#include "Db.h"
#include "VirtualRoDb.h"

// -------------------------------------------------------------------------
// Constructor
//
VirtualRoDb::VirtualRoDb( const char* name/*, int fstype */)
:   Db( name )
{
    MYMSG("VirtualRoDb::VirtualRoDb",3);
}

// Default constructor
//
VirtualRoDb::VirtualRoDb()
:   Db()
{
    throw MYRUNTIME_ERROR( "VirtualRoDb::VirtualRoDb: Default initialization is not allowed." );
}

// -------------------------------------------------------------------------
// Destructor
//
VirtualRoDb::~VirtualRoDb()
{
    MYMSG("VirtualRoDb::~VirtualRoDb",3);
}

// -------------------------------------------------------------------------
// Open: open the main database or input file and initialize a file 
// descriptor
//
void VirtualRoDb::Open()
{
    MYMSG("VirtualRoDb::Open",3);

    if( !GetDbName())
        throw MYRUNTIME_ERROR( "VirtualRoDb::Open: Unable to open input file." );

    if( GetDbDesc() != NULL )
        throw MYRUNTIME_ERROR( "VirtualRoDb::Open: Input file has been already opened." );

    SetDbDesc( fopen( GetDbMainName(), "r" ));
    if( GetDbDesc() == NULL ) {
        SetDbDesc( fopen( GetDbName(), "r" ));
        if( GetDbDesc() == NULL ) {
            throw MYRUNTIME_ERROR( "VirtualRoDb::Open: Failed to open input file." );
        }
    }

    try {
        SetNextSym( 0 );
        ReadTextHeader( GetDbDesc());

    } catch( myexception const& /*ex*/ )
    {
        //assume this is an input batch file of profiles and 
        //try reading profiles later
        if( fseek( GetDbDesc(), 0L, SEEK_SET ) != 0 ) {
            Close();
            throw MYRUNTIME_ERROR( 
            "VirtualRoDb::Open: Failed to fseek to the beginning of input file." );
        }
    }
}
