/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __VirtualRoDb_h__
#define __VirtualRoDb_h__

#include "liblib/mybase.h"
#include <stdio.h>
#include "Db.h"

// _________________________________________________________________________
// Class VirtualRoDb
//
// Virtual read-only database of profiles, which can also represent a 
// batch of input profiles
//
class VirtualRoDb: virtual public Db 
{
public:
    VirtualRoDb( const char* name );
    virtual ~VirtualRoDb();

    virtual void Open();

protected:
    explicit VirtualRoDb();

private:
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//

#endif//__VirtualRoDb_h__
