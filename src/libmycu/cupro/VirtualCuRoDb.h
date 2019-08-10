/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __VirtualCuRoDb_h__
#define __VirtualCuRoDb_h__

#include <stdio.h>
#include "liblib/mybase.h"
#include "libpro/srcpro/Db.h"
#include "libpro/srcpro/VirtualRoDb.h"
#include "libmycu/cupro/CuRoDb.h"

// _________________________________________________________________________
// Class VirtualCuRoDb
//
// Virtual read-only database of profiles, which can also represent a 
// batch of input profiles
//
class VirtualCuRoDb: public CuRoDb, public VirtualRoDb
{
public:
    VirtualCuRoDb( const char* name ): Db(name), CuRoDb(name), VirtualRoDb(name) {
        MYMSG("VirtualCuRoDb::VirtualCuRoDb",4);
    };
    virtual ~VirtualCuRoDb() {
        MYMSG("VirtualCuRoDb::~VirtualCuRoDb",4);
    };

	virtual void Open() {VirtualRoDb::Open();}

protected:
    explicit VirtualCuRoDb() {
        throw MYRUNTIME_ERROR("VirtualCuRoDb::VirtualCuRoDb: Default initialization prohibited.");
    };

private:
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//

#endif//__VirtualCuRoDb_h__
