/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuInputReader_h__
#define __CuInputReader_h__

#include <stdio.h>
#include "liblib/mybase.h"
#include "liblib/myfiler.h"
#include "libmycu/cupro/CuDbReader.h"

// _________________________________________________________________________
// Class CuInputReader
//
// Reader class for an input file which can be a database of profiles; 
// it can also represent a batch of input profiles
//
class CuInputReader: public CuDbReader
{
public:
    CuInputReader( const char* name, bool mapped ): CuDbReader(name, mapped) {
        MYMSG("CuInputReader::CuInputReader",4);
    };
    virtual ~CuInputReader() {
        MYMSG("CuInputReader::~CuInputReader",4);
    };

    virtual void Open();

protected:
    explicit CuInputReader() {
        throw MYRUNTIME_ERROR(
        "CuInputReader::CuInputReader: Default initialization prohibited.");
    };

private:
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
// Open: open the main database or input file and initialize a file 
// descriptor
inline
void CuInputReader::Open()
{
    MYMSG("CuInputReader::Open",3);

    try {
        OpenFile(Db::DbFlMAIN, GetDbMainName());
        SetOpenedFileType(Db::DbFlMAIN);
    } catch(myexception const& )
    {
        try {
            OpenFile(Db::DbFlMAIN, GetDbName());
            SetOpenedFileType(Db::DbFlMAIN);
        } catch(myexception const& )
        {
            throw MYRUNTIME_ERROR(
            "CuInputReader::Open: Failed to open input file.");
        }
    }

    ResetCurrentPageData();
    ResetProfileBufferData();

    SetDbSizeInBytes( ObtainDbSizeInBytes( GetOpenedFileType()));

    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "CuInputReader::Open: file_size %zu page_size %zu",
                GetDbSizeInBytes(), GetPageSize());
        MYMSG( msgbuf, 3 );
    MYMSGENDl

    if( GetMapped())
        MapFile(Db::DbFlMAIN);

    try {
        CacheProfile(true);//cache one page only

        ReadTextHeader( GetProfileBuffer());

    } catch(myexception const& )
    {
        //assume this is an input batch file of profiles and 
        //try reading profiles later
        ResetProfileBufferPosition();
    }
}

#endif//__CuInputReader_h__
