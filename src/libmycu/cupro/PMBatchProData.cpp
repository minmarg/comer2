/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdlib.h>
#include <string.h>

#include "liblib/msg.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "PMBatchProData.h"

// _________________________________________________________________________
// PMBatchProData constructor
//
PMBatchProData::PMBatchProData()
:   h_descs_( NULL ),
    h_files_( NULL ),
    h_descsend_( NULL ),
    h_filesend_( NULL )
{
    MYMSG("PMBatchProData::PMBatchProData",5);
    h_pmdata_ = (char**)malloc(pmv2DTotFlds * sizeof(void*));
    h_pmdatend_ = (char**)malloc(pmv2DTotFlds * sizeof(void*));
    szpmdata_ = (size_t*)malloc(pmv2DTotFlds * sizeof(size_t));

    if( h_pmdata_ == NULL || h_pmdatend_ == NULL || szpmdata_ == NULL )
        throw MYRUNTIME_ERROR("PMBatchProData::PMBatchProData: Not enough memory.");

    memset( h_pmdata_, 0, pmv2DTotFlds * sizeof(void*));
    memset( h_pmdatend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( szpmdata_, 0, pmv2DTotFlds * sizeof(size_t));
}

// -------------------------------------------------------------------------
// destructor
//
PMBatchProData::~PMBatchProData()
{
    MYMSG("PMBatchProData::~PMBatchProData",5);
    DestroyHostData();
    if( h_pmdata_ ) { free(h_pmdata_); h_pmdata_ = NULL; }
    if( h_pmdatend_ ) { free(h_pmdatend_); h_pmdatend_ = NULL; }
    if( szpmdata_ ) { free(szpmdata_); szpmdata_ = NULL; }
}

// -------------------------------------------------------------------------
// AllocNewSizeHostData: allocate new memory for profile data buffers;
// maxdatasize, maximum amount of memory that can be allocated for profile 
// data;
// initial call;
//
void PMBatchProData::AllocNewSizeHostDataInit( 
    const size_t maxdatasize, TAllocFunc allocfunc )
{
    MYMSG("PMBatchProData::AllocNewSizeHostData init",3);
    const mystring preamb = "PMBatchProData::AllocNewSizeHostData: ";
    const size_t initlen = 1000;//initial/predicted profile length
    //predicted occupied sizes
    size_t szpmdatafilled[pmv2DTotFlds];
    GetPMDataSize1( initlen, szpmdatafilled );
    //AllocHostDataCalc( maxdatasize, szpmdatafilled, &PMBatchProData::AllocHostData );
    AllocHostDataCalc( maxdatasize, szpmdatafilled, allocfunc );
}

// -------------------------------------------------------------------------
// ReallocNewSizeHostData: reallocate new memory for profile data buffers;
// newdatasize, amount of memory that can be allocated for profile data; 
// information stored currently in the buffers will be PRESERVED;
//
void PMBatchProData::ReallocNewSizeHostData( const size_t newdatasize )
{
    MYMSG("PMBatchProData::ReallocNewSizeHostData",3);
//     if( GetNoProsWritten() < 1 )
//         AllocNewSizeHostDataInit( newdatasize );
//     else
        AllocNewSizeHostDataHelper( newdatasize, &PMBatchProData::ReallocHostData );
}

// -------------------------------------------------------------------------
// AllocNewSizeHostData: allocate new memory for profile data buffers;
// maxdatasize, maximum amount of memory that can be allocated for profile 
// data; information stored currently in the buffers will be destroyed;
//
void PMBatchProData::AllocNewSizeHostData( const size_t maxdatasize )
{
    MYMSG("PMBatchProData::AllocNewSizeHostData",3);
//     if( GetNoProsWritten() < 1 )
//         AllocNewSizeHostDataInit( maxdatasize );
//     else
        AllocNewSizeHostDataHelper( maxdatasize, &PMBatchProData::AllocHostData );
}

// -------------------------------------------------------------------------
// AllocNewSizeHostDataHelper: allocate new memory for profile data buffers;
// maxdatasize, maximum amount of memory that can be allocated for profile 
// data;
//
void PMBatchProData::AllocNewSizeHostDataHelper( 
    const size_t maxdatasize, TAllocFunc allocfunc )
{
    MYMSG("PMBatchProData::AllocNewSizeHostDataHelper",3);
    const mystring preamb = "PMBatchProData::AllocNewSizeHostDataHelper: ";

    const size_t nflds = GetNoFields();
    size_t szpmdatafilled[pmv2DTotFlds/*nflds*/];//occupied sizes for each buffer
    int i;
    for( i = 0; i < (int)nflds; i++ ) {
#ifdef __DEBUG__
        if( h_pmdatend_[i] < h_pmdata_[i] )
            throw MYRUNTIME_ERROR( preamb + "Invalid buffer addresses." );
#endif
        szpmdatafilled[i] = (size_t)(h_pmdatend_[i]-h_pmdata_[i]);
        if( szpmdatafilled[i] < 1 )
            //zero allocation might have occurred if the last difference between the 
            // limit and cached profile data size is small
            return AllocNewSizeHostDataInit( maxdatasize, allocfunc );
            //throw MYRUNTIME_ERROR( preamb + "Invalid size of occupied memory." );
    }
    AllocHostDataCalc( maxdatasize, szpmdatafilled, allocfunc );
}

// -------------------------------------------------------------------------
// MemDeficit: verify whether a new portion of data overflows the allocated 
// memory;
// pmdataend points to the ends of filled buffers;
// 
bool PMBatchProData::MemDeficit( size_t prolen )
{
    MYMSG("PMBatchProData::MemDeficit",5);
    const mystring preamb = "PMBatchProData::MemDeficit: ";
    size_t szpmfieldfilled[pmv2DTotFlds];//occupied sizes for each buffer
    size_t szpmdata1[pmv2DTotFlds];
    size_t sztotal = 0;
    bool defc = false;
    int i, j;

    GetPMDataSize1( prolen, szpmdata1 );

    for( i = j = 0; i < pmv2DTotFlds; i++ ) {
#ifdef __DEBUG__
        if( h_pmdatend_[i] < h_pmdata_[i] )
            throw MYRUNTIME_ERROR( preamb + "Invalid size of occupied memory." );
#endif
        //occupied size of of a buffer 
        // plus the (upper bound of) size for a next profile
        szpmfieldfilled[i] = (size_t)(h_pmdatend_[i]-h_pmdata_[i]);
        szpmfieldfilled[i] += szpmdata1[i];
        
        sztotal += szpmfieldfilled[i];

        if( szpmdata_[i] < szpmfieldfilled[i] ) {
            if( !defc )
                j = i;
            defc = true;
        }
    }

    MYMSGBEG
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "%s %.1fM of filled data vs %.1fM of allocated", 
            preamb.c_str(), (float)(sztotal)/(float)ONEM,
            (float)GetTotalAllocSize()/(float)ONEM);
        MYMSG( msgbuf,5 );
        if( defc ) {
            sprintf( msgbuf, "%sBREAKPOINT: buf[%d]=%zu > alloc[%d]=%zu", 
                preamb.c_str(), j, szpmfieldfilled[j], j, szpmdata_[j]);
            MYMSG( msgbuf,5 );
        }
    MYMSGEND

    return defc;
}
