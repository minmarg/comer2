/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMBatchProData_h__
#define __PMBatchProData_h__

#include <stdlib.h>
#include <string.h>

#include <memory>

#include "liblib/mybase.h"
#include "libmycu/cupro/PM2DVectorFields.h"

// _________________________________________________________________________
// Class PMBatchProData
//
// Complete batch profile data for parallel processing
//
class PMBatchProData
{
    typedef void (PMBatchProData::*TAllocFunc)( size_t* );

public:
    PMBatchProData();
    ~PMBatchProData();

    mystring* GetDescs() { return h_descs_; }
    const mystring* GetDescs() const { return h_descs_; }
    mystring* GetFnames() { return h_files_; }
    const mystring* GetFnames() const { return h_files_; }
    char** GetPMDataEnd() { return h_pmdatend_; }
    char** GetPMData() { return h_pmdata_; }
    void GetPMData( char** pmdata );
    static void PMDataNextPro( char** pmdata );
    bool PMDataReachedEnd( char** pmdata );
    static size_t GetPMDataLen1At( char** pmdata );
    static float GetPMDataENO1At( char** pmdata );
    static size_t GetPMDataSize1At( char** pmdata );
    size_t GetNoProsAt( char** pmdata ) const;
//     size_t GetFilledSizeFromTo( char** pmdatbeg, char** pmdatend ) const;
    static size_t GetNoProsFromTo( char** pmdatbeg, char** pmdatend );
    static size_t GetNoPositsFromTo( char** pmdatbeg, char** pmdatend );
    static bool PMDataPresentFromTo( char** pmdatbeg, char** pmdatend );

    bool MemDeficit( size_t prolen );
    void AllocNewSizeHostDataInit( const size_t maxdatasize, TAllocFunc = &PMBatchProData::AllocHostData );
    void AllocNewSizeHostData( const size_t maxdatasize );
    void ReallocNewSizeHostData( const size_t newdatasize );
    void ReallocHostData();//for making occupation == allocation

    static constexpr size_t GetNoFields() { return pmv2DTotFlds; }
    static constexpr size_t GetNoHeaderFields() { return pps2DDist+1; }
    size_t GetNoProsWritten() const {//# profiles written in the buffers
        return (size_t)(h_pmdatend_[pps2DLen]-h_pmdata_[pps2DLen]) / SZINTYPE;
    }

    size_t GetTotalAllocSize() const;
    size_t GetTotalFilledSize() const;

    //data size for complete profile model data of one profile
    static void GetPMDataSize1( size_t prolen, size_t* szpmdata1 );// const;
    //upper bound of size for complete profile model data given number of positions
    static void GetPMDataSizeUB( size_t totprolen, size_t* szpmdata1 );
    static size_t GetPMDataSizeUBTotal( size_t totprolen );//size in total
    //lower bound of size for complete profile model data given number of positions
    static void GetPMDataSizeLB( size_t totprolen, size_t* szpmdata1 ) { 
        GetPMDataSize1( totprolen, szpmdata1 );
    };
    static size_t GetPMDataSizeLBTotal( size_t totprolen );//size in total

    //copy profile model header from the given addresses
    void CopyPMDataHeaderNext( mystring& desc, mystring& file, char** src );

protected:
    void DestroyHostData();
    void AllocHostData( size_t* szpmdata );
    void ReallocHostData( size_t* szpmdata );

    void ReallocMystrings( mystring** beg, mystring** end, size_t newsize );

    void AllocNewSizeHostDataHelper( const size_t maxdatasize, TAllocFunc );
    void AllocHostDataCalc( const size_t maxdatasize, size_t* szpmdatafilled, TAllocFunc );

private:
    //{{host data
    mystring* h_descs_;//profile descriptions
    mystring* h_files_;//profile filenames
    mystring* h_descsend_;//end pointer
    mystring* h_filesend_;//end pointer
    //profile-specific and position-specific profile model data for each profile in the buffer:
    char** h_pmdata_;//[pmv2DTotFlds];
    char** h_pmdatend_;//[pmv2DTotFlds];//pointers pointing to the ends of all h_pmdata_ buffers
    //}}
    size_t* szpmdata_;//[pmv2DTotFlds];//allocated size
};

////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// ReallocMystrings: reallocate descriptor or file string table
//
inline
void PMBatchProData::ReallocMystrings( mystring** beg, mystring** end, size_t newsize )
{
    MYMSG("PMBatchProData::ReallocMystrings",5);
    const mystring preamb = "PMBatchProData::ReallocMystrings: ";
#ifdef __DEBUG__
    if( !beg || !end || *end < *beg )
        throw MYRUNTIME_ERROR( preamb + "Invalid arguments.");
#endif
    int i;
    mystring* p;
    mystring* tmpstrs = new mystring[newsize];
    if( !tmpstrs )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
    for( i = 0, p = *beg; p < *end; p++, i++ )
        tmpstrs[i] = *p;
    *end = tmpstrs + (size_t)( *end - *beg );
    delete [] *beg;
    *beg = tmpstrs;
}

// -------------------------------------------------------------------------
// ReallocHostData: reallocate memory for host data
//
inline
void PMBatchProData::ReallocHostData( size_t* szpmdata )
{
    MYMSG("PMBatchProData::ReallocHostData",5);
    const mystring preamb = "PMBatchProData::ReallocHostData: ";
#ifdef __DEBUG__
    if( !szpmdata )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
#endif
    int i;
    size_t nopros = szpmdata[pps2DLen] / SZINTYPE + 1;//number of profiles to contain in the buffers
    MYMSGBEG
        char msgbuf[BUF_MAX];
        size_t sztotal = 0;
        for( i = 0; i < pmv2DTotFlds; i++ )
            sztotal += szpmdata[i];
        sprintf(msgbuf,"%sa total of %.1fM requested for %zu profiles (#fields=%d).",
            preamb.c_str(),(float)(sztotal)/(float)ONEM, nopros, pmv2DTotFlds );
        MYMSG( msgbuf, 5 );
    MYMSGEND

#ifdef __DEBUG__
    if( h_descsend_ < h_descs_ || h_filesend_ < h_files_ )
        throw MYRUNTIME_ERROR( preamb + "Invalid end addresses.");
#endif
    ReallocMystrings( &h_descs_, &h_descsend_, nopros );
    ReallocMystrings( &h_files_, &h_filesend_, nopros );

    for( i = 0; i < pmv2DTotFlds; i++ ) {
        if( szpmdata[i]) {
            char* tmp_pmdata_i = (char*)realloc( h_pmdata_[i], szpmdata[i]);
            if( !tmp_pmdata_i )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
            h_pmdatend_[i] = tmp_pmdata_i + (size_t)(h_pmdatend_[i]-h_pmdata_[i]);
            h_pmdata_[i] = tmp_pmdata_i;
        }
        else {
            free( h_pmdata_[i]);
            h_pmdata_[i] = NULL;
            h_pmdatend_[i] = NULL;
        }
        szpmdata_[i] = szpmdata[i];
    }
}

// -------------------------------------------------------------------------
// ReallocHostData: shrink allocated memory for host data to the 
// boundaries of occupied memory
//
inline
void PMBatchProData::ReallocHostData()
{
    MYMSG("PMBatchProData::ReallocHostData [free unused memory]",5);
    const mystring preamb = "PMBatchProData::ReallocHostData: ";

    const size_t nflds = GetNoFields();
    //size_t szpmdatafilled[nflds];//occupied sizes for each buffer
	std::unique_ptr<size_t[]> szpmdatafilled(new size_t[nflds]);//occupied sizes for each buffer
    int i;
    for( i = 0; i < (int)nflds; i++ ) {
#ifdef __DEBUG__
        if( h_pmdatend_[i] < h_pmdata_[i] )
            throw MYRUNTIME_ERROR( preamb + "Invalid buffer addresses." );
#endif
        szpmdatafilled[i] = (size_t)(h_pmdatend_[i]-h_pmdata_[i]);
#ifdef __DEBUG__
        //some fields can be unused
        //if( szpmdatafilled[i] < 1 )
        //    throw MYRUNTIME_ERROR( preamb + "Invalid size of occupied memory." );
#endif
    }
    ReallocHostData( szpmdatafilled.get());
}

// -------------------------------------------------------------------------
// DestroyHostData: destroy host data
//
inline
void PMBatchProData::DestroyHostData()
{
    MYMSG("PMBatchProData::DestroyHostData",5);
    if( h_descs_ ) { delete [] h_descs_; h_descs_ = NULL; }
    if( h_files_ ) { delete [] h_files_; h_files_ = NULL; }
    h_descsend_ = NULL;
    h_filesend_ = NULL;
    for( int i = 0; i < pmv2DTotFlds; i++ ) {
        if( h_pmdata_[i]) {
            free( h_pmdata_[i]);
            h_pmdata_[i] = NULL;
        }
        h_pmdatend_[i] = NULL;
        szpmdata_[i] = 0;
    }
}

// -------------------------------------------------------------------------
// AllocHostData: allocate memory for host data
//
inline
void PMBatchProData::AllocHostData( size_t* szpmdata )
{
    MYMSG("PMBatchProData::AllocHostData",5);
#ifdef __DEBUG__
    if( !szpmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::AllocHostData: Memory access error.");
#endif
    int i;
    size_t nopros = szpmdata[pps2DLen] / SZINTYPE + 1;//number of profiles to contain in the buffers
    MYMSGBEGl(5)
        char msgbuf[BUF_MAX];
        size_t sztotal = 0;
        for( i = 0; i < pmv2DTotFlds; i++ )
            sztotal += szpmdata[i];
        sprintf(msgbuf,"PMBatchProData::AllocHostData: "
            "a total of %.1fM requested for %zu profiles (#fields=%d).",
            (float)(sztotal)/(float)ONEM, nopros, pmv2DTotFlds );
        MYMSG( msgbuf, 5 );
    MYMSGEND

    DestroyHostData();

    h_descs_ = new mystring[nopros];
    if( !h_descs_)
        throw MYRUNTIME_ERROR("PMBatchProData::AllocHostData: Not enough memory.");
    h_descsend_ = h_descs_;
    h_files_ = new mystring[nopros];
    if( !h_files_ )
        throw MYRUNTIME_ERROR("PMBatchProData::AllocHostData: Not enough memory.");
    h_filesend_ = h_files_;
    for( i = 0; i < pmv2DTotFlds; i++ ) {
        if( szpmdata[i]) {
            h_pmdata_[i] = (char*)malloc(szpmdata[i]);
            if( !h_pmdata_[i] )
                throw MYRUNTIME_ERROR("PMBatchProData::AllocHostData: Not enough memory.");
            h_pmdatend_[i] = h_pmdata_[i];
        }
        szpmdata_[i] = szpmdata[i];
    }
}

// -------------------------------------------------------------------------
// AllocHostDataCalc: allocate memory for host data after calculating the 
// amount of memory for each buffer
//
inline
void PMBatchProData::AllocHostDataCalc( 
    const size_t maxdatasize, size_t* szpmdatafilled, TAllocFunc allocfunc )
{
    MYMSG("PMBatchProData::AllocHostDataCalc",3);
#ifdef __DEBUG__
    if( !szpmdatafilled )
        throw MYRUNTIME_ERROR("PMBatchProData::AllocHostDataCalc: Memory access error.");
#endif
    MYMSGBEGl(4)
        char msgbuf[BUF_MAX];
        sprintf(msgbuf,"PMBatchProData::AllocHostDataCalc: max size %zu (%.1fM).", 
                maxdatasize, (float)maxdatasize/(float)ONEM);
        MYMSG( msgbuf, 4 );
    MYMSGEND

    int i;
    size_t sztotal = 0;//total occupied size
    float proppmdat[pmv2DTotFlds];//proportions of how much data occupies each buffer
    size_t szpmdata[pmv2DTotFlds];//new sizes calculated wrt previous occupation

    for( i = 0; i < pmv2DTotFlds; i++ )
        sztotal += szpmdatafilled[i];

#ifdef __DEBUG__
    if( sztotal < 1 )
        throw MYRUNTIME_ERROR("PMBatchProData::AllocHostDataCalc: Invalid size of used memory.");
#endif

    for( i = 0; i < pmv2DTotFlds; i++ ) {
        proppmdat[i] = (float)szpmdatafilled[i]/(float)sztotal;
        szpmdata[i] = (size_t)(proppmdat[i]*maxdatasize);
    }

    //NOTE: increase the size of profile-specific buffers eightfold, as these are 
    // small compared to position-specific buffers and fill much faster recording an 
    // out-of-memory event;
    for( i = 0; i <= pps2DDist/*pps2DLen*/; i++ )
        szpmdata[i] *= 8;

    (this->*allocfunc)(szpmdata);
}

// -------------------------------------------------------------------------
// GetTotalAllocSize: get total size of memory allocated for profile model 
// data
//
inline
size_t PMBatchProData::GetTotalAllocSize() const
{
    MYMSG("PMBatchProData::GetTotalAllocSize",6);
    size_t sztotal = 0;//total allocated size
    for( int i = 0; i < pmv2DTotFlds; i++ )
        sztotal += szpmdata_[i];
    return sztotal;
}

// -------------------------------------------------------------------------
// GetTotalFilledSize: get total size of memory filled with profile model 
// data
//
inline
size_t PMBatchProData::GetTotalFilledSize() const
{
    MYMSG("PMBatchProData::GetTotalFilledSize",5);
    const mystring preamb = "PMBatchProData::GetTotalFilledSize: ";
    const size_t nflds = GetNoFields();
    size_t szfilled = 0;//total occupied size over all buffers
    int i;
    for( i = 0; i < (int)nflds; i++ ) {
#ifdef __DEBUG__
        if( h_pmdatend_[i] < h_pmdata_[i] )
            throw MYRUNTIME_ERROR( preamb + "Invalid buffer addresses." );
#endif
        szfilled += (size_t)(h_pmdatend_[i]-h_pmdata_[i]);
    }
    return szfilled;
}

// -------------------------------------------------------------------------
// GetPMData: get all addresses of the buffers containing profile data
//
inline
void PMBatchProData::GetPMData( char** pmdata )
{
    MYMSG("PMBatchProData::GetPMData",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMData: Null argument.");
#endif
    for( int i = 0; i < pmv2DTotFlds; i++ )
        pmdata[i] = h_pmdata_[i];
}

// -------------------------------------------------------------------------
// PMDataReachedEnd: return a flag of whether pointer of the buffers 
// reached the end of data
//
inline
bool PMBatchProData::PMDataReachedEnd( char** pmdata )
{
    MYMSG("PMBatchProData::PMDataReachedEnd",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::PMDataReachedEnd: Null argument.");
#endif
    for( int i = 0; i < pmv2DTotFlds; i++ )
        if( h_pmdatend_[i] <= pmdata[i] )
            return true;
    return false;
}

// -------------------------------------------------------------------------
// PMDataNextPro: move pointers to the buffers containing profile 
// data to the position of the next profile
//
inline
void PMBatchProData::PMDataNextPro( char** pmdata )
{
    MYMSG("PMBatchProData::PMDataNextPro",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::PMDataNextPro: Null argument.");
#endif
    int i;
    size_t prolen = (size_t)*(INTYPE*)pmdata[pps2DLen];

    for( i = 0; i < pmv2DNoElems; i++ )
        pmdata[pps2DBkgPrbs+i] += SZFPTYPE;
    pmdata[pps2DENO] += SZFPTYPE;
    pmdata[pps2DLen] += SZINTYPE;
    pmdata[pps2DDist] += SZLNTYPE;

    for( i = 0; i < ptr2DNoElems; i++ ) {
        pmdata[ptr2DTrnPrbs+i] += SZFPTYPE * (prolen+1);
//         pmdata[ptr2DTrnPrbsExp+i] += SZFPTYPE * (prolen+1);
    }

    for( i = 0; i < pmv2DNoElems; i++ )
        pmdata[pmv2DTrgFrqs+i] += SZFPTYPE * prolen;
    for( i = 0; i < pmv2DNoCVEls; i++ )
        pmdata[pmv2DCVentrs+i] += SZFPTYPE * prolen;
    pmdata[pmv2DCVprior] += SZFPTYPE * prolen;
    pmdata[pmv2DCVnorm2] += SZFPTYPE * prolen;
    for( i = 0; i < pmv2DNoSSSps; i++ )
        pmdata[pmv2DSSsprbs+i] += SZFPTYPE * prolen;
    pmdata[pmv2DHDP1prb] += SZFPTYPE * prolen;
    pmdata[pmv2DHDP1ind] += SZINTYPE * prolen;
    pmdata[pmv2DAddrPro] += SZINTYPE * prolen;
    pmdata[pmv2DAddrCV] += SZLNTYPE * prolen;
    pmdata[pmv2DAddrSS] += SZLNTYPE * prolen;
    pmdata[pmv2Daa] += SZCHTYPE * prolen;
    pmdata[pmv2DSSstate] += SZCHTYPE * prolen;
}

// -------------------------------------------------------------------------
// GetNoProsAt: get the number of profiles from the beginning of the 
// buffers
inline
size_t PMBatchProData::GetNoProsAt( char** pmdata ) const 
{
    MYMSG("PMBatchProData::GetNoProsAt",7);
#ifdef __DEBUG__
    if( !pmdata || pmdata[pps2DLen] < h_pmdata_[pps2DLen] )
        throw MYRUNTIME_ERROR("PMBatchProData::GetNoProsAt: Invalid argument.");
#endif
    return (size_t)(pmdata[pps2DLen]-h_pmdata_[pps2DLen]) / SZINTYPE;
}

// // -------------------------------------------------------------------------
// // GetFilledSizeAt: get the size of memory filled with profile model 
// // data from one position up to the other
// //
// inline
// size_t PMBatchProData::GetFilledSizeFromTo( char** pmdatbeg, char** pmdatend ) const
// {
//     MYMSG("PMBatchProData::GetFilledSizeFromTo",6);
//     const mystring preamb = "PMBatchProData::GetFilledSizeFromTo: ";
// #ifdef __DEBUG__
//     if( !pmdatbeg || !pmdatend )
//         throw MYRUNTIME_ERROR( preamb + "Null arguments.");
// #endif
//     const size_t nflds = GetNoFields();
//     size_t szfilled = 0;//total occupied size over all buffers
//     int i;
//     for( i = 0; i < (int)nflds; i++ ) {
// #ifdef __DEBUG__
//         if( pmdatend[i] < pmdatbeg[i] )
//             throw MYRUNTIME_ERROR( preamb + "Invalid argument buffer addresses." );
// #endif
//         szfilled += (size_t)(pmdatend[i]-pmdatbeg[i]);
//     }
//     return szfilled;
// }

// -------------------------------------------------------------------------
// GetNoProsFromTo: get the number of profiles present in the buffers;
// pmdatbeg and pmdatend are the beginning and terminal positions of buffers
//
inline
size_t PMBatchProData::GetNoProsFromTo( char** pmdatbeg, char** pmdatend )
{
    MYMSG("PMBatchProData::GetNoProsFromTo",6);
    const mystring preamb = "PMBatchProData::GetNoProsFromTo: ";
#ifdef __DEBUG__
    if( !pmdatbeg || !pmdatend )
        throw MYRUNTIME_ERROR( preamb + "Null arguments.");
    if( pmdatend[pps2DLen] < pmdatbeg[pps2DLen] )
        throw MYRUNTIME_ERROR( preamb + "Invalid argument addresses.");
#endif
    return (size_t)(pmdatend[pps2DLen]-pmdatbeg[pps2DLen]) / SZINTYPE;
}

// -------------------------------------------------------------------------
// GetNoPositsFromTo: get total number of positions calculated over all 
// profiles present in the buffers;
// pmdatbeg and pmdatend are the beginning and terminal positions of buffers
//
inline
size_t PMBatchProData::GetNoPositsFromTo( char** pmdatbeg, char** pmdatend )
{
    MYMSG("PMBatchProData::GetNoPositsFromTo",6);
    const mystring preamb = "PMBatchProData::GetNoPositsFromTo: ";
#ifdef __DEBUG__
    if( !pmdatbeg || !pmdatend )
        throw MYRUNTIME_ERROR( preamb + "Null arguments.");
    if( pmdatend[pmv2DTrgFrqs] < pmdatbeg[pmv2DTrgFrqs] )
        throw MYRUNTIME_ERROR( preamb + "Invalid argument addresses.");
#endif
    return (size_t)(pmdatend[pmv2DTrgFrqs+0]-pmdatbeg[pmv2DTrgFrqs+0])/SZFPTYPE;
}

// -------------------------------------------------------------------------
// PMDataPresentFromTo: get a flag of whether profile model data is preent 
// counting from the given address to the given terminal address;
// pmdatbeg and pmdatend are the beginning and terminal positions of buffers
//
inline
bool PMBatchProData::PMDataPresentFromTo( char** pmdatbeg, char** pmdatend )
{
    if( !pmdatbeg || !pmdatend)
        return false;
    for( int i = 0; i < pmv2DNoElems; i++ )
        if( pmdatend[pps2DBkgPrbs+i] <= pmdatbeg[pps2DBkgPrbs+i] )
            return false;
    if( pmdatend[pps2DENO] <= pmdatbeg[pps2DENO] )
        return false;
    if( pmdatend[pps2DLen] <= pmdatbeg[pps2DLen] )
        return false;
    if( pmdatend[pps2DDist] <= pmdatbeg[pps2DDist] )
        return false;

    //it's enough to check only the above addresses
//     for( int i = 0; i < ptr2DNoElems; i++ ) {
//         if( pmdatend[ptr2DTrnPrbs+i] <= pmdatbeg[ptr2DTrnPrbs+i] )
//             return false;
// //         if( pmdatend[ptr2DTrnPrbsExp+i] <= pmdatbeg[ptr2DTrnPrbsExp+i] )
// //             return false;
//     }
// 
//     for( int i = 0; i < pmv2DNoElems; i++ )
//         if( pmdatend[pmv2DTrgFrqs+i] <= pmdatbeg[pmv2DTrgFrqs+i] )
//             return false;
//     if( pmdatend[pmv2DHDP1prb] <= pmdatbeg[pmv2DHDP1prb] )
//         return false;
//     if( pmdatend[pmv2DHDP1ind] <= pmdatbeg[pmv2DHDP1ind] )
//         return false;
//     if( pmdatend[pmv2DAddrPro] <= pmdatbeg[pmv2DAddrPro] )
//         return false;
//     if( pmdatend[pmv2Daa] <= pmdatbeg[pmv2Daa] )
//         return false;
    return true;
}

// -------------------------------------------------------------------------
// GetPMDataLen1At: get the length of one profile saved at the given 
// position
inline
size_t PMBatchProData::GetPMDataLen1At( char** pmdata )
{
    MYMSG("PMBatchProData::GetPMDataLen1At",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMDataLen1At: Null argument.");
#endif
    size_t prolen = (size_t)*(INTYPE*)pmdata[pps2DLen];
    return prolen;
}

// -------------------------------------------------------------------------
// GetPMDataENO1At: get ENO of one profile saved at the given position
inline
float PMBatchProData::GetPMDataENO1At( char** pmdata )
{
    MYMSG("PMBatchProData::GetPMDataENO1At",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMDataENO1At: Null argument.");
#endif
    float proeno = (float)*(FPTYPE*)pmdata[pps2DENO];
    return proeno;
}

// -------------------------------------------------------------------------
// GetPMDataSize1At: get the size of complete profile model data of one 
// profile saved at the given position
inline
size_t PMBatchProData::GetPMDataSize1At( char** pmdata )
{
    MYMSG("PMBatchProData::GetPMDataSize1At",6);
#ifdef __DEBUG__
    if( !pmdata )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMDataSize1At: Null argument.");
#endif
    int i;
    size_t szpmdatafilled[pmv2DTotFlds];
    size_t szfilled = 0;//total occupied size over all buffers for one profile
    size_t prolen = (size_t)*(INTYPE*)pmdata[pps2DLen];
    GetPMDataSize1( prolen, szpmdatafilled );
    for( i = 0; i < pmv2DTotFlds; i++ ) {
        szfilled += szpmdatafilled[i];
    }
    return szfilled;
}

// -------------------------------------------------------------------------
// GetPMDataSize1: get the size of complete profile model data of one 
// profile
inline
void PMBatchProData::GetPMDataSize1( size_t prolen, size_t* szpmdata1 )// const
{
    MYMSG("PMBatchProData::GetPMDataSize1",6);
#ifdef __DEBUG__
    if( !szpmdata1 )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMDataSize1: Memory access error.");
#endif
    int i;
    memset( szpmdata1, 0, pmv2DTotFlds * sizeof(size_t));
    for( i = 0; i < pmv2DNoElems; i++ )
        szpmdata1[pps2DBkgPrbs+i] = SZFPTYPE;
    szpmdata1[pps2DENO] = SZFPTYPE;
    szpmdata1[pps2DLen] = SZINTYPE;
    szpmdata1[pps2DDist] = SZLNTYPE;

    for( i = 0; i < ptr2DNoElems; i++ ) {
        szpmdata1[ptr2DTrnPrbs+i] = SZFPTYPE * (prolen+1);
//         szpmdata1[ptr2DTrnPrbsExp+i] = SZFPTYPE * (prolen+1);
    }

    for( i = 0; i < pmv2DNoElems; i++ )
        szpmdata1[pmv2DTrgFrqs+i] = SZFPTYPE * prolen;
    for( i = 0; i < pmv2DNoCVEls; i++ )
        szpmdata1[pmv2DCVentrs+i] = SZFPTYPE * prolen;
    szpmdata1[pmv2DCVprior] = SZFPTYPE * prolen;
    szpmdata1[pmv2DCVnorm2] = SZFPTYPE * prolen;
    for( i = 0; i < pmv2DNoSSSps; i++ )
        szpmdata1[pmv2DSSsprbs+i] = SZFPTYPE * prolen;
    szpmdata1[pmv2DHDP1prb] = SZFPTYPE * prolen;
    szpmdata1[pmv2DHDP1ind] = SZINTYPE * prolen;
    szpmdata1[pmv2DAddrPro] = SZINTYPE * prolen;
    szpmdata1[pmv2DAddrCV] = SZLNTYPE * prolen;
    szpmdata1[pmv2DAddrSS] = SZLNTYPE * prolen;
    szpmdata1[pmv2Daa] = SZCHTYPE * prolen;
    szpmdata1[pmv2DSSstate] = SZCHTYPE * prolen;
}

// -------------------------------------------------------------------------
// GetPMDataSizeUB: get the upper bound of size of complete profile model 
// data given total number of profile positions;
// totprolen, total number of profile positions
inline
void PMBatchProData::GetPMDataSizeUB( size_t totprolen, size_t* szpmdata1 )
{
    MYMSG("PMBatchProData::GetPMDataSizeUB",6);
#ifdef __DEBUG__
    if( !szpmdata1 )
        throw MYRUNTIME_ERROR("PMBatchProData::GetPMDataSizeUB: Memory access error.");
#endif
    int i;
    memset( szpmdata1, 0, pmv2DTotFlds * sizeof(size_t));
    for( i = 0; i < pmv2DNoElems; i++ )
        szpmdata1[pps2DBkgPrbs+i] = SZFPTYPE * totprolen;
    szpmdata1[pps2DENO] = SZFPTYPE * totprolen;
    szpmdata1[pps2DLen] = SZINTYPE * totprolen;
    szpmdata1[pps2DDist] = SZLNTYPE * totprolen;

    for( i = 0; i < ptr2DNoElems; i++ ) {
        szpmdata1[ptr2DTrnPrbs+i] = SZFPTYPE * (totprolen+totprolen);
//         szpmdata1[ptr2DTrnPrbsExp+i] = SZFPTYPE * (totprolen+totprolen);
    }

    for( i = 0; i < pmv2DNoElems; i++ )
        szpmdata1[pmv2DTrgFrqs+i] = SZFPTYPE * totprolen;
    for( i = 0; i < pmv2DNoCVEls; i++ )
        szpmdata1[pmv2DCVentrs+i] = SZFPTYPE * totprolen;
    szpmdata1[pmv2DCVprior] = SZFPTYPE * totprolen;
    szpmdata1[pmv2DCVnorm2] = SZFPTYPE * totprolen;
    for( i = 0; i < pmv2DNoSSSps; i++ )
        szpmdata1[pmv2DSSsprbs+i] = SZFPTYPE * totprolen;
    szpmdata1[pmv2DHDP1prb] = SZFPTYPE * totprolen;
    szpmdata1[pmv2DHDP1ind] = SZINTYPE * totprolen;
    szpmdata1[pmv2DAddrPro] = SZINTYPE * totprolen;
    szpmdata1[pmv2DAddrCV] = SZLNTYPE * totprolen;
    szpmdata1[pmv2DAddrSS] = SZLNTYPE * totprolen;
    szpmdata1[pmv2Daa] = SZCHTYPE * totprolen;
    szpmdata1[pmv2DSSstate] = SZCHTYPE * totprolen;
}

// -------------------------------------------------------------------------
// GetPMDataSizeUBTotal: get the upper bound of size of complete profile 
// model data given total number of profile positions;
// totprolen, total number of profile positions;
// return the total size in bytes
inline
size_t PMBatchProData::GetPMDataSizeUBTotal( size_t totprolen )
{
    size_t szpmdata1[pmv2DTotFlds];
    size_t sztotal = 0;
    GetPMDataSizeUB( totprolen, szpmdata1 );
    for( int i = 0; i < pmv2DTotFlds; i++ )
        sztotal += szpmdata1[i];
    return sztotal;
}

// -------------------------------------------------------------------------
// GetPMDataSizeLBTotal: get the lower bound of size of complete profile 
// model data given total number of profile positions;
// totprolen, total number of profile positions;
// return the total size in bytes
inline
size_t PMBatchProData::GetPMDataSizeLBTotal( size_t totprolen )
{
    size_t szpmdata1[pmv2DTotFlds];
    size_t sztotal = 0;
    GetPMDataSizeLB( totprolen, szpmdata1 );
    for( int i = 0; i < pmv2DTotFlds; i++ )
        sztotal += szpmdata1[i];
    return sztotal;
}


// -------------------------------------------------------------------------
// CopyPMDataHeader: copy profile model header from the given addresses in 
// src and move dest pointers to the next profile
inline
void PMBatchProData::CopyPMDataHeaderNext( mystring& desc, mystring& file, char** src )
{
    MYMSG("PMBatchProData::CopyPMDataHeader",5);
#ifdef __DEBUG__
    if( !src )
        throw MYRUNTIME_ERROR("PMBatchProData::CopyPMDataHeader: Null argument.");
#endif
    int i;

#ifdef __DEBUG__
    if( !h_descs_ || !h_files_ )
        throw MYRUNTIME_ERROR("PMBatchProData::CopyPMDataHeader: Memory access error.");
    if( h_descsend_ < h_descs_ || h_filesend_ < h_files_ )
        throw MYRUNTIME_ERROR("PMBatchProData::CopyPMDataHeader: Invalid running end addresses.");
#endif
//     *h_descsend_++ = desc;
//     *h_filesend_++ = file;
    h_descsend_->move(desc); h_descsend_++;
    h_filesend_->move(file); h_filesend_++;

    for( i = 0; i < pmv2DNoElems; i++ ) {
#ifdef __DEBUG__
        if( h_pmdatend_[pps2DBkgPrbs+i] < h_pmdata_[pps2DBkgPrbs+i])
            throw MYRUNTIME_ERROR("PMBatchProData::CopyPMDataHeader: Invalid running end addresses.");
#endif
        *(FPTYPE*)h_pmdatend_[pps2DBkgPrbs+i] = *(FPTYPE*)src[pps2DBkgPrbs+i];
        h_pmdatend_[pps2DBkgPrbs+i] += SZFPTYPE;
    }
#ifdef __DEBUG__
    if( h_pmdatend_[pps2DENO] < h_pmdata_[pps2DENO] ||
        h_pmdatend_[pps2DLen] < h_pmdata_[pps2DLen] ||
        h_pmdatend_[pps2DDist] < h_pmdata_[pps2DDist] )
        throw MYRUNTIME_ERROR("PMBatchProData::CopyPMDataHeader: Invalid running end addresses.");
#endif
    *(FPTYPE*)h_pmdatend_[pps2DENO] = *(FPTYPE*)src[pps2DENO];
    h_pmdatend_[pps2DENO] += SZFPTYPE;

    *(INTYPE*)h_pmdatend_[pps2DLen] = *(INTYPE*)src[pps2DLen];
    h_pmdatend_[pps2DLen] += SZINTYPE;

    //NOTE:this field is updated while reading profile position-specific data
    //*(LNTYPE*)h_pmdatend_[pps2DDist] = *(LNTYPE*)src[pps2DDist];
    //h_pmdatend_[pps2DDist] += SZLNTYPE;
}

#endif//__PMBatchProData_h__
