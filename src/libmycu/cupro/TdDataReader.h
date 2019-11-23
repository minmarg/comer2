/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TdDataReader_h__
#define __TdDataReader_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <string.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <cuda_runtime_api.h>

#include "tsafety/TSCounterVar.h"

#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/CuDbReader.h"

//defined implies pinned memory allocation for data read
//#define PINNEDMEMORYDATA

#define TREADER_MSG_UNSET -1
#define TREADER_MSG_ERROR -2

struct DRDataDeleter {
    void operator()(void* p) const {
        if(p)
            std::free(p);
    };
};
struct DRHostDataDeleter {
    void operator()(void* p) const {
        if(p) {
            if(CLOptions::GetIO_UNPINNED())
                std::free(p);
            else
                cudaFreeHost(p);
        }
    };
};

// -------------------------------------------------------------------------
//
struct SbdbCData {
    SbdbCData()
    :   bdbCdata_(nullptr),
        bdbCdescs_(nullptr),
        bdbCprodescs_(nullptr),
        szbdbCdata_(0),
        szbdbCdescs_(0),
        nbdbCprodescs_(0)
    {
        memset( bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
        memset( bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
        memset( szpm2dvfields_, 0, pmv2DTotFlds * sizeof(size_t));
    }
    ~SbdbCData() {}
    void AllocateSpace(size_t chunkdatasize, size_t chunknpros);
    void AllocateSpaceForData(size_t requestedsize);
    bool AllocateSpaceForDescriptions(size_t requestedsize);
    void AllocateSpaceForDescPtrs(size_t requestedno);
    //
    std::unique_ptr<char,DRHostDataDeleter> bdbCdata_;
    char* bdbCpmbeg_[pmv2DTotFlds];//addresses of the beginnings of the fields
    char* bdbCpmend_[pmv2DTotFlds];//addresses of the endings of the fields
    size_t szpm2dvfields_[pmv2DTotFlds+1];//beginnings in bytes (sizes) of the fields written in bdbCdata_
    std::unique_ptr<char,DRDataDeleter> bdbCdescs_;
    std::unique_ptr<char*[]> bdbCprodescs_;//pointers to profile descriptions in bdbCdescs_
    size_t szbdbCdata_;//size allocated for bdbCdata_
    size_t szbdbCdescs_;//size allocated for bdbCdescs_
    size_t nbdbCprodescs_;//number of slots allocated for bdbCprodescs_
    TSCounterVar cnt_;//thread-safe counter of how many agents access the data
};

// _________________________________________________________________________
// Class TdDataReader
//
// thread class responsible fro reading profile data
//
class TdDataReader
{
    enum TDRLimits {
        tdrMAX_N_ENTRIES_READ = 67108864,//2^26
        tdrDEFDESCLEN = 1024//default description length
        //tdrSZREADARRAY = 4//size of array to contain read data
    };

public:
    enum TDRMsg {
        tdrmsgGetSize,
        tdrmsgGetData,
        tdrmsgTerminate
    };
    enum TDRResponseMsg {
        tdrrespmsgSize,
        tdrrespmsgDataReady,
        tdrrespmsgNoData,
        tdrrespmsgTerminating
    };

public:
    TdDataReader(
        const char* dbname,
        bool mapped,
        int nagents,
        Configuration*
    );
    TdDataReader();
    ~TdDataReader();

    std::mutex& GetPrivateMutex() {return mx_dataccess_;}

    //{{NOTE: messaging functions accessed from outside!
    void Notify(int msg) {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_dataccess_);
            req_msg_ = msg;
        }
        cv_msg_.notify_one();
    }
    int Wait(int rsp1, int rsp2 = TREADER_MSG_ERROR) {
        //wait until a response arrives
        std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
        cv_msg_.wait(lck_msg,
            [this,rsp1,rsp2]{
                return (rsp_msg_ == rsp1 || rsp_msg_ == rsp2 || 
                        rsp_msg_ == TREADER_MSG_ERROR);
            }
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        if( rsp_msg_!= TREADER_MSG_ERROR )
            rsp_msg_ = TREADER_MSG_UNSET;
        return rspmsg;
    }
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return rsp_msg_;
    }
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if( rsp_msg_!= TREADER_MSG_ERROR )
            rsp_msg_ = TREADER_MSG_UNSET;
    }
    //}}

    void GetbdbCdata(
        char**& bdbCprodescs,
        char**& bdbCpmbeg, char**& bdbCpmend,
        size_t*& szpm2dvfields,
        TSCounterVar*& tscnt,
        bool* lastchunk )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        bdbCprodescs = bdbCprodescs_;
        bdbCpmbeg = bdbCpmbeg_;
        bdbCpmend = bdbCpmend_;
        szpm2dvfields = szpm2dvfields_;
        tscnt = ccnt_;
        *lastchunk = lastchunk_;

//         if( bdbCpmbeg && bdbCpmend ) {
//             memcpy( bdbCpmbeg, bdbCpmbeg_, pmv2DTotFlds * sizeof(void*));
//             memcpy( bdbCpmend, bdbCpmend_, pmv2DTotFlds * sizeof(void*));
//         }
//         if( szpm2dvfields )
//             memcpy( szpm2dvfields, szpm2dvfields_, pmv2DTotFlds * sizeof(size_t));
    }

    void GetDbSize( size_t* nentries, size_t* dbsize ) const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        *nentries = ndbentries_;
        *dbsize = prodbsize_;
    }

    void SetChunkDataAttributes( 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        chunkdatasize_ = chunkdatasize;
        chunkdatalen_ = chunkdatalen;
        chunknpros_ = chunknpros;
    }

protected:
    void Execute( void* args );

    void OpenAndReadPreamble();
    bool ReadLengthsAndDescEndAddrs();
    bool GetData(size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros);
    void GetNextData(size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros);
    bool ReadDataChunk(SbdbCData*, size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros);
    void ReadDataChunkHelper( SbdbCData*, size_t lengthsndx, size_t npros, size_t totlen);

    void SetDbSize( size_t nentries, size_t dbsize ) {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        ndbentries_ = nentries;
        prodbsize_ = dbsize;
    }

    void ResetProfileCounters() {
        bdbCdata_from_ = (0UL);
        bdbCdata_to_ = (0UL);
        bdbCdata_poss_from_ = (0UL);
        bdbCdata_poss_to_ = (0UL);
    }

private:
    //static constants
    static const size_t s_szdatalign_;//data alignment size wrt cache line size
private:
    //thread section
    std::thread* tobj_;//thread object
private:
    //{{messaging
    std::condition_variable cv_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing data
    int req_msg_;//request message issued for thread
    int rsp_msg_;//private response message
    //}}
    //{{variables that determine the size of data chunks being read
    size_t chunkdatasize_;
    size_t chunkdatalen_;
    size_t chunknpros_;
    //}}
    //{{db name and configuration:
    CuDbReader dbobj_;
    Configuration* config_;
    //}}
    //{{database attributes:
    size_t prodbsize_;//profile database size in positions
    size_t ndbentries_;//number of database entries
    size_t addrtable_[pmv2DTotFlds];//table of section adresses in file
    size_t addrdesc_end_addrs_;//address of end addresses of descriptions in file
    size_t addrdescs_;//address of descriptions in file
    //}}
    //{{data:
    std::unique_ptr<int,DRDataDeleter> bdbClengths_;//profile lengths read
    size_t bdbClengths_from_, bdbClengths_to_;//number of lengths read (in the number of profiles)
    std::unique_ptr<size_t,DRDataDeleter> bdbCdesc_end_addrs_;//description end addresses read (numbers as above)
    //
    std::vector<SbdbCData> bdbCstruct_;//read data
    int pbdbCstruct_ndx_;//index of the currently read data block
    //
    size_t bdbCdata_from_, bdbCdata_to_;//profile data read in the number of profiles
    size_t bdbCdata_poss_from_, bdbCdata_poss_to_;//profile data read in the number of positions
    size_t addrdescproced_;//address of the last description processed
    bool endofdata_;//no data read on the last call
    //}}
    //{{addresses of read data to return
    char** bdbCprodescs_;//profile descriptions
    char** bdbCpmbeg_;//addresses of the beginnings of the fields
    char** bdbCpmend_;//addresses of the endings of the fields
    size_t* szpm2dvfields_;//beginnings in bytes (sizes) of the fields
    TSCounterVar* ccnt_;//counter associated with data to be return
    bool lastchunk_;//flasg of the last chunk to be returned to the master
    int nagents_;//number of agents accessing data
    //}}
    friend struct SbdbCData;
};

#endif//__TdDataReader_h__
