/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __DevCommThread_h__
#define __DevCommThread_h__

#include <stdio.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "liblib/mybase.h"
#include "libHDP/HDPbase.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cupro/CuBatchProcessing.cuh"

#define THREAD_MSG_UNSET -1
#define THREAD_MSG_ERROR -2
#define THREAD_MSG_ADDRESSEE_NONE -1

// _________________________________________________________________________
// Class DevCommThread
//
// thread class for communicating with a device (GPU)
//
class DevCommThread
{
public:
    enum TThreadMsg {
        tthreadmsgGetDataChunkSize,
        tthreadmsgProcessNewData,
        tthreadmsgProbe,
        tthreadmsgTerminate
    };
    enum TThreadResponseMsg {
        ttrespmsgDataReady,
        ttrespmsgInProgress,
        ttrespmsgProbed,
        ttrespmsgTerminating
    };

public:
    DevCommThread(
        int tid,
        DeviceProperties dprop,
        AlnWriter*,
        Configuration*,
        const mystring* queryfnames,
        const mystring* querydescs,
        char** querypmbeg,
        char** querypmend,
        const mystring* bdb1fnames,
        const mystring* bdb1descs,
        char** bdb1pmbeg,
        char** bdb1pmend,
        size_t prodbsize,
        size_t ndbentries 
    );
    DevCommThread();
    ~DevCommThread();

    std::mutex& GetPrivateMutex() {return mx_rsp_msg_;}
    std::condition_variable& GetPrivateCV() {return cv_rsp_msg_;}


    //{{NOTE: these functions should be called, and member variables accessed, 
    // only under locked mutex mx_rsp_msg_!
    void SetBcastMessageAndAddressee( int msg, int adr ) {
        req_msg_ = msg; msg_addressee_ = adr;
    }
    int GetMessage() const {return req_msg_;}
    int GetResponseMsg() const {return rsp_msg_;}
    void ResetResponseMsg() {
        if( rsp_msg_!= THREAD_MSG_ERROR ) 
            rsp_msg_ = THREAD_MSG_UNSET;
    }
    //}}


    void SetQueryAttributes( size_t nqyposs, float scorethld, float logevthld)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_nqyposs_ = nqyposs;
        mstr_set_scorethld_ = scorethld;
        mstr_set_logevthld_ = logevthld;
    }


    void SetMstrQueryBDbdata(
        int qrysernr,
        std::unique_ptr<PMBatchProData> bdbC,
        char** querypmbeg, char** querypmend,
        char** bdb1pmbeg, char** bdb1pmend,
        char** bdbCpmbeg, char** bdbCpmend )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_qrysernr_ = qrysernr;
        mstr_set_bdbC_ = std::move(bdbC);
        memcpy( mstr_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
        memcpy( mstr_set_querypmend_, querypmend, pmv2DTotFlds * sizeof(void*));
        if( bdb1pmbeg && bdb1pmend ) {
            memcpy( mstr_set_bdb1pmbeg_, bdb1pmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdb1pmend_, bdb1pmend, pmv2DTotFlds * sizeof(void*));
        }
        if( bdbCpmbeg && bdbCpmend ) {
            memcpy( mstr_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
        }
    }
    bool GetMstrDataEmpty() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return !mstr_set_bdbC_ &&
                mstr_set_querypmbeg_[0] == NULL && mstr_set_querypmend_[0] == NULL &&
                mstr_set_bdb1pmbeg_[0] == NULL && mstr_set_bdb1pmend_[0] == NULL &&
                mstr_set_bdbCpmbeg_[0] == NULL && mstr_set_bdbCpmend_[0] == NULL;
    }


    void GetChunkDataAttributes( 
        size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        *chunkdatasize = chunkdatasize_;
        *chunkdatalen = chunkdatalen_;
        *chunknpros = chunknpros_;
    }


protected:
    void Execute( void* args );


    // CopyDataOnMsgGetDataChunkSize: copy data set by the master on 
    // acceptance  of message tthreadmsgGetDataChunkSize
    void CopyDataOnMsgGetDataChunkSize( CuBatchProcessing& cbpc )
    {
        //safely read locked data
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        nqyposs_ = mstr_set_nqyposs_;
        cbpc.SetDbDetails( mstr_set_nqyposs_, mstr_set_prodbsize_, mstr_set_ndbentries_ );
        cbpc.SetScoreThreshold( mstr_set_scorethld_ );
        cbpc.SetLogEThreshold( mstr_set_logevthld_ );
    }
    void CalculateMaxDbDataChunkSize( CuBatchProcessing& );
    void SetChunkDataAttributes( 
        size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros)
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        chunkdatasize_ = chunkdatasize;
        chunkdatalen_ = chunkdatalen;
        chunknpros_ = chunknpros;
    }


    // CopyDataOnMsgProcessNewData: copy data set by the master on acceptance of
    // message tthreadmsgProcessNewData
    //
    void CopyDataOnMsgProcessNewData( CuBatchProcessing& )
    {
        //safely read and write addresses written by the master
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        qrysernr_ = mstr_set_qrysernr_;
        bdbC_ = std::move(mstr_set_bdbC_);
        memcpy( querypmbeg_, mstr_set_querypmbeg_, pmv2DTotFlds * sizeof(void*));
        memcpy( querypmend_, mstr_set_querypmend_, pmv2DTotFlds * sizeof(void*));
        memcpy( bdb1pmbeg_, mstr_set_bdb1pmbeg_, pmv2DTotFlds * sizeof(void*));
        memcpy( bdb1pmend_, mstr_set_bdb1pmend_, pmv2DTotFlds * sizeof(void*));
        memcpy( bdbCpmbeg_, mstr_set_bdbCpmbeg_, pmv2DTotFlds * sizeof(void*));
        memcpy( bdbCpmend_, mstr_set_bdbCpmend_, pmv2DTotFlds * sizeof(void*));
        //reset addresses so that master can write the addresses of following data
        mstr_set_qrysernr_ = -1;
        memset( mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
        memset( mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
        memset( mstr_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
        memset( mstr_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
        memset( mstr_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
        memset( mstr_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    }
    void ProcessScoreMatrix( CuBatchProcessing& );

private:
    //thread section
    int mytid_;//provate thread id
    std::thread* tobj_;//thread object
    //{{messaging
    // broadcasting attributes:
private:
    int req_msg_;//request message issued for thread
    int msg_addressee_;//addressee of a message seen by thread
    // response attributes:
    std::mutex mx_rsp_msg_;//mutex for messaging between the thread and the master
    std::condition_variable cv_rsp_msg_;//condition variable for messaging
    std::mutex mx_dataccess_;//mutex for accessing class data
    int rsp_msg_;//private response message to master
    //}}
    //
    //NOTE: these data are reinitialized on msg tthreadmsgGetDataChunkSize and
    // are supposed to be guarded by mutex mx_rsp_msg_
    size_t chunkdatasize_;
    size_t chunkdatalen_;
    size_t chunknpros_;
    //
    //properties of device the thread's communicating with:
    DeviceProperties dprop_;
    //results writer and configuration:
    AlnWriter* alnwriter_;
    Configuration* config_;
    //cached data: 
    const mystring* cached_queryfnames_;
    const mystring* cached_querydescs_;
    char** cached_querypmbeg_;
    char** cached_querypmend_;
    const mystring* cached_bdb1fnames_;
    const mystring* cached_bdb1descs_;
    char** cached_bdb1pmbeg_;
    char** cached_bdb1pmend_;
    //database arguments: 
    size_t mstr_set_prodbsize_;//profile database size in positions
    size_t mstr_set_ndbentries_;//number of database entries
    //query-specific arguments: 
    size_t nqyposs_;//saved query length guarded by mutex mx_msg_!
    size_t mstr_set_nqyposs_;//master-set query length guarded by mutex mx_msg_!
    float mstr_set_scorethld_;
    float mstr_set_logevthld_;
    //{{data arguments: 
    // master-set data/addresses:
    int mstr_set_qrysernr_;
    std::unique_ptr<PMBatchProData> mstr_set_bdbC_;
    char* mstr_set_querypmbeg_[pmv2DTotFlds];
    char* mstr_set_querypmend_[pmv2DTotFlds];
    char* mstr_set_bdb1pmbeg_[pmv2DTotFlds];
    char* mstr_set_bdb1pmend_[pmv2DTotFlds];
    char* mstr_set_bdbCpmbeg_[pmv2DTotFlds];
    char* mstr_set_bdbCpmend_[pmv2DTotFlds];
    // adresses saved by slave thread
    int qrysernr_;//query serial number
    std::unique_ptr<PMBatchProData> bdbC_;
    char* querypmbeg_[pmv2DTotFlds];
    char* querypmend_[pmv2DTotFlds];
    char* bdb1pmbeg_[pmv2DTotFlds];
    char* bdb1pmend_[pmv2DTotFlds];
    char* bdbCpmbeg_[pmv2DTotFlds];
    char* bdbCpmend_[pmv2DTotFlds];
    //}}
};

#endif//__DevCommThread_h__
