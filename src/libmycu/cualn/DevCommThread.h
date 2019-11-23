/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __DevCommThread_h__
#define __DevCommThread_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "tsafety/TSCounterVar.h"

#include "libHDP/HDPbase.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cupro/CuDeviceMemory.cuh"
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
        ttrespmsgChunkSizeReady,
        ttrespmsgInProgress,
        ttrespmsgProbed,
        ttrespmsgTerminating
    };

public:
    DevCommThread(
        int tid,
        CuDeviceMemory* dmem,
        int areano,
        Configuration*,
        AlnWriter*,
        const mystring* queryfnames,
        const mystring* querydescs,
        const char** bdb1descs,
        size_t prodbsize,
        size_t ndbentries 
    );
    DevCommThread();
    ~DevCommThread();

    std::mutex& GetPrivateMutex() {return mx_rsp_msg_;}
    std::condition_variable& GetPrivateCV() {return cv_rsp_msg_;}


    //{{NOTE: messaging functions accessed from outside!
    void Notify(int msg, int adr) {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_rsp_msg_);
            req_msg_ = msg;
            msg_addressee_ = adr;
        }
        cv_rsp_msg_.notify_one();
    }
    int waitForDataAccess() {
        int rsp = GetResponseAsync();
        if( rsp == THREAD_MSG_ERROR )
            return rsp;
        {   std::unique_lock<std::mutex> lck(mx_dataccess_);
            cv_dataccess_.wait(lck, [this]{return GetMstrDataEmpty();});
            //mutex release
        }
        return GetResponseAsync();
    }
    int Wait(int rsp) {
        //wait until a response arrives
        std::unique_lock<std::mutex> lck_msg(mx_rsp_msg_);
        cv_rsp_msg_.wait(lck_msg,
            [this,rsp]{return (rsp_msg_ == rsp || rsp_msg_ == THREAD_MSG_ERROR);}
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        if( rsp_msg_!= THREAD_MSG_ERROR )
            rsp_msg_ = THREAD_MSG_UNSET;
        return rspmsg;
    }
    bool IsIdle(bool* mstr_set_data_empty, bool wait = false) {
        std::unique_lock<std::mutex> lck_busy(mx_rsp_msg_, std::defer_lock);
        if( mstr_set_data_empty )
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            *mstr_set_data_empty = GetMstrDataEmpty();
        }
        //NOTE: when recycling the loop, lock may be acquired, 
        // although the thread is not idle;
        // data emptiness should be checked jointly!
        bool lck_acquired = true;
        if(wait)
            lck_busy.lock();
        else
            lck_acquired = lck_busy.try_lock();
        return lck_acquired;
        //release of the mutex if locked
    }
    int GetResponseAsync() const {
        //get a response if available
        std::unique_lock<std::mutex> lck_busy(mx_rsp_msg_, std::defer_lock);
        int rsp = THREAD_MSG_UNSET;
        if( lck_busy.try_lock())
            rsp = rsp_msg_;
        return rsp;
    }
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_rsp_msg_);
        return rsp_msg_;
    }
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_rsp_msg_);
        if( rsp_msg_!= THREAD_MSG_ERROR )
            rsp_msg_ = THREAD_MSG_UNSET;
    }
    //}}


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


    void SetQueryLen( size_t nqyposs )
    {
        //safely read locked data
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_nqyposs_ = nqyposs;
    }


    void SetMstrQueryBDbdata(
        bool lastchunk,
        int chunkno,
        int nqueries, int qrysernr, int qrystep,
        float* scorethlds, float* logevthlds,
        char** querypmbeg,
        char** bdb1pmbeg, char** bdb1pmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend, size_t* szpm2dvfields,
        TSCounterVar* tscnt )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_lastchunk_ = lastchunk;
        mstr_set_chunkno_ = chunkno;
        mstr_set_nqueries_ = nqueries;
        mstr_set_qrysernr_ = qrysernr;
        mstr_set_qrystep_ = qrystep;
        mstr_set_scorethlds_ = scorethlds;
        mstr_set_logevthlds_ = logevthlds;
        //
        memcpy( mstr_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
        mstr_set_bdbCdesc_ = bdbCdesc;
        mstr_set_cnt_ = tscnt;
        if( bdb1pmbeg && bdb1pmend ) {
            memcpy( mstr_set_bdb1pmbeg_, bdb1pmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdb1pmend_, bdb1pmend, pmv2DTotFlds * sizeof(void*));
        }
        if( bdbCpmbeg && bdbCpmend && szpm2dvfields ) {
            memcpy( mstr_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_szCpm2dvfields_, szpm2dvfields, (pmv2DTotFlds+1) * sizeof(size_t));
        }
    }

    void SetMstrQueryBDbdata_obs(
        int chunkno,
        int qrysernr,
        size_t nqyposs, float scorethld, float logevthld,
        char** querypmbeg, char** querypmend,
        char** bdb1pmbeg, char** bdb1pmend,
        const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend, size_t* szpm2dvfields,
        TSCounterVar* tscnt )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        mstr_set_chunkno_ = chunkno;
        mstr_set_qrysernr_ = qrysernr;
        mstr_set_nqyposs_ = nqyposs;
        mstr_set_scorethld_ = scorethld;
        mstr_set_logevthld_ = logevthld;
        //
        memcpy( mstr_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
        memcpy( mstr_set_querypmend_, querypmend, pmv2DTotFlds * sizeof(void*));
        mstr_set_bdbCdesc_ = bdbCdesc;
        mstr_set_cnt_ = tscnt;
        if( bdb1pmbeg && bdb1pmend ) {
            memcpy( mstr_set_bdb1pmbeg_, bdb1pmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdb1pmend_, bdb1pmend, pmv2DTotFlds * sizeof(void*));
        }
        if( bdbCpmbeg && bdbCpmend && szpm2dvfields ) {
            memcpy( mstr_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
            memcpy( mstr_set_szCpm2dvfields_, szpm2dvfields, (pmv2DTotFlds+1) * sizeof(size_t));
        }
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
    void GetArgsOnMsgGetDataChunkSize( CuBatchProcessing& cbpc )
    {
        //safely read locked data
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        nqyposs_ = mstr_set_nqyposs_;
        cbpc.SetQueryLen( nqyposs_ );
//         cbpc.SetScoreThreshold( mstr_set_scorethld_ );
//         cbpc.SetLogEThreshold( mstr_set_logevthld_ );
        mstr_set_nqyposs_ = 0;
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
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            lastchunk_ = mstr_set_lastchunk_;
            chunkno_ = mstr_set_chunkno_;
            nqueries_ = mstr_set_nqueries_;
            qrysernr_ = mstr_set_qrysernr_;
            qrystep_ = mstr_set_qrystep_;
            nqyposs_ = mstr_set_nqyposs_;/**/
            scorethld_ = mstr_set_scorethld_;/**/
            logevthld_ = mstr_set_logevthld_;/**/
            scorethlds_ = mstr_set_scorethlds_;
            logevthlds_ = mstr_set_logevthlds_;
            //
            bdbCdesc_ = mstr_set_bdbCdesc_;
            cnt_ = mstr_set_cnt_;
            memcpy( querypmbeg_, mstr_set_querypmbeg_, pmv2DTotFlds * sizeof(void*));
            memcpy( querypmend_, mstr_set_querypmend_, pmv2DTotFlds * sizeof(void*));/**/
            memcpy( bdb1pmbeg_, mstr_set_bdb1pmbeg_, pmv2DTotFlds * sizeof(void*));
            memcpy( bdb1pmend_, mstr_set_bdb1pmend_, pmv2DTotFlds * sizeof(void*));
            memcpy( bdbCpmbeg_, mstr_set_bdbCpmbeg_, pmv2DTotFlds * sizeof(void*));
            memcpy( bdbCpmend_, mstr_set_bdbCpmend_, pmv2DTotFlds * sizeof(void*));
            memcpy( szCpm2dvfields_, mstr_set_szCpm2dvfields_, (pmv2DTotFlds+1) * sizeof(size_t));
        }
        ResetMasterData();
    }

    void ResetMasterData()
    {
        {   std::lock_guard<std::mutex> lck(mx_dataccess_);
            //reset addresses so that master can write the addresses of next data
            mstr_set_lastchunk_ = false;
            mstr_set_chunkno_ = -1;
            mstr_set_nqueries_ = -1;
            mstr_set_qrysernr_ = -1;
            mstr_set_qrystep_ = 0;
            mstr_set_nqyposs_ = 0;
            mstr_set_scorethld_ = 0.0f;
            mstr_set_logevthld_ = 0.0f;
            mstr_set_scorethlds_ = NULL;
            mstr_set_logevthlds_ = NULL;
            mstr_set_bdbCdesc_ = NULL;
            mstr_set_cnt_ = NULL;
            memset( mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
            memset( mstr_set_szCpm2dvfields_, 0, (pmv2DTotFlds+1) * sizeof(size_t));
        }
        cv_dataccess_.notify_one();
    }

    bool GetMstrDataEmpty() {
        return  mstr_set_querypmbeg_[0] == NULL && mstr_set_querypmend_[0] == NULL &&
                mstr_set_bdb1pmbeg_[0] == NULL && mstr_set_bdb1pmend_[0] == NULL &&
            mstr_set_bdbCdesc_ == NULL && mstr_set_bdbCpmbeg_[0] == NULL && mstr_set_bdbCpmend_[0] == NULL
            /*&& mstr_set_szCpm2dvfields_[1] == 0*/;
    }

    void ProcessScoreMatrix( CuBatchProcessing& );

    void ProcessScoreMatrix_obs( CuBatchProcessing& );

private:
    //thread section
    int mytid_;//private thread id
    int myareano_;//device memory area no.
    std::thread* tobj_;//thread object
    //{{messaging
    // broadcasting attributes:
private:
    int req_msg_;//request message issued for thread
    int msg_addressee_;//addressee of a message seen by thread
    // response attributes:
    mutable std::mutex mx_rsp_msg_;//mutex for messaging between the thread and the master
    std::condition_variable cv_rsp_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing class data
    std::condition_variable cv_dataccess_;//condition variable for checking for data emptyness
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
    CuDeviceMemory* dmem_;
    //results writer and configuration:
    AlnWriter* alnwriter_;
    Configuration* config_;
    //cached data: 
    const mystring* cached_queryfnames_;
    const mystring* cached_querydescs_;
    const char** cached_bdb1descs_;
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
    bool mstr_set_lastchunk_;
    int mstr_set_chunkno_;
    int mstr_set_nqueries_;
    int mstr_set_qrysernr_;
    int mstr_set_qrystep_;
    float* mstr_set_scorethlds_;
    float* mstr_set_logevthlds_;
    char* mstr_set_querypmbeg_[pmv2DTotFlds];
    char* mstr_set_querypmend_[pmv2DTotFlds];
    char* mstr_set_bdb1pmbeg_[pmv2DTotFlds];
    char* mstr_set_bdb1pmend_[pmv2DTotFlds];
    const char** mstr_set_bdbCdesc_;
    char* mstr_set_bdbCpmbeg_[pmv2DTotFlds];
    char* mstr_set_bdbCpmend_[pmv2DTotFlds];
    size_t mstr_set_szCpm2dvfields_[pmv2DTotFlds+1];
    TSCounterVar* mstr_set_cnt_;
    // adresses saved by slave thread
    float scorethld_;//score threshold
    float logevthld_;//e-value threshold
    bool lastchunk_;
    int chunkno_;//data chunk serial number
    int nqueries_;//number of queries
    int qrysernr_;//query serial number
    int qrystep_;//step with which to process queries
    float* scorethlds_;
    float* logevthlds_;
    char* querypmbeg_[pmv2DTotFlds];
    char* querypmend_[pmv2DTotFlds];
    char* bdb1pmbeg_[pmv2DTotFlds];
    char* bdb1pmend_[pmv2DTotFlds];
    const char** bdbCdesc_;
    char* bdbCpmbeg_[pmv2DTotFlds];
    char* bdbCpmend_[pmv2DTotFlds];
    size_t szCpm2dvfields_[pmv2DTotFlds+1];
    TSCounterVar* cnt_;
    //}}
};

#endif//__DevCommThread_h__
