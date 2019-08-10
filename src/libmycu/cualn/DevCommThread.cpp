/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "liblib/mybase.h"
#include "libHDP/HDPbase.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/CuRoDb.h"
#include "libmycu/cupro/VirtualCuRoDb.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "libmycu/cupro/CuBatchProcessing.cuh"

#include "libHDP/HDPscores.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"

#include "DevCommThread.h"

// _________________________________________________________________________
// class statics
//

// _________________________________________________________________________
// Class DevCommThread
//
// Constructor
//
DevCommThread::DevCommThread(
    int tid,
    DeviceProperties dprop,
    AlnWriter* writer,
    Configuration* config,
    const mystring* queryfnames,
    const mystring* querydescs,
    char** querypmbeg,
    char** querypmend,
    const mystring* bdb1fnames,
    const mystring* bdb1descs,
    char** bdb1pmbeg,
    char** bdb1pmend,
    size_t prodbsize,
    size_t ndbentries )
:   mytid_(tid),
    tobj_(NULL),
    req_msg_(THREAD_MSG_UNSET),
    msg_addressee_(THREAD_MSG_ADDRESSEE_NONE),
    rsp_msg_(THREAD_MSG_UNSET),
    //
    chunkdatasize_(0UL),
    chunkdatalen_(0UL),
    chunknpros_(0UL),
    //
    dprop_(dprop),
    alnwriter_(writer),
    config_(config),
    cached_queryfnames_(queryfnames),
    cached_querydescs_(querydescs),
    cached_querypmbeg_(querypmbeg),
    cached_querypmend_(querypmend),
    cached_bdb1fnames_(bdb1fnames),
    cached_bdb1descs_(bdb1descs),
    cached_bdb1pmbeg_(bdb1pmbeg),
    cached_bdb1pmend_(bdb1pmend),
    mstr_set_prodbsize_(prodbsize),
    mstr_set_ndbentries_(ndbentries),
    //query-specific arguments:
    nqyposs_(0UL),
    mstr_set_nqyposs_(0UL),
    mstr_set_scorethld_(0.0f),
    mstr_set_logevthld_(0.0f),
    //data arguments:
    mstr_set_qrysernr_(-1),
    mstr_set_bdbC_(nullptr),
    qrysernr_(-1)
{
    MYMSG( "DevCommThread::DevCommThread", 3 );
    memset( mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    tobj_ = new std::thread( &DevCommThread::Execute, this, (void*)NULL );
}

// Default constructor
//
DevCommThread::DevCommThread()
:   tobj_(NULL),
    config_(NULL)
{
    throw MYRUNTIME_ERROR("DevCommThread::DevCommThread: Default initialization is prohibited.");
}

// Destructor
//
DevCommThread::~DevCommThread()
{
    MYMSG( "DevCommThread::~DevCommThread", 3 );
    if( tobj_ ) {
        tobj_->join();
        delete tobj_;
        tobj_ = NULL;
    }
}

// -------------------------------------------------------------------------
// Execute: thread's starting point and execution process
//
void DevCommThread::Execute( void* )
{
    MYMSG( "DevCommThread::Execute", 3 );
    myruntime_error mre;

    try {
        CuBatchProcessing cbpc( 
            alnwriter_,
            config_,
            dprop_,
            HDPSCORES.GetScores()!=NULL,
            MOptions::GetMINPP()<1.0f,
            dprop_.reqmem_,
            cached_queryfnames_,
            cached_querydescs_,
            cached_bdb1fnames_,
            cached_bdb1descs_
        );

        cbpc.CacheSSENNWeights();
        cbpc.CacheSSSScores( SSSSCORES );
        cbpc.CacheCVS2Scores( CVS2SCORES );
        cbpc.CacheHDPScores( HDPSCORES );

        cbpc.CacheData(
            cached_querypmbeg_, cached_querypmend_,
            cached_bdb1pmbeg_, cached_bdb1pmend_ );

        while(1) {
            //wait until the master bradcasts a message
            std::unique_lock<std::mutex> lck_msg(mx_rsp_msg_);

            cv_rsp_msg_.wait(lck_msg,
                [this]{return 
                    msg_addressee_ == mytid_
                    &&
                    ((0 <= req_msg_ && req_msg_ <= tthreadmsgTerminate) || 
                     req_msg_ == THREAD_MSG_ERROR
                    );}
            );

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "DevCommThread::Execute: Msg %d Adr %d",
                        req_msg_, msg_addressee_);
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //thread owns the lock after the wait;
            //read message req_msg_
            int reqmsg = req_msg_;

            //unset the message to avoid live cycle when starting over the loop
            req_msg_ = THREAD_MSG_UNSET;
            msg_addressee_ = THREAD_MSG_ADDRESSEE_NONE;

            //set response msg to error upon occurance of an exception
            rsp_msg_ = THREAD_MSG_ERROR;
            int rspmsg = rsp_msg_;

            //immediately read the master-set data and...
            switch(reqmsg) {
                case tthreadmsgGetDataChunkSize:
                        //master has set the query length, score and e-value thresholds;
                        CopyDataOnMsgGetDataChunkSize(cbpc);
                        break;
                case tthreadmsgProcessNewData:
                        //master has written data addresses
                        CopyDataOnMsgProcessNewData(cbpc);
                        break;
                case tthreadmsgProbe: break;
                case tthreadmsgTerminate: break;
                default: break;
            };

            switch(reqmsg) {
                case tthreadmsgGetDataChunkSize:
                        //master has set the query length, score and e-value thresholds;
                        CalculateMaxDbDataChunkSize( cbpc );
                        rspmsg = ttrespmsgDataReady;
                        break;
                case tthreadmsgProcessNewData:
                        //master has written data addresses
                        ProcessScoreMatrix( cbpc );
                        //master does not wait for a response nor requires data to read;
                        //unset response code
                        rspmsg = THREAD_MSG_UNSET;//ttrespmsgInProgress;
                        break;
                case tthreadmsgProbe:
                        cbpc.WaitForIdleChilds();
                        rspmsg = ttrespmsgProbed;
                        break;
                case tthreadmsgTerminate:
                        rspmsg = ttrespmsgTerminating;
                        break;
                default:
                        rspmsg = THREAD_MSG_UNSET;
                        break;
            };

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "DevCommThread::Execute: Msg %d Rsp %d",reqmsg, rspmsg );
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //save response code and ...
            rsp_msg_ = rspmsg;

            //send a message back to the master:
            //unlock the mutex to avoid to block the waiting master and 
            //  notify the master using the cv
            lck_msg.unlock();
            cv_rsp_msg_.notify_one();

            if( reqmsg < 0 || reqmsg == tthreadmsgTerminate)
                //terminate execution
                break;
        }

    } catch( myruntime_error const& ex ) {
        mre = ex;
    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Unknown exception caught.");
    }

    if( mre.isset()) {
        error( mre.pretty_format().c_str());
        {//notify the master
            std::lock_guard<std::mutex> lck_msg(mx_rsp_msg_);
            rsp_msg_ = THREAD_MSG_ERROR;
        }
        cv_rsp_msg_.notify_one();
        return;
    }
}

// -------------------------------------------------------------------------
// CalculateMaxDbDataChunkSize: calculate maximum database data chunk size 
// given query length previously set by the master thread
//
void DevCommThread::CalculateMaxDbDataChunkSize( CuBatchProcessing& cbpc )
{
    //safely write member variables as long as the master waiting for them is blocked
    size_t chunkdatasize = cbpc.CalcMaxDbDataChunkSize( nqyposs_ );
    size_t chunkdatalen = cbpc.GetCurrentMaxDbPos();
    size_t chunknpros = cbpc.GetCurrentMaxNDbPros();
    SetChunkDataAttributes( chunkdatasize, chunkdatalen, chunknpros);
}

// -------------------------------------------------------------------------
// ProcessScoreMatrix: batch processing of part of score matrix
//
void DevCommThread::ProcessScoreMatrix( CuBatchProcessing& cbpc )
{
    cbpc.ProcessScoreMatrix(
        qrysernr_,
        std::move(bdbC_),
        querypmbeg_, querypmend_, 
        bdb1pmbeg_[0]? bdb1pmbeg_: NULL, bdb1pmend_[0]? bdb1pmend_: NULL,
        bdbCpmbeg_[0]? bdbCpmbeg_: NULL, bdbCpmend_[0]? bdbCpmend_: NULL
    );
}
