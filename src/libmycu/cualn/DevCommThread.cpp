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
#include "libmycu/cupro/CuDeviceMemory.cuh"
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
    CuDeviceMemory* dmem,
    int areano,
    Configuration* config,
    AlnWriter* writer,
    const mystring* queryfnames,
    const mystring* querydescs,
    const char** bdb1descs,
    size_t prodbsize,
    size_t ndbentries )
:   
    mytid_(tid),
    myareano_(areano),
    tobj_(NULL),
    req_msg_(THREAD_MSG_UNSET),
    msg_addressee_(THREAD_MSG_ADDRESSEE_NONE),
    rsp_msg_(THREAD_MSG_UNSET),
    //
    chunkdatasize_(0UL),
    chunkdatalen_(0UL),
    chunknpros_(0UL),
    //
    dmem_(dmem),
    alnwriter_(writer),
    config_(config),
    cached_queryfnames_(queryfnames),
    cached_querydescs_(querydescs),
    cached_bdb1descs_(bdb1descs),
    mstr_set_prodbsize_(prodbsize),
    mstr_set_ndbentries_(ndbentries),
    //query-specific arguments:
    nqyposs_(0UL),
    mstr_set_nqyposs_(0UL),
    mstr_set_scorethld_(0.0f),
    mstr_set_logevthld_(0.0f),
    //data arguments:
    mstr_set_lastchunk_(false),
    mstr_set_chunkno_(-1),
    mstr_set_nqueries_(0UL),
    mstr_set_qrysernr_(-1),
    mstr_set_qrystep_(0UL),
    mstr_set_scorethlds_(NULL),
    mstr_set_logevthlds_(NULL),
    mstr_set_bdbCdesc_(NULL),
    mstr_set_cnt_(NULL),
    scorethld_(0.0f),
    logevthld_(0.0f),
    lastchunk_(false),
    chunkno_(-1),
    nqueries_(0UL),
    qrysernr_(-1),
    qrystep_(0UL),
    scorethlds_(NULL),
    logevthlds_(NULL),
    bdbCdesc_(NULL),
    cnt_(NULL)
{
    MYMSG( "DevCommThread::DevCommThread", 3 );
    if( dmem_ == NULL )
        throw MYRUNTIME_ERROR("DevCommThread::DevCommThread: Null device memory object.");
    memset( mstr_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( mstr_set_szCpm2dvfields_, 0, pmv2DTotFlds * sizeof(size_t));
    memset( querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( szCpm2dvfields_, 0, pmv2DTotFlds * sizeof(size_t));
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
            dmem_,
            myareano_,
            config_,
            alnwriter_,
            cached_queryfnames_,
            cached_querydescs_,
            cached_bdb1descs_
        );

        cbpc.SetDbDetails( mstr_set_prodbsize_, mstr_set_ndbentries_ );

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
                        GetArgsOnMsgGetDataChunkSize(cbpc);
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
                        rspmsg = ttrespmsgChunkSizeReady;
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
        ResetMasterData();
        return;
    }
}

// -------------------------------------------------------------------------
// CalculateMaxDbDataChunkSize: calculate maximum database data chunk size 
// given query length previously set by the master thread
//
void DevCommThread::CalculateMaxDbDataChunkSize( CuBatchProcessing& /*cbpc*/ )
{
    throw MYRUNTIME_ERROR("DevCommThread::CalculateMaxDbDataChunkSize: Currently should not be called!");
    //safely write member variables as long as the master waiting for them is blocked
//     size_t chunkdatasize = cbpc.CalcMaxDbDataChunkSize( nqyposs_ );
//     size_t chunkdatalen = cbpc.GetCurrentMaxDbPos();
//     size_t chunknpros = cbpc.GetCurrentMaxNDbPros();
//     SetChunkDataAttributes( chunkdatasize, chunkdatalen, chunknpros);
}

// -------------------------------------------------------------------------
// ProcessScoreMatrix: batch processing of part of score matrix
//
void DevCommThread::ProcessScoreMatrix( CuBatchProcessing& cbpc )
{
    //NOTE: lock so that other threads (if any) assigned to the same 
    // device wait for the transfer to complete
    {   std::lock_guard<std::mutex> lck(dmem_->GetDeviceProp().shdcnt_->get_mutex());
        int cntval = dmem_->GetDeviceProp().shdcnt_->get_under_lock();
        if( cntval != chunkno_ ) {
            dmem_->GetDeviceProp().shdcnt_->inc_under_lock();
            if( bdbCpmbeg_[0] && bdbCpmend_[0] /*&& szCpm2dvfields_[1]*/)
                cbpc.TransferCPMData(bdbCpmbeg_, bdbCpmend_, szCpm2dvfields_);
        }
    }

    char msgbuf[BUF_MAX];
    char* qrypmbeg[pmv2DTotFlds];
    int qrypercprcd = 0;

    memcpy( qrypmbeg, querypmbeg_, pmv2DTotFlds * sizeof(void*));
    PMBatchProData::PMDataSkipNPros( qrypmbeg, qrysernr_ );

    for(int q = qrysernr_, n = 0;
        q < nqueries_; 
        q += qrystep_, n++,
        PMBatchProData::PMDataSkipNPros( qrypmbeg, qrystep_ ))
    {
        MYMSGBEGl(3)
            sprintf(msgbuf, "Processing QUERY %d (%s)",q,cached_queryfnames_[q].c_str());
            MYMSG(msgbuf,3);
        MYMSGENDl
        //
        size_t qrylen = PMBatchProData::GetPMDataLen1At(qrypmbeg);
        //
        //increase # chunks being processed only if this is not the last chunk!
        //this ensures the processing of all chunks and informs the writer to begin 
        // writing once the last chunk has finished
        if( !lastchunk_ )
            alnwriter_->IncreaseQueryNParts(q);
        //
        cbpc.SetQueryLen(qrylen);
        cbpc.SetScoreThreshold(scorethlds_[q]);
        cbpc.SetLogEThreshold(logevthlds_[q]);
        //
        cbpc.ProcessScoreMatrix( q,
            qrylen,
            qrypmbeg,// querypmend_, 
            bdb1pmbeg_[0]? bdb1pmbeg_: NULL, bdb1pmend_[0]? bdb1pmend_: NULL,
            bdbCdesc_,
            bdbCpmbeg_[0]? bdbCpmbeg_: NULL, bdbCpmend_[0]? bdbCpmend_: NULL,
            szCpm2dvfields_,
            cnt_
        );
        //
        MYMSGBEGl(1)
            if( qrypercprcd !=(n+1)*10/nqueries_) {
                qrypercprcd = (n+1)*10/nqueries_;
                sprintf(msgbuf, "Worker %d processed %d%% of queries",mytid_,qrypercprcd*10);
                MYMSG(msgbuf,1);
            }
        MYMSGENDl
    }
}

// -------------------------------------------------------------------------
// ProcessScoreMatrix_obs: batch processing of part of score matrix
//
void DevCommThread::ProcessScoreMatrix_obs( CuBatchProcessing& cbpc )
{
    cbpc.SetQueryLen(nqyposs_);
    cbpc.SetScoreThreshold(scorethld_);
    cbpc.SetLogEThreshold(logevthld_);

    //NOTE: lock so that other threads (if any) assigned to the same 
    // device wait for the transfer to complete
    {   std::lock_guard<std::mutex> lck(dmem_->GetDeviceProp().shdcnt_->get_mutex());
        int cntval = dmem_->GetDeviceProp().shdcnt_->get_under_lock();
        if( cntval != chunkno_ ) {
            dmem_->GetDeviceProp().shdcnt_->inc_under_lock();
            if( bdbCpmbeg_[0] && bdbCpmend_[0] /*&& szCpm2dvfields_[1]*/)
                cbpc.TransferCPMData(bdbCpmbeg_, bdbCpmend_, szCpm2dvfields_);
        }
    }

    size_t qrylen = PMBatchProData::GetPMDataLen1At(querypmbeg_);

    cbpc.ProcessScoreMatrix(
        qrysernr_,
        qrylen,
        querypmbeg_,// querypmend_, 
        bdb1pmbeg_[0]? bdb1pmbeg_: NULL, bdb1pmend_[0]? bdb1pmend_: NULL,
        bdbCdesc_,
        bdbCpmbeg_[0]? bdbCpmbeg_: NULL, bdbCpmend_[0]? bdbCpmend_: NULL,
        szCpm2dvfields_,
        cnt_
    );
}
