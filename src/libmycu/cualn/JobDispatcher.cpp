/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdlib.h>
#include <string.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <vector>

#include "tsafety/TSCounterVar.h"

#include "libHDP/HDPbase.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
// #include "libmycu/cupro/CuRoDb.h"
// #include "libmycu/cupro/VirtualCuRoDb.h"
#include "libmycu/cupro/CuDbReader.h"
#include "libmycu/cupro/CuInputReader.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "libmycu/cupro/TdDataReader.h"
#include "libmycu/cupro/CuBatchProcessing.cuh"
#include "libmycu/cualn/DevCommThread.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cualn/Devices.h"

#include "libHDP/HDPscores.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"

#include "JobDispatcher.h"

// _________________________________________________________________________
// Class JobDispatcher
//
// Constructor
//
JobDispatcher::JobDispatcher(
    const char* configfile,
    const char* input,
    const char* database,
    const char* output,
    float eval_thld,
    int no_hits,
    int no_alns,
    bool showpars )
:
    configFile_( configfile ),
    input_( input ),
    database_( database ),
    output_( output ),
    writer_(NULL),

//     scoreSystem_( NULL ),

    eval_upper_( eval_thld ),
    max_no_hits_( no_hits ),
    max_no_alns_( no_alns ),
    show_pars_( showpars ),

    reader_(NULL),
    db_no_seqs_( 0 ),
    db_size_( 0ULL ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
//     HDPctbase_( NULL ),

    ssemodel_( 0 ),

    hsplength_( 0 ),
    hspscore_( 0 ),
    hspdistance_( 0 ),
    hspnohsps_( 0 )
{
}

// Default constructor
//
JobDispatcher::JobDispatcher()
:
    configFile_( NULL ),
    input_( NULL ),
    database_( NULL ),
    output_( NULL ),

    eval_upper_( 0.0f ),
    max_no_hits_( 0 ),
    max_no_alns_( 0 ),
    show_pars_( false ),

    reader_( NULL ),
    db_no_seqs_( 0 ),
    db_size_( 0ULL ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
//     HDPctbase_( NULL ),

    ssemodel_( 0 ),

    hsplength_( 0 ),
    hspscore_( 0 ),
    hspdistance_( 0 ),
    hspnohsps_( 0 )
{
    throw MYRUNTIME_ERROR("JobDispatcher::JobDispatcher: Default initialization is prohibited.");
}

// Destructor
//
JobDispatcher::~JobDispatcher()
{
    for(int tid = 0; tid < (int)hostworkers_.size(); tid++ ) {
        if( hostworkers_[tid]) {
            delete hostworkers_[tid];
            hostworkers_[tid] = NULL;
        }
    }
    for(int d = 0; d < (int)memdevs_.size(); d++ ) {
        if( memdevs_[d]) {
            delete memdevs_[d];
            memdevs_[d] = NULL;
        }
    }
    if( writer_ ) {
        writer_->Notify(AlnWriter::wrtthreadmsgTerminate);
        delete writer_;
        writer_ = NULL;
    }
    if( reader_ ) {
        reader_->Notify(TdDataReader::tdrmsgTerminate);
        delete reader_;
        reader_ = NULL;
    }
}

// -------------------------------------------------------------------------
// ReadConfiguration: read configurations
//
void JobDispatcher::ReadConfiguration()
{
    MYMSG( "JobDispatcher::ReadConfiguration", 3 );
    GetConfiguration(CTUngapped).SetFilename( GetConfigFile());
    GetConfiguration(CTGapped).SetFilename( GetConfigFile());

    //read parameters
    static const bool autogapcosts = true;
    static const int autogapopencostunused = 4;
    static const int autogapextncostunused = 1;
    GetConfiguration(CTUngapped).ReadUngapped();//read ungapped configuration
    GetConfiguration(CTGapped).SetAutoGapOpenCost( autogapcosts );
    GetConfiguration(CTGapped).SetGapOpenCost( autogapopencostunused );
    GetConfiguration(CTGapped).SetGapExtendCost( autogapextncostunused );
    GetConfiguration(CTGapped).Read();//read reference parameters
}



// -------------------------------------------------------------------------
// CreateReader: create a thread reading data from file
//
void JobDispatcher::CreateReader( 
    const char* dbname,
    bool mapped,
    int nqueries,
    Configuration* config )
{
    MYMSG( "JobDispatcher::CreateReader", 3 );
    reader_ = new TdDataReader( dbname, mapped, nqueries, config );
    if( reader_ == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::CreateReader: Failed to create the reading thread.");
}

// -------------------------------------------------------------------------
// GetDbSizeFromReader: wait for the reader to read database attributes and 
// get the db size and the number of entries
//
void JobDispatcher::GetDbSizeFromReader()
{
    MYMSG( "JobDispatcher::GetDbSizeFromReader", 3 );
    if( reader_ == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDbSizeFromReader: Database has not been opened.");
    reader_->Notify(TdDataReader::tdrmsgGetSize);
    if( reader_->Wait(TdDataReader::tdrrespmsgSize) != 
        TdDataReader::tdrrespmsgSize ) {
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDbSizeFromReader: Invalid response from the reader.");
    }
    reader_->GetDbSize( &db_no_seqs_, &db_size_ );
}

// -------------------------------------------------------------------------
// GetDbDataFromReader: instrcut the reader to get data and wait for data to
// be ready;
// return false if there are no data to read;
//
bool JobDispatcher::GetDbDataFromReader(
    char**& bdbCprodescs,
    char**& bdbCpmbeg, char**& bdbCpmend,
    size_t*& szpm2dvfields,
    TSCounterVar*& tscnt,
    bool* lastchunk )
{
    MYMSG( "JobDispatcher::GetDbDataFromReader", 3 );
    if( reader_ == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDbDataFromReader: Database has not been opened.");
    reader_->Notify(TdDataReader::tdrmsgGetData);
    int rsp = reader_->Wait(TdDataReader::tdrrespmsgDataReady, TdDataReader::tdrrespmsgNoData);
    if( rsp == TREADER_MSG_ERROR )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDbDataFromReader: The reader terminated with errors.");
    if( rsp != TdDataReader::tdrrespmsgDataReady && rsp != TdDataReader::tdrrespmsgNoData ) {
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDbDataFromReader: Invalid response from the reader.");
    }
    *lastchunk = true;
    rsp = (rsp != TdDataReader::tdrrespmsgNoData);
    if(rsp)
        reader_->GetbdbCdata(
            bdbCprodescs,
            bdbCpmbeg, bdbCpmend,
            szpm2dvfields,
            tscnt,
            lastchunk);
    return rsp;
}



// -------------------------------------------------------------------------
// CreateAlnWriter: create a thread for writing results to files
//
void JobDispatcher::CreateAlnWriter( 
    Configuration* config,
    const char* outdirname,
    const char* dbname,
    size_t prodbsize,
    size_t ndbentries,
    int nqueries )
{
    MYMSG( "JobDispatcher::CreateAlnWriter", 3 );
    writer_ = new AlnWriter( config, outdirname, dbname, prodbsize, ndbentries, nqueries );
    if( writer_ == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::CreateAlnWriter: Failed to create the writing thread.");
}

// -------------------------------------------------------------------------
// NotifyAlnWriter: notify the writer of the complete results for a query
//
void JobDispatcher::NotifyAlnWriter()
{
    MYMSG( "JobDispatcher::NotifyAlnWriter", 3 );
    if( writer_ ) {
        int rspcode = writer_->GetResponse();
        if( rspcode == WRITERTHREAD_MSG_ERROR || 
            rspcode == AlnWriter::wrttrespmsgTerminating)
            throw MYRUNTIME_ERROR(
            "JobDispatcher::NotifyAlnWriter: Results writer terminated with errors.");
        writer_->Notify(AlnWriter::wrtthreadmsgWrite);
    }
}

// -------------------------------------------------------------------------
// WaitForAlnWriterToFinish: notify the writer of the end of search and 
// wait for it to finish writings
//
void JobDispatcher::WaitForAlnWriterToFinish(bool error)
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::WaitForAlnWriterToFinish: Snd Msg %d",
                AlnWriter::wrtthreadmsgTerminate);
        MYMSG(msgbuf,3);
    MYMSGENDl
    if( writer_ == NULL ) {
        MYMSGBEGl(3)
            MYMSG("JobDispatcher::WaitForAlnWriterToFinish: Null Writer",3);
        MYMSGENDl
        return;
    }
    int rsp = WRITERTHREAD_MSG_UNSET;

    if(!error)
        rsp = writer_->WaitDone();

    if( rsp == WRITERTHREAD_MSG_ERROR ) {
        MYMSG("JobDispatcher::WaitForAlnWriterToFinish: Writer terminated with ERRORS.",1);
        return;
    }

    writer_->Notify(AlnWriter::wrtthreadmsgTerminate);
    //NOTE: do not wait, as the writer might have finished
//     int rsp = writer_->Wait(AlnWriter::wrttrespmsgTerminating);
//     if( rsp != AlnWriter::wrttrespmsgTerminating ) {
//         throw MYRUNTIME_ERROR(
//         "JobDispatcher::WaitForAlnWriterToFinish: Invalid response from the writer.");
//     }
//     MYMSGBEGl(3)
//         sprintf(msgbuf,"JobDispatcher::WaitForAlnWriterToFinish: Rcv Msg %d",rsp);
//         MYMSG(msgbuf,3);
//     MYMSGENDl
}



// -------------------------------------------------------------------------
// CreateDevMemoryConfigs: create memory configurations for all devices
//
void JobDispatcher::CreateDevMemoryConfigs(
    char** querypmbeg, char** querypmend,
    char** bdb1pmbeg, char** bdb1pmend,
    size_t nareasperdevice )
{
    MYMSG( "JobDispatcher::CreateDevMemoryConfigs", 3 );
    const mystring preamb = "JobDispatcher::CreateDevMemoryConfigs: ";

    for(int tid = 0; tid < DEVPROPs.GetNDevices(); tid++) {
        const DeviceProperties* dprop = DEVPROPs.GetDevicePropertiesAt(tid);
        if( dprop == NULL )
            continue;
        CuDeviceMemory* dmem = new CuDeviceMemory( 
            *dprop,
            HDPSCORES.GetScores()!=NULL,
            MOptions::GetMAPALN()==1,
            dprop->reqmem_,
            (int)nareasperdevice
        );
        if( dmem == NULL ) {
            warning((preamb+"Not enough memory to create all device memory configurations.").c_str());
            break;
        }
        memdevs_.push_back(dmem);
        dmem->CacheCompleteData(
            querypmbeg, querypmend,
            bdb1pmbeg, bdb1pmend
        );
    }
}



// -------------------------------------------------------------------------
// CreateWorkerThreads: create worker threads on the host side
//
void JobDispatcher::CreateWorkerThreads(
    const mystring* queryfnames,
    const mystring* querydescs,
    const char** bdb1descs,
    size_t prodbsize,
    size_t ndbentries )
{
    MYMSG( "JobDispatcher::CreateWorkerThreads", 3 );
    const mystring preamb = "JobDispatcher::CreateWorkerThreads: ";

    int tid = 0;

    for(int dmc = 0; dmc < (int)memdevs_.size(); dmc++) {
        if( memdevs_[dmc] == NULL )
            throw MYRUNTIME_ERROR(preamb + "Null device memory object.");
        for(int ano = 0; ano < memdevs_[dmc]->GetNAreas(); ano++, tid++) {
            DevCommThread* t = new DevCommThread( 
                tid,
                memdevs_[dmc],
                ano,
                GetConfiguration(),
                GetAlnWriter(),
                queryfnames, querydescs,
                bdb1descs,
                prodbsize, ndbentries
            );
            if( t == NULL ) {
                warning((preamb+"Not enough memory to create all required worker threads.").c_str());
                break;
            }
            hostworkers_.push_back(t);
        }
    }
}

// -------------------------------------------------------------------------
// GetDataChunkSize: get data chunk size the given worker can 
// process at a time
//
void JobDispatcher::GetDataChunkSize(
    int tid, size_t querylen,
    size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros)
{
    char msgbuf[BUF_MAX];
    bool waitfordata = (chunkdatasize && chunkdatalen && chunknpros);

    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::GetDataChunkSize: Snd Msg %d Adr %d %s",
                DevCommThread::tthreadmsgGetDataChunkSize, tid,
                waitfordata? "": "[no wait]");
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDataChunkSize: Worker thread is null.");

    tdc->SetQueryLen( querylen );
    tdc->Notify(DevCommThread::tthreadmsgGetDataChunkSize, tid);

    if( !waitfordata )
        return;

    int rsp = tdc->Wait(DevCommThread::ttrespmsgChunkSizeReady);
    if( rsp != DevCommThread::ttrespmsgChunkSizeReady ) {
        throw MYRUNTIME_ERROR(
        "JobDispatcher::GetDataChunkSize: Invalid response from the worker thread.");
    }
    MYMSGBEGl(3)
        sprintf(msgbuf,"JobDispatcher::GetDataChunkSize: Rcv Msg %d Adr %d",rsp,tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    tdc->GetChunkDataAttributes(chunkdatasize, chunkdatalen, chunknpros);
}

// -------------------------------------------------------------------------
// SubmitWorkerJob: submit a job for a given worker
//
void JobDispatcher::SubmitWorkerJob(
    int tid, bool lastchunk, int chunkno,
    int nqueries, int qrysernr, int qrystep,
    float* scorethlds, float* logevthlds,
    char** querypmbeg,
    char** bdb1pmbeg, char** bdb1pmend,
    const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend, size_t* szpm2dvfields,
    TSCounterVar* tscnt )
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::SubmitWorkerJob: Snd Msg %d Adr %d [no wait]",
                DevCommThread::tthreadmsgProcessNewData, tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::SubmitWorkerJob: Worker thread is null.");

    tdc->SetMstrQueryBDbdata(
        lastchunk,
        chunkno,
        nqueries, qrysernr, qrystep,
        scorethlds, logevthlds,
        querypmbeg,
        bdb1pmbeg, bdb1pmend,
        bdbCdesc, bdbCpmbeg, bdbCpmend, szpm2dvfields,
        tscnt );
    tdc->Notify(DevCommThread::tthreadmsgProcessNewData, tid);
}

// -------------------------------------------------------------------------
// SubmitWorkerJob_obs: submit a job for a given worker
//
void JobDispatcher::SubmitWorkerJob_obs(
    int tid,
    int chunkno,
    int qrysernr,
    size_t nqyposs, float scorethld, float logevthld,
    char** querypmbeg, char** querypmend,
    char** bdb1pmbeg, char** bdb1pmend,
    const char** bdbCdesc, char** bdbCpmbeg, char** bdbCpmend, size_t* szpm2dvfields,
    TSCounterVar* tscnt )
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::SubmitWorkerJob_obs: Snd Msg %d Adr %d [no wait]",
                DevCommThread::tthreadmsgProcessNewData, tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::SubmitWorkerJob_obs: Worker thread is null.");

    tdc->SetMstrQueryBDbdata_obs(
        chunkno,
        qrysernr,
        nqyposs, scorethld, logevthld,
        querypmbeg, querypmend,
        bdb1pmbeg, bdb1pmend,
        bdbCdesc, bdbCpmbeg, bdbCpmend, szpm2dvfields,
        tscnt );
    tdc->Notify(DevCommThread::tthreadmsgProcessNewData, tid);
}

// -------------------------------------------------------------------------
// WaitForWorker: wait until the worker can accept new part of data
//
void JobDispatcher::WaitForWorker( int tid )
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::WaitForWorker: Wrk %d",tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL ) {
        MYMSGBEGl(3)
            sprintf(msgbuf,"JobDispatcher::WaitForWorker: Null Worker %d",tid);
            MYMSG(msgbuf,3);
        MYMSGENDl
        return;
    }

    //wait until data can be accessed and get the last response if available
    int rsp =tdc->waitForDataAccess();

    if( rsp == THREAD_MSG_ERROR ) {
        sprintf(msgbuf,"JobDispatcher::WaitForWorker: Worker %d terminated with errors",tid);
        throw MYRUNTIME_ERROR(msgbuf);
    }

    //wait until the worker has finished processing
    tdc->IsIdle(NULL/*mstr_set_data_empty*/, true/*wait*/);
}

// -------------------------------------------------------------------------
// ProbeWorker: check for worker idle state
//
void JobDispatcher::ProbeWorker( int tid )
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::ProbeWorker: Snd Msg %d Adr %d",
                DevCommThread::tthreadmsgProbe, tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL ) {
        MYMSGBEGl(3)
            sprintf(msgbuf,"JobDispatcher::ProbeWorker: Null Worker %d",tid);
            MYMSG(msgbuf,3);
        MYMSGENDl
        return;
    }

    //first, make sure the worker is idle and there are no pending jobs
    bool idle = false, data_empty = false;
    for(; !idle || !data_empty;) idle = tdc->IsIdle(&data_empty, true/*wait*/);

    //then, send a probe message
    tdc->Notify(DevCommThread::tthreadmsgProbe, tid);
    int rsp = tdc->Wait(DevCommThread::ttrespmsgProbed);
    if( rsp != DevCommThread::ttrespmsgProbed ) {
        throw MYRUNTIME_ERROR(
        "JobDispatcher::ProbeWorker: Invalid response from the worker thread.");
    }
    MYMSGBEGl(3)
        sprintf(msgbuf,"JobDispatcher::ProbeWorker: Rcv Msg %d Adr %d",rsp,tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
}

// -------------------------------------------------------------------------
// WaitForAllWorkersToFinish: wait until all workers idle;
//
inline
void JobDispatcher::WaitForAllWorkersToFinish()
{
    for( int tid = 0; tid < (int)hostworkers_.size(); tid++ )
        ProbeWorker(tid);
}

// -------------------------------------------------------------------------
// TerminateWorker: terminate a worker
//
void JobDispatcher::TerminateWorker( int tid )
{
    char msgbuf[BUF_MAX];
    MYMSGBEGl(3)
        sprintf(msgbuf, "JobDispatcher::TerminateWorker: Snd Msg %d Adr %d [no wait]",
                DevCommThread::tthreadmsgTerminate, tid);
        MYMSG(msgbuf,3);
    MYMSGENDl
    DevCommThread* tdc = hostworkers_[tid];
    if( tdc == NULL ) {
        MYMSGBEGl(3)
            sprintf(msgbuf,"JobDispatcher::TerminateWorker: Null Worker %d",tid);
            MYMSG(msgbuf,3);
        MYMSGENDl
        return;
    }
    tdc->Notify(DevCommThread::tthreadmsgTerminate, tid);
    //do not wait for a response
}

// -------------------------------------------------------------------------
// TerminateAllWorkers: terminate all the workers
//
void JobDispatcher::TerminateAllWorkers()
{
    for( int tid = 0; tid < (int)hostworkers_.size(); tid++ )
        TerminateWorker(tid);
}



// -------------------------------------------------------------------------
// BcastMessage: broadcast a message to worker threads and 
// optionally wait for a response from each of them
//
template<typename TI, typename TO>
inline
void JobDispatcher::BcastMessage(
    bool waitforresponse,
    int msgcode, int addressee,
    TWriteDataForWorker1<TI> writefunc,
    const std::vector<TI>& datatowrite,
    TGetDataFromWorker1<TO> readfunc,
    std::vector<TO>& datatoread,
    std::unique_ptr<PMBatchProData> bdbC )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::BcastMessage: Msg %d Adr %d",
                msgcode, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
    //
    //iterate over all to send a message
    for(int tid = 0; tid < (int)hostworkers_.size(); tid++ ) {
        DevCommThread* tdc = hostworkers_[tid];
        if( tdc == NULL || (0<=addressee && addressee!=tid))
            continue;
        {//set the message under lock
            std::unique_lock<std::mutex> lck_brd_msg(tdc->GetPrivateMutex());
            //wait until the worker has no pending jobs
            tdc->GetPrivateCV().wait(lck_brd_msg,
                [tdc]{return 
                    tdc->GetMessage() < 0 || 
                    tdc->GetResponseMsg() == THREAD_MSG_ERROR;}
            );
            tdc->ResetResponseMsg();///reset response before beginning messaging
            tdc->SetBcastMessageAndAddressee(msgcode, tid);
            //write data before release of the mutex
            (this->*writefunc)(msgcode, tid, tdc, datatowrite, std::move(bdbC));
        }
        //lock released, notify the worker
        tdc->GetPrivateCV().notify_one();
    }

    //iterate over all to receive a private message
    for(int tid = 0; tid < (int)hostworkers_.size(); tid++ ) {
        DevCommThread* tdc = hostworkers_[tid];
        if( tdc == NULL || (0<=addressee && addressee!=tid))
            continue;
        int response = THREAD_MSG_UNSET;
        if(waitforresponse){
            //wait for a response
            std::unique_lock<std::mutex> lck_rsp_msg(tdc->GetPrivateMutex());
            tdc->GetPrivateCV().wait(lck_rsp_msg,
                [tdc]{return 
                    (0 <= tdc->GetResponseMsg() && 
                     tdc->GetResponseMsg() <= DevCommThread::ttrespmsgTerminating) || 
                    tdc->GetResponseMsg() == THREAD_MSG_ERROR
                    ;}
            );
            //lock on the mutex is back
            response = tdc->GetResponseMsg();
            tdc->ResetResponseMsg();
            //read data under the lock
            if( response != THREAD_MSG_ERROR )
                (this->*readfunc)(msgcode, response, tid, tdc, datatoread);
        }
        else {
            //check if not terminated
            std::unique_lock<std::mutex> 
                lck_chck_state(tdc->GetPrivateMutex(), std::defer_lock);
            if( lck_chck_state.try_lock()) {
                //get the state of the thread under the lock
                response = tdc->GetResponseMsg();
            }
        }
        //join threads that issued an error or terminated
        if( response == DevCommThread::ttrespmsgTerminating || 
            response == THREAD_MSG_ERROR ) {
            delete tdc;
            hostworkers_[tid] = NULL;
        }
    }
}
// -------------------------------------------------------------------------
// GetAvailableWorker: get a worker ready to accept new data for processing
//
void JobDispatcher::GetAvailableWorker( int* tid )
{
    MYMSG( "JobDispatcher::GetAvailableWorker", 7 );

    bool alltermed = true;//all workers terminated

    if( tid == NULL )
        return;

    *tid = JDSP_WORKER_BUSY;

    //first, try to find a waiting worker;
    //if not that, try to find a worker with an empty data slot
    for(int ptid = 0; ptid < (int)hostworkers_.size(); ptid++ ) {
        //iterate over all worker to check for availability
        DevCommThread* tdc = hostworkers_[ptid];
        if( tdc == NULL )
            continue;
        alltermed = false;
        bool idle, data_empty;
        idle = tdc->IsIdle(&data_empty);
        if( data_empty && (idle || *tid < 0)) {
            *tid = ptid;
            if( idle )
                //the worker is ready to accept new data and perform computation
                break;
        }
    }

    if( alltermed )
        *tid = JDSP_WORKER_NONE;

    MYMSGBEGl(5)
        char msgbuf[BUF_MAX];
        sprintf(msgbuf,"JobDispatcher::GetAvailableWorker: Wrk %d",*tid);
        MYMSG(msgbuf,5);
    MYMSGENDl
}
// -------------------------------------------------------------------------
// WaitForAvailableWorker: wait until a worker becomes ready to process data
// or accept new data for processing
inline
void JobDispatcher::WaitForAvailableWorker( int* tid )
{
    MYMSG( "JobDispatcher::WaitForAvailableWorker", 7 );
    if( tid == NULL )
        return;
    for( *tid = JDSP_WORKER_BUSY; *tid == JDSP_WORKER_BUSY; 
        GetAvailableWorker(tid));
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::WaitForAvailableWorker: Wrk %d",*tid);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
}
// -------------------------------------------------------------------------
// GetNextWorker: get the next worker irrespective of its busy or idle 
// status;
// NOTE: tid, the address of the current busy worker should be initialized
inline
void JobDispatcher::GetNextWorker( int* tid )
{
    MYMSG( "JobDispatcher::GetNextWorker", 7 );
    if( tid == NULL )
        return;

    int initid = (0 <= *tid && *tid < (int)hostworkers_.size())? *tid: 0;

    *tid = JDSP_WORKER_NONE;

    for(int ptid = initid+1; ; ptid++ ) {
        if((int)hostworkers_.size() <= ptid)
            ptid = 0;
        if((int)hostworkers_.size() <= ptid)
            break;
        DevCommThread* tdc = hostworkers_[ptid];
        if( tdc == NULL ) {
            if(ptid == initid)
                break;
            else
                continue;
        }
        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            sprintf( msgbuf, "JobDispatcher::GetNextWorker: Wrk %d",ptid);
            MYMSG( msgbuf, 3 );
        MYMSGENDl
        *tid = ptid;
        break;
    }
}

// =========================================================================
// =========================================================================



// -------------------------------------------------------------------------
// CalculateScoreThreshold: calculate score threshold given the e-value 
// threshold;
// scorethld, score threshold given e-value in the options;
// logevthld, log E-value threshold
//
void JobDispatcher::CalculateScoreThreshold(
    size_t nqyposs, size_t dbsize,
    float* scorethld, float* logevthld)
{
    MYMSG( "JobDispatcher::CalculateScoreThreshold", 5 );
    const mystring preamb = "JobDispatcher::CalculateScoreThreshold: ";
    size_t sspace = nqyposs * dbsize;
    float evaluethld = MOptions::GetEVAL();
    float explambda = GetConfiguration(CTGapped).GetLambda();
    float expK = GetConfiguration(CTGapped).GetK();
    float score = 0.0f;

    if( sspace < 1 ) {
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid search space (query or database size)." );
    }
    if( expK <= 0.0f || explambda <= 0.0f ) {
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid statistical parameters in the configuration file." );
    }

    if( evaluethld <= 0.0f ) {
        *scorethld = 999.0f;
        *logevthld = -99.0f;
        return;
    }
    else {
        *logevthld = logf(evaluethld);
        //
        if( evaluethld < 1.0f )
            evaluethld *= 10.0f;
        else
            evaluethld *= 3.0f;
    }

    evaluethld /= expK * sspace;

    if( evaluethld <= 0.0f ) {
        *scorethld = 999.0f;
        return;
    }

    evaluethld = logf(evaluethld);

    //score = (logf(expK * sspace) - evaluethld) / explambda;
    score = -evaluethld / explambda;

    if( score < 0.0f )
        score = 0.0f;

    *scorethld = score;
}

// -------------------------------------------------------------------------
// ReadProfiles: read profiles from a given database;
// db, (a database of) profiles;
// bdata, batch data of proffiles read from the database;
// maxprolen, maximum length of profiles read;
// maxdatasize, maximum amount of memory the profile data can occupy; it is 
//  ignored if readall is true;
// readall, if true, read all profiles of the database;
// stoppos, position of the datbase file where the read has stopped; 
//  ignored if readall is true;
//
void JobDispatcher::ReadProfiles( 
    CuDbReader& db, PMBatchProData& bdata, int* maxprolen,
    size_t maxdatasize, bool readall, TCharStream* stoppos )
{
    MYMSG( "JobDispatcher::ReadProfiles", 3 );
    const mystring preamb = "JobDispatcher::ReadProfiles: ";
    myruntime_error mre;
    char msgbuf[BUF_MAX];
    bool noteof = true;
    int prolen, scale, profnr = 0; 
    unsigned int distn = 0;//distance in positions to a profile
    mystring  pdesc, pfile;
    //*test*/const size_t nflds = bdata.GetNoFields();
    const size_t nheadflds = bdata.GetNoHeaderFields();
    char auxpsdat[nheadflds][SZFPTYPE];//auxiliary buffer for profile header, i.e., profile-specific data
    char* auxps[nheadflds];//pointers to the above structure
    //*test*/char* prevcopypmdata[nflds];

    for( int i = 0; i < (int)nheadflds; i++ )
        auxps[i] = auxpsdat[i];

    try {
        db.Open();

        if( stoppos )
            db.GetDbPos( stoppos );

        if( maxdatasize )
            bdata.AllocNewSizeHostDataInit( maxdatasize );

        for(; !db.Eof() && maxdatasize;) 
        {
            noteof = db.NextProfileHeader( &pdesc, &pfile, &prolen, &scale, auxps );

            if( maxprolen && *maxprolen < prolen )
                *maxprolen = prolen;
                
            if( !noteof ) {
                sprintf( msgbuf, "Unexpected end of file %d", profnr+1 );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            //avoid overflow
            if( INT_MAX <= (unsigned int)prolen ||
                (unsigned int)(INT_MAX - prolen) <= distn )
                break;

            if( bdata.MemDeficit( prolen )) {
                if( !readall )
                    break;
                maxdatasize *= 2;
                bdata.ReallocNewSizeHostData( maxdatasize );
            }

            //{{save data read from the profile header
            //save the beginning of a next profile
            //*test*/memcpy( prevcopypmdata, bdata.GetPMDataEnd(), nflds * sizeof(void*));
            bdata.CopyPMDataHeaderNext( pdesc, pfile, auxps );
            //}}

            noteof = db.Next( distn, profnr, prolen, scale, bdata.GetPMData(), bdata.GetPMDataEnd());

            distn += prolen;
            profnr++;

            if( stoppos )
                db.GetDbPos( stoppos );

            //test:
            //*test*/TextWriteProfileHD( stderr, pdesc, pfile, scale, bdata.GetPMData(), prevcopypmdata );

        }//for(;!myeof;)

    } catch( myexception const& ex ) {
        mre = ex;
    }

    db.Close();

    if( mre.isset())
        throw mre;

    if( readall )
        bdata.ReallocHostData();//free unused memory
}



// =========================================================================
// Run: launch worker threads for searching and aligning profiles
//
void JobDispatcher::Run()
{
    MYMSG( "JobDispatcher::Run", 3 );
    const mystring preamb = "JobDispatcher::Run: ";

    myruntime_error mre;
    char msgbuf[BUF_MAX];
    const size_t initquerydatasize = 100UL*ONEM;
    size_t chunkdatasize = 0UL;
    size_t chunkdatalen = 0UL;//maximum total length of db profile positions that can be processed at once
    size_t chunknpros = 0UL;//maximum number of db profiles allowed to be processed at once
    //access to input batch of profiles, which also can be a database:
    CuInputReader inpdb(GetInput(), CLOptions::GetIO_FILEMAP());
    PMBatchProData bquery;//batch of input profile data
    int maxquerylen = 0; 
    const size_t nflds = PMBatchProData::GetNoFields();
    char* querypmdata[nflds];//running pointers for queries

    std::vector<float> vec_scthld;
    std::vector<float> vec_evthld;

    char** bdbCprodescs;//profile descriptions read
    char** bdbCpmbeg, **bdbCpmend;//profile data encapsulated in a data chunk
    size_t* szpm2dvfields;//sizes of profile data fields
    TSCounterVar* tscnt;//counter of how many agents read the data
    bool lastchunk = false;//indicator of the last chank
    float scorethld, logevthld;//score and e-value thresholds that depend on the query length
    int tid = JDSP_WORKER_BUSY;//worker id
    int didmm = -1;//id of a device with the minimum requested memory

    ReadConfiguration();


    if( DEVPROPs.GetNDevices() < 1 ) {
        warning("There is no available device to run the program on.");
        message("Please use a version to run on CPU.");
        return;
    }


    try {
        //INITIALIZATION...
        MYMSG( "Reading queries", 1 );
        ReadProfiles( inpdb, bquery, &maxquerylen, initquerydatasize, true/*readall*/, NULL/*stoppos*/);
        inpdb.Destroy();//explicitly destroy unused resources here

        size_t nqueries = PMBatchProData::GetNoProsFromTo( bquery.GetPMData(), bquery.GetPMDataEnd());
        bool querinf = PMBatchProData::PMDataPresentFromTo( bquery.GetPMData(), bquery.GetPMDataEnd());

        if( !querinf || nqueries < 1 )
            throw MYRUNTIME_ERROR( preamb + "Query data has not been cached." );

        CreateReader( 
            GetDatabase(),
            CLOptions::GetIO_FILEMAP(),
            (int)nqueries,
            GetConfiguration());

        vec_scthld.resize(nqueries, -1.0f);
        vec_evthld.resize(nqueries, -1.0f);

        if(( didmm = DEVPROPs.GetDevIdWithMinRequestedMem()) < 0 )
            throw MYRUNTIME_ERROR( preamb + "Failed to obtain required device id." );

        GetDbSizeFromReader();

        //create the results writing thread
        CreateAlnWriter( 
            GetConfiguration(),
            output_,
            GetDatabase(),
            GetDbSize(),
            GetNoSequences(),
            (int)nqueries );

        //create memory configuration for each device
        CreateDevMemoryConfigs(
            bquery.GetPMData(), bquery.GetPMDataEnd(),
            NULL/*bdb1pmbeg*/, NULL/*bdb1pmend*/,
            CLOptions::GetDEV_NGRIDS()/*nareasperdevice*/
        );

        if( memdevs_.size() < 1 )
            throw MYRUNTIME_ERROR(preamb+"Failed to create device memory configuration.");

        //initialize the memory sections of all devices
        // according to the given maximum query length
        for(int dmc = 0; dmc < (int)memdevs_.size(); dmc++ ) {
            size_t szchdata;
            if( memdevs_[dmc] == NULL)
                continue;
            szchdata = memdevs_[dmc]->CalcMaxDbDataChunkSize(maxquerylen);
            //get globally valid data chunk size
            if( didmm == dmc) {
                chunkdatasize = szchdata;
                chunkdatalen = memdevs_[dmc]->GetCurrentMaxDbPos();
                chunknpros = memdevs_[dmc]->GetCurrentMaxNDbPros();
                //make sure max # profiles does not exceed technical specifications
                chunknpros = SLC_MIN(chunknpros,(size_t)memdevs_[dmc]->GetDeviceProp().maxGridSize_[1]);
            }
        }

        //create worker threads
        CreateWorkerThreads(
                bquery.GetFnames(), bquery.GetDescs(),
                NULL/*bdb1desc*/,
                GetDbSize(), GetNoSequences());

        if( hostworkers_.size() < 1 )
            throw MYRUNTIME_ERROR(preamb+"Failed to create any worker threads.");

        MYMSGBEGl(1)
            char strbuf[BUF_MAX];
            sprintf(strbuf,"Processing in chunks of: size %zu length %zu profiles %zu",
                chunkdatasize, chunkdatalen, chunknpros);
            MYMSG(strbuf, 1);
        MYMSGENDl

        if( chunkdatasize < 1 || chunkdatalen < 1 || chunknpros < 1)
            throw MYRUNTIME_ERROR(preamb+"Invalid calculated data chunk size.");

        GetReader()->SetChunkDataAttributes(chunkdatasize, chunkdatalen, chunknpros);

        bquery.GetPMData(querypmdata);

        //calculate thresholds for each query in advance
        MYMSG((preamb + "Calculating thresholds").c_str(),3);
        char* plen = querypmdata[pps2DLen];
        for( size_t i = 0; i < nqueries; i++, plen += PMBatchProData::GetLengthFieldSize()) {
            size_t nqyposs = PMBatchProData::GetLengthValue(plen);
            CalculateScoreThreshold( nqyposs, GetDbSize(), &scorethld, &logevthld);
            vec_scthld[i] = scorethld;
            vec_evthld[i] = logevthld;
        }

        //CHUNKS..
        for( int chkno = 0; tid != JDSP_WORKER_NONE; chkno++ )
        {
            if( ! GetDbDataFromReader(
                bdbCprodescs,  bdbCpmbeg, bdbCpmend,  szpm2dvfields,
                tscnt,  &lastchunk))
                //if no data, finish
                break;

            //NOTE: wait for all workers to finish; new chunk implies new data!
            for( tid = 0; tid < (int)hostworkers_.size(); tid++ )
                WaitForWorker(tid);

            MYMSGBEGl(1)
                size_t nCpros = PMBatchProData::GetNoProsFromTo( bdbCpmbeg, bdbCpmend );
                size_t nposits = PMBatchProData::GetNoPositsFromTo(bdbCpmbeg, bdbCpmend);
                sprintf(msgbuf,"%s[=====] Processing database CHUNK No. %d: %.1f%% db positions (%zu pros.)",
                        NL, chkno, (float)nposits*100.0f/(float)GetDbSize(),nCpros);
                MYMSG(msgbuf,1);
            MYMSGENDl

            //QUERIES passed to WORKERS...
            for( tid = 0; tid < (int)hostworkers_.size(); tid++ ) {
                //
                ProcessPart(tid, lastchunk, chkno+1,
                    (int)nqueries, tid/*start-with query no.*/, (int)hostworkers_.size()/*step*/,
                    vec_scthld.data(), vec_evthld.data(),
                    querypmdata,
                    NULL/*bdb1pmbeg*/, NULL/*bdb1pmend*/,
                    (const char**)bdbCprodescs, bdbCpmbeg, bdbCpmend, szpm2dvfields,
                    tscnt);
            }
        }

        //wait for the workers to finish
        for( tid = 0; tid < (int)hostworkers_.size(); tid++ )
            WaitForWorker(tid);

        WaitForAllWorkersToFinish();

    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Exception caught.");
    }

    //workers' children may still have job; wait for the writer to finish first
    WaitForAlnWriterToFinish( mre.isset());

    TerminateAllWorkers();

    if( mre.isset())
        throw mre;

    MYMSG("Done.",1);
}

// -------------------------------------------------------------------------
// ProcessPart: process a part of data, a number of databases 
// profiles for a number of queries;
// tid, worker id;
// chunkno, chunk serial number;
// nqueries, total number of queries;
// qrysernr, query number for the wroker to start computation with;
// qrystep, step in the number of queries;
// scorethld, score thresholds for each query;
// logevthld, e-value thresholds for each query;
// querypmbeg, beginnings of the query fields;
// bdb1pmbeg, beginnings of the cached db profile fields;
// bdb1pmend, ends of the cached db profile fields;
// bdbCprodescs, profile descriptions for each db profile in a data chunk;
// bdbCpmbeg, beginnings of the db profile fields of a chunk;
// bdbCpmend, ends of the db profile fields of a chunk;
// szpm2dvfields, sizes of the db profile fields of a chunk;
// tscnt, thread-safe counter of processing agents;
//
void JobDispatcher::ProcessPart(
    int tid, bool lastchunk, int chunkno,
    int nqueries, int qrysernr, int qrystep,
    float* scorethlds, float* logevthlds,
    char** querypmbeg,
    char** bdb1pmbeg, char** bdb1pmend,
    const char** bdbCprodescs, char** bdbCpmbeg, char** bdbCpmend, size_t* szpm2dvfields,
    TSCounterVar* tscnt )
{
    MYMSG( "JobDispatcher::ProcessPart", 4 );
    const mystring preamb = "JobDispatcher::ProcessPart: ";

#ifdef __DEBUG__
    if( !querypmbeg )
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if( !bdb1pmbeg && !bdb1pmend && !bdbCpmbeg && !bdbCpmend )
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if((bdbCpmbeg || bdbCpmend) && (!bdbCprodescs || !szpm2dvfields || !tscnt))
        throw MYRUNTIME_ERROR( preamb + "Inconsistent arguments." );
#endif

    bool bdb1inf = PMBatchProData::PMDataPresentFromTo( bdb1pmbeg, bdb1pmend );
    bool bdbCinf = PMBatchProData::PMDataPresentFromTo( bdbCpmbeg, bdbCpmend );

    if( !bdb1inf && !bdbCinf )
        throw MYRUNTIME_ERROR( preamb + "Db profile data is empty." );

    SubmitWorkerJob(
        tid, lastchunk, chunkno,
        nqueries, qrysernr, qrystep,
        scorethlds, logevthlds,
        querypmbeg,
        bdb1pmbeg, bdb1pmend,
        bdbCprodescs, bdbCpmbeg, bdbCpmend, szpm2dvfields,
        tscnt);
}
