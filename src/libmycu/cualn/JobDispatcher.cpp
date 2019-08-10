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
    scorethld_(0.0f),
    logevthld_(0.0f),

    eval_upper_( eval_thld ),
    max_no_hits_( no_hits ),
    max_no_alns_( no_alns ),
    show_pars_( showpars ),

    prodb_( database, !CLOptions::GetIO_NOFILEMAP()),
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
    scorethld_(0.0f),
    logevthld_(0.0f),

    eval_upper_( 0.0f ),
    max_no_hits_( 0 ),
    max_no_alns_( 0 ),
    show_pars_( false ),

    prodb_( NULL ),
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
    if( writer_ ) {
        writer_->Notify(AlnWriter::wrtthreadmsgTerminate);
        delete writer_;
        writer_ = NULL;
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
// CreateWorkerThreads: create worker threads on the host side
//
void JobDispatcher::CreateWorkerThreads(
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
{
    MYMSG( "JobDispatcher::CreateWorkerThreads", 3 );
    const mystring preamb = "JobDispatcher::CreateWorkerThreads: ";

    for(int tid = 0; tid < DEVPROPs.GetNDevices(); tid++) {
        if( DEVPROPs.GetDevicePropertiesAt(tid) == NULL ) {
            hostworkers_.push_back(NULL);
            continue;
        }
        DevCommThread* t = new DevCommThread( 
            tid,
            *DEVPROPs.GetDevicePropertiesAt(tid),
            GetAlnWriter(),
            GetConfiguration(),
            queryfnames, querydescs,
            querypmbeg, querypmend,
            //
            bdb1fnames, bdb1descs,
            bdb1pmbeg, bdb1pmend,
            prodbsize, ndbentries
        );
        if( t == NULL ) {
            warning((preamb+"Not enough memory to create all required worker threads.").c_str());
            break;
        }
        hostworkers_.push_back(t);
    }
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
        if( tdc == NULL )
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
            tdc->SetBcastMessageAndAddressee(
                msgcode, (addressee<0 || addressee==tid)? tid: -1);
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
void JobDispatcher::GetAvailableWorker( 
    int* tid,
    size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros )
{
    MYMSG( "JobDispatcher::GetAvailableWorker", 7 );

    bool alltermed = true;//all workers terminated

    if( tid == NULL )
        return;

    *tid = JDSP_WORKER_BUSY;

    if( chunkdatasize ) *chunkdatasize = 0UL;
    if( chunkdatalen ) *chunkdatalen = 0UL;
    if( chunknpros ) *chunknpros = 0UL;

    //first, try to find a waiting worker;
    //if not that, try to find a worker with an empty data slot
    for(int busy = 0; *tid < 0 && busy < 2; busy++ ) {
        //iterate over all worker to check for availability
        for(int ptid = 0; *tid < 0 && ptid < (int)hostworkers_.size(); ptid++ ) {
            DevCommThread* tdc = hostworkers_[ptid];
            if( tdc == NULL )
                continue;
            alltermed = false;
            //verify if a worker is ready
            {
                std::unique_lock<std::mutex> 
                    lck_busy(tdc->GetPrivateMutex(), std::defer_lock);
                if( busy? 1: lck_busy.try_lock()) {
                    bool empty = tdc->GetMstrDataEmpty();
                    if( empty ) {
                        //the worker is ready to accept new data and perform computation;
                        MYMSGBEGl(3)
                            char msgbuf[BUF_MAX];
                            sprintf( msgbuf, "JobDispatcher::GetAvailableWorker: Wrk %d",ptid);
                            MYMSG( msgbuf, 3 );
                        MYMSGENDl
                        *tid = ptid;
                        //read data under lock
                        tdc->GetChunkDataAttributes( chunkdatasize, chunkdatalen, chunknpros);
                    }
                }
                //release of the mutex if locked
            }
        }
    }

    if( alltermed )
        *tid = JDSP_WORKER_NONE;
}
// -------------------------------------------------------------------------
// WaitForAvailableWorker: wait until a worker becomes ready to process data
// or accept new data for processing
inline
void JobDispatcher::WaitForAvailableWorker( 
    int* tid,
    size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros )
{
    MYMSG( "JobDispatcher::WaitForAvailableWorker", 7 );
    if( tid == NULL )
        return;
    for( *tid = JDSP_WORKER_BUSY; *tid == JDSP_WORKER_BUSY; 
        GetAvailableWorker(tid, chunkdatasize, chunkdatalen, chunknpros));
}
// -------------------------------------------------------------------------
// GetNextWorker: get the next worker irrespective of its busy or idle 
// status;
// NOTE: tid, the address of the current busy worker should be initialized
inline
void JobDispatcher::GetNextWorker( 
    int* tid,
    size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros )
{
    MYMSG( "JobDispatcher::GetNextWorker", 7 );
    if( tid == NULL )
        return;

    *tid = JDSP_WORKER_NONE;

    int initid = (0 <= *tid && *tid < (int)hostworkers_.size())? *tid: 0;

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
        //read data under lock
        tdc->GetChunkDataAttributes( chunkdatasize, chunkdatalen, chunknpros);
        break;
    }
}
// -------------------------------------------------------------------------
// WaitForAllWorkersToFinish: wait until all workers idle;
//
inline
void JobDispatcher::WaitForAllWorkersToFinish()
{
    std::vector<int> dummy;
    BcastMessage(
        true,//waitforresponse
        DevCommThread::tthreadmsgProbe, 
        -1,//addressee: broadcast to all threads
        &JobDispatcher::VoidWriteBeforeSending, dummy,
        &JobDispatcher::ReadOnResponseOfMsgProbe, dummy );
}

// =========================================================================
// WriteForMsgGetDataChunkSize: write data under lock before sending message
// MsgGetDataChunkSize
inline
void JobDispatcher::WriteForMsgGetDataChunkSize( 
    int reqmsg, int addressee, DevCommThread* tdc, 
    const std::vector<size_t>& params,
    std::unique_ptr<PMBatchProData> )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::WriteForMsgGetDataChunkSize: "
                "Msg %d Adr %d", reqmsg, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
    if( reqmsg != DevCommThread::tthreadmsgGetDataChunkSize || tdc == NULL )
        throw MYRUNTIME_ERROR("JobDispatcher::WriteForMsgGetDataChunkSize: Invalid arguments.");
    tdc->SetQueryAttributes( params[0], GetScoreThreshold(), GetLogEThreshold());
}
// ReadOnResponseOfMsgGetDataChunkSize: read data under lock after 
// receiving a response from a worker
inline
void JobDispatcher::ReadOnResponseOfMsgGetDataChunkSize(
    int reqmsg, int rspmsg, int addressee, DevCommThread* tdc, std::vector<size_t>& )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::ReadOnResponseOfMsgGetDataChunkSize: "
                "Msg %d Rsp %d Adr %d", reqmsg, rspmsg, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
    if( reqmsg != DevCommThread::tthreadmsgGetDataChunkSize || 
        rspmsg != DevCommThread::ttrespmsgDataReady || tdc == NULL )
        throw MYRUNTIME_ERROR("JobDispatcher::ReadOnResponseOfMsgGetDataChunkSize: Invalid arguments.");
}
// -------------------------------------------------------------------------
// WriteForMsgProcessNewData: write data under lock before sending message
// MsgProcessNewData
inline
void JobDispatcher::WriteForMsgProcessNewData( 
    int reqmsg, int addressee, DevCommThread* tdc, 
    const std::vector<char**>& params,
    std::unique_ptr<PMBatchProData> bdbC )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::WriteForMsgProcessNewData: "
                "Msg %d Adr %d", reqmsg, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
    if( reqmsg != DevCommThread::tthreadmsgProcessNewData || tdc == NULL )
        throw MYRUNTIME_ERROR("JobDispatcher::WriteForMsgProcessNewData: Invalid arguments.");
    tdc->SetMstrQueryBDbdata(
        (int)(size_t)params[0]/*qrysernr*/,
        std::move(bdbC),
        params[1]/*querypmbeg*/, params[2]/*querypmend*/,
        params[3]/*bdb1pmbeg*/, params[4]/*bdb1pmend*/,
        params[5]/*bdbCpmbeg*/, params[6]/*bdbCpmend*/
    );
}
// -------------------------------------------------------------------------
// ReadOnResponseOfMsgProbe: read a response after sending a probe message;
template<typename T>
inline
void JobDispatcher::ReadOnResponseOfMsgProbe(
    int reqmsg, int rspmsg, int addressee, DevCommThread* tdc, std::vector<T>& )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::ReadOnResponseOfMsgProbe: "
                "Msg %d Rsp %d Adr %d", reqmsg, rspmsg, addressee);
        MYMSG( msgbuf, 4 );
    MYMSGENDl
    if( tdc == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::ReadOnResponseOfMsgProbe: Invalid arguments.");
    if( reqmsg != DevCommThread::tthreadmsgProbe || 
        rspmsg != DevCommThread::ttrespmsgProbed || tdc == NULL )
        throw MYRUNTIME_ERROR(
        "JobDispatcher::ReadOnResponseOfMsgProbe: Invalid response from a worker.");
}
// -------------------------------------------------------------------------
// DummyWriteForMsgTerminate: void write function before sending a message 
// to workers
template<typename T>
inline
void JobDispatcher::VoidWriteBeforeSending( 
    int reqmsg, int addressee, DevCommThread*, const std::vector<T>&,
    std::unique_ptr<PMBatchProData> )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::VoidWriteBeforeSending: "
                "Msg %d Adr %d", reqmsg, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
}
// VoidReadOnReceive: void read function after receiving a response from a 
// worker
template<typename T>
inline
void JobDispatcher::VoidReadOnReceive(
    int reqmsg, int rspmsg, int addressee, DevCommThread*, std::vector<T>& )
{
    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "JobDispatcher::VoidReadOnReceive: "
                "Msg %d Rsp %d Adr %d", reqmsg, rspmsg, addressee);
        MYMSG( msgbuf, 3 );
    MYMSGENDl
}
// =========================================================================



// -------------------------------------------------------------------------
// CalculateScoreThreshold: calculate score threshold given the e-value 
// threshold
//
void JobDispatcher::CalculateScoreThreshold( size_t nqyposs, size_t dbsize )
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
        SetScoreThreshold( 999.0f );
        SetLogEThreshold( -99.0f );
        return;
    }
    else {
        SetLogEThreshold( logf(evaluethld));
        //
        if( evaluethld < 1.0f )
            evaluethld *= 10.0f;
        else
            evaluethld *= 3.0f;
    }

    evaluethld /= expK * sspace;

    if( evaluethld <= 0.0f ) {
        SetScoreThreshold( 999.0f );
        return;
    }

    evaluethld = logf(evaluethld);

    //score = (logf(expK * sspace) - evaluethld) / explambda;
    score = -evaluethld / explambda;

    if( score < 0.0f )
        score = 0.0f;

    SetScoreThreshold( score );
}

// -------------------------------------------------------------------------
// ReadProfiles: read profiles from a given database;
// db, (a database of) profiles;
// bdata, batch data of proffiles read from the database;
// maxdatasize, maximum amount of memory the profile data can occupy; it is 
//  ignored if readall is true;
// readall, if true, read all profiles of the database;
// stoppos, position of the datbase file where the read has stopped; 
//  ignored if readall is true;
//
void JobDispatcher::ReadProfiles( 
    CuDbReader& db, PMBatchProData& bdata, 
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



// -------------------------------------------------------------------------
// Run: launch worker threads for searching and aligning profiles
//
void JobDispatcher::Run()
{
    MYMSG( "JobDispatcher::Run", 3 );
    const mystring preamb = "JobDispatcher::Run: ";

    const int cpufreemem = CLOptions::GetCPU_FREEMEM();

    myruntime_error mre;
    char msgbuf[BUF_MAX];
    const size_t initquerydatasize = 100UL*ONEM;
    size_t maxdatasize = GetMaxCacheSize();
    size_t chunkdatasize = 0UL;
    size_t chunkdatalen = 0UL;//maximum total length of db profile positions that can be processed at once
    size_t chunknpros = 0UL;//maximum number of db profiles allowed to be processed at once
    size_t lenpm1, lenpmcum, npros;//runnning length and the total length of profile data and #profiles
    size_t szpm1, szpmcum;//runnning size and the accumulated size of profile data
    //access to input batch of profiles, which also can be a database:
    CuInputReader inpdb(GetInput(), !CLOptions::GetIO_NOFILEMAP());
    PMBatchProData bquery;//batch of input profile data
    PMBatchProData bdb1;//the first and ...
    std::unique_ptr<PMBatchProData> bdbC;//a running chunk read from the database 
    TCharStream stpat;
    bool bdb1complete = false;//whether bdb1 contains complete db profile data
    bool noteof = true;
    int prolen = 0, scale, profnr = 0; 
    unsigned int distn = 0;//distance in positions to a profile
    mystring  pdesc, pfile;
    const size_t nflds = bdb1.GetNoFields();
    const size_t nheadflds = bdb1.GetNoHeaderFields();
    char auxpsdat[nheadflds][SZFPTYPE];//auxiliary buffer for profile header, i.e., profile-specific data
    char* auxps[nheadflds];//pointers to the above structure
    char* querypmdata[nflds];//running pointers for queries
    char* querypmprcd[nflds];//pointers at the boundary of the last processed query data
    char* db1pmdata[nflds];//running pointers of db1 data
    char* db1pmprcd[nflds];//pointers at the boundary of the last processed db1 data
    //*test*/char* prevcopypmdata[nflds];
    int tid = JDSP_WORKER_BUSY;//worker id

    ReadConfiguration();

    for( int i = 0; i < (int)nheadflds; i++ )
        auxps[i] = auxpsdat[i];


    if( DEVPROPs.GetNDevices() < 1 ) {
        warning("There is no available device to run the program on.");
        message("Please use a version to run on CPU.");
        return;
    }


    try {

        MYMSG( "Reading queries", 1 );
        ReadProfiles( inpdb, bquery, initquerydatasize, true/*readall*/, NULL/*stoppos*/);
        inpdb.Destroy();//explicitly destroy unused resources here

        MYMSG( "Reading database profiles to be cached", 1 );
        ReadProfiles( prodb_, bdb1, maxdatasize, false/*readall*/, &stpat );

        size_t nqueries = PMBatchProData::GetNoProsFromTo( bquery.GetPMData(), bquery.GetPMDataEnd());
        bool querinf = PMBatchProData::PMDataPresentFromTo( bquery.GetPMData(), bquery.GetPMDataEnd());
        bool bdb1inf = PMBatchProData::PMDataPresentFromTo( bdb1.GetPMData(), bdb1.GetPMDataEnd());


        if( !querinf || nqueries < 1 )
            throw MYRUNTIME_ERROR( preamb + "Query data has not been cached." );


        //create the results writing thread
        CreateAlnWriter( 
            GetConfiguration(),
            output_,
            prodb_.GetDbName(),
            prodb_.GetDbSize(),
            prodb_.GetNoSequences(),
            (int)nqueries );

        //create worker threads
        CreateWorkerThreads(
                bquery.GetFnames(), bquery.GetDescs(),
                bquery.GetPMData(), bquery.GetPMDataEnd(),
                bdb1inf? bdb1.GetFnames(): NULL,
                bdb1inf? bdb1.GetDescs(): NULL,
                bdb1inf? bdb1.GetPMData(): NULL,
                bdb1inf? bdb1.GetPMDataEnd(): NULL,
                prodb_.GetDbSize(), prodb_.GetNoSequences());

        if( hostworkers_.size() < 1 )
            throw MYRUNTIME_ERROR(preamb+"Failed to create any worker threads.");


        if( prodb_.GetNoSequences() <= bdb1.GetNoProsWritten())
            bdb1complete = true;

        if( !bdb1complete )
            prodb_.Open();

        bquery.GetPMData(querypmdata);
        memcpy( querypmprcd, querypmdata, nflds * sizeof(void*));


        for( int qrysernr = 0;
            !bquery.PMDataReachedEnd(querypmprcd) && tid != JDSP_WORKER_NONE;
            GetAlnWriter()->DereaseNPartsAndTrigger(qrysernr),
            qrysernr++,
            memcpy( querypmprcd, querypmdata, nflds * sizeof(void*)))
        {
            MYMSGBEGl(1)
                char msgbuf[BUF_MAX];
                size_t qbndx = bquery.GetNoProsAt(querypmprcd);
                sprintf(msgbuf, "Processing QUERY %d (%s)",
                        qrysernr, bquery.GetFnames()[qbndx].c_str());
                MYMSG( msgbuf, 1 );
            MYMSGENDl

            bquery.PMDataNextPro(querypmdata);

            size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmprcd, querypmdata );
            CalculateScoreThreshold( nqyposs, prodb_.GetDbSize());


            std::vector<size_t> params {nqyposs};

            //broadcast message to workers and get the size of data chunk
            BcastMessage(
                true,//waitforresponse
                DevCommThread::tthreadmsgGetDataChunkSize, 
                -1,//addressee: broadcast to all threads
                &JobDispatcher::WriteForMsgGetDataChunkSize, params,
                &JobDispatcher::ReadOnResponseOfMsgGetDataChunkSize, params );

            WaitForAvailableWorker(&tid, &chunkdatasize, &chunkdatalen, &chunknpros);


            szpmcum = 0;
            lenpmcum = 0;
            npros = 0;
            bdb1.GetPMData(db1pmdata);
            memcpy( db1pmprcd, db1pmdata, nflds * sizeof(void*));

            for(; !bdb1.PMDataReachedEnd(db1pmdata) && tid != JDSP_WORKER_NONE; 
                bdb1.PMDataNextPro(db1pmdata))
            {
                szpm1 = bdb1.GetPMDataSize1At(db1pmdata);
                lenpm1 = bdb1.GetPMDataLen1At(db1pmdata);
                if( chunkdatasize < szpmcum + szpm1 || 
                    chunkdatalen < lenpmcum + lenpm1 || chunknpros < npros + 1 ) {
                    //NOTE: COMPUTATION HERE!
                    ProcessPart(
                        qrysernr, prodb_.GetDbSize(),
                        &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                        &bquery, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata,
                        std::move(bdbC));
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                    memcpy( db1pmprcd, db1pmdata, nflds * sizeof(void*));
                }
                szpmcum += szpm1;
                lenpmcum += lenpm1;
                npros++;
            }

            if(tid == JDSP_WORKER_NONE)
                break;

            if( bdb1complete ) {
                if( szpmcum ) {
                    //NOTE: COMPUTATION HERE!
                    ProcessPart(
                        qrysernr, prodb_.GetDbSize(),
                        &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                        &bquery, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata,
                        std::move(bdbC));
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                }
                continue;
            }

            MYMSG( "Resuming database scan", 1 );

            prodb_.SetDbPos( &stpat );

            distn = 0;
            profnr = 0;


            //NOTE: the usual case is that space is allocated for two data 
            // chunks at a time: one is being processed while the other is 
            // put in a queue;
            // data chunks grow in size as the database is sorted by profile 
            // length;
            // therefore, the requirement for memory resources consistently, 
            // although inconsiderably, grow too;
            bdbC.reset(new PMBatchProData);
            if( !bdbC )
                throw MYRUNTIME_ERROR("JobDispatcher::Run: Not enough memory");
            bdbC->AllocNewSizeHostDataInit( chunkdatasize - szpmcum );

            for(; !prodb_.Eof() && tid != JDSP_WORKER_NONE;) 
            {
                noteof = prodb_.NextProfileHeader( &pdesc, &pfile, &prolen, &scale, auxps );
                lenpm1 = prolen;

                if( !noteof ) {
                    sprintf( msgbuf, "Unexpected end of file %d", profnr+1 );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                //avoid overflow
                bool overflow = ( INT_MAX <= (unsigned int)prolen ||
                        (unsigned int)(INT_MAX - prolen) <= distn );

                if( overflow ||
                    bdbC->MemDeficit( prolen ) ||
                    chunkdatalen < lenpmcum + lenpm1 || chunknpros < npros + 1 ) {
                    //NOTE: free allocated but unused memory: use the required size:
                    if( cpufreemem )
                        bdbC->ReallocHostData();
                    //NOTE: COMPUTATION HERE!
                    if( szpmcum )
                        ProcessPart(
                            qrysernr, prodb_.GetDbSize(),
                            &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                            &bquery, querypmprcd, querypmdata, 
                            &bdb1, db1pmprcd, db1pmdata,
                            std::move(bdbC), bdbC->GetPMData(), bdbC->GetPMDataEnd());
                    else
                        ProcessPart(
                            qrysernr, prodb_.GetDbSize(),
                            &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                            &bquery, querypmprcd, querypmdata,
                            NULL, NULL, NULL,
                            std::move(bdbC), bdbC->GetPMData(), bdbC->GetPMDataEnd());
                    distn = 0;/**/
                    profnr = 0;/**/
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                    bdbC.reset(new PMBatchProData);
                    if( !bdbC )
                        throw MYRUNTIME_ERROR("JobDispatcher::Run: Not enough memory");
                    bdbC->AllocNewSizeHostData( chunkdatasize );
                }

                lenpmcum += lenpm1;
                npros++;

                //{{save data read from the profile header
                //save the beginning of a next profile
                //*test*/memcpy( prevcopypmdata, bdbC.GetPMDataEnd(), nflds * sizeof(void*));
                bdbC->CopyPMDataHeaderNext( pdesc, pfile, auxps );
                //}}

                noteof = prodb_.Next( distn, profnr, prolen, scale, bdbC->GetPMData(), bdbC->GetPMDataEnd());

                distn += prolen;
                profnr++;

                //test:
                //*test*/TextWriteProfileHD( stderr, pdesc, pfile, scale, bdbC.GetPMData(), prevcopypmdata );

            }//for(;!myeof;)

            if( bdbC->GetNoProsWritten() && tid != JDSP_WORKER_NONE) {
                //NOTE: free allocated but unused memory: use the required size:
                if( cpufreemem )
                    bdbC->ReallocHostData();
                //NOTE: COMPUTATION HERE!
                if( szpmcum )
                    ProcessPart(
                        qrysernr, prodb_.GetDbSize(),
                        &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                        &bquery, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata,
                        std::move(bdbC), bdbC->GetPMData(), bdbC->GetPMDataEnd());
                else
                    ProcessPart(
                        qrysernr, prodb_.GetDbSize(),
                        &tid, &chunkdatasize, &chunkdatalen, &chunknpros,
                        &bquery, querypmprcd, querypmdata, 
                        NULL, NULL, NULL,
                        std::move(bdbC), bdbC->GetPMData(), bdbC->GetPMDataEnd());
                szpmcum = 0;
                lenpmcum = 0;
                npros = 0;
            }

        }//for(;!bquery.PMDataReachedEnd(querypmprcd);)

        WaitForAllWorkersToFinish();

    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Exception caught.");
    }

    prodb_.Close();

    std::vector<int> dummy;
    BcastMessage(
        true,//waitforresponse
        DevCommThread::tthreadmsgTerminate, 
        -1,//addressee: broadcast to all threads
        &JobDispatcher::VoidWriteBeforeSending, dummy,
        &JobDispatcher::VoidReadOnReceive, dummy );

    if( mre.isset())
        throw mre;

    MYMSG( "Done.", 1 );
}

// -------------------------------------------------------------------------
// ProcessPart: process a part of data, a number of databases profiles for a
// number of queries
//
void JobDispatcher::ProcessPart(
    int qrysernr, size_t dbsize,
    int* tid, size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros,
    const PMBatchProData* query, char** querypmbeg, char** querypmend, 
    const PMBatchProData* bdb1, char** bdb1pmbeg, char** bdb1pmend,
    std::unique_ptr<PMBatchProData> bdbC, char** bdbCpmbeg, char** bdbCpmend )
{
    MYMSG( "JobDispatcher::ProcessPart", 4 );
    const mystring preamb = "JobDispatcher::ProcessPart: ";

#ifdef __DEBUG__
    if( !query || !querypmbeg || !querypmend )
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if( bdb1 && ( !bdb1pmbeg || !bdb1pmend ))
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if( bdbC && ( !bdbCpmbeg || !bdbCpmend ))
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
#endif

    bool querinf = PMBatchProData::PMDataPresentFromTo( querypmbeg, querypmend );
    bool bdb1inf = PMBatchProData::PMDataPresentFromTo( bdb1pmbeg, bdb1pmend );
    bool bdbCinf = PMBatchProData::PMDataPresentFromTo( bdbCpmbeg, bdbCpmend );

    if( !querinf )
        throw MYRUNTIME_ERROR( preamb + "Query data is empty." );
    if( !bdb1inf && !bdbCinf )
        throw MYRUNTIME_ERROR( preamb + "Db profile data is empty." );

    MYMSGBEGl(1)
        char msgbuf[BUF_MAX];
        size_t nposits = 0;
        if( bdb1 && bdb1pmbeg && bdb1pmend )
            nposits += bdb1->GetNoPositsFromTo(bdb1pmbeg, bdb1pmend);
        if( bdbC && bdbCpmbeg && bdbCpmend )
            nposits += bdbC->GetNoPositsFromTo(bdbCpmbeg, bdbCpmend);
        sprintf(msgbuf, "   %d%% db positions submitted for processing",
                (int)rintf((float)nposits*100.0f/(float)dbsize));
        MYMSG( msgbuf, 1 );
    MYMSGENDl

    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        mystring strbuf = preamb;
        sprintf( msgbuf, "%s[=====] Processing QUERY(-ies) ",NL);
        strbuf += msgbuf;
        size_t qbndx = query->GetNoProsAt( querypmbeg );
        size_t qendx = query->GetNoProsAt( querypmend );
        if( qbndx + 1 >= qendx ) {
            sprintf( msgbuf, "%zu (", qbndx );
            strbuf += msgbuf;
            strbuf += query->GetFnames()[qbndx] + "): ";
        } else {
            sprintf( msgbuf, "%zu-%zu (", qbndx, qendx-1 );
            strbuf += msgbuf;
            strbuf += query->GetFnames()[qbndx] + "-" + query->GetFnames()[qendx-1] + "): ";
        }
        if( bdb1 ) {
            size_t b1bndx = bdb1->GetNoProsAt( bdb1pmbeg );
            size_t b1endx = bdb1->GetNoProsAt( bdb1pmend );
            if( b1bndx < b1endx ) {
                strbuf += NL;
                strbuf += " cached db pros ";
                sprintf( msgbuf, "%zu-%zu (", b1bndx, b1endx-1 );
                strbuf += msgbuf;
                strbuf += bdb1->GetFnames()[b1bndx] + "-" + bdb1->GetFnames()[b1endx-1] + ")";
            }
        }
        if( bdbC ) {
            size_t bCbndx = bdbC->GetNoProsAt( bdbCpmbeg );
            size_t bCendx = bdbC->GetNoProsAt( bdbCpmend );
            if( bCbndx < bCendx ) {
                sprintf( msgbuf, "%s + %zu resumed db pros (", NL, bCendx );
                strbuf += msgbuf;
                strbuf += bdbC->GetFnames()[bCbndx] + "-" + bdbC->GetFnames()[bCendx-1] + ")";
            }
        }
        MYMSG( strbuf.c_str(), 3 );
    MYMSGENDl


    GetAlnWriter()->IncreaseQueryNParts(qrysernr);

    std::vector<char**> params {
        (char**)(size_t)qrysernr,
        //the responsibility of bdbC to be delegated to the worker
        querypmbeg, querypmend,
        bdb1inf? bdb1pmbeg: NULL, bdb1inf? bdb1pmend: NULL,
        bdbCinf? bdbCpmbeg: NULL, bdbCinf? bdbCpmend: NULL
    };
    //send message to a dedicated worker...
    BcastMessage(
        false,//waitforresponse
        DevCommThread::tthreadmsgProcessNewData, 
        *tid,//one addressee
        &JobDispatcher::WriteForMsgProcessNewData, params,
        &JobDispatcher::VoidReadOnReceive, params/*unused*/,
        std::move(bdbC)
    );
    //and get the index of a next available worker
    //WaitForAvailableWorker(tid, chunkdatasize, chunkdatalen, chunknpros);
    //NOTE: do not block and select the next worker to 
    // proceed this thread with preparing next data portion
    GetNextWorker(tid, chunkdatasize, chunkdatalen, chunknpros);
}





#if 0
// -------------------------------------------------------------------------
// Run: launch worker threads for searching and aligning profiles
//
void JobDispatcher::Run_obs()
{
    MYMSG( "JobDispatcher::Run", 3 );
    const mystring preamb = "JobDispatcher::Run: ";

    myruntime_error mre;
    char msgbuf[BUF_MAX];
    const size_t initdatasize = 100L*ONEM;
    size_t maxdatasize = 220UL*ONEM;//220L*ONEM;
    size_t chunkdatasize = 0UL;
    size_t chunkdatalen = 0UL;//maximum total length of db profile positions that can be processed at once
    size_t chunknpros = 0UL;//maximum number of db profiles allowed to be processed at once
    size_t lenpm1, lenpmcum, npros;//runnning length and the total length of profile data and #profiles
    size_t szpm1, szpmcum;//runnning size and the accumulated size of profile data
    PMBatchProData bdb1, bdbC;//the first and a running chunk read from the database 
    fpos_t stpat;
    bool bdb1complete = false;//whether bdb1 contains complete db profile data
    bool noteof = true;
    int prolen = 0, scale, profnr = 0; 
    unsigned int distn = 0;//distance in positions to a profile
    mystring  pdesc, pfile;
    const size_t nflds = bdbC.GetNoFields();
    const size_t nheadflds = bdbC.GetNoHeaderFields();
    char auxpsdat[nheadflds][SZFPTYPE];//auxiliary buffer for profile header, i.e., profile-specific data
    char* auxps[nheadflds];//pointers to the above structure
    char* querypmdata[nflds];//running pointers for queries
    char* querypmprcd[nflds];//pointers at the boundary of the last processed query data
    char* db1pmdata[nflds];//running pointers of db1 data
    char* db1pmprcd[nflds];//pointers at the boundary of the last processed db1 data
    //*test*/char* prevcopypmdata[nflds];

    ReadConfiguration();

    //nitialize object after creating configuration
    CuBatchProcessing cbpc( 
        GetConfiguration(), 
        HDPSCORES.GetScores()!=NULL,
        MOptions::GetMINPP()<1.0f,
        4000UL*ONEM
    );

    cbpc.CacheSSENNWeights();
    cbpc.CacheSSSScores( SSSSCORES );
    cbpc.CacheCVS2Scores( CVS2SCORES );
    cbpc.CacheHDPScores( HDPSCORES );

    for( int i = 0; i < (int)nheadflds; i++ )
        auxps[i] = auxpsdat[i];

    try {

        ReadProfiles( inpdb_, bquery_, initdatasize, true/*readall*/, NULL/*stoppos*/);

        ReadProfiles( prodb_, bdb1, maxdatasize, false/*readall*/, &stpat );

        bool querinf = PMBatchProData::PMDataPresentFromTo( bquery_.GetPMData(), bquery_.GetPMDataEnd());
        bool bdb1inf = PMBatchProData::PMDataPresentFromTo( bdb1.GetPMData(), bdb1.GetPMDataEnd());

        if( !querinf )
            throw MYRUNTIME_ERROR( preamb + "Query data is not cached. Increase memory limits." );

        cbpc.CacheData(
            bquery_.GetPMData(), bquery_.GetPMDataEnd(),
            bdb1inf? bdb1.GetPMData(): NULL, bdb1inf? bdb1.GetPMDataEnd(): NULL );

        if( prodb_.GetNoSequences() <= bdb1.GetNoProsWritten())
            bdb1complete = true;

        if( !bdb1complete )
            prodb_.Open();

        bquery_.GetPMData(querypmdata);
        memcpy( querypmprcd, querypmdata, nflds * sizeof(void*));

// int kk=0;
        for( ; !bquery_.PMDataReachedEnd(querypmprcd); 
            memcpy( querypmprcd, querypmdata, nflds * sizeof(void*)))
        {
            MYMSG( "JobDispatcher::Run: *** New query ***", 3 );

            bquery_.PMDataNextPro(querypmdata);
            //test query:
            //*test*/TextWriteProfileHD( stderr, pdesc="", pfile="", 
            //            pmodel::PMProfileModel::ReadWriteScale, 
            //            bquery_.GetPMData(), querypmdata );

            size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmprcd, querypmdata );
            CalculateScoreThreshold( nqyposs, prodb_.GetDbSize());

            cbpc.SetDbDetails( nqyposs, prodb_.GetDbSize(), prodb_.GetNoSequences());

            chunkdatasize = cbpc.CalcMaxDbDataChunkSize( nqyposs );
            chunkdatalen = cbpc.GetCurrentMaxDbPos();
            chunknpros = cbpc.GetCurrentMaxNDbPros();

            szpmcum = 0;
            lenpmcum = 0;
            npros = 0;
            bdb1.GetPMData(db1pmdata);
            memcpy( db1pmprcd, db1pmdata, nflds * sizeof(void*));
    
            for( ; !bdb1.PMDataReachedEnd(db1pmdata); bdb1.PMDataNextPro(db1pmdata))
            {
                szpm1 = bdb1.GetPMDataSize1At(db1pmdata);
                lenpm1 = bdb1.GetPMDataLen1At(db1pmdata);
                if( chunkdatasize < szpmcum + szpm1 || 
                    chunkdatalen < lenpmcum + lenpm1 || chunknpros < npros + 1 ) {
                    //NOTE: COMPUTATION HERE!
                    ProcessPart_obs( 
                        cbpc,
                        &bquery_, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata );
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                    memcpy( db1pmprcd, db1pmdata, nflds * sizeof(void*));
// break;
                }
                szpmcum += szpm1;
                lenpmcum += lenpm1;
                npros++;
            }
// if(++kk>1/*180*/) break;

            if( bdb1complete ) {
                if( szpmcum ) {
                    //NOTE: COMPUTATION HERE!
                    ProcessPart_obs( 
                        cbpc,
                        &bquery_, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata );
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                }
                continue;
            }

            MYMSG( "JobDispatcher::Run: Resuming db scan", 3 );

            prodb_.SetDbPos( &stpat );

            distn = 0;
            profnr = 0;

            bdbC.AllocNewSizeHostDataInit( chunkdatasize - szpmcum );

            for(; !prodb_.Eof();) 
            {
                noteof = prodb_.NextProfileHeader( &pdesc, &pfile, &prolen, &scale, auxps );
                lenpm1 = prolen;

                if( !noteof ) {
                    sprintf( msgbuf, "Unexpected end of file %d", profnr+1 );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if( bdbC.MemDeficit( prolen ) || 
                    chunkdatalen < lenpmcum + lenpm1 || chunknpros < npros + 1 ) {
                    //NOTE: COMPUTATION HERE!
                    if( szpmcum )
                        ProcessPart_obs(
                            cbpc,
                            &bquery_, querypmprcd, querypmdata, 
                            &bdb1, db1pmprcd, db1pmdata,
                            &bdbC, bdbC.GetPMData(), bdbC.GetPMDataEnd());
                    else
                        ProcessPart_obs(
                            cbpc,
                            &bquery_, querypmprcd, querypmdata,
                            NULL, NULL, NULL,
                            &bdbC, bdbC.GetPMData(), bdbC.GetPMDataEnd());
                    distn = 0;/**/
                    profnr = 0;/**/
                    szpmcum = 0;
                    lenpmcum = 0;
                    npros = 0;
                    bdbC.AllocNewSizeHostData( chunkdatasize );
                }

                lenpmcum += lenpm1;
                npros++;

                //{{save data read from the profile header
                //save the beginning of a next profile
                //*test*/memcpy( prevcopypmdata, bdbC.GetPMDataEnd(), nflds * sizeof(void*));
                bdbC.CopyPMDataHeaderNext( pdesc, pfile, auxps );
                //}}

                noteof = prodb_.Next( distn, profnr, prolen, scale, bdbC.GetPMData(), bdbC.GetPMDataEnd());

                distn += prolen;
                profnr++;

                //test:
                //*test*/TextWriteProfileHD( stderr, pdesc, pfile, scale, bdbC.GetPMData(), prevcopypmdata );

            }//for(;!myeof;)

            if( bdbC.GetNoProsWritten()) {
                //NOTE: COMPUTATION HERE!
                if( szpmcum )
                    ProcessPart_obs(
                        cbpc,
                        &bquery_, querypmprcd, querypmdata, 
                        &bdb1, db1pmprcd, db1pmdata,
                        &bdbC, bdbC.GetPMData(), bdbC.GetPMDataEnd());
                else
                    ProcessPart_obs(
                        cbpc,
                        &bquery_, querypmprcd, querypmdata, 
                        NULL, NULL, NULL,
                        &bdbC, bdbC.GetPMData(), bdbC.GetPMDataEnd());
                szpmcum = 0;
                lenpmcum = 0;
                npros = 0;
            }

        }//for(;!bquery_.PMDataReachedEnd(querypmdata);)

    } catch( myexception const& ex ) {
        mre = ex;
    }

    prodb_.Close();

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// ProcessPart: process a part of data, a number of databases profiles for a
// number of queries
//
void JobDispatcher::ProcessPart_obs( 
    CuBatchProcessing& cbpc,
    const PMBatchProData* query, char** querypmbeg, char** querypmend, 
    const PMBatchProData* bdb1, char** bdb1pmbeg, char** bdb1pmend,
    const PMBatchProData* bdbC, char** bdbCpmbeg, char** bdbCpmend )
{
    MYMSG( "JobDispatcher::ProcessPart", 4 );
    const mystring preamb = "JobDispatcher::ProcessPart: ";

#ifdef __DEBUG__
    if( !query || !querypmbeg || !querypmend )
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if( bdb1 && ( !bdb1pmbeg || !bdb1pmend ))
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
    if( bdbC && ( !bdbCpmbeg || !bdbCpmend ))
        throw MYRUNTIME_ERROR( preamb + "Null arguments." );
#endif

    bool querinf = PMBatchProData::PMDataPresentFromTo( querypmbeg, querypmend );
    bool bdb1inf = PMBatchProData::PMDataPresentFromTo( bdb1pmbeg, bdb1pmend );
    bool bdbCinf = PMBatchProData::PMDataPresentFromTo( bdbCpmbeg, bdbCpmend );

    if( !querinf )
        throw MYRUNTIME_ERROR( preamb + "Query data is empty." );
    if( !bdb1inf && !bdbCinf )
        throw MYRUNTIME_ERROR( preamb + "Db profile data is empty." );

    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        mystring strbuf = preamb + "Processing query(-ies) ";
        size_t qbndx = query->GetNoProsAt( querypmbeg );
        size_t qendx = query->GetNoProsAt( querypmend );
        if( qbndx + 1 >= qendx ) {
            sprintf( msgbuf, "%lu (", qbndx );
            strbuf += msgbuf;
            strbuf += query->GetFnames()[qbndx] + "): ";
        } else {
            sprintf( msgbuf, "%lu-%lu (", qbndx, qendx-1 );
            strbuf += msgbuf;
            strbuf += query->GetFnames()[qbndx] + "-" + query->GetFnames()[qendx-1] + "): ";
        }
        if( bdb1 ) {
            size_t b1bndx = bdb1->GetNoProsAt( bdb1pmbeg );
            size_t b1endx = bdb1->GetNoProsAt( bdb1pmend );
            if( b1bndx < b1endx ) {
                strbuf += "cached db pros ";
                sprintf( msgbuf, "%lu-%lu (", b1bndx, b1endx-1 );
                strbuf += msgbuf;
                strbuf += bdb1->GetFnames()[b1bndx] + "-" + bdb1->GetFnames()[b1endx-1] + ")";
            }
        }
        if( bdbC ) {
            size_t bCbndx = bdbC->GetNoProsAt( bdbCpmbeg );
            size_t bCendx = bdbC->GetNoProsAt( bdbCpmend );
            if( bCbndx < bCendx ) {
                strbuf += " + # resumed db pros ";
                sprintf( msgbuf, "%lu (", bCendx );
                strbuf += msgbuf;
                strbuf += bdbC->GetFnames()[bCbndx] + "-" + bdbC->GetFnames()[bCendx-1] + ")";
            }
        }
        MYMSG( strbuf.c_str(), 3 );
    MYMSGENDl

    cbpc.SetScoreThreshold( GetScoreThreshold());
    cbpc.SetLogEThreshold( GetLogEThreshold());
    cbpc.ProcessScoreMatrix(
        querypmbeg, querypmend, 
        bdb1inf? bdb1pmbeg: NULL, bdb1inf? bdb1pmend: NULL,
        bdbCinf? bdbCpmbeg: NULL, bdbCinf? bdbCpmend: NULL
    );
}
#endif
