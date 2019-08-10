/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __JobDispatcher_h__
#define __JobDispatcher_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include <memory>
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
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cualn/DevCommThread.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cualn/Devices.h"

#define JDSP_WORKER_BUSY -1
#define JDSP_WORKER_NONE -2

class CuBatchProcessing;

class JobDispatcher;

template<typename T>
using TWriteDataForWorker1 = 
        void (JobDispatcher::*)(int req, int addr, DevCommThread*, const std::vector<T>&,
            std::unique_ptr<PMBatchProData> );
template<typename T>
using TGetDataFromWorker1 = 
        void (JobDispatcher::*)(int req, int rsp, int addr, DevCommThread*, std::vector<T>& );

// _________________________________________________________________________
// Class JobDispatcher
//
// Implementation of distributing jobs over CPU and GPU threads
//
class JobDispatcher
{
public:
    JobDispatcher(
            const char* configfile,
            const char* input,
            const char* database,
            const char* output,
            float eval_thld,
            int no_hits,
            int no_alns,
            bool showpars
    );
    ~JobDispatcher();

    void            Run();
    void            Run_obs();

    const char*     GetConfigFile() const { return configFile_; }
    const char*     GetInput() const            { return input_; }
    const char*     GetDatabase() const         { return database_; }
    const char*     GetOutput() const           { return output_; }

    double          GetEvalueUpper() const  { return eval_upper_; }
    int             GetMaxNoHits() const  { return max_no_hits_; }
    int             GetMaxNoAlns() const  { return max_no_alns_; }
    bool            GetShowPars() const          { return show_pars_; }

    void            SetTarFrMix( int value ) { tfrmix_ = value; }
    bool            GetTarFrMixHDPCtx() const { return tfrmix_ == tfrmixHDPCtx; }

    void            SetScoAdjment( int value ) { scoadj_ = value; }
    bool            GetScoAdjmentHDPCtx() const { return scoadj_ == scoadjHDPCtx || scoadj_ == scoadjHDPsco; }

    void            SetHDPbase( const HDPbase* value ) { HDPbase_ = value; }
    const HDPbase*  GetHDPbase() const { return HDPbase_; }

//     void            SetHDPctbase( const HDPbase* value ) { HDPctbase_ = value; }
//     const HDPbase*  GetHDPctbase() const { return HDPctbase_; }


    float           GetScoreThreshold() const { return scorethld_; }
    void            SetScoreThreshold( float value ) { scorethld_ = value; }

    float           GetLogEThreshold() const { return logevthld_; }
    void            SetLogEThreshold( float value ) { logevthld_ = value; }


    void            SetSSEModel( int modndx ) { ssemodel_ = modndx; }
    int             GetSSEModel() const { return ssemodel_; }


    int             GetHSPLength() const                        { return hsplength_; }
    void            SetHSPLength( int value )                   { hsplength_ = value; }

    int             GetHSPScore() const                         { return hspscore_; }
    void            SetHSPScore( int value )                    { hspscore_ = value; }

    int             GetHSPDistance() const                      { return hspdistance_; }
    void            SetHSPDistance( int value )                 { hspdistance_ = value; }

    int             GetHSPNoHSPs() const                        { return hspnohsps_; }
    void            SetHSPNoHSPs( int value )                   { hspnohsps_ = value; }

    void            PrintMethodName( FILE* fp ) const;//scoring method name
    void            PrintParameterTable( FILE* ) const;

protected:
    explicit JobDispatcher();


    size_t GetMaxCacheSize();

    void ReadConfiguration();

    void CreateAlnWriter( 
        Configuration* config,
        const char* outdirname,
        const char* dbname,
        size_t prodbsize,
        size_t ndbentries,
        int nqueries
    );
    void NotifyAlnWriter();

    void CreateWorkerThreads(
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

    void BlockAndNotifyWriter() {
        WaitForAllWorkersToFinish();
        NotifyAlnWriter();
    }

    template<typename TI, typename TO>
    void BcastMessage(
        bool waitforresponse,
        int msgcode, int addressee,
        TWriteDataForWorker1<TI> writefunc,
        const std::vector<TI>&,
        TGetDataFromWorker1<TO> readfunc,
        std::vector<TO>&,
        std::unique_ptr<PMBatchProData> = nullptr
    );
    void GetAvailableWorker( 
        int* tid, size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros );
    void GetNextWorker( 
        int* tid, size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros );
    void WaitForAvailableWorker( 
        int* tid,
        size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros );

    void WaitForAllWorkersToFinish();

    void WriteForMsgGetDataChunkSize(
        int reqmsg, int addressee, DevCommThread* tdc, const std::vector<size_t>& params,
        std::unique_ptr<PMBatchProData> );
    void ReadOnResponseOfMsgGetDataChunkSize(
        int reqmsg, int rspmsg, int addressee, DevCommThread* tdc, std::vector<size_t>& data );

    void WriteForMsgProcessNewData(
        int reqmsg, int addressee, DevCommThread* tdc, const std::vector<char**>& params,
        std::unique_ptr<PMBatchProData> );

    template<typename T>
    void ReadOnResponseOfMsgProbe(
        int reqmsg, int rspmsg, int addressee, DevCommThread* tdc, std::vector<T>& );

    template<typename T>
    void VoidWriteBeforeSending( int reqmsg, int addressee, DevCommThread*, const std::vector<T>&,
        std::unique_ptr<PMBatchProData> );
    template<typename T>
    void VoidReadOnReceive( int reqmsg, int rspmsg, int addressee, DevCommThread*, std::vector<T>& );



    void    CalculateScoreThreshold( size_t nqyposs, size_t dbsize );
    void    ReadProfiles( 
        CuDbReader& db, PMBatchProData& bdata, 
        size_t maxdatasize, bool readall, TCharStream* stoppos );

    void    ProcessPart(
        int qrysernr, size_t dbsize,
        int* tid, size_t* chunkdatasize, size_t* chunkdatalen, size_t* chunknpros,
        const PMBatchProData*, char** querypmbeg, char** querypmend, 
        const PMBatchProData*, char** bdb1pmbeg, char** bdb1pmend,
        std::unique_ptr<PMBatchProData>, char** bdbCpmbeg = NULL, char** bdbCpmend = NULL );
    void    ProcessPart_obs(
        CuBatchProcessing& cbpc,
        const PMBatchProData*, char** querypmbeg, char** querypmend, 
        const PMBatchProData*, char** bdb1pmbeg, char** bdb1pmend,
        const PMBatchProData* = NULL, char** bdbCpmbeg = NULL, char** bdbCpmend = NULL );



    const AlnWriter* GetAlnWriter() const {return writer_;}
    AlnWriter* GetAlnWriter() {return writer_;}

    const Configuration&    GetConfiguration( TConfigType ) const;
    Configuration&          GetConfiguration( TConfigType );
    Configuration*          GetConfiguration() { return config_; }



    size_t      GetNoSequences() const { return db_no_seqs_; }
    void        SetNoSequences( size_t value ) { db_no_seqs_ = value; }

    uint64_mt   GetDbSize() const { return db_size_; }
    void        SetDbSize( uint64_mt value ) { db_size_ = value; }

private:
    const char*             configFile_;//configuration file of statistical parameters
    const char*             input_;//input's filename
    const char*             database_;//profile database
    const char*             output_;//pattern for output file (null=standard output)
//     AbstractScoreMatrix<>*    scoreSystem_;//profile scoring system
    Configuration           config_[NoCTypes];//configurations
    std::vector<DevCommThread*> hostworkers_;//host worker threads
    AlnWriter*              writer_;//alignment results writer

    float                   scorethld_;//score threshold given e-value in the options
    float                   logevthld_;//log E-value threshold

    float                   eval_upper_;//e-value upper threshold
    int                     max_no_hits_;//maximum number of hits to show in the result list
    int                     max_no_alns_;//maximum number of alignments to show in the output
    bool                    show_pars_;//whether to show statistical parameters below alignments

    CuDbReader              prodb_;//profile database
    size_t                  db_no_seqs_;//number of profiles in the database
    uint64_mt               db_size_;//size of the database in positions

    int                     tfrmix_;//mixing of target frequencies
    int                     scoadj_;////score adjustment
    const HDPbase*          HDPbase_;//HDP base data structure
//     const HDPbase*          HDPctbase_;//HDP ctx base structure

    int                     ssemodel_;//index of a model for the estimation of stat. sign.

    int                     hsplength_;//minimum length for HSPs
    int                     hspscore_;//minimum HSP score
    int                     hspdistance_;//maximum distance between the HSPs
    int                     hspnohsps_;//minimum number of HSPs in a diagonal
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// GetConfiguration: returns reference to parameter configuration object
//
inline
const Configuration& JobDispatcher::GetConfiguration( TConfigType ct ) const
{
#ifdef __DEBUG__
    if( NoCTypes <= ct )
        throw MYRUNTIME_ERROR( "JobDispatcher: Wrong configuration index." );
#endif
    return config_[ct];
}
inline
Configuration& JobDispatcher::GetConfiguration( TConfigType ct )
{
#ifdef __DEBUG__
    if( NoCTypes <= ct )
        throw MYRUNTIME_ERROR( "JobDispatcher: Wrong configuration index." );
#endif
    return config_[ct];
}

// -------------------------------------------------------------------------
// GetMaxCacheSize: get the cache size for caching db profile data to GPU
inline
size_t JobDispatcher::GetMaxCacheSize()
{
    MYMSG( "JobDispatcher::GetMaxCacheSize", 5 );
    const mystring preamb = "JobDispatcher::GetMaxCacheSize: ";

    const float cacheperc = CLOptions::GetDEV_CACHEP();
    size_t szcache = 0;

    for(int tid = 0; tid < DEVPROPs.GetNDevices(); tid++) {
        const DeviceProperties* d1prop = DEVPROPs.GetDevicePropertiesAt(tid);
        if( d1prop== NULL )
            continue;
        szcache = (size_t)((float)d1prop->reqmem_ * cacheperc);
        break;
    }

    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        mystring strbuf = preamb;
        sprintf( msgbuf, " %zu (%zuMB)",szcache,szcache/ONEM);
        strbuf += msgbuf;
        MYMSG( strbuf.c_str(), 3 );
    MYMSGENDl

    return szcache;
}

#endif//__JobDispatcher_h__
