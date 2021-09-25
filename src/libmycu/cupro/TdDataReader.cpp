/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <algorithm>
#include <thread>

#include <cuda_runtime_api.h>

#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/CuRoDb.h"
#include "libmycu/cupro/VirtualCuRoDb.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"

#include "TdDataReader.h"

const std::vector<std::string> dummyvec;

// _________________________________________________________________________
// statics
//
const size_t TdDataReader::s_szdatalign_ = CUL2CLINESIZE;

// _________________________________________________________________________
// Class TdDataReader
//
// Constructor
// 
// dbname, database name;
// mapped, flag of whether the database is mapped;
// nagents, number of agents processing each read chunk of data;
// config, configuration;
//
TdDataReader::TdDataReader(
    const std::vector<std::string>& dbnamelist,
    bool mapped,
    int nagents,
    Configuration* config )
:   tobj_(NULL),
    req_msg_(TREADER_MSG_UNSET),
    rsp_msg_(TREADER_MSG_UNSET),
    //
    chunkdatasize_(0UL),
    chunkdatalen_(0UL),
    chunknpros_(0UL),
    //
    dbnamelist_(dbnamelist),
    dbobj_(nullptr),
    config_(config),
    mapped_(mapped),
    totprodbsize_(0UL),
    totndbentries_(0UL),
    //
    prodbsize_(0UL),
    ndbentries_(0UL),
    addrdesc_end_addrs_(0UL),
    addrdescs_(0UL),
    //data:
    bdbClengths_(nullptr),
    bdbClengths_from_(0UL),
    bdbClengths_to_(0UL),
    bdbCdesc_end_addrs_(nullptr),
    bdbCstruct_(CLOptions::GetIO_NBUFFERS()),
    pbdbCstruct_ndx_(-1),
    bdbCdata_from_(0UL),
    bdbCdata_to_(0UL),
    bdbCdata_poss_from_(0UL),
    bdbCdata_poss_to_(0UL),
    addrdescproced_(0),
    endofdata_(false),
    bdbCprodescs_(NULL),
    bdbCpmbeg_(NULL),
    bdbCpmend_(NULL),
    szpm2dvfields_(NULL),
    ccnt_(NULL),
    lastchunk_(true),
    nagents_(nagents)
{
    MYMSG( "TdDataReader::TdDataReader", 3 );
    memset( addrtable_, 0, pmv2DTotFlds * sizeof(size_t));
    tobj_ = new std::thread( &TdDataReader::Execute, this, (void*)NULL );
}

// Default constructor
//
TdDataReader::TdDataReader()
:   tobj_(NULL),
    dbnamelist_(dummyvec),
    config_(NULL),
    bdbCprodescs_(NULL),
    bdbCpmbeg_(NULL),
    bdbCpmend_(NULL),
    szpm2dvfields_(NULL),
    ccnt_(NULL)
{
    throw MYRUNTIME_ERROR("TdDataReader::TdDataReader: Default initialization is prohibited.");
}

// Destructor
//
TdDataReader::~TdDataReader()
{
    MYMSG( "TdDataReader::~TdDataReader", 3 );
    if( tobj_ ) {
        Notify(tdrmsgTerminate);
        tobj_->join();
        delete tobj_;
        tobj_ = NULL;
    }
}

// -------------------------------------------------------------------------
// Execute: thread's starting point and execution process
//
void TdDataReader::Execute( void* )
{
    MYMSG( "TdDataReader::Execute", 3 );
    myruntime_error mre;
    size_t dbndx = 0;//database index


    //get-data functional with control over multiple databases
    std::function<int(size_t,size_t,size_t)> lfGetDataMeta = 
    [this, &dbndx](size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros)
    {
        if(dbnamelistsrtd_.size() <= dbndx)
            return tdrrespmsgNoData;
        while( !GetData(chunkdatasize, chunkdatalen, chunknpros)) {
            if(dbnamelistsrtd_.size() <= ++dbndx)
                return tdrrespmsgNoData;
            OpenAndReadPreamble(dbnamelistsrtd_[dbndx]);
        }
        if(dbndx+1 < dbnamelistsrtd_.size())
            lastchunk_ = false;
        return tdrrespmsgDataReady;
    };


    try {
        CalculateTotalDbSize();

        if(dbnamelistsrtd_.size() < 1)
            throw MYRUNTIME_ERROR("No valid databases.");

        OpenAndReadPreamble(dbnamelistsrtd_[dbndx]);

        while(1) {
            //wait until the master bradcasts a message
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);

            cv_msg_.wait(lck_msg,
                [this]{return 
                    ((0 <= req_msg_ && req_msg_ <= tdrmsgTerminate) || 
                    req_msg_ == TREADER_MSG_ERROR
                    );}
            );

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "TdDataReader::Execute: Msg %d",req_msg_);
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //thread owns the lock after the wait;
            //read message req_msg_
            int reqmsg = req_msg_;

            //unset the message to avoid live cycle when starting over the loop
            req_msg_ = TREADER_MSG_UNSET;

            //set response msg to error upon occurance of an exception
            rsp_msg_ = TREADER_MSG_ERROR;
            int rspmsg = rsp_msg_;

            size_t chunkdatasize = chunkdatasize_;
            size_t chunkdatalen = chunkdatalen_;
            size_t chunknpros = chunknpros_;

            //immediately read the master-set data and...
            switch(reqmsg) {
                case tdrmsgGetSize:
                        rspmsg = tdrrespmsgSize;
                        break;
                case tdrmsgGetData:
                        ;;
                        rspmsg = lfGetDataMeta(chunkdatasize, chunkdatalen, chunknpros);
                        ;;
                        break;
                case tdrmsgTerminate:
                        rspmsg = tdrrespmsgTerminating;
                        break;
                default:
                        rspmsg = TREADER_MSG_UNSET;
                        break;
            };

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "TdDataReader::Execute: Msg %d Rsp %d",reqmsg, rspmsg );
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //save response code and ...
            rsp_msg_ = rspmsg;

            //send a message back to the master:
            //unlock the mutex to avoid to block the waiting master and 
            //  notify the master using the cv
            lck_msg.unlock();
            cv_msg_.notify_one();

            if( reqmsg < 0 || reqmsg == tdrmsgTerminate)
                //terminate execution
                break;

            //if not end of file, read the next portion of data in advance;
            //NOTE: the lock has been released
            if( !endofdata_ && reqmsg == tdrmsgGetData )
                GetNextData(chunkdatasize, chunkdatalen, chunknpros);
        }

    } catch( myruntime_error const& ex ) {
        mre = ex;
    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Unknown exception caught.");
    }

    if( mre.isset())
        error( mre.pretty_format().c_str());

    if(dbobj_)
        dbobj_->Close();//the thread closes a db

    if( mre.isset()) {
        {//notify the master
            std::lock_guard<std::mutex> lck_msg(mx_dataccess_);
            rsp_msg_ = TREADER_MSG_ERROR;
        }
        cv_msg_.notify_one();
        return;
    }
    {//set the exit flag
        std::lock_guard<std::mutex> lck_msg(mx_dataccess_);
        rsp_msg_ = TREADER_MSG_EXIT;
    }
    cv_msg_.notify_one();
}

// =========================================================================
// CalculateTotalDbSize: sort database names by size and calculate total 
// database size over all databases
// 
void TdDataReader::CalculateTotalDbSize()
{
    std::lock_guard<std::mutex> lck_msg(mx_dataccess_);

    //set the error flag at first
    rsp_msg_ = TREADER_MSG_ERROR;

    MYMSG( "TdDataReader::CalculateTotalDbSize", 3 );
    const mystring preamb = "TdDataReader::CalculateTotalDbSize: ";

    //index-size pairs:
    std::vector<std::pair<size_t,size_t>> dbsizes;
    size_t totnentries = 0UL, totdbsize = 0UL;

    dbsizes.reserve(dbnamelist_.size());

    for(size_t i = 0; i < dbnamelist_.size(); i++) {
        CuDbReader db(dbnamelist_[i].c_str());
        try {
            db.Open();//will be closed on destruction
        } catch( myruntime_error const& ex ) {
            error( ex.pretty_format().c_str());
            mystring strbuf = "Database skipped: ";
            strbuf += dbnamelist_[i].c_str();
            warning(strbuf.c_str());
            continue;
        }
        const size_t szdbverstr = strlen(Db::patstrDBBINVER);
        const size_t dzdbver = strlen(pmodel::PMProfileModel::GetBinaryDataVersion());
        size_t szdbattr = szdbverstr + dzdbver + sizeof(size_t) + sizeof(size_t);
        size_t sztable = pmodel::PMProfileModel::v2_2_NSECTIONS * sizeof(size_t);
        char* p = NULL;

        std::unique_ptr<char,DRDataDeleter>
            locbuf((char*)std::malloc(szdbattr+sztable), DRDataDeleter());

        db.ReadData( 0, p = locbuf.get(), szdbattr+sztable );//READ

        if( strncmp(p,Db::patstrDBBINVER,szdbverstr)) {
            mystring strbuf = "Database skipped due to wrong binary format: ";
            strbuf += dbnamelist_[i].c_str();
            warning(strbuf.c_str());
            continue;
        }
        p += szdbverstr;

        if( strncmp(p,pmodel::PMProfileModel::GetBinaryDataVersion(),dzdbver)) {
            mystring strbuf = "Database skipped due to inconsistent version number: ";
            strbuf += dbnamelist_[i].c_str();
            warning(strbuf.c_str());
            continue;
        }
        p += dzdbver;

        //number of profiles
        size_t nentries = *(size_t*)p;
        p += sizeof(size_t);
        if( nentries < 1 ) {
            mystring strbuf = "Database skipped due to invalid number of profiles: ";
            strbuf += dbnamelist_[i].c_str();
            warning(strbuf.c_str());
            continue;
        }

        //database size
        size_t dbsize = *(size_t*)p;
        p += sizeof(size_t);
        if( dbsize < 1 ) {
            mystring strbuf = "Database skipped due to invalid size: ";
            strbuf += dbnamelist_[i].c_str();
            warning(strbuf.c_str());
            continue;
        }

        totnentries += nentries;
        totdbsize += dbsize;
        dbsizes.push_back(std::make_pair(i,dbsize));
    }

    SetDbSize( totnentries, totdbsize );

    //sort by db size
    std::sort(dbsizes.begin(), dbsizes.end(),
        [](const std::pair<size_t,size_t>& p1, const std::pair<size_t,size_t>& p2) {
            return p1.second > p2.second;
        });

    //save valid sorted database names
    dbnamelistsrtd_.reserve(dbsizes.size());
    for(const std::pair<size_t,size_t>& pr: dbsizes) {
        dbnamelistsrtd_.push_back(dbnamelist_[pr.first]);
    }

    rsp_msg_ = TREADER_MSG_UNSET;
}



// =========================================================================
//
void GetNextValue( size_t* value, char*& pfrom, const size_t filesize )
{
    static const mystring errortable = "GetNextValue: Invalid address table in file.";
    *value = *(size_t*)pfrom;
    if( filesize <= *value )
        throw MYRUNTIME_ERROR(errortable);
    pfrom += sizeof(size_t);
}

// -------------------------------------------------------------------------
// OpenAndReadPreamble: open a database and read preamble, including 
// database attributes, address table and profile lengths
// 
void TdDataReader::OpenAndReadPreamble(const std::string& dbname)
{
    MYMSG( "TdDataReader::OpenAndReadPreamble", 3 );
    const mystring preamb = "TdDataReader::OpenAndReadPreamble: ";

    //close will be called on destruction:
    dbobj_.reset( new CuDbReader(dbname.c_str(), mapped_));
    if(!dbobj_)
        throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

    dbobj_->Open();

    const size_t filesize = dbobj_->GetDbSizeInBytes();
    const size_t szdbverstr = strlen(Db::patstrDBBINVER);
    const size_t dzdbver = strlen(pmodel::PMProfileModel::GetBinaryDataVersion());
    size_t szdbattr = szdbverstr + dzdbver + sizeof(size_t) + sizeof(size_t);
    size_t sztable = pmodel::PMProfileModel::v2_2_NSECTIONS * sizeof(size_t);

    const int noress = NUMAA;//# residues
    const int vsize = NUMAA-1;//context vector length
    int r, t;
    char* p = NULL;

    std::unique_ptr<char,DRDataDeleter>
        locbuf((char*)std::malloc(szdbattr+sztable), DRDataDeleter());

    bdbClengths_from_ = (0UL);
    bdbClengths_to_ = (0UL);

    endofdata_ = false;
    pbdbCstruct_ndx_ = -1;

    ResetProfileCounters();

    dbobj_->ReadData( 0, p = locbuf.get(), szdbattr+sztable );//READ

    if( strncmp(p,Db::patstrDBBINVER,szdbverstr))
        throw MYRUNTIME_ERROR( preamb + "Wrong database binary format.");
    p += szdbverstr;

    if( strncmp(p,pmodel::PMProfileModel::GetBinaryDataVersion(),dzdbver))
        throw MYRUNTIME_ERROR( preamb + "Inconsistent database version number.");
    p += dzdbver;

    //number of profiles
    size_t nentries = *(size_t*)p;
    p += sizeof(size_t);
    if( nentries < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles read.");

    //database size
    size_t dbsize = *(size_t*)p;
    p += sizeof(size_t);
    if( dbsize < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid database size read.");

    prodbsize_ = dbsize;
    ndbentries_ = nentries;

    //fill in address table
    GetNextValue( addrtable_ + pps2DLen, p, filesize );//lengths
    GetNextValue( addrtable_ + pps2DENO, p, filesize );//ENOs
    for( r = 0; r < noress; r++ )
        GetNextValue( addrtable_ + pps2DBkgPrbs+r, p, filesize );//bkg probabilities for each r
    for( t = 0; t < P_NSTATES-2; t++ )
        GetNextValue( addrtable_ + ptr2DTrnPrbs+t, p, filesize );//trn probabilities for each t
    for( r = 0; r < noress; r++ )
        GetNextValue( addrtable_ + pmv2DTrgFrqs+r, p, filesize );//target probabilities for each r
    for( t = 0; t < vsize; t++ )
        GetNextValue( addrtable_ + pmv2DCVentrs+t, p, filesize );//context vector entries (for each t)
    GetNextValue( addrtable_ + pmv2DCVprior, p, filesize );//context vector probabilities
    GetNextValue( addrtable_ + pmv2DCVnorm2, p, filesize );//squared norms of context vector
    for( t = 0; t < SS_NSTATES; t++ )
        GetNextValue( addrtable_ + pmv2DSSsprbs+t, p, filesize );//secondary structure state probabilities for each t
    GetNextValue( addrtable_ + pmv2DHDP1prb, p, filesize );//HDP1 cluster membership probabilities
    GetNextValue( addrtable_ + pmv2DHDP1ind, p, filesize );//HDP1 cluster indices
    GetNextValue( addrtable_ + pmv2Daa, p, filesize );//amino acid sequences
    GetNextValue( addrtable_ + pmv2DSSstate, p, filesize );//secondary structure state sequences

    GetNextValue(&addrdesc_end_addrs_, p, filesize );//end addresses of profile descriptions
    GetNextValue(&addrdescs_, p, filesize );//profile descriptions

    //address of the first description
    addrdescproced_ = addrdescs_;

    //read lengths and end addresses
    if( !ReadLengthsAndDescEndAddrs())
        throw MYRUNTIME_ERROR( preamb + "Failed to read profile lengths.");
}

// -------------------------------------------------------------------------
// ReadLengthsAndDescEndAddrs: read a number of or all lengths and end 
// addresses of profile descriptions at once;
// return false if there are no data to read;
// 
bool TdDataReader::ReadLengthsAndDescEndAddrs()
{
    const mystring preamb = "TdDataReader::ReadLengthsAndDescEndAddrs: ";

    if( ndbentries_ <= bdbClengths_to_ )
        //all data has been processed
        return false;

    if(!dbobj_)
        throw MYRUNTIME_ERROR(preamb + "NULL Db object.");

    //new number of profiles:
    bdbClengths_from_ = bdbClengths_to_;
    bdbClengths_to_ = bdbClengths_from_ + 
        SLC_MIN(ndbentries_ - bdbClengths_from_, (size_t)tdrMAX_N_ENTRIES_READ);
    size_t diffpron = bdbClengths_to_ - bdbClengths_from_;//difference in the number of profiles

    //read lengths
    size_t startaddr = addrtable_[pps2DLen] + bdbClengths_from_ * sizeof(int);
    size_t allocsize = diffpron * sizeof(int);
    bdbClengths_.reset((int*)std::malloc(allocsize));

    dbobj_->ReadData(startaddr, bdbClengths_.get(), allocsize);//READ

    //read end addresses of profile descriptions
    startaddr = addrdesc_end_addrs_ + bdbClengths_from_ * sizeof(size_t);
    allocsize = diffpron * sizeof(size_t);
    bdbCdesc_end_addrs_.reset((size_t*)std::malloc(allocsize));

    dbobj_->ReadData( startaddr, bdbCdesc_end_addrs_.get(), allocsize);//READ

    return true;
}

// -------------------------------------------------------------------------
// GetData: get portion of data given maximum data chunk size and profile 
// values;
// chunkdatasize, data chunk size;
// chunkdatalen, total profile length limit;
// chunknpros, limit number of profiles;
// 
bool TdDataReader::GetData(
    size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros )
{
    MYMSG( "TdDataReader::GetData", 3 );
    const mystring preamb = "TdDataReader::GetData: ";
    bool dataread = !endofdata_;

    if( pbdbCstruct_ndx_ < 0 ) {
        //starting reading the file
        pbdbCstruct_ndx_ = 0;
        dataread = ReadDataChunk(&bdbCstruct_[pbdbCstruct_ndx_], chunkdatasize, chunkdatalen, chunknpros);
    }
    bdbCprodescs_ = bdbCstruct_[pbdbCstruct_ndx_].bdbCprodescs_.get();
    bdbCpmbeg_ = bdbCstruct_[pbdbCstruct_ndx_].bdbCpmbeg_;
    bdbCpmend_ = bdbCstruct_[pbdbCstruct_ndx_].bdbCpmend_;
    szpm2dvfields_ = bdbCstruct_[pbdbCstruct_ndx_].szpm2dvfields_;
    ccnt_ = &bdbCstruct_[pbdbCstruct_ndx_].cnt_;
    lastchunk_ = (ndbentries_ <= bdbCdata_to_);
    return dataread;
}

// -------------------------------------------------------------------------
// GetNextData: get next portion of data given maximum data chunk size and 
// profile values;
// chunkdatasize, data chunk size;
// chunkdatalen, total profile length limit;
// chunknpros, limit number of profiles;
// 
void TdDataReader::GetNextData(
    size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros)
{
    MYMSG( "TdDataReader::GetNextData", 3 );
    const mystring preamb = "TdDataReader::GetNextData: ";

    if( pbdbCstruct_ndx_ < 0 )
        throw MYRUNTIME_ERROR(preamb + "Memory access error.");

    pbdbCstruct_ndx_++;
    if( pbdbCstruct_ndx_ >= (int)bdbCstruct_.size())
        pbdbCstruct_ndx_ = 0;

    endofdata_ = ! ReadDataChunk(&bdbCstruct_[pbdbCstruct_ndx_], chunkdatasize, chunkdatalen, chunknpros);

    MYMSGBEGl(3)
        char strbuf[BUF_MAX];
        sprintf(strbuf,"TdDataReader::GetNextData: eof %d",endofdata_);
        MYMSG(strbuf,3);
    MYMSGENDl
}

// -------------------------------------------------------------------------
// ReadDataChunk: read a data chunk from file;
// return false if there are no data to read;
// 
bool TdDataReader::ReadDataChunk( 
    SbdbCData* pbdbC, 
    size_t chunkdatasize, size_t chunkdatalen, size_t chunknpros )
{
    MYMSG( "TdDataReader::ReadDataChunk", 3 );
    const mystring preamb = "TdDataReader::ReadDataChunk: ";

#ifdef __DEBUG__
    if( pbdbC == NULL )
        throw MYRUNTIME_ERROR(preamb + "Memory access error.");
#endif

    MYMSG("TdDataReader::ReadDataChunk: Waiting until release of data...", 5 );
    //NOTE: wait for all data-using agents to finish their work
    pbdbC->cnt_.wait0();
    MYMSG("TdDataReader::ReadDataChunk: Wait finished", 5 );

    if( chunkdatasize < 1 || chunkdatalen < 1 || chunknpros < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid requested data sizes.");

    pbdbC->AllocateSpace(chunkdatasize, chunknpros);

    int* plens = bdbClengths_.get();
    char strbuf[BUF_MAX];
    int prolen;
    size_t nskipped = 0;//#skipped profiles
    size_t iinit, i, szpm1, szpmtotal = 0, totlen = 0, npros = 0;

    //print warning if there profiles have been skipped
    std::function<void(size_t)> lfWarnifSkipped = [&strbuf](size_t nskpd) {
        if(nskpd) {
            sprintf(strbuf, 
              "IMPOSSIBLE ALIGNMENT with %zu Db profile(s) due to "
              "memory limit (GPU or --dev-mem): SKIPPED.", nskpd);
            warning(strbuf);
        }
    };

    //{{NOTE: skip all profiles (sorted by length) that do not fit into memory
    for(;;) {
        if( bdbClengths_to_ <= bdbCdata_to_ ) {
            //current buffer has been processed, read another block
            if( !ReadLengthsAndDescEndAddrs()) {
                lfWarnifSkipped(nskipped);
                return false;
            }
            plens = bdbClengths_.get();
        }
#ifdef __DEBUG__
        if( plens == NULL )
            throw MYRUNTIME_ERROR( preamb + "Memory access error.");
#endif
        if( bdbCdata_to_ < bdbClengths_from_ )
            throw MYRUNTIME_ERROR( preamb + "Invalid data index.");

        prolen = plens[bdbCdata_to_ - bdbClengths_from_];
        szpm1 = PMBatchProData::GetPMDataSize1((size_t)prolen);
        if(szpm1 < chunkdatasize && (size_t)prolen < chunkdatalen && 1 < chunknpros) {
            lfWarnifSkipped(nskipped);
            break;
        }
        nskipped++;
        bdbCdata_to_++;
        bdbCdata_poss_to_ += prolen;
    }
    //}}

    bdbCdata_from_ = bdbCdata_to_;
    bdbCdata_poss_from_ = bdbCdata_poss_to_;

    if( bdbCdata_from_ < bdbClengths_from_ )
        throw MYRUNTIME_ERROR( preamb + "Invalid data index.");

    for( iinit = i = bdbCdata_from_ - bdbClengths_from_;
        bdbCdata_to_ < bdbClengths_to_; 
        bdbCdata_to_++, i++ )
    {
        prolen = plens[i];
        if( prolen < 1 ) {
            sprintf(strbuf, "Invalid profile length: %d; profile no. %zu",
                prolen, bdbCdata_to_);
            throw MYRUNTIME_ERROR( preamb + strbuf);
        }
        szpm1 = PMBatchProData::GetPMDataSize1((size_t)prolen);
        if( chunkdatasize < szpmtotal + szpm1 || 
            chunkdatalen < totlen + prolen ||
            chunknpros < npros+1 )
            break;
        totlen += prolen;
        szpmtotal += szpm1;
        npros++;
        bdbCdata_poss_to_ += prolen;
        //
        if((size_t)INT_MAX <= totlen )
            throw MYRUNTIME_ERROR( preamb + "Overflow detected. Data chunk size must be reduced.");
    }

    if( szpmtotal == 0 || totlen == 0 || npros == 0 )
        throw MYRUNTIME_ERROR( preamb + "Data not read plausibly due to memory shortage.");

    ReadDataChunkHelper( pbdbC, iinit, npros, totlen);

    //NOTE: set the number of agents expected to access the read data
    pbdbC->cnt_.set(nagents_);

    return true;
}

// -------------------------------------------------------------------------
// ReadDataChunkHelper: read a data chunk from file, format and pack it for
// processing;
// lengthsndx, starting index of the lengths structure;
// npros, number of profiles to read;
// totlen, total number of positions the profiles contain;
// 
void TdDataReader::ReadDataChunkHelper(
    SbdbCData* pbdbC, size_t lengthsndx, size_t npros, size_t totlen)
{
    MYMSGBEGl(3)
        char strbuf[BUF_MAX];
        sprintf(strbuf,"TdDataReader::ReadDataChunkHelper: %zu pros of length %zu",npros,totlen);
        MYMSG(strbuf,3);
    MYMSGENDl
    const mystring preamb = "TdDataReader::ReadDataChunkHelper: ";

    if(!dbobj_)
        throw MYRUNTIME_ERROR(preamb + "NULL Db object.");

#ifdef __DEBUG__
    if( pbdbC == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
#endif

    if( bdbCdata_to_ <= bdbCdata_from_ ||
        bdbCdata_poss_to_ <= bdbCdata_poss_from_ )
        throw MYRUNTIME_ERROR( preamb + "Invalid data pointers.");

    int* plens = bdbClengths_.get();
    char* pbdbCdata = pbdbC->bdbCdata_.get();
    char* pbdbCdesc = pbdbC->bdbCdescs_.get();
    char** pbdbCprodesc = pbdbC->bdbCprodescs_.get();
    char* pbdbCdist = NULL, *pbdbCaddrPro = NULL, *pbdbCaddrSS = NULL;
    char* pbdbCSSsprbs[pmv2DNoSSSps] = {NULL,NULL,NULL};
    char* pendaddrs = NULL;
    //
    char** bdbCpmbeg = pbdbC->bdbCpmbeg_;
    char** bdbCpmend = pbdbC->bdbCpmend_;
    size_t* szpm2dvfields = pbdbC->szpm2dvfields_;
    //
    int prolen;
    size_t szdat, szval, n, nlen;

#ifdef __DEBUG__
    if( plens == NULL || pbdbCdata == NULL || pbdbCdesc == NULL || pbdbCprodesc == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
#endif

    szpm2dvfields[0] = 0;

    //pack data
    //READ//bkg probabilities for each n//[pps2DBkgPrbs]
    for( n = 0; n < pmv2DNoElems; n++ ) {
        bdbCpmbeg[pps2DBkgPrbs+n] = pbdbCdata;//
        szdat = SZFPTYPE * npros;//size to read
        dbobj_->ReadData( 
            addrtable_[pps2DBkgPrbs+n] + SZFPTYPE * bdbCdata_from_,
            pbdbCdata, szdat );
        szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
        bdbCpmend[pps2DBkgPrbs+n] = pbdbCdata + szdat;//
        pbdbCdata += szval;
        szpm2dvfields[pps2DBkgPrbs+n+1] = szpm2dvfields[pps2DBkgPrbs+n] + szval;//
    }
    //READ//ENOs//[pps2DENO]
    bdbCpmbeg[pps2DENO] = pbdbCdata;//
    szdat = SZFPTYPE * npros;
    dbobj_->ReadData( 
        addrtable_[pps2DENO] + SZFPTYPE * bdbCdata_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pps2DENO] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pps2DENO+1] = szpm2dvfields[pps2DENO] + szval;//
    //mem copy//lengths//[pps2DLen]
    bdbCpmbeg[pps2DLen] = pbdbCdata;//
    szdat = SZINTYPE * npros;
    memcpy( pbdbCdata, plens + lengthsndx, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pps2DLen] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pps2DLen+1] = szpm2dvfields[pps2DLen] + szval;//
    //NOTE: omit distances to be filled in later in the final loop//[pps2DDist]
    bdbCpmbeg[pps2DDist] = pbdbCdata;//
    pbdbCdist = pbdbCdata;//save the address of distances
    szdat = SZLNTYPE * npros;
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pps2DDist] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pps2DDist+1] = szpm2dvfields[pps2DDist] + szval;//
    //READ//trn probabilities for each n//[ptr2DTrnPrbs]
    for( n = 0; n < ptr2DNoElems; n++ ) {
        bdbCpmbeg[ptr2DTrnPrbs+n] = pbdbCdata;//
        szdat = SZFPTYPE * (npros + totlen);//size to read
        dbobj_->ReadData(
            addrtable_[ptr2DTrnPrbs+n] + SZFPTYPE * (bdbCdata_from_ + bdbCdata_poss_from_),
            pbdbCdata, szdat );
        szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
        bdbCpmend[ptr2DTrnPrbs+n] = pbdbCdata + szdat;//
        pbdbCdata += szval;
        szpm2dvfields[ptr2DTrnPrbs+n+1] = szpm2dvfields[ptr2DTrnPrbs+n] + szval;//
    }
    //READ//target probabilities for each n//[pmv2DTrgFrqs]
    for( n = 0; n < pmv2DNoElems; n++ ) {
        bdbCpmbeg[pmv2DTrgFrqs+n] = pbdbCdata;//
        szdat = SZFPTYPE * totlen;
        dbobj_->ReadData(
            addrtable_[pmv2DTrgFrqs+n] + SZFPTYPE * bdbCdata_poss_from_,
            pbdbCdata, szdat );
        szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
        bdbCpmend[pmv2DTrgFrqs+n] = pbdbCdata + szdat;//
        pbdbCdata += szval;
        szpm2dvfields[pmv2DTrgFrqs+n+1] = szpm2dvfields[pmv2DTrgFrqs+n] + szval;//
    }
    //READ//context vector entries for each n//[pmv2DCVentrs]
    for( n = 0; n < pmv2DNoCVEls; n++ ) {
        bdbCpmbeg[pmv2DCVentrs+n] = pbdbCdata;//
        szdat = SZFPTYPE * totlen;
        dbobj_->ReadData(
            addrtable_[pmv2DCVentrs+n] + SZFPTYPE * bdbCdata_poss_from_,
            pbdbCdata, szdat );
        szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
        bdbCpmend[pmv2DCVentrs+n] = pbdbCdata + szdat;//
        pbdbCdata += szval;
        szpm2dvfields[pmv2DCVentrs+n+1] = szpm2dvfields[pmv2DCVentrs+n] + szval;//
    }
    //READ//context vector probabilities//[pmv2DCVprior]
    bdbCpmbeg[pmv2DCVprior] = pbdbCdata;//
    szdat = SZFPTYPE * totlen;
    dbobj_->ReadData( 
        addrtable_[pmv2DCVprior] + SZFPTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DCVprior] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DCVprior+1] = szpm2dvfields[pmv2DCVprior] + szval;//
    //READ//squared norms of context vectors//[pmv2DCVnorm2]
    bdbCpmbeg[pmv2DCVnorm2] = pbdbCdata;//
    szdat = SZFPTYPE * totlen;
    dbobj_->ReadData( 
        addrtable_[pmv2DCVnorm2] + SZFPTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DCVnorm2] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DCVnorm2+1] = szpm2dvfields[pmv2DCVnorm2] + szval;//
    //READ//secondary structure state probabilities for each n//[pmv2DSSsprbs]
    for( n = 0; n < pmv2DNoSSSps; n++ ) {
        bdbCpmbeg[pmv2DSSsprbs+n] = pbdbCdata;//
        pbdbCSSsprbs[n] = pbdbCdata;//save the beginning addresses of SS state probabilities
        szdat = SZFPTYPE * totlen;
        dbobj_->ReadData(
            addrtable_[pmv2DSSsprbs+n] + SZFPTYPE * bdbCdata_poss_from_,
            pbdbCdata, szdat );
        szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
        bdbCpmend[pmv2DSSsprbs+n] = pbdbCdata + szdat;//
        pbdbCdata += szval;
        szpm2dvfields[pmv2DSSsprbs+n+1] = szpm2dvfields[pmv2DSSsprbs+n] + szval;//
    }
    //READ//HDP1 cluster membership probabilities//[pmv2DHDP1prb]
    bdbCpmbeg[pmv2DHDP1prb] = pbdbCdata;//
    szdat = SZFPTYPE * totlen;
    dbobj_->ReadData( 
        addrtable_[pmv2DHDP1prb] + SZFPTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DHDP1prb] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DHDP1prb+1] = szpm2dvfields[pmv2DHDP1prb] + szval;//
    //READ//HDP1 cluster indices//[pmv2DHDP1ind]
    bdbCpmbeg[pmv2DHDP1ind] = pbdbCdata;//
    szdat = SZINTYPE * totlen;
    dbobj_->ReadData( 
        addrtable_[pmv2DHDP1ind] + SZINTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DHDP1ind] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DHDP1ind+1] = szpm2dvfields[pmv2DHDP1ind] + szval;//
    //NOTE: profile local addresses to be filled in later in the final loop//[pmv2DAddrPro]
    bdbCpmbeg[pmv2DAddrPro] = pbdbCdata;//
    pbdbCaddrPro = pbdbCdata;//save the address for later reference
    szdat = SZINTYPE * totlen;
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DAddrPro] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DAddrPro+1] = szpm2dvfields[pmv2DAddrPro] + szval;//
    //NOTE: context vector local addresses WON'T BE FILLED BUT
    // THE SPACE IS CURRENTLY RESERVED!//[pmv2DAddrCV]
    bdbCpmbeg[pmv2DAddrCV] = pbdbCdata;//
    szdat = SZLNTYPE * totlen;
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DAddrCV] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DAddrCV+1] = szpm2dvfields[pmv2DAddrCV] + szval;//
    //NOTE: SS info local addresses to be filled in later in the final loop//[pmv2DAddrSS]
    bdbCpmbeg[pmv2DAddrSS] = pbdbCdata;//
    pbdbCaddrSS = pbdbCdata;//save the address for later reference
    szdat = SZLNTYPE * totlen;
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DAddrSS] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DAddrSS+1] = szpm2dvfields[pmv2DAddrSS] + szval;//
    //READ//amino acid sequences//[pmv2Daa]
    bdbCpmbeg[pmv2Daa] = pbdbCdata;//
    szdat = SZCHTYPE * totlen;
    dbobj_->ReadData( 
        addrtable_[pmv2Daa] + SZCHTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2Daa] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2Daa+1] = szpm2dvfields[pmv2Daa] + szval;//
    //READ//secondary structure state sequences//[pmv2DSSstate]
    bdbCpmbeg[pmv2DSSstate] = pbdbCdata;//
    szdat = SZCHTYPE * totlen;
    dbobj_->ReadData(
        addrtable_[pmv2DSSstate] + SZCHTYPE * bdbCdata_poss_from_,
        pbdbCdata, szdat );
    szval = ALIGN_UP(szdat, TdDataReader::s_szdatalign_);
    bdbCpmend[pmv2DSSstate] = pbdbCdata + szdat;//
    pbdbCdata += szval;
    szpm2dvfields[pmv2DSSstate+1] = szpm2dvfields[pmv2DSSstate] + szval;//

    //READ//end addresses of profile descriptions
    szdat = sizeof(size_t) * npros;
    dbobj_->ReadData(
        addrdesc_end_addrs_ + sizeof(size_t) * bdbCdata_from_,
        pbdbCdata, szdat );
    pendaddrs = pbdbCdata;//save the address of description end addresses
    pbdbCdata += szdat;
    //READ//profile descriptions
    //address of the last description to be processed:
    size_t addrnextlastdesc = *(((size_t*)pbdbCdata)-1);
    if( addrnextlastdesc <= addrdescproced_ )
        throw MYRUNTIME_ERROR( preamb + "Invalid address of descriptions read.");
    szdat = addrnextlastdesc - addrdescproced_;
    if( pbdbC->szbdbCdescs_ < szdat ) {
        pbdbC->AllocateSpaceForDescriptions(szdat);
        pbdbCdesc = pbdbC->bdbCdescs_.get();
    }
    dbobj_->ReadData(addrdescproced_, pbdbCdesc, szdat);

    size_t addrprevdesc = addrdescproced_;

    addrdescproced_ = addrnextlastdesc;

    //Now fill in deferred structures
    for( nlen = lengthsndx, n = 0, szdat = 0; n < npros; n++, nlen++ ) {
        prolen = plens[nlen];
        //
        //write distances//[pps2DDist]
        *(LNTYPE*)pbdbCdist = (LNTYPE)szdat;
        pbdbCdist += SZLNTYPE;
        //write profile local addresses//[pmv2DAddrPro]
        //NOTE: expensive loop over all positions
        std::fill((INTYPE*)pbdbCaddrPro, (INTYPE*)(pbdbCaddrPro + SZINTYPE * prolen), (INTYPE)n);
        pbdbCaddrPro += SZINTYPE * prolen;
        //write SS info local addresses;
        //NOTE: only the first positions are filled with the rest being UNINITIALIZED!
        //NOTE: THE SPACE IS CURRENTLY RESERVED!//[pmv2DAddrSS]
        if( pbdbCSSsprbs[0]||pbdbCSSsprbs[1]||pbdbCSSsprbs[2])
            *(LNTYPE*)pbdbCaddrSS = (LNTYPE)szdat;
        else
            *(LNTYPE*)pbdbCaddrSS = (LNTYPE)-1;//SS information is not used
        pbdbCaddrSS += SZLNTYPE * prolen;
        pbdbCSSsprbs[0] += SZFPTYPE * prolen;
        pbdbCSSsprbs[1] += SZFPTYPE * prolen;
        pbdbCSSsprbs[2] += SZFPTYPE * prolen;
        //initialize description pointers for each of npros profiles
        pbdbCprodesc[n] = pbdbCdesc;
        //pbdbCdesc += strlen(pbdbCprodesc[n])+1;//plus 1 to pass over 0
        pbdbCdesc += *(size_t*)pendaddrs - addrprevdesc;
        addrprevdesc = *(size_t*)pendaddrs;
        pendaddrs += sizeof(size_t);
        //
        szdat += (size_t)(LNTYPE)prolen;
        if((size_t)INT_MAX <= szdat )
            throw MYRUNTIME_ERROR( preamb + "Overflow detected. Data chunk size must be reduced.");
    }
}



// =========================================================================
// AllocateSpace: allocate space for profile data and descriptions given the 
// limits of data chunk size and the number of profiles
// 
void SbdbCData::AllocateSpace(size_t chunkdatasize, size_t chunknpros)
{
    //allocate space for profile data plus space for the end addresses of profile descriptions
    AllocateSpaceForData(
        chunkdatasize + TdDataReader::s_szdatalign_ * pmv2DTotFlds/*for data alignment*/ +
        chunknpros * sizeof(size_t));
    //
    AllocateSpaceForDescriptions(chunknpros * TdDataReader::tdrDEFDESCLEN);
    AllocateSpaceForDescPtrs(chunknpros);
}

// -------------------------------------------------------------------------
// AllocateSpaceForData: allocate space for profile data
// 
void SbdbCData::AllocateSpaceForData( size_t requestedsize )
{
    if( requestedsize <= szbdbCdata_ )
        return;
    char strbuf[BUF_MAX];
    sprintf(strbuf, "SbdbCData::AllocateSpaceForData: "
        "New allocation for profile data, %zuB", requestedsize);
    MYMSG( strbuf, 3 );
    szbdbCdata_ = requestedsize;

    if(CLOptions::GetIO_UNPINNED()) {
        bdbCdata_.reset((char*)std::malloc(szbdbCdata_));
    }
    else {
        char* h_mpinned;
        MYCUDACHECK( cudaMallocHost((void**)&h_mpinned, szbdbCdata_) );
        MYCUDACHECKLAST;
        bdbCdata_.reset(h_mpinned);
    }
    if( bdbCdata_.get() == NULL )
        throw MYRUNTIME_ERROR("SbdbCData::AllocateSpaceForData: Memory allocation failed.");
}

// -------------------------------------------------------------------------
// AllocateSpaceForDescriptions: allocate space for profile descriptions
// 
bool SbdbCData::AllocateSpaceForDescriptions( size_t requestedsize )
{
    if( requestedsize <= szbdbCdescs_ )
        return false;
    char strbuf[BUF_MAX];
    sprintf(strbuf, "SbdbCData::AllocateSpaceForDescriptions: "
        "New allocation for descriptions, %zuB", requestedsize);
    MYMSG( strbuf, 3 );
    szbdbCdescs_ = requestedsize;
    bdbCdescs_.reset((char*)std::malloc(szbdbCdescs_));
    if( bdbCdescs_.get() == NULL )
        throw MYRUNTIME_ERROR("SbdbCData::AllocateSpaceForDescriptions: Not enough memory.");
    return true;
}

// -------------------------------------------------------------------------
// AllocateSpaceForDescPtrs: allocate space for pointers to the 
// descriptions of each profile
// 
void SbdbCData::AllocateSpaceForDescPtrs(size_t requestedno)
{
    if( requestedno <= nbdbCprodescs_ )
        return;
    char strbuf[BUF_MAX];
    sprintf(strbuf, "SbdbCData::AllocateSpaceForDescPtrs: "
        "New number of description pointers, %zu", requestedno);
    MYMSG( strbuf, 3 );
    nbdbCprodescs_ = requestedno;
    bdbCprodescs_.reset(new char*[nbdbCprodescs_+1]);//NOTE +1
    if( bdbCprodescs_.get() == NULL )
        throw MYRUNTIME_ERROR("SbdbCData::AllocateSpaceForDescPtrs: Not enough memory.");
}
