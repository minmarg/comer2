/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <math.h>
#include <cmath>

#include <memory>
#include <utility>
#include <functional>
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "tsafety/TSCounterVar.h"

#include "liblib/fmtdescription.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/btckcoords.h"
#include "libmycu/cualn/Devices.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "libmycu/cuss/CuBatchSS_com.h"

#include "CuBatchProcessingFinalizer.h"

static cudaStream_t streamdummy;

// _________________________________________________________________________
// class statics
//

// _________________________________________________________________________
// Class CuBatchProcessingFinalizer
//
// Constructor
//
CuBatchProcessingFinalizer::CuBatchProcessingFinalizer(
    cudaStream_t& strcopyres,
    DeviceProperties dprop, 
    AlnWriter* writer,
    Configuration* config,
    const mystring* queryfnames,
    const mystring* querydescs,
    const char** bdb1descs )
:   tobj_(NULL),
    req_msg_(CUBPTHREAD_MSG_UNSET),
    rsp_msg_(CUBPTHREAD_MSG_UNSET),
    //
    strcopyres_(strcopyres),
    dprop_(dprop),
    alnwriter_(writer),
    config_(config),
    //cached data:
    cached_queryfnames_(queryfnames),
    cached_querydescs_(querydescs),
    cached_bdb1descs_(bdb1descs),
    //data arguments:
    cubp_set_qrysernr_(-1),
    cubp_set_nqyposs_(0),
    cubp_set_qyeno_(0.0f),
    cubp_set_deltalen_(0U),
    cubp_set_sspace_(0.0f),
    cubp_set_logevthld_(0.0f),
    cubp_set_bdbCdesc_(NULL),
    cubp_set_cnt_(NULL),
    cubp_set_ndb1pros_(0UL),
    cubp_set_ndbCpros_(0UL),
    cubp_set_querprosOmtd_(0UL),
    cubp_set_ndb1prosOmtd_(0UL),
    cubp_set_ndbCprosOmtd_(0UL),
    cubp_set_querposoffset_(0UL),
    cubp_set_bdb1posoffset_(0UL),
    cubp_set_bdbCposoffset_(0UL),
    cubp_set_nposits_(0U),
    cubp_set_npros_(0U),
    cubp_set_h_results_(NULL),
    cubp_set_sz_alndata_(0UL),
    cubp_set_sz_alns_(0UL),
    cubp_set_sz_dbalnlen2_(0),
    annotations_(nullptr),
    alignments_(nullptr)
{
    MYMSG( "CuBatchProcessingFinalizer::CuBatchProcessingFinalizer", 3 );
    memset( cubp_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    tobj_ = new std::thread( &CuBatchProcessingFinalizer::Execute, this, (void*)NULL );
}

// Default constructor
//
CuBatchProcessingFinalizer::CuBatchProcessingFinalizer()
:   tobj_(NULL),
    strcopyres_(streamdummy),
    config_(NULL)
{
    throw MYRUNTIME_ERROR(
    "CuBatchProcessingFinalizer::CuBatchProcessingFinalizer: "
    "Default initialization is prohibited.");
}

// Destructor
//
CuBatchProcessingFinalizer::~CuBatchProcessingFinalizer()
{
    MYMSG( "CuBatchProcessingFinalizer::~CuBatchProcessingFinalizer", 3 );
    if( tobj_ ) {
        tobj_->join();
        delete tobj_;
        tobj_ = NULL;
    }
}

// -------------------------------------------------------------------------
// Execute: thread's starting point and execution process
//
void CuBatchProcessingFinalizer::Execute( void* )
{
    MYMSG( "CuBatchProcessingFinalizer::Execute", 3 );
    myruntime_error mre;

    try {
        const int outfmt = CLOptions::GetB_FMT();

        MYCUDACHECK( cudaSetDevice( dprop_.devid_ ));
        MYCUDACHECKLAST;

        while(1) {
            //wait until a message arrives
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);

            cv_msg_.wait(lck_msg,
                [this]{return 
                    ((0 <= req_msg_ && req_msg_ <= cubpthreadmsgTerminate) || 
                     req_msg_ == CUBPTHREAD_MSG_ERROR
                    );}
            );

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "CuBatchProcessingFinalizer::Execute: Msg %d",req_msg_);
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //thread owns the lock after the wait;
            //read message req_msg_
            int reqmsg = req_msg_;

            //unset the message to avoid live cycle when starting over the loop
            req_msg_ = CUBPTHREAD_MSG_UNSET;

            //set response msg to error upon occurance of an exception
            rsp_msg_ = CUBPTHREAD_MSG_ERROR;
            int rspmsg = rsp_msg_;

            switch(reqmsg) {
                case cubpthreadmsgFinalize:
                        //data addresses have been written already;
                        //make sure data transfer has finished
                        MYCUDACHECK( cudaStreamSynchronize(strcopyres_));
                        MYCUDACHECKLAST;
                        ;;
                        if(outfmt==CLOptions::ofJSON)
                            CompressResultsJSON();
                        else
                            CompressResultsPlain();
                        //results have been processed and the buffers are no longer needed;
                        //decrease the counter of data chunk processed by an agent
                        if( cubp_set_cnt_)
                            cubp_set_cnt_->dec();
                        SortCompressedResults();
                        PassResultsToWriter();
                        ;;
                        //parent does not wait for a response nor requires data to read;
                        //unset response code
                        rspmsg = CUBPTHREAD_MSG_UNSET;//cubptrespmsgFinalizing;
                        break;
                case cubpthreadmsgTerminate:
                        rspmsg = cubptrespmsgTerminating;
                        break;
                default:
                        //rspmsg = CUBPTHREAD_MSG_UNSET;
                        break;
            };

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "CuBatchProcessingFinalizer::Execute: Msg %d Rsp %d",reqmsg, rspmsg );
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //save response code
            rsp_msg_ = rspmsg;

            //unlock the mutex and notify the parent using the cv
            lck_msg.unlock();
            cv_msg_.notify_one();

            if( reqmsg < 0 || reqmsg == cubpthreadmsgTerminate)
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
        SetResponseError();
        cv_msg_.notify_one();
        return;
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// PassResultsToWriter: transfer the addresses of sorted results to the 
// alignment writer;
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::PassResultsToWriter()
{
    MYMSG( "CuBatchProcessingFinalizer::PassResultsToWriter", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::PassResultsToWriter: ";

    if( alnwriter_ ) {
        alnwriter_->PushPartOfResults(
            cubp_set_qrysernr_,
            cubp_set_nqyposs_,
            cached_queryfnames_? cached_queryfnames_[cubp_set_querprosOmtd_].c_str(): NULL,
            cached_querydescs_? cached_querydescs_[cubp_set_querprosOmtd_].c_str(): NULL,
            cubp_set_deltalen_,
            cubp_set_sspace_,
            cubp_set_logevthld_,
            std::move(annotations_),
            std::move(alignments_),
            std::move(srtindxs_),
            std::move(logevalues_),
            std::move(alnptrs_),
            std::move(annotptrs_)
        );
    }
}

// -------------------------------------------------------------------------
// SortCompressedResults: sort formatted results;
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::SortCompressedResults()
{
    MYMSG( "CuBatchProcessingFinalizer::SortCompressedResults", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::SortCompressedResults: ";

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(preamb + "Null compressed results.");

    if( srtindxs_->size() !=logevalues_->size() ||
        logevalues_->size() != alnptrs_->size() ||
        logevalues_->size() != annotptrs_->size())
        throw MYRUNTIME_ERROR(preamb + "Inconsistent result sizes.");

    std::sort(srtindxs_->begin(), srtindxs_->end(),
        [this](size_t n1, size_t n2) {
            return (*logevalues_)[n1] < (*logevalues_)[n2];
        });
}



// -------------------------------------------------------------------------
// PrintCompressedResults: print formatted alignments
//
void CuBatchProcessingFinalizer::PrintCompressedResults() const
{
    MYMSG( "CuBatchProcessingFinalizer::PrintCompressedResults", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::PrintCompressedResults: ";

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(preamb + "Null compressed results.");

    if( srtindxs_->size() !=logevalues_->size() ||
        logevalues_->size() != alnptrs_->size() ||
        logevalues_->size() != annotptrs_->size())
        throw MYRUNTIME_ERROR(preamb + "Inconsistent result sizes.");

    for(size_t i = 0; i < srtindxs_->size(); i++ ) {
        fprintf(stdout,"%s",(*annotptrs_)[(*srtindxs_)[i]]);
    }

    fprintf(stdout,"%s",NL);

    for(size_t i = 0; i < srtindxs_->size(); i++ ) {
        fprintf(stdout,"%f%s",(*logevalues_)[(*srtindxs_)[i]],NL);
        fprintf(stdout,"%s",(*alnptrs_)[(*srtindxs_)[i]]);
    }
}
