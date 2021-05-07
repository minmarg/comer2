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
#include <functional>
#include <algorithm>
#include <numeric>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "extsp/psl.h"
#include "liblib/fmtdescription.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"

#include "AlnWriter.h"

// _________________________________________________________________________
// Class AlnWriter
//
// Constructor
//
AlnWriter::AlnWriter( 
    Configuration* config,
    const char* outdirname,
    const std::vector<std::string>& dbnamelist,
    size_t prodbsize,
    size_t ndbentries,
    int nqueries )
:   tobj_(NULL),
    req_msg_(WRITERTHREAD_MSG_UNSET),
    rsp_msg_(WRITERTHREAD_MSG_UNSET),
    refLambda_(-1.0f),
    refK_(-1.0f),
    expGappedLambda_(-1.0f),
    expGappedK_(-1.0f),
    //
    mstr_set_outdirname_(outdirname),
    //mstr_set_dbnamelist_(dbnamelist),
    mstr_set_prodbsize_(prodbsize),
    mstr_set_ndbentries_(ndbentries),
    //
    vec_nqyposs_(nqueries,0),
    vec_qryname_(nqueries,NULL),
    vec_qrydesc_(nqueries,NULL),
    vec_deltalen_(nqueries,0U),
    vec_sspace_(nqueries,0.0f),
    vec_logevthld_(nqueries,0.0f),
    vec_annotations_(nqueries),
    vec_alignments_(nqueries),
    vec_srtindxs_(nqueries),
    vec_logevalues_(nqueries),
    vec_alnptrs_(nqueries),
    vec_annotptrs_(nqueries),
    //
    p_finalindxs_(NULL),
    //
    //let the grand master decide on triggering the write by 
    // initializing query parts to 1
    queryparts_(nqueries,1),
    qrysernr_(-1)
{
    MYMSG( "AlnWriter::AlnWriter", 3 );

    if( config ) {
        const Configuration& ungapped_config = config[CTUngapped];
        const Configuration& gapped_config = config[CTGapped];
        //
        refLambda_ = ungapped_config.GetLambda();
        refK_ = ungapped_config.GetK();
        expGappedLambda_ = gapped_config.GetLambda();
        expGappedK_ = gapped_config.GetK();
    }

    for(const std::string& dbn: dbnamelist) {
        //std::string::npos+1==0
        std::string::size_type pos = dbn.rfind(DIRSEP) + 1;
        mstr_set_dbnamelist_.push_back(dbn.substr(pos));
    }

    tobj_ = new std::thread( &AlnWriter::Execute, this, (void*)NULL );
}

// Default constructor
//
AlnWriter::AlnWriter()
:   tobj_(NULL)
{
    throw MYRUNTIME_ERROR(
    "AlnWriter::AlnWriter: Default initialization is prohibited.");
}

// Destructor
//
AlnWriter::~AlnWriter()
{
    MYMSG( "AlnWriter::~AlnWriter", 3 );
    if( tobj_ ) {
        tobj_->join();
        delete tobj_;
        tobj_ = NULL;
    }
}

// -------------------------------------------------------------------------
// Execute: thread's starting point and execution process
//
void AlnWriter::Execute( void* )
{
    MYMSG( "AlnWriter::Execute", 3 );
    myruntime_error mre;

    try {
        while(1) {
            //wait until a message arrives
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);

            cv_msg_.wait(lck_msg,
                [this]{return 
                    ((0 <= req_msg_ && req_msg_ <= wrtthreadmsgTerminate) || 
                     req_msg_ == WRITERTHREAD_MSG_ERROR
                    );}
            );

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "AlnWriter::Execute: Msg %d",req_msg_);
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //thread owns the lock after the wait;
            //read message req_msg_
            int reqmsg = req_msg_;

            //unset the message to avoid live cycle when starting over the loop
            req_msg_ = WRITERTHREAD_MSG_UNSET;

            //set response msg to error upon exception
            rsp_msg_ = WRITERTHREAD_MSG_ERROR;
            int rspmsg = rsp_msg_;

            switch(reqmsg) {
                case wrtthreadmsgWrite:
                        ;;
                        if((int)queryparts_.size() <= qrysernr_ || qrysernr_ < 0 )
                            throw MYRUNTIME_ERROR(
                            "AlnWriter::Execute: Invalid query serial number.");
                        MergeResults();
                        WriteResults();
                        ReleaseAllocations();
                        qrysernr_ = -1;
                        ;;
                        //parent does not wait for a response nor requires data to read;
                        //unset response code
                        rspmsg = WRITERTHREAD_MSG_UNSET;
                        break;
                case wrtthreadmsgTerminate:
                        rspmsg = wrttrespmsgTerminating;
                        break;
                default:
                        rspmsg = WRITERTHREAD_MSG_UNSET;
                        break;
            };

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "AlnWriter::Execute: Msg %d Rsp %d",reqmsg, rspmsg );
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //save response code
            rsp_msg_ = rspmsg;

            //unlock the mutex and notify waiting threads
            lck_msg.unlock();
            cv_msg_.notify_all();

            if( reqmsg < 0 || reqmsg == wrtthreadmsgTerminate)
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
        cv_msg_.notify_all();
        return;
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MergeResults: merge parts of results
//
void AlnWriter::MergeResults()
{
    MYMSG( "AlnWriter::MergeResults", 4 );
    mystring preamb = "AlnWriter::MergeResults: ";

    InitializeVectors();

    if( vec_srtindxs_[qrysernr_].size() != vec_annotations_[qrysernr_].size() ||
        vec_srtindxs_[qrysernr_].size() != vec_alignments_[qrysernr_].size() ||
        vec_srtindxs_[qrysernr_].size() != vec_logevalues_[qrysernr_].size() ||
        vec_srtindxs_[qrysernr_].size() != vec_alnptrs_[qrysernr_].size() ||
        vec_srtindxs_[qrysernr_].size() != vec_annotptrs_[qrysernr_].size())
        throw MYRUNTIME_ERROR(preamb + "Inconsistent result sizes.");

    //partially merged and new merged vectors:
    std::vector<int>* p_prtmerged = &finalsrtindxs_;
    std::vector<int>* p_newmerged = &finalsrtindxs_dup_;
    p_finalindxs_ = NULL;

    //number of entries merged already
    size_t nmerged = 0UL;

    for(size_t i = 0; i < vec_srtindxs_[qrysernr_].size(); i++) {
        if( !vec_srtindxs_[qrysernr_][i] || !vec_logevalues_[qrysernr_][i])
            continue;
        if( vec_srtindxs_[qrysernr_][i]->size() < 1 )
            continue;
        std::merge(
            p_prtmerged->begin(), p_prtmerged->begin() + nmerged,
            p_prtmerged->begin() + nmerged, 
                    p_prtmerged->begin() + nmerged + vec_srtindxs_[qrysernr_][i]->size(),
            p_newmerged->begin(),
            [this](size_t n1, size_t n2) {
                return
                    ( *vec_logevalues_[qrysernr_][allsrtvecs_[n1]] )[allsrtindxs_[n1]] < 
                    ( *vec_logevalues_[qrysernr_][allsrtvecs_[n2]] )[allsrtindxs_[n2]];
            }
        );
        //save the pointer pointing to the results
        p_finalindxs_ = p_newmerged;
        //swap pointers:
        p_newmerged = p_prtmerged;
        p_prtmerged = p_finalindxs_;
        //number of merged entries has increased:
        nmerged += vec_srtindxs_[qrysernr_][i]->size();
    }
}

// -------------------------------------------------------------------------
// InitializeVectors: initialize index vectors
inline
void AlnWriter::InitializeVectors()
{
    MYMSG( "AlnWriter::InitializeVectors", 5 );
    int nrecords = GetTotalNumberOfRecords();

    allsrtindxs_.reserve(nrecords);
    allsrtvecs_.reserve(nrecords);
    finalsrtindxs_.reserve(nrecords);
    finalsrtindxs_dup_.reserve(nrecords);

    allsrtindxs_.clear();
    allsrtvecs_.clear();
    finalsrtindxs_.clear();
    finalsrtindxs_dup_.clear();

    for(size_t i = 0; i < vec_srtindxs_[qrysernr_].size(); i++) {
        if( !vec_srtindxs_[qrysernr_][i] || vec_srtindxs_[qrysernr_][i]->size() < 1 )
            continue;
        allsrtindxs_.insert(allsrtindxs_.end(),
            vec_srtindxs_[qrysernr_][i]->begin(), vec_srtindxs_[qrysernr_][i]->end());
        allsrtvecs_.insert(allsrtvecs_.end(), vec_srtindxs_[qrysernr_][i]->size(), (int)i);
        //sequentially enumerate all records:
        //std::vector<int>::iterator finalprevend = finalsrtindxs_.end();
        int szfinal = (int)finalsrtindxs_.size();
        finalsrtindxs_.resize(szfinal + vec_srtindxs_[qrysernr_][i]->size());
        //std::iota(finalprevend, finalsrtindxs_.end(), szfinal);
        std::iota(finalsrtindxs_.begin() + szfinal, finalsrtindxs_.end(), szfinal);
    }

    finalsrtindxs_dup_ = finalsrtindxs_;
}



// =========================================================================
// BufferData: buffer data and write the buffer contents to a file when 
// it is full;
// fp, file pointer;
// buffer, buffer to store data in;
// szbuffer, size of the buffer;
// outptr, varying address of the pointer pointing to a location in the 
//  buffer;
// offset, outptr offset from the beginning of the buffer (filled size);
// data, data to store in the buffer;
// szdata, size of the data;
//
void AlnWriter::BufferData( 
    FILE* fp, 
    char* const buffer, const int szbuffer, char*& outptr, int& offset, 
    const char* data, int szdata )
{
    MYMSG( "AlnWriter::BufferData", 9);
    while( szdata > 0 ) {
        if( szbuffer <= offset + szdata ) {
            int left = szbuffer - offset;
            if( left >= 0 ) {
                if( left)
                    strncpy(outptr, data, left);
                WriteToFile( fp, buffer, szbuffer );//WRITE TO FILE
                data += left;
                szdata -= left;
                outptr = buffer;
                offset = 0;
            }
        }
        if( offset + szdata <= szbuffer && szdata > 0 ) {
            strncpy(outptr, data, szdata);
            outptr += szdata;
            offset += szdata;
            szdata = 0;
        }
    }
}

// -------------------------------------------------------------------------
// WriteToFile: write data to file
//
void AlnWriter::WriteToFile( FILE* fp, char* data, int szdata )
{
    if( fwrite(data, sizeof(char), szdata, fp) != (size_t)szdata )
        throw MYRUNTIME_ERROR("AlnWriter::WriteToFile: write to a file failed.");
}

// -------------------------------------------------------------------------
// GetOutputFilename: make a name for the output file of alignments;
// outfilename, filename to be constructed;
// outdirname, output directory name given;
// qryname, query name;
// qrynr, query serial number;
//
void AlnWriter::GetOutputFilename( 
    mystring& outfilename,
    const char* outdirname,
    const char* qryname,
    const int qrynr )
{
    size_t pos;
    int outfmt = CLOptions::GetB_FMT();
    const char* ext = (outfmt==CLOptions::ofJSON)? "json": "out";
    char tail[20];
    sprintf( tail, "__%d.%s", qrynr,ext);
    outfilename = outdirname;
    if( !outfilename.empty() && outfilename[outfilename.length()-1] != DIRSEP )
        outfilename += DIRSEP;
    mystring qn = my_basename(qryname);
    if((pos = qn.rfind('.')) != mystring::npos && qn.length() - pos <= 4 )
        qn = qn.substr(0,pos);
    outfilename += qn + tail;
}

