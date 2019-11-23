/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __AlnWriter_h__
#define __AlnWriter_h__

#include "liblib/mybase.h"

#include <stdio.h>
// #include <math.h>
#include <cmath>

#include <memory>
#include <utility>
#include <functional>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "libpro/srcpro/Configuration.h"

#define WRITERTHREAD_MSG_UNSET -1
#define WRITERTHREAD_MSG_ERROR -2

struct MyDataDestroyer {
    void operator()(char* p) const {
        std::free(p);
    };
};

// _________________________________________________________________________
// Class AlnWriter
//
// alignment writer thread
//
class AlnWriter
{
    enum {
        szTmpBuffer = KBYTE,
        szWriterBuffer = TIMES4(KBYTE)
    };

public:
    enum TWriterThreadMsg {
        wrtthreadmsgWrite,
        wrtthreadmsgTerminate
    };
    enum TWriterThreadResponse {
        wrttrespmsgWriting,
        wrttrespmsgTerminating
    };

public:
    AlnWriter( 
        Configuration*,
        const char* outdirname,
        const char* dbname,
        size_t prodbsize,
        size_t ndbentries,
        int nqueries
    );
    AlnWriter();
    ~AlnWriter();

    std::mutex& GetPrivateMutex() {return mx_dataccess_;}

    //{{NOTE: messaging functions accessed from outside!
    void Notify(int msg) {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_dataccess_);
            req_msg_ = msg;
        }
        cv_msg_.notify_all();
    }
    int Wait(int rsp) {
        //wait until a response
        std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
        cv_msg_.wait(lck_msg,
            [this,rsp]{return (rsp_msg_ == rsp || rsp_msg_ == WRITERTHREAD_MSG_ERROR);}
        );
        //lock is back; unset the response
        int rspmsg = rsp_msg_;
        if( rsp_msg_!= WRITERTHREAD_MSG_ERROR )
            rsp_msg_ = WRITERTHREAD_MSG_UNSET;
        return rspmsg;
    }
    int WaitDone() {
        for( size_t i = 0; i < queryparts_.size(); )
        {
            //wait until all queries have been processed
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);
            cv_msg_.wait(lck_msg,
                [this,i]{return 
                    queryparts_[i] <= 0 ||
                    rsp_msg_ == WRITERTHREAD_MSG_ERROR;}
            );
            if( rsp_msg_ == WRITERTHREAD_MSG_ERROR )
                return rsp_msg_;
            if( req_msg_ != WRITERTHREAD_MSG_UNSET )
                continue;
            i++;
        }
        return rsp_msg_;
    }
    void IncreaseQueryNParts( int qrysernr ) {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
            queryparts_[qrysernr]++;
    }
    void DereaseNPartsAndTrigger( int qrysernr ) {
        std::unique_lock<std::mutex> lck(mx_dataccess_);
        cv_msg_.wait(lck,
            [this]{return req_msg_ == WRITERTHREAD_MSG_UNSET;}
        );
        if( --queryparts_[qrysernr] <= 0 ) {
            //this is the last part for the given query:
            //trigger write to a file
            qrysernr_ = qrysernr;
            req_msg_ = wrtthreadmsgWrite;
            lck.unlock();
            cv_msg_.notify_all();
        }
    }
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return rsp_msg_;
    }
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if( rsp_msg_!= WRITERTHREAD_MSG_ERROR )
            rsp_msg_ = WRITERTHREAD_MSG_UNSET;
    }
    //}}


    void PushPartOfResults( 
        int qrysernr,
        int nqyposs,
        const char* qryname,
        const char* qrydesc,
        unsigned int deltalen,
        float sspace,
        float logevthld,
        std::unique_ptr<char,MyDataDestroyer> annotations,
        std::unique_ptr<char,MyDataDestroyer> alignments,
        std::unique_ptr<std::vector<int>> srtindxs,
        std::unique_ptr<std::vector<float>> logevalues,
        std::unique_ptr<std::vector<char*>> alnptrs,
        std::unique_ptr<std::vector<char*>> annotptrs )
    {
        std::unique_lock<std::mutex> lck(mx_dataccess_);
        //{{NOTE: [inserted]
        cv_msg_.wait(lck,
            [this]{return req_msg_ == WRITERTHREAD_MSG_UNSET;}
        );
        //}}
        if((int)queryparts_.size() <= qrysernr || qrysernr < 0 )
            throw MYRUNTIME_ERROR(
            "AlnWriter::PushPartOfResults: Invalid query serial number.");
        vec_nqyposs_[qrysernr] = nqyposs;
        vec_qryname_[qrysernr] = qryname;
        vec_qrydesc_[qrysernr] = qrydesc;
        vec_deltalen_[qrysernr] = deltalen;
        vec_sspace_[qrysernr] = sspace;
        vec_logevthld_[qrysernr] = logevthld;
        //
        vec_annotations_[qrysernr].push_back(std::move(annotations));
        vec_alignments_[qrysernr].push_back(std::move(alignments));
        vec_srtindxs_[qrysernr].push_back(std::move(srtindxs));
        vec_logevalues_[qrysernr].push_back(std::move(logevalues));
        vec_alnptrs_[qrysernr].push_back(std::move(alnptrs));
        vec_annotptrs_[qrysernr].push_back(std::move(annotptrs));
        //{{NOTE: [commented out] the following statement must be the last
        //lck.unlock();
        //DereaseNPartsAndTrigger(qrysernr);
        //}}
        if( --queryparts_[qrysernr] <= 0 ) {
            //this is the last part for the given query:
            //trigger write to a file
            qrysernr_ = qrysernr;
            req_msg_ = wrtthreadmsgWrite;
            lck.unlock();
            cv_msg_.notify_all();
        }
    }


protected:
    void Execute( void* args );

    void SetResponseError() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        rsp_msg_ = WRITERTHREAD_MSG_ERROR;
    }

    void MergeResults();
    void WriteResults();
    void BufferData( 
        FILE* fp, 
        char* const buffer, const int szbuffer, char*& outptr, int& offset, 
        const char* data, int szdata );
    void WriteToFile( FILE* fp, char* data, int szdata );
    void GetOutputFilename( 
        mystring& outfilename,
        const char* outdirname,
        const char* qryname,
        const int qrynr );
    int WriteProgname( char*& outptr );
    int WriteQueryDescription( 
        char*& outptr,
        int maxsize,
        bool printname,
        const int qrylen,
        const char* name,
        const char* desc,
        const int width );
    int WriteSearchInformation( 
        char*& outptr,
        int maxsize,
        const char* dbname,
        const size_t dbsize,
        const size_t ndbentries,
        const float logevthrld,
        const int indent,
        const int annotlen,
        const bool found );
    int WriteSummary( 
        char*& outptr,
        const float refLambda, const float refK,
        const float expGappedLambda, const float expGappedK,
        const int qrylen,
        const size_t dbsize,
        const size_t ndbentries,
        const float sspace,
        const int deltalen );

    int GetTotalNumberOfRecords() const
    {
        if((int)queryparts_.size() <= qrysernr_ || qrysernr_ < 0 )
            throw MYRUNTIME_ERROR(
            "AlnWriter::GetTotoalNumberOfRecords: Invalid query serial number.");
        int ntot = 0;
        for(size_t i = 0; i < vec_srtindxs_[qrysernr_].size(); i++) {
            if( vec_srtindxs_[qrysernr_][i])
                ntot += (int)vec_srtindxs_[qrysernr_][i]->size();
        }
        return ntot;
    }

    void InitializeVectors();

    void ReleaseAllocations() {
        if((int)queryparts_.size() <= qrysernr_ || qrysernr_ < 0 )
            throw MYRUNTIME_ERROR(
            "AlnWriter::ReleaseAllocations: Invalid query serial number.");
        vec_annotations_[qrysernr_].clear();
        vec_alignments_[qrysernr_].clear();
        vec_srtindxs_[qrysernr_].clear();
        vec_logevalues_[qrysernr_].clear();
        vec_alnptrs_[qrysernr_].clear();
        vec_annotptrs_[qrysernr_].clear();
    }

private:
    //thread section
    std::thread* tobj_;//thread object
private:
    //{{messaging
    std::condition_variable cv_msg_;//condition variable for messaging
    mutable std::mutex mx_dataccess_;//mutex for accessing class data
    int req_msg_;//request message issued for thread
    int rsp_msg_;//private response message
    //}}
    //{{data from configuration:
    float refLambda_;//reference lambda parameter
    float refK_;//reference parameter K
    float expGappedLambda_;//experimental gapped lambda
    float expGappedK_;//experimental gapped K
    //}}
    const char* mstr_set_outdirname_;//output directory name
    const char* mstr_set_dbname_;//database name
    size_t mstr_set_prodbsize_;//profile database size in positions
    size_t mstr_set_ndbentries_;//number of database entries
    //{{query and summary data:
    std::vector<int> vec_nqyposs_;//query length
    std::vector<const char*> vec_qryname_;//query name
    std::vector<const char*> vec_qrydesc_;//query description
    std::vector<unsigned int> vec_deltalen_;//length adjustment
    std::vector<float> vec_sspace_;//search space size
    std::vector<float> vec_logevthld_;//log e-value threshold
    //}}
    //{{vectors of formatted results for QUERIES:
    std::vector<std::vector< std::unique_ptr<char,MyDataDestroyer> >> vec_annotations_;
    std::vector<std::vector< std::unique_ptr<char,MyDataDestroyer> >> vec_alignments_;
    std::vector<std::vector< std::unique_ptr<std::vector<int>> >> vec_srtindxs_;//index vectors of sorted log e-values
    std::vector<std::vector< std::unique_ptr<std::vector<float>> >> vec_logevalues_;//2D vector of log e-values for queries
    std::vector<std::vector< std::unique_ptr<std::vector<char*>> >> vec_alnptrs_;//2D vector of alignments for queries
    std::vector<std::vector< std::unique_ptr<std::vector<char*>> >> vec_annotptrs_;//2D vector of annotations for queries
    //}}
    //{{sorted indices over all parts (vectors) of results for a query:
    std::vector<int> allsrtindxs_;//indices along all vectors
    std::vector<int> allsrtvecs_;//corresponding vector indices (part numbers)
    std::vector<int> finalsrtindxs_;//globally (over all parts) sorted indices 
    std::vector<int> finalsrtindxs_dup_;//duplicate of globally sorted indices (for efficient memory management)
    std::vector<int>* p_finalindxs_;//pointer to the final vector of sorted indices
    //}}
    //buffer for writing to file:
    char buffer_[szWriterBuffer];
    std::vector<int> queryparts_;//vector of the number of parts for each query serial number
    int qrysernr_;//query serial number
};

// -------------------------------------------------------------------------

#endif//__AlnWriter_h__
