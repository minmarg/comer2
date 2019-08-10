/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchProcessingFinalizer_h__
#define __CuBatchProcessingFinalizer_h__

#include "liblib/mybase.h"

#include <stdio.h>
#include <math.h>

#include <memory>
#include <functional>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
// #include "libmycu/cupro/CuBatchProcessing.cuh"

#define CUBPTHREAD_MSG_UNSET -1
#define CUBPTHREAD_MSG_ERROR -2

// _________________________________________________________________________
// Class CuBatchProcessingFinalizer
//
// thread class for finalizing results calculated on a device (GPU)
//
class CuBatchProcessingFinalizer
{
public:
    enum TCUBPThreadMsg {
        cubpthreadmsgFinalize,
        cubpthreadmsgTerminate
    };
    enum TCUBPThreadResponse {
        cubptrespmsgFinalizing,
        cubptrespmsgTerminating
    };

public:
    CuBatchProcessingFinalizer(
        cudaStream_t& strcopyres,
        DeviceProperties dprop, 
        AlnWriter*,
        Configuration*,
        const mystring* queryfnames,
        const mystring* querydescs,
        const mystring* bdb1fnames,
        const mystring* bdb1descs
    );
    CuBatchProcessingFinalizer();
    ~CuBatchProcessingFinalizer();

    std::mutex& GetPrivateMutex() {return mx_dataccess_;}

    //{{NOTE: messaging functions accessed from outside!
    void Notify(int msg) {
        {//mutex must be unlocked before notifying
            std::lock_guard<std::mutex> lck(mx_dataccess_);
            req_msg_ = msg;
        }
        cv_msg_.notify_one();
    }
    int GetResponse() const {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        return rsp_msg_;
    }
    void ResetResponse() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        if( rsp_msg_!= CUBPTHREAD_MSG_ERROR )
            rsp_msg_ = CUBPTHREAD_MSG_UNSET;
    }
    //}}


    void SetCuBPBDbdata( 
        int qrysernr,
        size_t nqyposs, float qyeno,
        unsigned int deltalen,
        float sspace,
        float logevthld,
        std::unique_ptr<PMBatchProData> bdbC,
        char** querypmbeg, char** querypmend, 
        char** bdb1pmbeg, char** bdb1pmend,
        char** bdbCpmbeg, char** bdbCpmend,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t querprosOmtd,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        unsigned int nposits,
        unsigned int npros,
        const char* h_results,
        size_t sz_alndata, size_t sz_alns,
        int dbalnlen2 )
    {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        cubp_set_qrysernr_ = qrysernr;
        cubp_set_nqyposs_ = (int)nqyposs;
        cubp_set_qyeno_ = qyeno;
        cubp_set_deltalen_ = deltalen;
        cubp_set_sspace_ = sspace;
        cubp_set_logevthld_ = logevthld;
        cubp_set_bdbC_ = std::move(bdbC);
        memcpy( cubp_set_querypmbeg_, querypmbeg, pmv2DTotFlds * sizeof(void*));
        memcpy( cubp_set_querypmend_, querypmend, pmv2DTotFlds * sizeof(void*));
        if( bdb1pmbeg && bdb1pmend ) {
            memcpy( cubp_set_bdb1pmbeg_, bdb1pmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( cubp_set_bdb1pmend_, bdb1pmend, pmv2DTotFlds * sizeof(void*));
        }
        if( bdbCpmbeg && bdbCpmend ) {
            memcpy( cubp_set_bdbCpmbeg_, bdbCpmbeg, pmv2DTotFlds * sizeof(void*));
            memcpy( cubp_set_bdbCpmend_, bdbCpmend, pmv2DTotFlds * sizeof(void*));
        }
        cubp_set_ndb1pros_ = ndb1pros;
        cubp_set_ndbCpros_ = ndbCpros;
        cubp_set_querprosOmtd_ = querprosOmtd;
        cubp_set_ndb1prosOmtd_ = ndb1prosOmtd;
        cubp_set_ndbCprosOmtd_ = ndbCprosOmtd;
        cubp_set_querposoffset_ = querposoffset;
        cubp_set_bdb1posoffset_ = bdb1posoffset;
        cubp_set_bdbCposoffset_ = bdbCposoffset;
        cubp_set_nposits_ = nposits;
        cubp_set_npros_ = npros;
        cubp_set_h_results_ = h_results;
        cubp_set_sz_alndata_ = sz_alndata;
        cubp_set_sz_alns_ = sz_alns;
        cubp_set_sz_dbalnlen2_ = dbalnlen2;
    }


protected:
    void Execute( void* args );

    void SetResponseError() {
        std::lock_guard<std::mutex> lck(mx_dataccess_);
        rsp_msg_ = CUBPTHREAD_MSG_ERROR;
    }

    void CompressResults();
    void SortCompressedResults();
    void PassResultsToWriter();
    void GetSizeOfCompressedResults(
        size_t* szannot, size_t* szalns, size_t* szalnswodesc) const;
    void PrintCompressedResults() const;


    float GetBitScore( float logeval ) {
        return cubp_set_sspace_ <= 0.0f? 0.0f:
            (logf(cubp_set_sspace_)-logeval)/SLC_LN2;
    }
    double GetEvalue( float logeval ) {
        return SLC_LOG_DP_MAX < logeval? 1.0e+308: exp((double)logeval);
    }
    double GetPvalue( double evalue ) {
        return evalue < 0.01? evalue: 1.0 - exp(-evalue);
    }


    void MakeAnnotation( 
        char*& outptr,
        const mystring* name,
        const mystring* desc,
        const bool printname,
        const int maxoutlen,
        const int width,
        const float score,
        const double evalue) const;

    void FormatScores(
        char*& outptr,
        unsigned int prondx,
        unsigned int orgprondx,
        unsigned int alnlen,
        float score,
        float logeval,
        double evalue,
        int dbprolen );

    void FormatAlignment(
        char*& outptr,
        unsigned int prondx,
        unsigned int orgprondx,
        unsigned int dbpro2dst,
        int alnlen,
        int dbprolen,
        const int width,
        bool printsss,
        const bool qrysssinuse );

    void FormatFooter(
        char*& outptr, unsigned int prondx );

    
    const char* GetBegOfAlns() {
        return cubp_set_h_results_ + cubp_set_sz_alndata_;
    }

    const char* GetAlnSectionAt( const char* ptr, const int sctndx ) const {
        return ptr + sctndx * cubp_set_sz_dbalnlen2_;
    }

    template<typename T>
    T GetOutputAlnDataField(unsigned int prondx, unsigned int field) const
        {
            return *(T*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+field);
        }

    template<typename T>
    T GetProfileField( char* const pmbeg[pmv2DTotFlds],
        unsigned int ndx, unsigned int field) const
        {
            return ((T*)(pmbeg[field]))[ndx];
        }

    template<typename T>
    T GetQueryFieldPos(unsigned int field, unsigned int pos) const;
    template<typename T>
    T GetDbProfileField(unsigned int orgprondx, unsigned int field) const;
    template<typename T>
    T GetDbProfileFieldPos(unsigned int orgprondx, unsigned int field, unsigned int pos) const;

    void GetDbProfileNameDesc(
        const mystring*& name, const mystring*& desc, unsigned int orgprondx) const;


    void ReserveVectors( int capacity ) {
        srtindxs_.reset( new std::vector<int>);
        logevalues_.reset( new std::vector<float>);
        alnptrs_.reset( new std::vector<char*>);
        annotptrs_.reset( new std::vector<char*>);
        if( capacity > 0 ) {
            if(srtindxs_) srtindxs_->reserve(capacity);
            if(logevalues_) logevalues_->reserve(capacity);
            if(alnptrs_) alnptrs_->reserve(capacity);
            if(annotptrs_) annotptrs_->reserve(capacity);
        }
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
    //
    //properties of device the thread is associated with:
    cudaStream_t& strcopyres_;
    DeviceProperties dprop_;
    //results writer and configuration object:
    AlnWriter* alnwriter_;
    Configuration* config_;
    //{{cached host data
    const mystring* cached_queryfnames_;
    const mystring* cached_querydescs_;
    const mystring* cached_bdb1fnames_;
    const mystring* cached_bdb1descs_;
    //}}
    //{{data arguments: 
    // cubp-set data/addresses:
    int cubp_set_qrysernr_;//query serial number
    int cubp_set_nqyposs_;//query length
    float cubp_set_qyeno_;//query ENO
    unsigned int cubp_set_deltalen_;//length adjustment
    float cubp_set_sspace_;//search space size
    float cubp_set_logevthld_;//log e-value threshold
    std::unique_ptr<PMBatchProData> cubp_set_bdbC_;
    char* cubp_set_querypmbeg_[pmv2DTotFlds];
    char* cubp_set_querypmend_[pmv2DTotFlds];
    char* cubp_set_bdb1pmbeg_[pmv2DTotFlds];
    char* cubp_set_bdb1pmend_[pmv2DTotFlds];
    char* cubp_set_bdbCpmbeg_[pmv2DTotFlds];
    char* cubp_set_bdbCpmend_[pmv2DTotFlds];
    size_t cubp_set_ndb1pros_;
    size_t cubp_set_ndbCpros_;
    size_t cubp_set_querprosOmtd_;
    size_t cubp_set_ndb1prosOmtd_;
    size_t cubp_set_ndbCprosOmtd_;
    size_t cubp_set_querposoffset_;
    size_t cubp_set_bdb1posoffset_;
    size_t cubp_set_bdbCposoffset_;
    unsigned int cubp_set_nposits_;//total number of positions in the transfered results
    unsigned int cubp_set_npros_;//number of profiles in the transfered results
    const char* cubp_set_h_results_;//cubp-set host-side results
    size_t cubp_set_sz_alndata_;//size of alignment data of results
    size_t cubp_set_sz_alns_;//size of alignments of results
    int cubp_set_sz_dbalnlen2_;
    //}}
    //{{formatted results:
    std::unique_ptr<char,MyDataDestroyer> annotations_;
    std::unique_ptr<char,MyDataDestroyer> alignments_;
    std::unique_ptr<std::vector<int>> srtindxs_;//index vector of sorted log e-values
    std::unique_ptr<std::vector<float>> logevalues_;//vector of log e-values
    std::unique_ptr<std::vector<char*>> alnptrs_;//vector of alignments
    std::unique_ptr<std::vector<char*>> annotptrs_;//vector of annotations
    //}}
};

// -------------------------------------------------------------------------
// GetQueryFieldPos: get a field of the query profile at a given position;
// NOTE: the pointer of the query structure under process is accessed!
template<typename T>
inline
T CuBatchProcessingFinalizer::GetQueryFieldPos(
    unsigned int field, unsigned int pos) const
{
    return GetProfileField<T>(cubp_set_querypmbeg_, pos, field);
}
// GetDbProfileField: get a field of the given Db profile;
// orgprondx, index of a profile over all pm data structures;
template<typename T>
inline
T CuBatchProcessingFinalizer::GetDbProfileField(
    unsigned int orgprondx, unsigned int field) const
{
    if( orgprondx < cubp_set_ndb1pros_) {
        //NOTE: cubp_set_bdb1pmbeg_ is already shifted:
        // do not add # omitted profiles!
        //orgprondx += cubp_set_ndb1prosOmtd_;
        return GetProfileField<T>(cubp_set_bdb1pmbeg_, orgprondx, field);
    } else {
        //NOTE: same as above!
        //orgprondx += /*cubp_set_ndbCprosOmtd_ */- cubp_set_ndb1pros_;
		orgprondx -= (unsigned int)cubp_set_ndb1pros_;
		return GetProfileField<T>(cubp_set_bdbCpmbeg_, orgprondx, field);
    }
}
// GetDbProfileFieldPos: get a field of the given Db profile at a given position;
// orgprondx, index of a profile over all pm data structures;
template<typename T>
inline
T CuBatchProcessingFinalizer::GetDbProfileFieldPos(
    unsigned int orgprondx, unsigned int field, unsigned int pos) const
{
    if( orgprondx < cubp_set_ndb1pros_) {
        return GetProfileField<T>(cubp_set_bdb1pmbeg_, pos, field);
    } else {
        return GetProfileField<T>(cubp_set_bdbCpmbeg_, pos, field);
    }
}

// GetDbProfileNameDesc: get the db profile name and description;
// orgprondx, index of a profile over all pm data structures;
inline
void CuBatchProcessingFinalizer::GetDbProfileNameDesc(
    const mystring*& name, const mystring*& desc,
    unsigned int orgprondx) const
{
    if( orgprondx < cubp_set_ndb1pros_) {
        orgprondx += (unsigned int)cubp_set_ndb1prosOmtd_;
        name = cached_bdb1fnames_ + orgprondx;
        desc = cached_bdb1descs_ + orgprondx;
    } else {
        orgprondx += (unsigned int)(cubp_set_ndbCprosOmtd_ - cubp_set_ndb1pros_);
        name = cubp_set_bdbC_->GetFnames() + orgprondx;
        desc = cubp_set_bdbC_->GetDescs() + orgprondx;
    }
}

#endif//__CuBatchProcessingFinalizer_h__
