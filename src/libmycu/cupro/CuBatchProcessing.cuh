/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchProcessing_h__
#define __CuBatchProcessing_h__

#include "liblib/mybase.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <memory>
#include <mutex>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "tsafety/TSCounterVar.h"

#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cusco/CuBatchScoreMatrix.cuh"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "libmycu/cuss/CuBatchSS.cuh"
#include "libmycu/cumapdp/CuBatchMAPDP_com.h"
#include "libmycu/cumapdp/CuBatchMAPDP.cuh"
#include "CuBatchProcessingFinalizer.h"
#include "CuBatchProcessing_com.cuh"
#include "CuDeviceMemory.cuh"

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchProcessing
// Batch computation of profile-profile substitution scores using the CUDA 
// architecture or CPU threads
//
class CuBatchProcessing
{
public:
    CuBatchProcessing(
        CuDeviceMemory* dmem,
        int dareano,
        Configuration config[NoCTypes], 
        AlnWriter* writer,
        const mystring* queryfnames,
        const mystring* querydescs,
        const char** bdb1descs,
        int precscale = 1 );

    virtual ~CuBatchProcessing();

    void TransferCPMData(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields )
    {
        dmem_->TransferCPMData(
            bdbCpmbeg,
            bdbCpmend,
            szCpm2dvfields);
    }

    void SetDbDetails( uint64_mt db_len, size_t no_sequences ) {
        dblength_ = db_len;
        ndbsequences_ = no_sequences;
    }

    size_t GetQueryLen() const { return querylength_; }
    void SetQueryLen( size_t value ) { querylength_ = value; }

    float GetScoreThreshold() const { return scorethld_; }
    void SetScoreThreshold( float value ) { scorethld_ = value; }

    float GetLogEThreshold() const { return logevthld_; }
    void SetLogEThreshold( float value ) { logevthld_ = value; }

    void ProcessScoreMatrix(
        int qrysernr,
        size_t nqyposs,
        char** querypmbeg,
        //char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        const char** bdbCdesc,
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields,
        TSCounterVar* cnt
    );

    void WaitForIdleChilds() {
        if( cbpfin_ )
            std::lock_guard<std::mutex> lck(cbpfin_->GetPrivateMutex());
    }

    size_t GetCurrentMaxDbPos() const { return dmem_->GetCurrentMaxDbPos(); }
    size_t GetCurrentMaxNDbPros() const { return dmem_->GetCurrentMaxNDbPros(); }
    unsigned int GetCurrentDbxPadding() const { return curdbxpad_; }

    size_t GetCurrentMaxDbPosPass2() const { return dmem_->GetCurrentMaxDbPosPass2(); }
    size_t GetCurrentMaxNDbProsPass2() const { return dmem_->GetCurrentMaxNDbProsPass2(); }
    unsigned int GetCurrentDbxPaddingPass2() const { return dbxpadphase2_; }

protected:
    explicit CuBatchProcessing();

    void ResetDbDetails() {
        querylength_ = 0UL;
        dblength_ = 0UL;
        ndbsequences_ = 0UL;
    }

    void BatchProcessScoreMatrixDevice(
        int qrysernr,
        size_t nqyposs,
        char** querypmbeg,
        //char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        const char** bdbCdesc,
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields,
        TSCounterVar* cnt
    );

    void DestroyCBSMObject(){ if(cbsm_) delete cbsm_; }
    void DestroyCBDPObject(){ if(cbdp_) delete cbdp_; }
    void DestroyCBSSObject(){ if(cbss_) delete cbss_; }
    void DestroyCBMAPDPObject(){ if(cbmapdp_) delete cbmapdp_; }

    void SetCurrentDbxPadding( unsigned int value ) { curdbxpad_ = value; }

    void SetCurrentDbxPaddingPass2( unsigned int value ) { dbxpadphase2_ = value; }
    void SetCurrentDbAlnLengthWithPaddingPass2( unsigned int value ) { dbalnlenphase2_ = value; }
    unsigned int GetCurrentDbAlnLengthWithPaddingPass2() const { return dbalnlenphase2_; }

    size_t GetOffsetOfHeapSection(int sec) const {return dmem_->GetOffsetOfHeapSection(devareano_,sec);}
    char* GetAddrOfHeapSection(int sec) const {
        return dmem_->GetHeap()+dmem_->GetOffsetOfHeapSection(devareano_,sec);
    }
    char* GetAddrOfHeapSS2Data() const {
        return dmem_->GetHeap()+dmem_->GetEndOfSS2Data(devareano_);
    }

    void TESTPrintProProScores1(
        char** querypmbeg, char** querypmend,
        char** bdb1pmbeg, char** bdb1pmend, char** bdbCpmbeg, char** bdbCpmend );


    //{{ ===== results section =====
    void TransferResultsDevice(
        int qrysernr,
        size_t nqyposs, float qyeno,
        char** querypmbeg,
        //char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        const char** bdbCdesc,
        char** bdbCpmbeg,
        char** bdbCpmend,
        TSCounterVar* cnt,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t querprosOmtd, size_t ndb1prosOmtd, size_t ndbCprosOmtd,
        size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
        unsigned int nposits,
        unsigned int npros
    );
    void CheckHostResultsSync();
    void HostFreeResults()
    {
        MYMSG("CuBatchProcessing::HostFreeResults",4);
        if(h_results_) {
            if( lockedresmem_ ) {
                MYCUDACHECK( cudaFreeHost(h_results_));
                MYCUDACHECKLAST;
            }
            else
                free(h_results_);
        }
        h_results_ = NULL;
        lockedresmem_ = false;
    }
    void HostAllocResults( size_t szresults )
    {
        HostFreeResults();
        lockedresmem_ = true;
        if( cudaSuccess !=
            cudaHostAlloc((void**)&h_results_, szresults, cudaHostAllocDefault)) {
            cudaGetLastError();
            h_results_ = NULL;
            lockedresmem_ = false;
            h_results_ = (char*)malloc(szresults);
            if( h_results_ == NULL )
                throw MYRUNTIME_ERROR(
                "CuBatchProcessing::HostAllocResults: Not enough memory.");
        }
        MYMSGBEGl(4)
            char msgbuf[BUF_MAX];
            sprintf(msgbuf, "CuBatchProcessing::HostAllocResults: sz %zu lckd %d",
                szresults,lockedresmem_);
            MYMSG( msgbuf,4 );
        MYMSGENDl
    }
    //}}


private:
    CuDeviceMemory* dmem_;//device memory configuration
    const int devareano_;//device area number
    CuBatchScoreMatrix* cbsm_;//batch score matrix
    CuBatchDP* cbdp_;//batch dynamic programming object
    CuBatchSS* cbss_;//object for the batch calculation of statistical significance
    CuBatchMAPDP* cbmapdp_;//batch MAP dynamic programming object
    CuBatchProcessingFinalizer* cbpfin_;//results finalizer
    float scorethld_;//score threshold for pass-2 alignments
    float logevthld_;//log E-value threshold
    //{{db details
    size_t querylength_;
    uint64_mt dblength_;
    size_t ndbsequences_;
    //}}
    unsigned int curdbxpad_;//current padding in positions along the x (db) profile positions
    unsigned int dbxpadphase2_;//padding in positions along the x (db) axis in phase 2
    unsigned int dbalnlenphase2_;//total alignment length with padding in phase 2
    //{{host allocated pointers for results
    cudaStream_t streamcopyres_;//stream for copying results
    char* h_results_;//results received from device
    bool lockedresmem_;//page-locked memory for results; also specifies that a stream is initialized
    size_t sz_mem_results_;//size of allocated memory for results
    size_t limit_beg_results_;//the boundary where the results sections begin
    //}}
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__CuBatchProcessing_h__
