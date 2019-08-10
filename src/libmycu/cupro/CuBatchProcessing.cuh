/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchProcessing_h__
#define __CuBatchProcessing_h__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <memory>
#include <mutex>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"
#include "libHDP/HDPscores.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cusco/CuBatchScoreMatrix.cuh"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "libmycu/cuss/CuBatchSS.cuh"
#include "libmycu/cumapdp/CuBatchMAPDP_com.h"
#include "libmycu/cumapdp/CuBatchMAPDP.cuh"
#include "CuBatchProcessingFinalizer.h"
#include "CuBatchProcessing_com.cuh"

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchProcessing
// Batch computation of profile-profile substitution scores using the CUDA 
// architecture or CPU threads
//
class CuBatchProcessing
{
public:
    //relevant positions of sections
    enum TDevDataSections {
        ddsEndOfPadding,//end of padding to align data
        ddsEndOfSSSTable,//end of score table for secondary structure predictions
        ddsEndOfCVS2SMap,//end of a map between context vector scores and translated scores
        ddsEndOfHDPscores,//end of scores of HDP cluster membership predictions
        ddsEndOfConstantData=ddsEndOfHDPscores,
        //
        ddsEndOfCached,//end of cached profile model data in the allocated block of memory
        ddsEndOfDbChunk,//end of db profile model data chunk in the allocated block of memory
        ddsEndOfOrgScores,//end of the section of originally calculated scores
        ddsEndOfModScores,//end of the section of modular scores w/o certain terms
        //
        // --- phase 1 division ---
        //end of the section of temporary diagonal scores for dynamic programming matrices:
        ddsEndOfDPDiagScores,
        ddsEndOfDPBottomScores,//end of the section of bottom-line scores for dynamic programming matrices
        ddsEndOfDPMaxCoords,//end of the section of the coordinates of maximum alignment scores
        ddsEndOfDPBackTckData,//end of the section of DP backtracking data
        //
        nDevDataSectionsDiv1,
        //
        // --- phase 2 division ---
        ddsEndOfDP2DiagScores=nDevDataSectionsDiv1,//end of phase 2 diagonal scores
        ddsEndOfDP2BottomScores,//end of phase 2 bottom-line scores
//**/ddsEndOfDP2DiagBlkProbScales,//block-specific scale factors of forward probabilitie
//**/ddsEndOfDP2PrbScales,//final scale factors for forward/backward probabilities
        ddsEndOfDP2MaxCoords,//end of phase 2 section of coordinates of alignment scores
        ddsEndOfDP2BackTckData,//end of phase 2 section of backtracking data
        ddsEndOfDP2FwdMtx,//end of phase 2 section of forward probability matrix
//ddsEndOfDP2BwdMtx,//end of phase 2 section of backward probability matrix
        ddsEndOfSS2Data,//end of phase 2 data for calculation of statistical significance [profile-specific]
        ddsEndOfDP2AlnData,//end of phase 2 section of alignment data [profile-specific, output for host]
        ddsEndOfDP2Alns,//end of phase 2 section of drawn alignments [profile-specific, output for host]
        //
        nDevDataSectionsDiv2
    };
    enum TDevMatrixSections {
        dmsOrgScores,//originally calculated scores
        dmsModScores,//modular scores w/o certain terms
        nDevMatrixSections
    };

    CuBatchProcessing(
        AlnWriter* writer,
        Configuration config[NoCTypes], 
        DeviceProperties dprop, 
        bool hdp1scoresinuse, 
        bool mapalninuse, 
        size_t deviceallocsize, 
        const mystring* queryfnames,
        const mystring* querydescs,
        const mystring* bdb1fnames,
        const mystring* bdb1descs,
        int precscale = 1 );
    virtual ~CuBatchProcessing();

    void CacheSSENNWeights();
    void CacheSSSScores( const SSSScores& );
    void CacheCVS2Scores( const CVS2Scores& );
    void CacheHDPScores( const HDPscores& );
    //
    void CacheData(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend
    );

    void SetDbDetails( size_t query_len, uint64_mt db_len, size_t no_sequences ) {
        querylength_ = query_len;
        dblength_ = db_len;
        ndbsequences_ = no_sequences;
    }

    size_t CalcMaxDbDataChunkSize( size_t nqyposs );
    size_t CalcMaxDbDataChunkSizeObs(
        char** querypmbeg,
        char** querypmend
    );

    float GetScoreThreshold() const { return scorethld_; }
    void SetScoreThreshold( float value ) { scorethld_ = value; }

    float GetLogEThreshold() const { return logevthld_; }
    void SetLogEThreshold( float value ) { logevthld_ = value; }

    void ProcessScoreMatrix(
        int qrysernr,
        std::unique_ptr<PMBatchProData> bdbC,
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        char** bdbCpmbeg,
        char** bdbCpmend
    );

    void WaitForIdleChilds() {
        if( cbpfin_ )
            std::lock_guard<std::mutex> lck(cbpfin_->GetPrivateMutex());
    }

    bool GetHDP1ScoresInUse() const { return hdp1scoresinuse_; }
    bool GetMAPAlnInUse() const { return mapalninuse_; }
    bool GetModScoreMatrixInUse() const { 
        return MOptions::GetSSSWGT() > 0.0f || MOptions::GetCVSWGT() > 0.0f;
    }

    size_t GetCurrentMaxDbPos() const { return curmaxdbpos_; }
    size_t GetCurrentMaxNDbPros() const { return curmaxndbpros_; }
    unsigned int GetCurrentDbxPadding() const { return curdbxpad_; }

    size_t GetCurrentMaxDbPosPass2() const { return curmaxdbpospass2_; }
    size_t GetCurrentMaxNDbProsPass2() const { return curmaxndbprospass2_; }
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
        std::unique_ptr<PMBatchProData> bdbC,
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        char** bdbCpmbeg,
        char** bdbCpmend
    );

    void CacheSSENNWeightsDevice();
    void CacheSSSScoresDevice( const SSSScores& );
    void CacheCVS2ScoresDevice( const CVS2Scores& );
    void CacheHDPScoresDevice( const HDPscores& );

    void CacheDataDevice(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend
    );
    void CopyPMDataToDevice(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        char** bdbCpmbeg,
        char** bdbCpmend,
        char*& dev_pckdpm,
        size_t szmaxsize,
        size_t cmndx
    );
    void PackPMDataForDevice( 
        size_t szpm2dvfields[], char*& tmpbuf, size_t& szalloc,
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        char** bdbCpmbeg,
        char** bdbCpmend
    );
    void DestroyCBSMObject(){ if(cbsm_) delete cbsm_; }
    void DestroyCBDPObject(){ if(cbdp_) delete cbdp_; }
    void DestroyCBSSObject(){ if(cbss_) delete cbss_; }
    void DestroyCBMAPDPObject(){ if(cbmapdp_) delete cbmapdp_; }
    void DestroyTextureObject( cudaTextureObject_t& texObj );
    void FreeDevicePtr( char*& d_ptr );

    void SetCurrentMaxDbPos( size_t value ) { curmaxdbpos_ = value; }
    void SetCurrentMaxNDbPros( size_t value ) { curmaxndbpros_ = value; }
    void SetCurrentDbxPadding( unsigned int value ) { curdbxpad_ = value; }

    void SetCurrentMaxDbPosPass2( size_t value ) { curmaxdbpospass2_ = value; }
    void SetCurrentMaxNDbProsPass2( size_t value ) { curmaxndbprospass2_ = value; }
    void SetCurrentDbxPaddingPass2( unsigned int value ) { dbxpadphase2_ = value; }
    void SetCurrentDbAlnLengthWithPaddingPass2( unsigned int value ) { dbalnlenphase2_ = value; }
    unsigned int GetCurrentDbAlnLengthWithPaddingPass2() const { return dbalnlenphase2_; }

    void TESTPrintProProScores1(
        char** querypmbeg, char** querypmend,
        char** bdb1pmbeg, char** bdb1pmend, char** bdbCpmbeg, char** bdbCpmend );

    size_t GetMemAlignment() const {
        size_t cszalnment = 256UL;
        cszalnment = SLC_MAX( cszalnment, deviceProp_.textureAlignment_ );
        return cszalnment;
    }

    size_t GetSizeOfDPDiagScores( size_t maxdbposs ) const;
    size_t GetSizeOfDPBottomScores( size_t maxdbposs ) const;
    size_t GetSizeOfDPMaxCoords( size_t maxdbposs ) const;
    size_t GetSizeOfDPBackTckData( size_t queryposs, size_t maxdbposs ) const;

    size_t GetSizeOfDP2DiagScores( size_t dbposs ) const;
    size_t GetSizeOfDP2BottomScores( size_t dbposs ) const;
//size_t GetSizeOfDP2DiagBlkProbScales( size_t queryposs, size_t dbposs, size_t dbpros ) const;
//size_t GetSizeOfDP2DP2PrbScales( size_t queryposs, size_t dbposs, size_t dbpros ) const;
    size_t GetSizeOfDP2MaxCoords( size_t dbposs ) const;
    size_t GetSizeOfDP2BackTckData( size_t queryposs, size_t dbposs ) const;
    size_t GetSizeOfDP2FwdMtx( size_t queryposs, size_t dbposs ) const;
//size_t GetSizeOfDP2BwdMtx( size_t queryposs, size_t dbposs ) const;
    size_t GetSizeOfSS2Data( size_t dbposs, size_t dbpros ) const;
    size_t GetSizeOfDP2AlnData( size_t dbposs, size_t dbpros ) const;
    size_t GetSizeOfDP2Alns( size_t queryposs, size_t dbposs, size_t dbpros, bool sssinuse = true ) const;

    size_t GetMaxAllowedNumDbProfilePositions( size_t queryposs, bool sssinuse = true ) const;

    void GetTotalMemoryReqs( 
        size_t nqyposs, size_t maxdbposs, size_t maxamountfordbposs,
        //{{db profile positions
        size_t* maxsizedbposs,
        //}}{{scores
        size_t* szsmatrix,
        //}}{{phase1
        size_t* szdpdiag, size_t* szdpbottom, size_t* szdpmaxcoords, size_t* szdpbtckdat,
        //}}{{phase2
        size_t* maxdbposspass2, size_t* maxndbprospass2, 
        size_t* szdp2diag, size_t* szdp2bottom, 
//**/size_t* szdp2diagblkprobscales, size_t* szdp2probscales,
        size_t* szdp2maxcoords, size_t* szdp2btckdat,
        size_t* szdp2fwdmtx,// size_t* szdp2bwdmtx, 
        size_t* szss2data, size_t* szdp2alndata, size_t* szdp2alns,
        //}}{{overall size
        size_t* szovlpos
        //}}
    ) const;

    void MsgAddressTable( mystring preamb, int level ) const;

    static size_t GetAbsDiff( size_t val1, size_t val2 ) { 
        return (val1 < val2)? val2 - val1: val1 - val2; 
    }


    //{{ ===== results section =====
    void TransferResultsDevice(
        int qrysernr,
        size_t nqyposs, float qyeno,
        std::unique_ptr<PMBatchProData> bdbC,
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend,
        char** bdbCpmbeg,
        char** bdbCpmend,
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
    CuBatchScoreMatrix* cbsm_;//batch score matrix
    CuBatchDP* cbdp_;//batch dynamic programming object
    CuBatchSS* cbss_;//object for the batch calculation of statistical significance
    CuBatchMAPDP* cbmapdp_;//batch MAP dynamic programming object
    CuBatchProcessingFinalizer* cbpfin_;//results finalizer
    bool hdp1scoresinuse_;//whether HDP1 scores are being used
    bool mapalninuse_;//whether MAP alignment calculation is being used
    float scorethld_;//score threshold for pass-2 alignments
    float logevthld_;//log E-value threshold
    //{{db details
    size_t querylength_;
    uint64_mt dblength_;
    size_t ndbsequences_;
    //}}
    //{{attributes required for a device to process serialized data
    SerializedScoresAttr ssssattr_;//attributes of sereliazed sss scores
    SerializedCVS2ScoresAttr cvs2sattr_;//attributes of sereliazed cvs2s scores
    SerializedScoresAttr hdpsattr_;//attributes of sereliazed HDP scores
    //}}
    size_t deviceallocsize_;//allocated size for a device (max limit of memory used in a device)
    size_t curmaxdbpos_;//current maximum number of db profile positions for calculations
    size_t curmaxndbpros_;//current maximum number of db profiles to be processed
    unsigned int curdbxpad_;//current padding in positions along the x (db) profile positions
    size_t curmaxdbpospass2_;//maximum number of db profile positions for calculations in phase 2
    size_t curmaxndbprospass2_;//maximum number of db profiles in phase 2
    unsigned int dbxpadphase2_;//padding in positions along the x (db) axis in phase 2
    unsigned int dbalnlenphase2_;//total alignment length with padding in phase 2
    //{{host allocated pointers
    char** h_querypmbeg_;//beginning address of queries
    char** h_querypmend_;//terminal address of queries
    char** h_bdb1pmbeg_;//beginning address of cached database profiles
    char** h_bdb1pmend_;//terminal address of cached database profiles
    char** h_bdbCpmbeg_;//beginning address of new profiles read from the database 
    char** h_bdbCpmend_;//terminal address of new profiles read from the database 
    //}}
    //{{host allocated pointers for results
    cudaStream_t streamcopyres_;//stream for copying results
    char* h_results_;//results received from device
    bool lockedresmem_;//page-locked memory for results; also specifies that a stream is initialized
    size_t sz_mem_results_;//size of allocated memory for results
    size_t limit_beg_results_;//the boundary where the results sections begin
    //}}
    //{{device data
    DeviceProperties deviceProp_;//the properties of device deviceno_
    cudaTextureObject_t hdp1sTexObj_;//texture object for the HDP scores
    char* d_heap_;//global heap containing all data written, generated, and read
    //device pointers to the beginnings of sections in the heap (allocated block of memory):
    size_t sz_heapsections_[nDevDataSectionsDiv2];
    //}}
};

// -------------------------------------------------------------------------
// INLINES ...
//
// -------------------------------------------------------------------------
// MsgAddressTable: print as message the address table of sections 
inline
void CuBatchProcessing::MsgAddressTable( const mystring preamb, const int level ) const
{
    MYMSGBEGl(level)
        char msgbuf[BUF_MAX];
        sprintf(msgbuf, "Address table: %p", d_heap_ );
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfPadding", sz_heapsections_[ddsEndOfPadding]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfSSSTable", sz_heapsections_[ddsEndOfSSSTable]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfCVS2SMap", sz_heapsections_[ddsEndOfCVS2SMap]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfHDPscores [texture]", sz_heapsections_[ddsEndOfHDPscores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfCached", sz_heapsections_[ddsEndOfCached]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDbChunk", sz_heapsections_[ddsEndOfDbChunk]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfOrgScores", sz_heapsections_[ddsEndOfOrgScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfModScores", sz_heapsections_[ddsEndOfModScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        MYMSGnonl((preamb + "------- phase 1 division").c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPDiagScores", sz_heapsections_[ddsEndOfDPDiagScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPBottomScores", sz_heapsections_[ddsEndOfDPBottomScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPMaxCoords", sz_heapsections_[ddsEndOfDPMaxCoords]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPBackTckData", sz_heapsections_[ddsEndOfDPBackTckData]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "------- phase 2 division; %.2f of positions", CLOptions::GetDEV_PASS2MEMP());
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2DiagScores", sz_heapsections_[ddsEndOfDP2DiagScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2BottomScores", sz_heapsections_[ddsEndOfDP2BottomScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
//**/sprintf(msgbuf, "        +%zu ddsEndOfDP2DiagBlkProbScales", sz_heapsections_[ddsEndOfDP2DiagBlkProbScales]);
//MYMSGnonl((preamb + msgbuf).c_str(),3);
//**/sprintf(msgbuf, "        +%zu ddsEndOfDP2PrbScales", sz_heapsections_[ddsEndOfDP2PrbScales]);
//MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2MaxCoords", sz_heapsections_[ddsEndOfDP2MaxCoords]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2BackTckData", sz_heapsections_[ddsEndOfDP2BackTckData]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2FwdMtx", sz_heapsections_[ddsEndOfDP2FwdMtx]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
//sprintf(msgbuf, "        +%zu ddsEndOfDP2BwdMtx", sz_heapsections_[ddsEndOfDP2BwdMtx]);
//MYMSG((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfSS2Data", sz_heapsections_[ddsEndOfSS2Data]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2AlnData", sz_heapsections_[ddsEndOfDP2AlnData]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2Alns", sz_heapsections_[ddsEndOfDP2Alns]);
        MYMSG((preamb + msgbuf).c_str(),3);
    MYMSGENDl

    if( deviceallocsize_ < sz_heapsections_[nDevDataSectionsDiv1-1] ||
        deviceallocsize_ < sz_heapsections_[nDevDataSectionsDiv2-1] )
        throw MYRUNTIME_ERROR( preamb + "Out of range of allocated memory.");
}

// -------------------------------------------------------------------------
// GetSizeOfDPDiagScores: get size of the buffers of diagonal scores 
// used for dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDPDiagScores( size_t maxdbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDPDiagScores", 8 );
    const size_t cszalnment = GetMemAlignment();
    //NOTE:two diagonal buffers (each comprising a number of DP states) are 
    // required for the execution of DP and one for tracking max score while 
    // performing DP;
    //NOTE:the worst case scenario is all profiles of unit length, 
    // which in practice will not almost appear;
    //NOTE:add diagonal width to both ends of the db positions;
    //NOTE:use this size of memory for buffers (1st line), as a series of 
    // block diagonals are launched in the device, which in theory may require 
    // this much memory;
    size_t szdpdiag = 
            (maxdbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            ( dpdssDiagM * nTDPDiagScoreSubsections + 1/*dpdssDiagM*/
#ifdef CUDP_CALC_CORRSCORES_INLINE
            + 
              1/*dpdssCorrDiagM*/ + 
              CUDP_CORR_NSCORES/*dpdssCorrDiag1*/ + CUDP_CORR_NSCORES/*dpdssCorrDiag2*/ +
              1/*dpdssCorrDiagS1*/ + 1/*dpdssCorrDiagS2*/ 
#endif
            ) * sizeof(CUBDP_TYPE);
    szdpdiag = ALIGN_UP( szdpdiag, cszalnment );
    return szdpdiag;
}

// -------------------------------------------------------------------------
// GetSizeOfDPBottomScores: get size of the buffer of bottom scores 
// used for dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDPBottomScores( size_t maxdbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDPBottomScores", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbottom = 
            (maxdbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            ( nTDPDiagScoreSubsections/*dpbssBottm*/ 
#ifdef CUDP_CALC_CORRSCORES_INLINE
            +
              CUDP_CORR_NSCORES/*dpbssCorrBottm*/ + 1/*dpbssCorrBottmS*/
#endif
            ) * sizeof(CUBDP_TYPE);
    szdpbottom = ALIGN_UP( szdpbottom, cszalnment );
    return szdpbottom;
}

// -------------------------------------------------------------------------
// GetSizeOfDPDPMaxCoords: get size of the buffer of the coordinates of 
// maximum alignment scores calculated by dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDPMaxCoords( size_t maxdbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDPMaxCoords", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpdiagcoords = 
            (maxdbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            sizeof(uint);
    szdpdiagcoords = ALIGN_UP( szdpdiagcoords, cszalnment );
    return szdpdiagcoords;
}

// -------------------------------------------------------------------------
// GetSizeOfDPBackTckData: get size of the buffer of backtracking 
// information obtained from dynamic programming
inline
size_t CuBatchProcessing::GetSizeOfDPBackTckData( size_t queryposs, size_t maxdbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDPBackTckData", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbtckdat = 
            (maxdbposs * queryposs) * 
            sizeof(char);
    szdpbtckdat = ALIGN_UP( szdpbtckdat, cszalnment );
    return szdpbtckdat;
}


// phase 2 sizes -----------------------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// GetSizeOfDP2DiagScores: get size of the buffers of diagonal scores 
// used for phase2 dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDP2DiagScores( size_t dbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2DiagScores", 8 );
    const size_t cszalnment = GetMemAlignment();
    // see GetSizeOfDPDiagScores for comments
    size_t szdpdiag = 
            (dbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            ( dpdssDiagM * nTDPDiagScoreSubsections + 1/*dpdssDiagM*/) * 
            sizeof(float);
    szdpdiag = ALIGN_UP( szdpdiag, cszalnment );
    return szdpdiag;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2BottomScores: get size of the buffer of bottom scores 
// used for phase2 dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDP2BottomScores( size_t dbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2BottomScores", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbottom = 
            (dbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            ( nTDPDiagScoreSubsections/*dpbssBottm*/) * 
            sizeof(float);
    szdpbottom = ALIGN_UP( szdpbottom, cszalnment );
    return szdpbottom;
}


// // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// // GetSizeOfDP2DiagBlkProbScales: get size of the buffer of diagonal 
// // block-specific scale factors for forward probabilities
// inline
// size_t CuBatchProcessing::GetSizeOfDP2DiagBlkProbScales( 
//     size_t queryposs, size_t dbposs, size_t dbpros ) const
// {
//     MYMSG( "CuBatchProcessing::GetSizeOfDP2DiagBlkProbScales", 8 );
//     const size_t cszalnment = GetMemAlignment();
//     //approx. memory requirements for one profile: 
//     // l_query * (l_target/CUMAPDP_2DCACHE_DIM_D)
//     size_t szdp2blkprobscales = 
//             ((dbposs + TIMES2(CUMAPDP_2DCACHE_DIM_D) * dbpros + CUMAPDP_2DCACHE_DIM_D-1) / 
//             CUMAPDP_2DCACHE_DIM_D * queryposs) * 
//             sizeof(float);
//     szdp2blkprobscales = ALIGN_UP( szdp2blkprobscales, cszalnment );
//     return szdp2blkprobscales;
// }
// // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// // GetSizeOfDP2DP2PrbScales: get size of the buffer of final (rescaled) 
// // scale factors for forward/backward probabilities
// inline
// size_t CuBatchProcessing::GetSizeOfDP2DP2PrbScales( 
//     size_t queryposs, size_t /*dbposs*/, size_t dbpros ) const
// {
//     MYMSG( "CuBatchProcessing::GetSizeOfDP2DP2PrbScales", 8 );
//     const size_t cszalnment = GetMemAlignment();
//     //two values (nTMAPDPProbScales) are saved for each query position per db profile
//     size_t szdp2probscales = 
//             (dbpros * nTMAPDPProbScales * queryposs) * 
//             sizeof(float);
//     szdp2probscales = ALIGN_UP( szdp2probscales, cszalnment );
//     return szdp2probscales;
// }


// -------------------------------------------------------------------------
// GetSizeOfDP2MaxCoords: get size of the buffer of the coordinates of 
// maximum alignment scores calculated by phase2 dynamic programming
//
inline
size_t CuBatchProcessing::GetSizeOfDP2MaxCoords( size_t dbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2MaxCoords", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpdiagcoords = 
            (dbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            sizeof(uint);
    szdpdiagcoords = ALIGN_UP( szdpdiagcoords, cszalnment );
    return szdpdiagcoords;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2BackTckData: get size of the buffer of backtracking 
// information obtained from phase2 dynamic programming
inline
size_t CuBatchProcessing::GetSizeOfDP2BackTckData( size_t queryposs, size_t dbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2BackTckData", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbtckdat = 
            (dbposs * queryposs) * 
            sizeof(char);
    szdpbtckdat = ALIGN_UP( szdpbtckdat, cszalnment );
    return szdpbtckdat;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2FwdMtx: get size of the matrix of forward probabilities 
// obtained from phase2 dynamic programming
inline
size_t CuBatchProcessing::GetSizeOfDP2FwdMtx( size_t queryposs, size_t dbposs ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2FwdMtx", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbtckdat = 
            (dbposs * queryposs) * 
            sizeof(float);
    szdpbtckdat = ALIGN_UP( szdpbtckdat, cszalnment );
    return szdpbtckdat;
}

// // -------------------------------------------------------------------------
// // GetSizeOfDP2BwdMtx: get size of the matrix of backward probabilities 
// // obtained from phase2 dynamic programming
// inline
// size_t CuBatchProcessing::GetSizeOfDP2BwdMtx( size_t queryposs, size_t dbposs) const
// {
//     MYMSG( "CuBatchProcessing::GetSizeOfDP2BwdMtx", 8 );
//     const size_t cszalnment = GetMemAlignment();
//     size_t szdpbtckdat = 
//             (dbposs * queryposs) * 
//             sizeof(float);
//     szdpbtckdat = ALIGN_UP( szdpbtckdat, cszalnment );
//     return szdpbtckdat;
// }

// -------------------------------------------------------------------------
// GetSizeOfSS2Data: get size of the data associated with statistical 
// significance calculations in phase 2
inline
size_t CuBatchProcessing::GetSizeOfSS2Data( size_t /*dbposs*/, size_t dbpros ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfSS2Data", 8 );
    const size_t cszalnment = GetMemAlignment();
    //const int nscores = CUSS_N_DIFF_SCORES;
    //const int nscaledscores = CUSS_N_DIFF_SCALED_SCORES;
    const int alignedntotscores = CUSS_ALIGNED_N_DIFF_TOTAL_SCORES;
    size_t szss2dat = 
            //(dbpros * (nscores + nscaledscores)) * 
            (dbpros * alignedntotscores) * 
            sizeof(float);
    szss2dat = ALIGN_UP( szss2dat, cszalnment );
    return szss2dat;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2AlnData: get size of the data where all relative alignment 
// information of phase 2 calculations is to be placed
inline
size_t CuBatchProcessing::GetSizeOfDP2AlnData( size_t /*dbposs*/, size_t dbpros ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2AlnData", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdp2alndat = 
            (dbpros * nTDP2OutputAlnData) * 
            sizeof(float);
    szdp2alndat = ALIGN_UP( szdp2alndat, cszalnment );
    return szdp2alndat;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2Alns: get size of the buffer dedicated to alignments 
// themselves 
inline
size_t CuBatchProcessing::GetSizeOfDP2Alns( 
    size_t queryposs, size_t dbposs, size_t dbpros, bool sssinuse ) const
{
    MYMSG( "CuBatchProcessing::GetSizeOfDP2Alns", 8 );
    const size_t cszalnment = GetMemAlignment();
    //maximum alignment length: l_query + l_target
    size_t szdp2alns = 
            ((dbposs + dbpros * (queryposs+1)) *
             (sssinuse? nTDP2OutputAlignmentSSS: nTDP2OutputAlignment)) *
            sizeof(char);
    szdp2alns = ALIGN_UP( szdp2alns, cszalnment );
    return szdp2alns;
}



// =========================================================================
// GetMaxAllowedNumDbProfilePositions: get the maximum allowed number of Db 
// profile positions given query length to avoid aoverflow
inline
size_t CuBatchProcessing::GetMaxAllowedNumDbProfilePositions( 
    size_t queryposs, bool sssinuse ) const
{
    MYMSG( "CuBatchProcessing::GetMaxAllowedNumDbProfilePositions", 7 );
    size_t deplen1 = //dependent maximum length
        (size_t)INT_MAX / queryposs;
    size_t deplen2 = (size_t)(//maximum length dependent on the number of alignments
        (float)INT_MAX / (sssinuse? nTDP2OutputAlignmentSSS: nTDP2OutputAlignment) /
        (1.0f + (float)(queryposs+1+1/*compensate rounding*/) /
            (float)CLOptions::GetDEV_EXPCT_DBPROLEN())
    );
    size_t indeplen1 = //independent maximum length
        (size_t)INT_MAX / ( dpdssDiagM * nTDPDiagScoreSubsections + 1/*dpdssDiagM*/
#ifdef CUDP_CALC_CORRSCORES_INLINE
        + 
            1/*dpdssCorrDiagM*/ + 
            CUDP_CORR_NSCORES/*dpdssCorrDiag1*/ + CUDP_CORR_NSCORES/*dpdssCorrDiag2*/ +
            1/*dpdssCorrDiagS1*/ + 1/*dpdssCorrDiagS2*/ 
#endif
        ) - TIMES2(CUDP_2DCACHE_DIM_D);
    //maximum length independent on the query length but 
    // dependent on the number of db profiles:
    size_t indeplen2 = 
        (size_t)INT_MAX / CUSS_ALIGNED_N_DIFF_TOTAL_SCORES *
            CLOptions::GetDEV_EXPCT_DBPROLEN();

    //take minimum:
    size_t maxallwdposs = SLC_MIN(deplen1, deplen2);
    maxallwdposs = SLC_MIN(maxallwdposs, indeplen1);
    maxallwdposs = SLC_MIN(maxallwdposs, indeplen2);
    return maxallwdposs;
}



// overall memory requirements ---------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// GetTotalMemoryReqs: get total memeory requirements given the number of 
// query and db profile positions
// NOTE: addresses are assumed to be valid and they are not verified
inline
void CuBatchProcessing::GetTotalMemoryReqs( 
    size_t nqyposs, size_t maxdbposs, size_t maxamountfordbposs,
    //{{db profile positions
    size_t* maxsizedbposs,
    //}}{{scores
    size_t* szsmatrix,
    //}}{{phase1
    size_t* szdpdiag, size_t* szdpbottom, size_t* szdpmaxcoords, size_t* szdpbtckdat,
    //}}{{phase2
    size_t* maxdbposspass2, size_t* maxndbprospass2, 
    size_t* szdp2diag, size_t* szdp2bottom, 
//**/size_t* szdp2diagblkprobscales, size_t* szdp2probscales,
    size_t* szdp2maxcoords, size_t* szdp2btckdat,
    size_t* szdp2fwdmtx,// size_t* szdp2bwdmtx, 
    size_t* szss2data, size_t* szdp2alndata, size_t* szdp2alns,
    //}}{{overall size
    size_t* szovlpos
    //}}
    ) const
{
    MYMSG( "CuBatchProcessing::GetTotalMemoryReqs", 7 );
    const size_t cszalnment = GetMemAlignment();
    const float pass2memperc = CLOptions::GetDEV_PASS2MEMP();
    size_t sum = 0UL;//temporary sum variable 
    //{{memory for db profile data itself
    *maxsizedbposs = PMBatchProData::GetPMDataSizeUBTotal( maxdbposs );
    *maxsizedbposs = SLC_MIN( *maxsizedbposs, maxamountfordbposs );
    *maxsizedbposs = ALIGN_UP( *maxsizedbposs, cszalnment );
    //}}{{memory for scores
    *szsmatrix = nqyposs * maxdbposs * sizeof(CUBSM_TYPE);
    //}}{{phase-1 memory requirements
    *szdpdiag = GetSizeOfDPDiagScores( maxdbposs );
    *szdpbottom = GetSizeOfDPBottomScores( maxdbposs );
    *szdpmaxcoords = GetSizeOfDPMaxCoords( maxdbposs );
    *szdpbtckdat = GetSizeOfDPBackTckData( nqyposs, maxdbposs );
    //}}{{phase-2 memory requirements
    *maxdbposspass2 = (size_t)(maxdbposs * pass2memperc);
    *maxdbposspass2 = ALIGN_DOWN( *maxdbposspass2, CUL2CLINESIZE );
    *maxndbprospass2 = *maxdbposspass2 / CLOptions::GetDEV_EXPCT_DBPROLEN();
    *szdp2diag = GetSizeOfDP2DiagScores( *maxdbposspass2 );
    *szdp2bottom = GetSizeOfDP2BottomScores( *maxdbposspass2 );
//**/*szdp2diagblkprobscales = GetSizeOfDP2DiagBlkProbScales( nqyposs, *maxdbposspass2, *maxndbprospass2 );
//**/*szdp2probscales = GetSizeOfDP2DP2PrbScales( nqyposs, *maxdbposspass2, *maxndbprospass2 );
    *szdp2maxcoords = GetSizeOfDP2MaxCoords( *maxdbposspass2 );
    *szdp2btckdat = GetSizeOfDP2BackTckData( nqyposs, *maxdbposspass2 );
    *szdp2fwdmtx = GetSizeOfDP2FwdMtx( nqyposs, *maxdbposspass2 );
//*szdp2bwdmtx = GetSizeOfDP2BwdMtx( nqyposs, *maxdbposspass2 );
    *szss2data = GetSizeOfSS2Data( *maxdbposspass2, *maxndbprospass2 );
    *szdp2alndata = GetSizeOfDP2AlnData( *maxdbposspass2, *maxndbprospass2 );
    *szdp2alns = GetSizeOfDP2Alns( nqyposs, *maxdbposspass2, *maxndbprospass2 );
    sum = *szdp2diag + *szdp2bottom + 
//*szdp2diagblkprobscales + *szdp2probscales + 
          *szdp2maxcoords + *szdp2btckdat + 
          *szdp2fwdmtx + 
//*szdp2bwdmtx + 
          *szss2data + *szdp2alndata + *szdp2alns;
    //}}
    *szovlpos = *szdpdiag + *szdpbottom  + *szdpmaxcoords + *szdpbtckdat;
    if( *szovlpos < sum )
        *szovlpos = sum;
    *szovlpos += *maxsizedbposs + *szsmatrix;
    if( GetModScoreMatrixInUse())
        *szovlpos += *szsmatrix;
}

#endif//__CuBatchProcessing_h__
