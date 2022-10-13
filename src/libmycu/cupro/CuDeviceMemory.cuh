/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuDeviceMemory_h__
#define __CuDeviceMemory_h__

#include "liblib/mybase.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <memory>
#include <mutex>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libHDP/HDPscores.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "libmycu/cumapdp/CuBatchMAPDP_com.h"
#include "CuBatchProcessing_com.cuh"

////////////////////////////////////////////////////////////////////////////
// CLASS CuDeviceMemory
// Device memory arrangement for profile-profile computation
//
class CuDeviceMemory
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
        //
        ddsBegOfOrgScores,//beginning of the section of originally calculated scores
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
        ddsEndOfDP2MaxCoords,//end of phase 2 section of coordinates of alignment scores
        ddsEndOfDP2BackTckData,//end of phase 2 section of backtracking data
        ddsEndOfDP2FwdMtx,//end of phase 2 section of forward probability matrix
//ddsEndOfDP2BwdMtx,//end of phase 2 section of backward probability matrix
        ddsEndOfSS2Data,//end of phase 2 data for calculation of statistical significance [profile-specific]
        ddsEndOfDP2AlnData,//end of phase 2 section of alignment data [profile-specific, output for host]
        ddsEndOfDP2Alns,//end of phase 2 section of drawn alignments [profile-specific, output for host]
        //
        // --- variables ---
        ddsBegOfGlobVars = ddsEndOfDP2Alns,
        ddsEndOfGlobVars,//end of the section of global variables
        //
        nDevDataSectionsDiv2
    };
    enum TDevMatrixSections {
        dmsOrgScores,//originally calculated scores
        dmsModScores,//modular scores w/o certain terms
        nDevMatrixSections
    };
    //attributes of profiles passed through the DP phase
    enum TDevGlobVariables {
        dgvNPassedPros,//number of profiles passed to a further stage of processing
        dgvNPosits,//total number of positions of passed profiles
        dgvMaxProLen,//maximum profile length over passed profiles
        dgvMutex,//mutex for registering profiles passing to the next stage
        nDevGlobVariables
    };


    CuDeviceMemory(
        DeviceProperties dprop, 
        bool hdp1scoresinuse, 
        bool mapalninuse,
        bool Xuninf,
        size_t deviceallocsize,
        int nareas );

    virtual ~CuDeviceMemory();

    int GetNAreas() const {return nareas_;}
    const DeviceProperties& GetDeviceProp() const { return deviceProp_;}
    char* GetHeap() const {return d_heap_;}
    char** GetHstQueryPMBeg() const {return h_querypmbeg_;}
    char** GetHstQueryPMEnd() const {return h_querypmend_;}
    char** GetHstBdb1PMBeg() const {return h_bdb1pmbeg_;}
    char** GetHstBdb1PMEnd() const {return h_bdb1pmend_;}

    SerializedScoresAttr& GetSSSSAttr() { return ssssattr_;}
    SerializedCVS2ScoresAttr& GetCVS2SAttr() {return cvs2sattr_;}
    SerializedScoresAttr& GetHDPSAttr() {return hdpsattr_;}
    cudaTextureObject_t& GetHDP1STexObj() {return hdp1sTexObj_;}

    void CacheCompleteData(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend
    );

    void TransferCPMData(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields
    );

    size_t CalcMaxDbDataChunkSize( size_t nqyposs );

    bool GetHDP1ScoresInUse() const { return hdp1scoresinuse_; }
    bool GetMAPAlnInUse() const { return mapalninuse_; }
    bool GetXUninformative() const { return Xuninf_; }
    bool GetModScoreMatrixInUse() const { 
        return MOptions::GetSSSWGT() > 0.0f || MOptions::GetCVSWGT() > 0.0f;
    }

    size_t GetCurrentMaxDbPos() const { return curmaxdbpos_; }
    size_t GetCurrentMaxNDbPros() const { return curmaxndbpros_; }

    size_t GetCurrentMaxDbPosPass2() const { return curmaxdbpospass2_; }
    size_t GetCurrentMaxNDbProsPass2() const { return curmaxndbprospass2_; }

    size_t GetOffsetOfHeapSection(int ano, int s) const;

    //{{addressses determined at runtime:
    size_t GetEndOfSS2Data(int ano) const {
        return GetMAPAlnInUse()?
            sz_heapsections_[ano][ddsEndOfDP2FwdMtx]
        :   sz_heapsections_[ano][ddsEndOfDPBackTckData];
    }
    //}}

protected:
    explicit CuDeviceMemory();

    size_t CalcMaxDbDataChunkSizeHelper( size_t nqyposs, size_t leftoversize );

    void CacheSSENNWeights();
    void CacheSSSScores( const SSSScores& );
    void CacheCVS2Scores( const CVS2Scores& );
    void CacheHDPScores( const HDPscores& );
    void CacheData(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend
    );

    void CacheSSENNWeightsDevice();
    void CacheSSSScoresDevice( const SSSScores& );
    void CacheCVS2ScoresDevice( const CVS2Scores& );
    void CacheHDPScoresDevice( const HDPscores& );

    void TransferCPMDataDevice(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields
    );
    void CacheDataDevice(
        char** querypmbeg,
        char** querypmend,
        char** bdb1pmbeg,
        char** bdb1pmend
    );
    void CopyCPMDataToDevice(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields,
        char* dev_pckdpm,
        size_t szmaxsize,
        size_t cmndx 
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
    void DestroyTextureObject( cudaTextureObject_t& texObj );
    void FreeDevicePtr( char*& d_ptr );

    void SetCurrentMaxDbPos( size_t value ) { curmaxdbpos_ = value; }
    void SetCurrentMaxNDbPros( size_t value ) { curmaxndbpros_ = value; }

    void SetCurrentMaxDbPosPass2( size_t value ) { curmaxdbpospass2_ = value; }
    void SetCurrentMaxNDbProsPass2( size_t value ) { curmaxndbprospass2_ = value; }

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
    size_t GetSizeOfDP2MaxCoords( size_t dbposs ) const;
    size_t GetSizeOfDP2BackTckData( size_t queryposs, size_t dbposs ) const;
    size_t GetSizeOfDP2FwdMtx( size_t queryposs, size_t dbposs ) const;
//size_t GetSizeOfDP2BwdMtx( size_t queryposs, size_t dbposs ) const;
    size_t GetSizeOfSS2Data( size_t dbposs, size_t dbpros ) const;
    size_t GetSizeOfDP2AlnData( size_t dbposs, size_t dbpros ) const;
    size_t GetSizeOfDP2Alns( size_t queryposs, size_t dbposs, size_t dbpros, bool sssinuse = true ) const;
    size_t GetSizeOfGlobVariables() const;


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
        //}}{{variables
        size_t* szglbvars,
        //}}{{overall size
        size_t* szovlpos
        //}}
    ) const;

    void MsgAddressTable( int areano, mystring preamb, int level ) const;

    static size_t GetAbsDiff( size_t val1, size_t val2 ) { 
        return (val1 < val2)? val2 - val1: val1 - val2; 
    }

private:
    bool hdp1scoresinuse_;//whether HDP1 scores are being used
    bool mapalninuse_;//whether MAP alignment calculation is being used
    bool Xuninf_;//whether positions masked with Xs are to be treated as uninformative
    size_t querylength_;//length to adjust memory configuration to
    //{{attributes required for a device to process serialized data
    SerializedScoresAttr ssssattr_;//attributes of sereliazed sss scores
    SerializedCVS2ScoresAttr cvs2sattr_;//attributes of sereliazed cvs2s scores
    SerializedScoresAttr hdpsattr_;//attributes of sereliazed HDP scores
    //}}
    size_t deviceallocsize_;//allocated size for a device (max limit of memory used in a device)
    int nareas_;//number of distinct computation areas in device memory
    size_t curmaxdbpos_;//current maximum number of db profile positions for calculations
    size_t curmaxndbpros_;//current maximum number of db profiles to be processed
    size_t curmaxdbpospass2_;//maximum number of db profile positions for calculations in phase 2
    size_t curmaxndbprospass2_;//maximum number of db profiles in phase 2
    //{{host allocated pointers
    char** h_querypmbeg_;//beginning address of queries
    char** h_querypmend_;//terminal address of queries
    char** h_bdb1pmbeg_;//beginning address of cached database profiles
    char** h_bdb1pmend_;//terminal address of cached database profiles
    //}}
    //{{device data
    DeviceProperties deviceProp_;//the properties of device deviceno_
    cudaTextureObject_t hdp1sTexObj_;//texture object for the HDP scores
    char* d_heap_;//global heap containing all data written, generated, and read
    //device pointers to the beginnings of sections in the heap (allocated block of memory)
    // for each area
    size_t (*sz_heapsections_)[nDevDataSectionsDiv2];
    //}}
};

// -------------------------------------------------------------------------
// INLINES ...
//
inline
size_t CuDeviceMemory::GetOffsetOfHeapSection(int ano, int s) const
{
#ifdef __DEBUG__
    if( ano < 0 || ano >= nareas_ || s < 0 || s >= nDevDataSectionsDiv2)
        throw MYRUNTIME_ERROR("CuDeviceMemory::GetSizeOfHeapSection: Memory access error.");
#endif
    return sz_heapsections_[ano][s];
}
// -------------------------------------------------------------------------
// MsgAddressTable: print as message the address table of sections 
// ano, the number of area 
inline
void CuDeviceMemory::MsgAddressTable( int ano, const mystring preamb, const int level ) const
{
    MYMSGBEGl(level)
        char msgbuf[BUF_MAX];
        sprintf(msgbuf, "Address table: %p [area no. %d]", d_heap_, ano );
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfPadding", sz_heapsections_[ano][ddsEndOfPadding]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfSSSTable", sz_heapsections_[ano][ddsEndOfSSSTable]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfCVS2SMap", sz_heapsections_[ano][ddsEndOfCVS2SMap]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfHDPscores [texture]", sz_heapsections_[ano][ddsEndOfHDPscores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfCached", sz_heapsections_[ano][ddsEndOfCached]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDbChunk", sz_heapsections_[ano][ddsEndOfDbChunk]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsBegOfOrgScores", sz_heapsections_[ano][ddsBegOfOrgScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfOrgScores", sz_heapsections_[ano][ddsEndOfOrgScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfModScores", sz_heapsections_[ano][ddsEndOfModScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        MYMSGnonl((preamb + "------- phase 1 division").c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPDiagScores", sz_heapsections_[ano][ddsEndOfDPDiagScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPBottomScores", sz_heapsections_[ano][ddsEndOfDPBottomScores]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPMaxCoords", sz_heapsections_[ano][ddsEndOfDPMaxCoords]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDPBackTckData", sz_heapsections_[ano][ddsEndOfDPBackTckData]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        if( GetMAPAlnInUse()) {
            sprintf(msgbuf, "------- phase 2 division; %.2f of positions", CLOptions::GetDEV_PASS2MEMP());
            MYMSGnonl((preamb + msgbuf).c_str(),3);
            sprintf(msgbuf, "        +%zu ddsEndOfDP2DiagScores", sz_heapsections_[ano][ddsEndOfDP2DiagScores]);
            MYMSGnonl((preamb + msgbuf).c_str(),3);
            sprintf(msgbuf, "        +%zu ddsEndOfDP2BottomScores", sz_heapsections_[ano][ddsEndOfDP2BottomScores]);
            MYMSGnonl((preamb + msgbuf).c_str(),3);
            sprintf(msgbuf, "        +%zu ddsEndOfDP2MaxCoords", sz_heapsections_[ano][ddsEndOfDP2MaxCoords]);
            MYMSGnonl((preamb + msgbuf).c_str(),3);
            sprintf(msgbuf, "        +%zu ddsEndOfDP2BackTckData", sz_heapsections_[ano][ddsEndOfDP2BackTckData]);
            MYMSGnonl((preamb + msgbuf).c_str(),3);
            sprintf(msgbuf, "        +%zu ddsEndOfDP2FwdMtx", sz_heapsections_[ano][ddsEndOfDP2FwdMtx]);
            MYMSGnonl((preamb + msgbuf).c_str(),3);
//sprintf(msgbuf, "        +%zu ddsEndOfDP2BwdMtx", sz_heapsections_[ano][ddsEndOfDP2BwdMtx]);
//MYMSG((preamb + msgbuf).c_str(),3);
        }
        MYMSGnonl((preamb + "------- SS and alignment data division").c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfSS2Data", sz_heapsections_[ano][ddsEndOfSS2Data]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2AlnData", sz_heapsections_[ano][ddsEndOfDP2AlnData]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfDP2Alns", sz_heapsections_[ano][ddsEndOfDP2Alns]);
        MYMSGnonl((preamb + msgbuf).c_str(),3);
        sprintf(msgbuf, "        +%zu ddsEndOfGlobVars", sz_heapsections_[ano][ddsEndOfGlobVars]);
        MYMSG((preamb + msgbuf).c_str(),3);
    MYMSGENDl

    if( deviceallocsize_ < sz_heapsections_[ano][nDevDataSectionsDiv1-1] ||
        deviceallocsize_ < sz_heapsections_[ano][nDevDataSectionsDiv2-1] )
        throw MYRUNTIME_ERROR( preamb + "Out of range of allocated memory.");
}

// -------------------------------------------------------------------------
// GetSizeOfDPDiagScores: get size of the buffers of diagonal scores 
// used for dynamic programming
//
inline
size_t CuDeviceMemory::GetSizeOfDPDiagScores( size_t maxdbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDPDiagScores", 8 );
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
size_t CuDeviceMemory::GetSizeOfDPBottomScores( size_t maxdbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDPBottomScores", 8 );
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
size_t CuDeviceMemory::GetSizeOfDPMaxCoords( size_t maxdbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDPMaxCoords", 8 );
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
size_t CuDeviceMemory::GetSizeOfDPBackTckData( size_t queryposs, size_t maxdbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDPBackTckData", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbtckdat = 
            ((maxdbposs/**/ + TIMES2(CUDP_2DCACHE_DIM_D)/**/) * queryposs) * 
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
size_t CuDeviceMemory::GetSizeOfDP2DiagScores( size_t dbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2DiagScores", 8 );
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
size_t CuDeviceMemory::GetSizeOfDP2BottomScores( size_t dbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2BottomScores", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbottom = 
            (dbposs + TIMES2(CUDP_2DCACHE_DIM_D)) * 
            ( nTDPDiagScoreSubsections/*dpbssBottm*/) * 
            sizeof(float);
    szdpbottom = ALIGN_UP( szdpbottom, cszalnment );
    return szdpbottom;
}


// -------------------------------------------------------------------------
// GetSizeOfDP2MaxCoords: get size of the buffer of the coordinates of 
// maximum alignment scores calculated by phase2 dynamic programming
//
inline
size_t CuDeviceMemory::GetSizeOfDP2MaxCoords( size_t dbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2MaxCoords", 8 );
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
size_t CuDeviceMemory::GetSizeOfDP2BackTckData( size_t queryposs, size_t dbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2BackTckData", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szdpbtckdat = 
            ((dbposs/**/ + TIMES2(CUDP_2DCACHE_DIM_D)/**/) * queryposs) * 
            sizeof(char);
    szdpbtckdat = ALIGN_UP( szdpbtckdat, cszalnment );
    return szdpbtckdat;
}

// -------------------------------------------------------------------------
// GetSizeOfDP2FwdMtx: get size of the matrix of forward probabilities 
// obtained from phase2 dynamic programming
inline
size_t CuDeviceMemory::GetSizeOfDP2FwdMtx( size_t queryposs, size_t dbposs ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2FwdMtx", 8 );
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
// size_t CuDeviceMemory::GetSizeOfDP2BwdMtx( size_t queryposs, size_t dbposs) const
// {
//     MYMSG( "CuDeviceMemory::GetSizeOfDP2BwdMtx", 8 );
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
size_t CuDeviceMemory::GetSizeOfSS2Data( size_t /*dbposs*/, size_t dbpros ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfSS2Data", 8 );
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
size_t CuDeviceMemory::GetSizeOfDP2AlnData( size_t /*dbposs*/, size_t dbpros ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2AlnData", 8 );
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
size_t CuDeviceMemory::GetSizeOfDP2Alns( 
    size_t queryposs, size_t dbposs, size_t dbpros, bool sssinuse ) const
{
    MYMSG( "CuDeviceMemory::GetSizeOfDP2Alns", 8 );
    const size_t cszalnment = GetMemAlignment();
    //maximum alignment length: l_query + l_target
    size_t szdp2alns = 
            ((dbposs + dbpros * (queryposs+1)) *
             (sssinuse? nTDP2OutputAlignmentSSS: nTDP2OutputAlignment)) *
            sizeof(char);
    szdp2alns = ALIGN_UP( szdp2alns, cszalnment );
    return szdp2alns;
}



// -------------------------------------------------------------------------
// GetSizeOfGlobVariables: get size of the section of global variables
inline
size_t CuDeviceMemory::GetSizeOfGlobVariables() const
{
    MYMSG( "CuDeviceMemory::GetSizeOfGlobVariables", 8 );
    const size_t cszalnment = GetMemAlignment();
    size_t szglbvars = nDevGlobVariables * sizeof(int);
    szglbvars = ALIGN_UP( szglbvars, cszalnment );
    return szglbvars;
}



// =========================================================================
// GetMaxAllowedNumDbProfilePositions: get the maximum allowed number of Db 
// profile positions given query length to avoid aoverflow
inline
size_t CuDeviceMemory::GetMaxAllowedNumDbProfilePositions( 
    size_t queryposs, bool sssinuse ) const
{
    MYMSG( "CuDeviceMemory::GetMaxAllowedNumDbProfilePositions", 7 );
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
void CuDeviceMemory::GetTotalMemoryReqs( 
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
    size_t* szdp2maxcoords, size_t* szdp2btckdat,
    size_t* szdp2fwdmtx,// size_t* szdp2bwdmtx, 
    size_t* szss2data, size_t* szdp2alndata, size_t* szdp2alns,
    //}}{{variables
    size_t* szglbvars,
    //}}{{overall size
    size_t* szovlpos
    //}}
    ) const
{
    MYMSG( "CuDeviceMemory::GetTotalMemoryReqs", 7 );
    const size_t cszalnment = GetMemAlignment();
    const float pass2memperc = CLOptions::GetDEV_PASS2MEMP();
    size_t sum = 0UL;//temporary sum variable 
    //{{memory for db profile data itself
    *maxsizedbposs = PMBatchProData::GetPMDataSizeUBTotal( maxdbposs );
    *maxsizedbposs = SLC_MIN( *maxsizedbposs, maxamountfordbposs );
    *maxsizedbposs += pmv2DTotFlds * cszalnment;//all fields aligned;
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
    *szdp2maxcoords = GetSizeOfDP2MaxCoords( *maxdbposspass2 );
    *szdp2btckdat = GetSizeOfDP2BackTckData( nqyposs, *maxdbposspass2 );
    *szdp2fwdmtx = GetSizeOfDP2FwdMtx( nqyposs, *maxdbposspass2 );
//*szdp2bwdmtx = GetSizeOfDP2BwdMtx( nqyposs, *maxdbposspass2 );
    *szss2data = GetSizeOfSS2Data( *maxdbposspass2, *maxndbprospass2 );
    *szdp2alndata = GetSizeOfDP2AlnData( *maxdbposspass2, *maxndbprospass2 );
    *szdp2alns = GetSizeOfDP2Alns( nqyposs, *maxdbposspass2, *maxndbprospass2 );
    *szglbvars = GetSizeOfGlobVariables();
    size_t ssoutdat = *szss2data + *szdp2alndata + *szdp2alns;
    sum = *szdp2diag + *szdp2bottom + 
          *szdp2maxcoords + *szdp2btckdat + 
          *szdp2fwdmtx + 
//*szdp2bwdmtx + 
          ssoutdat;
    //}}
    *szovlpos = *szdpdiag + *szdpbottom  + *szdpmaxcoords + *szdpbtckdat;
    if( !GetMAPAlnInUse())
        *szovlpos += ssoutdat;
    else if( *szovlpos < sum )
        *szovlpos = sum;
    *szovlpos += *maxsizedbposs + *szsmatrix;
    *szovlpos += *szglbvars;
    if( GetModScoreMatrixInUse())
        *szovlpos += *szsmatrix;
}

#endif//__CuDeviceMemory_h__
