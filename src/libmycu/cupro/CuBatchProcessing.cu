/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory>
#include <mutex>

#include "tsafety/TSCounterVar.h"

#include "extsp/psl.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/TCTXVECT.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cualn/AlnWriter.h"
#include "libmycu/cupro/SerializedScores.cuh"
#include "libmycu/cupro/SerializedScoresCtor.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresCtor.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cusco/CuBatchScoreMatrix.cuh"
#include "libmycu/cudp/CuBatchDP.cuh"
#include "libmycu/cudp/CuBatchDP_final.cuh"
#include "libmycu/cuss/CuBatchSS.cuh"
#include "libmycu/cumapdp/CuBatchMAPDP.cuh"
#include "CuBatchProcessingFinalizer.h"
#include "CuBatchProcessing_com.cuh"
#include "CuBatchProcessing.cuh"
#include "CuDeviceMemory.cuh"

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// constructor
//
CuBatchProcessing::CuBatchProcessing(
        CuDeviceMemory* dmem,
        int dareano,
        Configuration config[NoCTypes], 
        AlnWriter* writer,
        const mystring* queryfnames,
        const mystring* querydescs,
        const char** bdb1descs,
        int precscale )
:
    dmem_(dmem),
    devareano_(dareano),
    scorethld_(0.0f),
    logevthld_(0.0f),
    querylength_(0UL),
    dblength_(0UL),
    ndbsequences_(0UL),
    curdbxpad_( 0 ),
    dbxpadphase2_( 0 ),
    dbalnlenphase2_( 0 ),
    //
    h_results_(NULL),
    lockedresmem_(false),
    sz_mem_results_(0UL),
    limit_beg_results_(0UL)
{
    MYMSG( "CuBatchProcessing::CuBatchProcessing", 4 );

    if( dmem_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Null device memory object.");

    if( devareano_ < 0 || dmem_->GetNAreas() <= devareano_ )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Invalid device memory area number.");

    cbsm_ = new CuBatchScoreMatrix( config, precscale );
    if( cbsm_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");

    cbdp_ = new CuBatchDP();
    if( cbdp_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");

    cbss_ = new CuBatchSS();
    if( cbss_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");

    cbmapdp_ = new CuBatchMAPDP();
    if( cbmapdp_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");

    MYCUDACHECK( cudaSetDevice( dmem_->GetDeviceProp().devid_ ));
    MYCUDACHECKLAST;

    cbpfin_ = NULL;

    MYCUDACHECK( cudaStreamCreate(&streamcopyres_));
    MYCUDACHECKLAST;

    cbpfin_ = new CuBatchProcessingFinalizer(
            streamcopyres_, 
            dmem_->GetDeviceProp(),
            writer,
            config,
            queryfnames,
            querydescs,
            bdb1descs
    );

    if( cbpfin_ == NULL )
        throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");
}

// -------------------------------------------------------------------------
// default constructor
//
CuBatchProcessing::CuBatchProcessing()
:   devareano_(-1)
{
    throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: "
                "Default initialization is not allowed.");
}

// -------------------------------------------------------------------------
// destructor
//
CuBatchProcessing::~CuBatchProcessing()
{
    MYMSG( "CuBatchProcessing::~CuBatchProcessing", 4 );
    DestroyCBSMObject();
    DestroyCBMAPDPObject();
    DestroyCBSSObject();
    DestroyCBDPObject();
    //
    if(cbpfin_) {
        cbpfin_->Notify(CuBatchProcessingFinalizer::cubpthreadmsgTerminate);
        delete cbpfin_;
        cbpfin_ = NULL;
    }

    MYCUDACHECK( cudaStreamDestroy( streamcopyres_ ));
    MYCUDACHECKLAST;

    //free memory after the finalizer has finished
    HostFreeResults();

    //device can still be in use by Devices
//     MYCUDACHECK( cudaDeviceReset());
//     MYCUDACHECKLAST;
}



// -------------------------------------------------------------------------
// CheckHostResultsSync: verify whether, and perform if needed, 
// synchronization of the results transfer is required
inline
void CuBatchProcessing::CheckHostResultsSync()
{
//     size_t szresults = //size allocted for results
//         sz_heapsections_[ddsEndOfDP2Alns] - 
//         sz_heapsections_[ddsEndOfSS2Data];

    //synchronize only if the current end of phase-1 section
    // overlaps with the beginning of the space allocated previously for 
    // results (OR the required size has increased --- moved to finalization)
    if( limit_beg_results_ < GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfDPBackTckData)
        /*|| sz_mem_results_ < szresults*/) 
    {
        //make sure of no results data transfer:
        MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
        MYCUDACHECKLAST;
    }

    limit_beg_results_ = GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data);
}

// -------------------------------------------------------------------------
// ComputeScoreMatrix: compute score matrices for part of query and database 
// profiles
//
void CuBatchProcessing::ProcessScoreMatrix(
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
    TSCounterVar* cnt )
{
    MYMSG( "CuBatchProcessing::ComputeScoreMatrix", 6 );
    BatchProcessScoreMatrixDevice(
        qrysernr,
        nqyposs,
        querypmbeg,// querypmend,
        bdb1pmbeg, bdb1pmend,
        bdbCdesc, bdbCpmbeg, bdbCpmend, szCpm2dvfields,
        cnt );
}





// =========================================================================
// -------------------------------------------------------------------------
// BatchProcessScoreMatrixDevice: compute score matrices for part of query 
// and database profiles on device;
// querypmbeg, relative beginning address of queries;
// querypmend, relative terminal address of queries;
// bdb1pmbeg, relative beginning address of cached database profiles;
// bdb1pmend, relative terminal address of cached database profiles;
// bdbCpmbeg, beginning address of new profiles read from the database; 
// bdbCpmend, terminal address of new profiles read from the database;
// NOTE: query and cached db data addresses are relative to the beginning 
// addresses stored in host memebers h_querypmbeg_, h_querypmend_, 
// h_bdb1pmbeg_, and h_bdb1pmend_;
//
void CuBatchProcessing::BatchProcessScoreMatrixDevice(
    int qrysernr,
    size_t nqyposs,
    char** querypmbeg,
    //char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend,
    const char** bdbCdesc,
    char** bdbCpmbeg,
    char** bdbCpmend,
    size_t* /*szCpm2dvfields*/,
    TSCounterVar* cnt )
{
    MYMSG( "CuBatchProcessing::BatchProcessScoreMatrixDevice", 4 );
    const mystring preamb = "CuBatchProcessing::BatchProcessScoreMatrixDevice: ";
    myruntime_error mre;
    char msgbuf[KBYTE];

#ifdef __DEBUG__
    if( !querypmbeg || /*!querypmend || */!dmem_->GetHstQueryPMBeg())
        throw MYRUNTIME_ERROR( preamb + "Null query addresses.");
    if( cbsm_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch score matrix object.");
    if( cbdp_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch DP object.");
    if( cbss_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch SS object.");
    if( cbmapdp_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch MAP DP object.");
#endif

    unsigned int attrpassed[CuDeviceMemory::dgvMutex] = {0,0,0};
    //attributes for the allowed number of passed profiles:
    unsigned int attrpassedallwd[CuDeviceMemory::dgvMutex] = {0,0,0};

    size_t ndb1pros = 0;//total number of profiles in the first buffer
    size_t ndbCpros = 0;//total number of profiles in the complementary buffer
    size_t dbpro1len = 0;//the length of the first profile, which is the largest length

    size_t querprosOmtd = PMBatchProData::GetNoProsFromTo( dmem_->GetHstQueryPMBeg(), querypmbeg );
    size_t ndb1prosOmtd = 0;//number of profiles omitted in the first buffer
    size_t ndbCprosOmtd = 0;//number of profiles omitted in the complementary buffer

//     uint querposoffset = (uint)PMBatchProData::GetNoProsFromTo( dmem_->GetHstQueryPMBeg(), querypmbeg );
    size_t querposoffset = PMBatchProData::GetNoPositsFromTo( dmem_->GetHstQueryPMBeg(), querypmbeg );
    size_t bdb1posoffset = 0;
    size_t bdbCposoffset = 0;

    //const size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmbeg, querypmend );
    const float qyeno = PMBatchProData::GetPMDataENO1At( querypmbeg );
    size_t ndb1poss = 0UL;
    size_t ndbCposs = 0UL;

    if( bdbCpmbeg && bdbCpmend ) {
        ndbCposs = PMBatchProData::GetNoPositsFromTo( bdbCpmbeg, bdbCpmend );
        ndbCpros = PMBatchProData::GetNoProsFromTo( bdbCpmbeg, bdbCpmend );
        dbpro1len = PMBatchProData::GetPMDataLen1At( bdbCpmbeg );
    }
    if( bdb1pmbeg && bdb1pmend ) {
        ndb1poss = PMBatchProData::GetNoPositsFromTo( bdb1pmbeg, bdb1pmend );
        ndb1pros = PMBatchProData::GetNoProsFromTo( bdb1pmbeg, bdb1pmend );
#ifdef __DEBUG__
        if( !dmem_->GetHstBdb1PMBeg())
            throw MYRUNTIME_ERROR( preamb + "Null cached data address.");
#endif
        bdb1posoffset = PMBatchProData::GetNoPositsFromTo( dmem_->GetHstBdb1PMBeg(), bdb1pmbeg );
        ndb1prosOmtd = PMBatchProData::GetNoProsFromTo( dmem_->GetHstBdb1PMBeg(), bdb1pmbeg );
        //bdb1pmbeg..bdb1pmend come in front of bdbCpmbeg..bdbCpmend, 
        //the first profile length here is hence the largest
        dbpro1len = PMBatchProData::GetPMDataLen1At( bdb1pmbeg );
    }

    size_t ndbxposs = ndb1poss + ndbCposs;

    //padding along the x axis
    const size_t dbxpad = ALIGN_UP( ndbxposs, CUL2CLINESIZE ) - ndbxposs;
    SetCurrentDbxPadding((unsigned int)dbxpad);

// #ifdef __DEBUG__
    size_t szsspace = nqyposs * ndbxposs;//search space
    if( ndb1pros + ndbCpros < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles of part of db.");
    if( szsspace < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of search space.");
    if( GetCurrentMaxDbPos() < ndb1poss + ndbCposs + dbxpad ) {
        sprintf(msgbuf, "Invalid accumulated number of Db profile positions: %zu < %zu+%zu+%zu.",
            GetCurrentMaxDbPos(), ndb1poss, ndbCposs, dbxpad );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }
// #endif

    cudaStream_t streamproc;
    MYCUDACHECK( cudaStreamCreate(&streamproc));

    try {
        cbsm_->ComputeScoreMatrixDevice(
            streamproc,
            dmem_->GetMAPAlnInUse(),
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfPadding)/*sssscores*/,
            dmem_->GetSSSSAttr(),
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSSSTable)/*cvs2scores*/,
            dmem_->GetCVS2SAttr(),
            dmem_->GetHDP1STexObj(),
            dmem_->GetHDPSAttr(),
            dmem_->GetHDP1ScoresInUse(),
            nqyposs,
            ndb1poss,
            ndbCposs,
            dbxpad,
            querposoffset,
            bdb1posoffset,
            bdbCposoffset,
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)/*outscores*/,
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfOrgScores)/*outmodscores*/
        );

        if( !dmem_->GetMAPAlnInUse()) {
            //make sure results are not being copied before overwriting them:
            MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
            MYCUDACHECKLAST;
        }

        cbdp_->PerformCompleteDynProgDevice(
            streamproc,
            GetScoreThreshold(),
            ndb1pros,
            ndbCpros,
            querprosOmtd,
            ndb1prosOmtd,
            ndbCprosOmtd,
            dbpro1len,
            nqyposs,
            ndb1poss,
            ndbCposs,
            dbxpad,
            querposoffset,
            bdb1posoffset,
            bdbCposoffset,
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)/*scores [in]*/,
            (CUBDP_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfModScores)/*tmpdpdiagbuffers [out]*/,
            (CUBDP_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPDiagScores)/*tmpdpbotbuffer [out]*/,
            (unsigned int*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPBottomScores)/*maxcoordsbuf [out]*/,
            (char*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPMaxCoords)/*btckdata [out]*/,
            //[out:]
            (unsigned int*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfGlobVars)/*globvarsbuf [out]*/
        );

        if( querylength_ && dblength_ && ndbsequences_ )
            cbsm_->ComputeLengthAdjustment( querylength_, dblength_, ndbsequences_ );

        MYCUDACHECK( cudaStreamSynchronize(streamproc));
        MYCUDACHECKLAST;

        //make sure results have been copied before overwriting them:
        MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
        MYCUDACHECKLAST;

        MYCUDACHECK( cudaMemcpy( attrpassed,
            GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfGlobVars),
            sizeof(attrpassed[0]) * CuDeviceMemory::dgvMutex,
            cudaMemcpyDeviceToHost));
        MYCUDACHECKLAST;

        MYMSGBEGl(1)
            sprintf(msgbuf,"%6c%u profiles passed to the realignment phase",' ',
                    attrpassed[CuDeviceMemory::dgvNPassedPros]);
            MYMSG(msgbuf,1);
        MYMSGENDl

        MYMSGBEGl(3)
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sPhase-1 results summary: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf,
                "# db_pros_passed_to_phase_2= %u (thld= %.1f), total positions= %u, max len= %u",
                attrpassed[CuDeviceMemory::dgvNPassedPros], GetScoreThreshold(), 
                attrpassed[CuDeviceMemory::dgvNPosits], attrpassed[CuDeviceMemory::dgvMaxProLen]);
            strbuf += msgbuf;
            //MYMSG(strbuf.c_str(),3);
            sprintf(msgbuf,
                "%s(phase-2 limit values: # profiles= %zu, # positions= %zu)",
                NL, GetCurrentMaxNDbProsPass2(), GetCurrentMaxDbPosPass2());
            //strbuf = preamb + msgbuf;
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        for( int i = 0; i < CuDeviceMemory::dgvMutex; i++)
            attrpassedallwd[i] = attrpassed[i];

        if( GetCurrentMaxNDbProsPass2() < attrpassed[CuDeviceMemory::dgvNPassedPros] ||
            GetCurrentMaxDbPosPass2() < attrpassed[CuDeviceMemory::dgvNPosits]) 
        {
            warning("RESULTS WILL NOT BE COMPLETE. "
                    "Increase the value of option --dev-pass2memp and rerun.");
            if( GetCurrentMaxNDbProsPass2() < attrpassed[CuDeviceMemory::dgvNPassedPros])
                attrpassedallwd[CuDeviceMemory::dgvNPassedPros] = (unsigned int)GetCurrentMaxNDbProsPass2();
            if( GetCurrentMaxDbPosPass2() < attrpassed[CuDeviceMemory::dgvNPosits])
                attrpassedallwd[CuDeviceMemory::dgvNPosits] = (unsigned int)GetCurrentMaxDbPosPass2();
        }

        cbss_->CalculateAlnStatisticsDevice(
            streamproc,
            dmem_->GetXUninformative(),
            ndb1pros,
            ndbCpros,
            ndb1prosOmtd,
            ndbCprosOmtd,
            nqyposs,
            ndb1poss,
            ndbCposs,
            dbxpad,
            querposoffset,
            bdb1posoffset,
            bdbCposoffset,
            //
            MOptions::GetSSEMODEL(),
            qyeno,
            (float)cbsm_->GetSearchSpace(),
            cbsm_->GetRefLambda(), cbsm_->GetRefK(),
            cbsm_->GetExpGappedLambda(), cbsm_->GetExpGappedK(),
            //
            (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)/*scores [in]*/,
            (CUBDP_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfModScores)/*tmpdpdiagbuffers [in]*/,
            //(CUBDP_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPDiagScores)/*tmpdpbotbuffer [in]*/,
            (float*)GetAddrOfHeapSS2Data()/*tmpss2datbuffers [out]*/,
            (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data)/*dp2alndatbuffers [out]*/,
            attrpassedallwd/*[in]*/
        );

        //padding along the x axis for db profile data in the phase-2 processing
        const size_t dbxpad2 =
            ALIGN_UP( attrpassedallwd[CuDeviceMemory::dgvNPosits], CUL2CLINESIZE ) - 
                    attrpassedallwd[CuDeviceMemory::dgvNPosits];
        SetCurrentDbxPaddingPass2((unsigned int)dbxpad2);
        if( GetCurrentMaxDbPosPass2() < (size_t)attrpassedallwd[CuDeviceMemory::dgvNPosits] + dbxpad2 )
        {
            sprintf(msgbuf, "Too large number of Db profile positions in pass 2: %zu < %u+%zu.",
                GetCurrentMaxDbPosPass2(), attrpassedallwd[CuDeviceMemory::dgvNPosits], dbxpad2);
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }
        //set db length for alignments, including padding
        size_t cbdbalnlen2 = 
            attrpassedallwd[CuDeviceMemory::dgvNPosits] +
            attrpassedallwd[CuDeviceMemory::dgvNPassedPros] * (nqyposs+1);
        cbdbalnlen2 = ALIGN_UP( cbdbalnlen2, CUL2CLINESIZE );
        SetCurrentDbAlnLengthWithPaddingPass2((unsigned int)cbdbalnlen2);

        if( dmem_->GetMAPAlnInUse())
            cbmapdp_->PerformCompleteMAPDynProgDevice(
                streamproc,
                GetLogEThreshold(),
                dmem_->GetModScoreMatrixInUse(),
                ndb1pros,
                ndbCpros,
                querprosOmtd,
                ndb1prosOmtd,
                ndbCprosOmtd,
                attrpassed[CuDeviceMemory::dgvMaxProLen],
                nqyposs,
                ndb1poss,
                ndbCposs,
                dbxpad,
                querposoffset,
                bdb1posoffset,
                bdbCposoffset,
                (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)/*scores [in]*/,
                (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfOrgScores)/*mod scores [in]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfModScores)/*tmpdpdiagbuffers [tmp]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2DiagScores)/*tmpdpbotbuffer [tmp]*/,
                (unsigned int*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2BottomScores)/*maxcoordsbuf [tmp]*/,
                (char*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2MaxCoords)/*btckdata [tmp]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2BackTckData)/*tmpfwdprobs [tmp]*/,
//(float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2FwdMtx)/*tmpbwdprobs [tmp]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2FwdMtx)/*tmpss2datbuffers [tmp]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data)/*dp2alndatbuffers [in]*/,
                (char*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2AlnData)/*output alignments [out]*/,
                attrpassedallwd,//[in]
                dbxpad2,
                cbdbalnlen2
            );
        else//(GetMAPAlnInUse()==false)
            cbdp_->FinalizeALNDynProgDevice(
                streamproc,
                GetLogEThreshold(),
                ndb1pros,
                querprosOmtd,
                ndb1prosOmtd,
                ndbCprosOmtd,
                nqyposs,
                ndb1poss,
                ndbCposs,
                dbxpad,
                querposoffset,
                bdb1posoffset,
                bdbCposoffset,
                (CUBSM_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)/*scores [in]*/,
                (CUBDP_TYPE*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfModScores)/*tmpdpdiagbuffers [in]*/,
                (float*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data)/*dp2alndatbuffers [in/out]*/,
                attrpassedallwd,//[in]
                dbxpad2,
                cbdbalnlen2,
                (unsigned int*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPBottomScores)/*maxcoordsbuf [in]*/,
                (char*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDPMaxCoords)/*btckdata [in]*/,
                (char*)GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfDP2AlnData)/*output alignments [out]*/
            );

        MYCUDACHECK( cudaStreamSynchronize(streamproc));
        MYCUDACHECKLAST;

        TransferResultsDevice(
            qrysernr,
            nqyposs, qyeno,
            querypmbeg,// querypmend,
            bdb1pmbeg, bdb1pmend,
            bdbCdesc, bdbCpmbeg, bdbCpmend,
            cnt,
            ndb1pros,
            ndbCpros,
            querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
            querposoffset, bdb1posoffset, bdbCposoffset,
            attrpassedallwd[CuDeviceMemory::dgvNPosits],
            attrpassedallwd[CuDeviceMemory::dgvNPassedPros]);

// if( nqyposs == /*26*/99999/*394*//*1074*//*967*//*1122*/ )
// cbsm_->TESTPrintProProScores1(
//     dmem_->GetHstQueryPMBeg(), dmem_->GetHstQueryPMEnd(),
//     querypmbeg, querypmend,
//     bdb1pmbeg, bdb1pmend, bdbCpmbeg, bdbCpmend, GetCurrentDbxPadding(),
//     GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfOrgScores) - GetOffsetOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores)
//     ,(CUBSM_TYPE*)(dmem_->GetHeap() + GetOffsetOfHeapSection(CuDeviceMemory::ddsBegOfOrgScores))/*outscores*/
//     // ,(CUBSM_TYPE*)(dmem_->GetHeap() + GetOffsetOfHeapSection(ddsEndOfOrgScores))/*outmodscores*/
// );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    MYCUDACHECK( cudaStreamDestroy( streamproc ));
    MYCUDACHECKLAST;

    if( mre.isset())
        throw mre;
}





// =========================================================================
// TransferResultsDevice: finalize processing performed on device, 
// copy results and format them;
// nqyposs, query length;
// qyeno, query ENO;
// querypmbeg, relative beginning address of queries;
// querypmend, relative terminal address of queries;
// bdb1pmbeg, relative beginning address of cached database profiles;
// bdb1pmend, relative terminal address of cached database profiles;
// bdbCdesc, descriptions of new profiles read from the database;
// bdbCpmbeg, beginning address of new profiles read from the database; 
// bdbCpmend, terminal address of new profiles read from the database;
// cnt, counter associated with how many agents access (iterations 
// performed on) the data read from the database;
// NOTE: query and cached db data addresses are relative to the beginning 
// addresses stored in host memebers h_querypmbeg_, h_querypmend_, 
// h_bdb1pmbeg_, and h_bdb1pmend_;
//
void CuBatchProcessing::TransferResultsDevice(
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
    unsigned int npros )
{
    MYMSG( "CuBatchProcessing::TransferResultsDevice", 3 );
    mystring preamb = "CuBatchProcessing::TransferResultsDevice: ";
    myruntime_error mre;

    size_t szresults = //size allocated for results
        GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfDP2Alns) - 
        GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data);

    //make sure of no pending transfer in this stream
    MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
    MYCUDACHECKLAST;

#ifdef __DEBUG__
//     if( !h_results_)
//         throw MYRUNTIME_ERROR( preamb + "Null results address.");
    if( !cbpfin_)
        throw MYRUNTIME_ERROR( preamb + "Null finalizer object.");
#endif

    try {
        int rspcode = cbpfin_->GetResponse();
        if( rspcode == CUBPTHREAD_MSG_ERROR )
            throw MYRUNTIME_ERROR(preamb + "Finalizer terminated with errors.");

        if( nposits && npros )
        {//lock this section as the data is about to change
            std::lock_guard<std::mutex> lck(cbpfin_->GetPrivateMutex());
            //
            //allocate memory if the required size has increased
            if( sz_mem_results_ < szresults ) {
                //allocate more memory to reduce the frequency of allocation operations
                size_t szalloc = szresults + szresults/2;
                HostAllocResults(szalloc);
                sz_mem_results_ = szalloc;
            }
            if( lockedresmem_ ) {
                MYCUDACHECK( cudaMemcpyAsync(
                    h_results_, GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data),
                    szresults, cudaMemcpyDeviceToHost, streamcopyres_));
                MYCUDACHECKLAST;
            }
            else {
                MYCUDACHECK( cudaMemcpy(
                    h_results_, GetAddrOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data),
                    szresults, cudaMemcpyDeviceToHost));
                MYCUDACHECKLAST;
            }
        }

        cbpfin_->SetCuBPBDbdata(
            qrysernr,
            nqyposs, qyeno,
            (unsigned int)cbsm_->GetDeltaLength(),
            (float)cbsm_->GetSearchSpace(),
            GetLogEThreshold(),
            querypmbeg,// querypmend, 
            bdb1pmbeg, bdb1pmend,
            bdbCdesc, bdbCpmbeg, bdbCpmend,
            cnt,
            ndb1pros,
            ndbCpros,
            querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
            querposoffset, bdb1posoffset, bdbCposoffset,
            nposits,
            npros,
            h_results_,
            GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfDP2AlnData)-
                GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfSS2Data),
            GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfDP2Alns)-
                GetOffsetOfHeapSection(CuDeviceMemory::ddsEndOfDP2AlnData),
            GetCurrentDbAlnLengthWithPaddingPass2()
        );
        cbpfin_->Notify(CuBatchProcessingFinalizer::cubpthreadmsgFinalize);

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}
