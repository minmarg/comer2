/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory>
#include <mutex>

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

// -------------------------------------------------------------------------
// static data to be transfered to device CMEM
// 
// *** *************************************************** ***
// *** RE MODEL: genpro-fit2-rnd3-fr9-r0.03-RE-altseed     ***
// *** *************************************************** ***
// *** ONE LAYER                                           ***
// *** *************************************************** ***
static const float SSE_NN_WEIGHTS_[] = {
    //
    //WEIGHTS_MU: rows -- units, columns -- inputs (0th -- biases)
     3.1179029f,     0.2732032f,    -1.7865388f,    -2.7370650f,    -2.4113938f,
     2.8802990f,    -0.0845794f,     0.7335855f,    -0.4155596f,    -3.5943765f,
    -2.3533863f,    -0.1646748f,     2.6470928f,     1.6489798f,    -0.7481717f,
    -0.1763256f,     0.3270687f,     0.5918245f,     0.4900585f,     0.0253580f,

     2.0547640f,     1.4051396f,    -1.7874640f,     2.0048461f,     0.6947985f,
    //
    //WEIGHTS_L_MU: rows -- units, columns -- inputs (0th -- biases)
    // *** EVD of scores grouped by Lmb
    // *** Exclusively-Lmb-based
    // *** ONE LAYER
    //(pad the last column of the first layer with zeros for parallel processing)
    -0.0605839f,     0.0533123f,    -1.0046299f,     0.9682557f,     0.0f,
    -0.0173019f,    -0.0051053f,     0.6215832f,    -0.5038284f,     0.0f,
     0.0600411f,    -0.0527767f,     1.0036328f,    -0.9667071f,     0.0f,
    -5.9555713f,    -3.1437676f,     1.8754493f,     3.4776044f,     0.0f,

     2.8086869f,    -0.5250554f,     0.3085814f,     0.5244239f,     2.4166220f,
    //
    //WEIGHTS_L_SCALE: rows -- units, columns -- inputs (0th -- biases)
    // *** regularized
    //(pad the last column of the first layer with zeros for parallel processing)
    -4.4935289f,   -11.5416831f,     0.4601384f,     4.1707612f,     0.0f,
     1.4710020f,     0.7746900f,     0.5248285f,    -2.2621661f,     0.0f,
    -1.4710009f,    -0.7746948f,    -0.5248293f,     2.2621667f,     0.0f,
     4.8137710f,     1.0000337f,     2.0837171f,    -5.3104535f,     0.0f,

     2.5137293f,     4.2316344f,    -0.8620832f,     0.8620837f,     3.1901578f,
    //
    //
    //Adjustemnt parameters for model's SS18 variant 1
     0.1f, 12.0f,  1.0f, 4.0f,  0.05f, 1.0f,   0.1f,  0.35f,  0.65f,
    //Adjustemnt parameters for model's SS18 variant 2
     0.6f, 16.0f,  1.0f,12.0f,  0.3f, -1.7f,   0.1f,  0.45f,  0.65f
};



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// constructor
//
CuBatchProcessing::CuBatchProcessing(
        AlnWriter* writer,
        Configuration config[],
        DeviceProperties dprop, 
        bool hdp1scoresinuse,
        bool mapalninuse,
        size_t deviceallocsize,
        const mystring* queryfnames,
        const mystring* querydescs,
        const mystring* bdb1fnames,
        const mystring* bdb1descs,
        int precscale )
:
    hdp1scoresinuse_( hdp1scoresinuse ),
    mapalninuse_( mapalninuse ),
    scorethld_(0.0f),
    logevthld_(0.0f),
    querylength_(0UL),
    dblength_(0UL),
    ndbsequences_(0UL),
    deviceallocsize_( deviceallocsize ),
    curmaxdbpos_( 0UL ),
    curmaxndbpros_( 0UL ),
    curdbxpad_( 0 ),
    curmaxdbpospass2_( 0UL ),
    curmaxndbprospass2_( 0UL ),
    dbxpadphase2_( 0 ),
    dbalnlenphase2_( 0 ),
    h_querypmbeg_(NULL),
    h_querypmend_(NULL),
    h_bdb1pmbeg_(NULL),
    h_bdb1pmend_(NULL),
    h_bdbCpmbeg_(NULL),
    h_bdbCpmend_(NULL),
    //
    h_results_(NULL),
    lockedresmem_(false),
    sz_mem_results_(0UL),
    limit_beg_results_(0UL),
    //
    deviceProp_(dprop),
    hdp1sTexObj_(0),
    d_heap_(NULL)
{
    MYMSG( "CuBatchProcessing::CuBatchProcessing", 4 );

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

    MYCUDACHECK( cudaSetDevice( deviceProp_.devid_ ));
    MYCUDACHECKLAST;

    cbpfin_ = NULL;

    for( int i = 0; i < nDevDataSectionsDiv2; i++ )
        sz_heapsections_[i] = 0UL;

    if( deviceallocsize_ ) {

        const size_t cszalnment = GetMemAlignment();

        MYCUDACHECK( cudaStreamCreate(&streamcopyres_));
        MYCUDACHECKLAST;

        MYCUDACHECK( cudaMalloc((void**)&d_heap_, deviceallocsize_ ));
        MYCUDACHECKLAST;

        size_t heapaligned = ALIGN_UP((size_t)d_heap_, cszalnment );
        sz_heapsections_[ddsEndOfPadding] = (size_t)d_heap_ - heapaligned;

        cbpfin_ = new CuBatchProcessingFinalizer(
                streamcopyres_, 
                deviceProp_,
                writer,
                config,
                queryfnames,
                querydescs,
                bdb1fnames,
                bdb1descs
        );
        if( cbpfin_ == NULL )
            throw MYRUNTIME_ERROR("CuBatchProcessing::CuBatchProcessing: Not enough memory.");
    }
}

// -------------------------------------------------------------------------
// default constructor
//
CuBatchProcessing::CuBatchProcessing()
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
    DestroyTextureObject( hdp1sTexObj_ );
    FreeDevicePtr( d_heap_ );
    //
    if(cbpfin_) {
        cbpfin_->Notify(CuBatchProcessingFinalizer::cubpthreadmsgTerminate);
        delete cbpfin_;
        cbpfin_ = NULL;
    }
    if( deviceallocsize_ ) {
        MYCUDACHECK( cudaStreamDestroy( streamcopyres_ ));
        MYCUDACHECKLAST;
    }
    //free memory after the finalizer has finished
    HostFreeResults();

    MYCUDACHECK( cudaDeviceReset());
    MYCUDACHECKLAST;
}





// -------------------------------------------------------------------------
// CacheSSENNWeights: cache NN weights for significance estimation
// 
void CuBatchProcessing::CacheSSENNWeights()
{
    MYMSG( "CuBatchProcessing::CacheSSENNWeights", 6 );
    if( MOptions::GetSSEMODEL() == 0 )
        return;
    CacheSSENNWeightsDevice();
}

// -------------------------------------------------------------------------
// CacheSSSScores: cache a score table of secondary structure predictions
// 
void CuBatchProcessing::CacheSSSScores( const SSSScores& sssscores )
{
    MYMSG( "CuBatchProcessing::CacheSSSScores", 6 );
    if( MOptions::GetSSSWGT() <= 0.0f )
        return;
    CacheSSSScoresDevice( sssscores );
}

// -------------------------------------------------------------------------
// CacheCVS2Scores: cache a map between cv scores and translated scores
// 
void CuBatchProcessing::CacheCVS2Scores( const CVS2Scores& cvs2scores )
{
    MYMSG( "CuBatchProcessing::CacheCVS2Scores", 6 );
    if( MOptions::GetCVSWGT() <= 0.0f )
        return;
    CacheCVS2ScoresDevice( cvs2scores );
}

// -------------------------------------------------------------------------
// CacheCVS2Scores: cache a map between cv scores and translated scores
// 
void CuBatchProcessing::CacheHDPScores( const HDPscores& hdpscores )
{
    MYMSG( "CuBatchProcessing::CacheHDPScores", 6 );
    if( !GetHDP1ScoresInUse() || MOptions::GetADJWGT() <= 0.0f )
        return;
    CacheHDPScoresDevice( hdpscores );
}


// -------------------------------------------------------------------------
// CacheData: cache query and db profile model data; host adresses are 
// saved for counting of positions later
// 
void CuBatchProcessing::CacheData(
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend )
{
    MYMSG( "CuBatchProcessing::CacheData", 6 );
    CacheDataDevice(
        h_querypmbeg_ = querypmbeg, 
        h_querypmend_ = querypmend, 
        h_bdb1pmbeg_ = bdb1pmbeg, 
        h_bdb1pmend_ = bdb1pmend );
}



// -------------------------------------------------------------------------
// GetMaxDbDataChunkSize: recalculate and return a maximum allowed number of 
// positions and memory size for db profile model data (profiles) so that 
// the implied search space does not exceed allocated memory (representing 
// the maximum allowed limit);
// querypmbeg, relative beginning address of queries;
// querypmend, relative terminal address of queries;
// 
size_t CuBatchProcessing::CalcMaxDbDataChunkSize( size_t nqyposs )
{
    MYMSG( "CuBatchProcessing::GetMaxDbDataChunkSize", 4 );
    const mystring preamb = "CuBatchProcessing::GetMaxDbDataChunkSize: ";
    const size_t cszalnment = GetMemAlignment();
    if( nqyposs < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of query positions.");
    //
    size_t tmpmaxndbpros, maxndbpros = 0UL;//maximum number of db profiles
    size_t tmpmaxndbprospass2, maxndbprospass2 = 0UL;//maximum number of db profiles in phase 2
    size_t tmpmaxdbposspass2, maxdbposspass2 = 0UL;//maximum number of db profile positions in phase 2
    //
    size_t tmpmaxdbposs, maxdbposs = 0UL;//maximum n. of db profile positions
    size_t tmpmaxsizedbposs, maxsizedbposs = 0UL;//maximum size (in bytes) for db profile positions
    size_t tmpszsmatrix, szsmatrix = 0UL;//maximum size for scores
    size_t tmpszdpdiag, szdpdiag = 0UL;//maximum size for DP diagonal score buffers
    size_t tmpszdpbottom, szdpbottom = 0UL;//maximum size for DP bottom score buffers
    size_t tmpszdpmaxcoords, szdpmaxcoords = 0UL;//maximum size for coordinates of max alignment scores
    size_t tmpszdpbtckdat, szdpbtckdat = 0UL;//maximum size for DP backtracking data
    //
    size_t tmpszdp2diag, szdp2diag = 0UL;//maximum size for phase-2 DP diagonal score buffers
    size_t tmpszdp2bottom, szdp2bottom = 0UL;//maximum size for phase-2 DP bottom score buffers
//**/size_t tmpszdp2diagblkprobscales, szdp2diagblkprobscales = 0UL;//maximum size for block-specific scale factors
//**/size_t tmpszdp2probscales, szdp2probscales = 0UL;//maximum size for final scale factors
    size_t tmpszdp2maxcoords, szdp2maxcoords = 0UL;//maximum size for coordinates of max alignment scores in phase 2 
    size_t tmpszdp2btckdat, szdp2btckdat = 0UL;//maximum size for phase-2 DP backtracking data
    size_t tmpszdp2fwdmtx, szdp2fwdmtx = 0UL;//maximum size for phase-2 forward probability matrix
//size_t tmpszdp2bwdmtx, szdp2bwdmtx = 0UL;//maximum size for phase-2 backward probability matrix
    size_t tmpszss2data, szss2data = 0UL;//maximum size for phase-2 statistical significance data
    size_t tmpszdp2alndata, szdp2alndata = 0UL;//maximum size for phase-2 alignment data
    size_t tmpszdp2alns, szdp2alns = 0UL;//maximum size for phase-2 alignments themselves
    //
    size_t tmpszovlpos, szovlpos = 0UL;//overall size of position-dependent data excluding maxsizedbposs
    size_t tmpdiff;
#ifdef __DEBUG__
    if( deviceallocsize_ < sz_heapsections_[ddsEndOfCached])
    {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: %zu.", deviceallocsize_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }
#endif

    size_t nmatrices = nDevMatrixSections;
    if( !GetModScoreMatrixInUse())
        nmatrices--;

    //seqch space divided by the n. of db profile positions
    //note: assuming matrices of type CUBSM_TYPE
    const size_t qsspace = nqyposs * nmatrices * sizeof(CUBSM_TYPE);
    const size_t leftoversize = deviceallocsize_ - sz_heapsections_[ddsEndOfCached];
    //maximum number of binary search iterations
    const int maxit = 10;
    const size_t desireddiff = ONEM;
    //maximum allowed number of Db profile positions to avoid overflow
    const size_t locMAXALLWDPOSS = GetMaxAllowedNumDbProfilePositions(nqyposs);
    //percentage of leftoversize corresponding to the max value for maxsizedbposs:
    // to keep the balance of (maximize) the number of db profile positions
    const float maxperc = 0.3f;
    size_t maxamountofleftover = (size_t)(leftoversize * maxperc);
    maxamountofleftover = ALIGN_UP( maxamountofleftover, cszalnment );
    //initial guess:
    tmpmaxdbposs = leftoversize / qsspace;
    tmpmaxdbposs = ALIGN_DOWN( tmpmaxdbposs, CUL2CLINESIZE );
    tmpmaxndbpros = tmpmaxdbposs / CLOptions::GetDEV_EXPCT_DBPROLEN();

    MYMSGBEGl(3)
        char msgbuf[BUF_MAX];
        mystring strbuf = preamb;
        sprintf(msgbuf,"Max allowed # db chunk positions: %zu",locMAXALLWDPOSS);
        strbuf += msgbuf;
        MYMSG(strbuf.c_str(),3);
    MYMSGENDl

    GetTotalMemoryReqs( 
        nqyposs, tmpmaxdbposs, maxamountofleftover,
        //db profile positions
        &tmpmaxsizedbposs,
        //scores
        &tmpszsmatrix,
        //phase1
        &tmpszdpdiag, &tmpszdpbottom, &tmpszdpmaxcoords, &tmpszdpbtckdat,
        //phase2
        &tmpmaxdbposspass2, &tmpmaxndbprospass2, 
        &tmpszdp2diag, &tmpszdp2bottom, 
//&tmpszdp2diagblkprobscales, &tmpszdp2probscales,
        &tmpszdp2maxcoords, &tmpszdp2btckdat,
        &tmpszdp2fwdmtx,// &tmpszdp2bwdmtx, 
        &tmpszss2data, &tmpszdp2alndata, &tmpszdp2alns,
        //overall size
        &tmpszovlpos
    );

    tmpdiff = tmpmaxdbposs;
    //binary search for finding maximum number of db profile positions
    for( int i = 0; i < maxit; i++ ) {
        MYMSGBEGl(5)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sAdjustment of total size for db profiles: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf, "it %d: %zu; maxdbpos= %zu; szdpdiag= %zu, szdpbottom= %zu, "
                            "szdpmaxcoords= %zu, szdpbtckdat= %zu; szovlpos= %zu; "
                            "leftoversize= %zu", i, 
                    tmpmaxsizedbposs, tmpmaxdbposs, tmpszdpdiag, tmpszdpbottom, 
                    tmpszdpmaxcoords, tmpszdpbtckdat, tmpszovlpos, leftoversize );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),5);
        MYMSGENDl
        tmpdiff /= 2;
        tmpdiff = ALIGN_DOWN( tmpdiff, CUL2CLINESIZE );
        if( tmpdiff == 0 )
            break;
        if(tmpszovlpos < leftoversize && 
           tmpmaxdbposs < locMAXALLWDPOSS ) 
        {
            //save this configuration
            maxdbposs = tmpmaxdbposs; 
            maxndbpros = tmpmaxndbpros;
            maxsizedbposs = tmpmaxsizedbposs;
            szsmatrix = tmpszsmatrix;
            szdpdiag = tmpszdpdiag; 
            szdpbottom = tmpszdpbottom;
            szdpmaxcoords = tmpszdpmaxcoords; 
            szdpbtckdat = tmpszdpbtckdat;
            maxdbposspass2 = tmpmaxdbposspass2; 
            maxndbprospass2 = tmpmaxndbprospass2; 
            szdp2diag = tmpszdp2diag; 
            szdp2bottom = tmpszdp2bottom; 
//szdp2diagblkprobscales = tmpszdp2diagblkprobscales; 
//szdp2probscales = tmpszdp2probscales;
            szdp2maxcoords = tmpszdp2maxcoords; 
            szdp2btckdat = tmpszdp2btckdat;
            szdp2fwdmtx = tmpszdp2fwdmtx; 
//szdp2bwdmtx = tmpszdp2bwdmtx; 
            szss2data = tmpszss2data; 
            szdp2alndata = tmpszdp2alndata; 
            szdp2alns = tmpszdp2alns;
            szovlpos = tmpszovlpos;
            //
            tmpmaxdbposs += tmpdiff;
        }
        else
            tmpmaxdbposs -= tmpdiff;

        tmpmaxndbpros = tmpmaxdbposs / CLOptions::GetDEV_EXPCT_DBPROLEN();

        if( maxdbposs && 
            GetAbsDiff(tmpszovlpos,leftoversize) < desireddiff )
            break;

        GetTotalMemoryReqs( 
            nqyposs, tmpmaxdbposs, maxamountofleftover,
            //db profile positions
            &tmpmaxsizedbposs,
            //scores
            &tmpszsmatrix,
            //phase1
            &tmpszdpdiag, &tmpszdpbottom, &tmpszdpmaxcoords, &tmpszdpbtckdat,
            //phase2
            &tmpmaxdbposspass2, &tmpmaxndbprospass2, 
            &tmpszdp2diag, &tmpszdp2bottom, 
//&tmpszdp2diagblkprobscales, &tmpszdp2probscales,
            &tmpszdp2maxcoords, &tmpszdp2btckdat,
            &tmpszdp2fwdmtx,// &tmpszdp2bwdmtx, 
            &tmpszss2data, &tmpszdp2alndata, &tmpszdp2alns,
            //overall size
            &tmpszovlpos
        );
    }


    if( maxdbposs < 1 || maxsizedbposs < 1 || 
        szovlpos < 1 || leftoversize <= szovlpos) {
        char msgbuf[KBYTE];
        sprintf( msgbuf, "Insufficient amount of leftover device memory: %zu%s"
                    "(maxszdbpos= %zu, szdpdiag= %zu, szdpbottom= %zu, "
                    "szdpmaxcoords= %zu, szdpbtckdat= %zu; ndbpos= %zu; "
                    "szovlpos= %zu).", 
                 leftoversize, NL, maxsizedbposs, szdpdiag, szdpbottom, 
                 szdpmaxcoords, szdpbtckdat, maxdbposs, szovlpos );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }


    MYMSGBEGl(3)
        char msgbuf[KBYTE];
        mystring strbuf = preamb;
        sprintf(msgbuf, "%sTotal size for db profiles "
                "(max %.1f of leftover): %zu; maxdbpos= %zu [maxndbpros= %zu]%s",
                NL, maxperc, maxsizedbposs, maxdbposs, maxndbpros, NL );
        strbuf += msgbuf;
        sprintf(msgbuf, " (dev alloc= %zu; dev leftover= %zu;%s "
                "szdpdiag= %zu, szdpbottom= %zu, "
                "szdpmaxcoords= %zu, szdpbtckdat= %zu; szovlpos= %zu)", 
                deviceallocsize_, leftoversize, NL, szdpdiag, szdpbottom, 
                szdpmaxcoords, szdpbtckdat, szovlpos );
        strbuf += msgbuf;
        MYMSG(strbuf.c_str(),3);
    MYMSGENDl

    SetCurrentMaxDbPos(maxdbposs);
    SetCurrentMaxNDbPros(maxndbpros);
    //{{phase 2
    SetCurrentMaxDbPosPass2( maxdbposspass2 );
    SetCurrentMaxNDbProsPass2( maxndbprospass2 );
    //}}

    sz_heapsections_[ddsEndOfDbChunk] = sz_heapsections_[ddsEndOfCached] + maxsizedbposs;
    sz_heapsections_[ddsEndOfOrgScores] = sz_heapsections_[ddsEndOfDbChunk] + szsmatrix;
    sz_heapsections_[ddsEndOfModScores] = sz_heapsections_[ddsEndOfOrgScores];
    if( GetModScoreMatrixInUse())
        sz_heapsections_[ddsEndOfModScores] += szsmatrix;
    //phase-1 sections
    sz_heapsections_[ddsEndOfDPDiagScores] = sz_heapsections_[ddsEndOfModScores] + szdpdiag;
    sz_heapsections_[ddsEndOfDPBottomScores] = sz_heapsections_[ddsEndOfDPDiagScores] + szdpbottom;
    sz_heapsections_[ddsEndOfDPMaxCoords] = sz_heapsections_[ddsEndOfDPBottomScores] + szdpmaxcoords;
    sz_heapsections_[ddsEndOfDPBackTckData] = sz_heapsections_[ddsEndOfDPMaxCoords] + szdpbtckdat;
    //phase-2 sections: reused previous occupations
    sz_heapsections_[ddsEndOfDP2DiagScores] = sz_heapsections_[ddsEndOfModScores] + szdp2diag;
    sz_heapsections_[ddsEndOfDP2BottomScores] = sz_heapsections_[ddsEndOfDP2DiagScores] + szdp2bottom;
//**/sz_heapsections_[ddsEndOfDP2DiagBlkProbScales] = sz_heapsections_[ddsEndOfDP2BottomScores] + szdp2diagblkprobscales;
//**/sz_heapsections_[ddsEndOfDP2PrbScales] = sz_heapsections_[ddsEndOfDP2DiagBlkProbScales] + szdp2probscales;
//sz_heapsections_[ddsEndOfDP2MaxCoords] = sz_heapsections_[ddsEndOfDP2PrbScales] + szdp2maxcoords;
    sz_heapsections_[ddsEndOfDP2MaxCoords] = sz_heapsections_[ddsEndOfDP2BottomScores] + szdp2maxcoords;
    sz_heapsections_[ddsEndOfDP2BackTckData] = sz_heapsections_[ddsEndOfDP2MaxCoords] + szdp2btckdat;
    sz_heapsections_[ddsEndOfDP2FwdMtx] = sz_heapsections_[ddsEndOfDP2BackTckData] + szdp2fwdmtx;
//sz_heapsections_[ddsEndOfDP2BwdMtx] = sz_heapsections_[ddsEndOfDP2FwdMtx] + szdp2bwdmtx;
//sz_heapsections_[ddsEndOfSS2Data] = sz_heapsections_[ddsEndOfDP2BwdMtx] + szss2data;
    sz_heapsections_[ddsEndOfSS2Data] = sz_heapsections_[ddsEndOfDP2FwdMtx] + szss2data;
    sz_heapsections_[ddsEndOfDP2AlnData] = sz_heapsections_[ddsEndOfSS2Data] + szdp2alndata;
    sz_heapsections_[ddsEndOfDP2Alns] = sz_heapsections_[ddsEndOfDP2AlnData] + szdp2alns;

    MsgAddressTable( preamb, 3 );


    //if needed, allocate host memeory here once
    // memory configuration has changed
    CheckHostResultsSync();


    return maxsizedbposs;
}

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
    if( limit_beg_results_ < sz_heapsections_[ddsEndOfDPBackTckData] 
        /*|| sz_mem_results_ < szresults*/) 
    {
        //make sure of no results data transfer:
        MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
        MYCUDACHECKLAST;
    }
    limit_beg_results_ = sz_heapsections_[ddsEndOfSS2Data];
}

// -------------------------------------------------------------------------
// GetMaxDbDataChunkSizeObs: recalculate and return a maximum allowed 
// number of positions and memory size for db profile model data (profiles) 
// so that the implied search space does not exceed allocated memory 
// (representing the maximum allowed limit);
// querypmbeg, relative beginning address of queries;
// querypmend, relative terminal address of queries;
// 
size_t CuBatchProcessing::CalcMaxDbDataChunkSizeObs(
    char** querypmbeg,
    char** querypmend )
{
    MYMSG( "CuBatchProcessing::GetMaxDbDataChunkSizeObs", 4 );
    const mystring preamb = "CuBatchProcessing::GetMaxDbDataChunkSizeObs: ";
    const size_t cszalnment = GetMemAlignment();
    const size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmbeg, querypmend );
    if( nqyposs < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of query positions.");
    //
    const float pass2memperc = CLOptions::GetDEV_PASS2MEMP();
    size_t maxndbpros = 0UL;//maximum number of db profiles
    size_t maxndbprospass2 = 0UL;//maximum number of db profiles in phase 2
    size_t maxdbposspass2 = 0UL;//maximum number of db profile positions in phase 2
    //
    size_t maxdbposs = 0UL;//maximum n. of db profile positions
    size_t maxsizedbposs = 0UL;//maximum size (in bytes) for db profile positions
    size_t szdpdiag = 0UL;//maximum size for DP diagonal score buffers
    size_t szdpbottom = 0UL;//maximum size for DP bottom score buffers
    size_t szdpmaxcoords = 0UL;//maximum size for coordinates of max alignment scores
    size_t szdpbtckdat = 0UL;//maximum size for DP backtracking data
    //
    size_t szdp2diag = 0UL;//maximum size for phase-2 DP diagonal score buffers
    size_t szdp2bottom = 0UL;//maximum size for phase-2 DP bottom score buffers
    size_t szdp2maxcoords = 0UL;//maximum size for coordinates of max alignment scores in phase 2 
    size_t szdp2btckdat = 0UL;//maximum size for phase-2 DP backtracking data
    size_t szdp2fwdmtx = 0UL;//maximum size for phase-2 forward probability matrix
    size_t szss2data = 0UL;//maximum size for phase-2 statistical significance data
    size_t szdp2alndata = 0UL;//maximum size for phase-2 alignment data
    size_t szdp2alns = 0UL;//maximum size for phase-2 alignments themselves
    //
    size_t sum = 0UL;//temporary sum variable 
    size_t szovlpos = 0UL;//overall size of position-dependent data excluding maxsizedbposs
#ifdef __DEBUG__
    if( deviceallocsize_ < sz_heapsections_[ddsEndOfCached])
    {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: %zu.", deviceallocsize_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }
#endif

    size_t nmatrices = nDevMatrixSections;
    if( !GetModScoreMatrixInUse())
        nmatrices--;

    //seqch space divided by the n. of db profile positions
    //note: assuming matrices of type CUBSM_TYPE
    const size_t qsspace = nqyposs * nmatrices * sizeof(CUBSM_TYPE);
    size_t leftoversize = deviceallocsize_ - sz_heapsections_[ddsEndOfCached];
    //percentage of leftoversize corresponding to the max value for maxsizedbposs:
    // to keep the balance of (maximize) the number of db profile positions
    const float maxperc = 0.3f;
    size_t maxamountofleftover = (size_t)(leftoversize * maxperc);
    maxamountofleftover = ALIGN_UP( maxamountofleftover, cszalnment );
    maxdbposs = leftoversize / qsspace;
    maxdbposs = ALIGN_DOWN( maxdbposs, CUL2CLINESIZE );
    maxndbpros = maxdbposs / CLOptions::GetDEV_EXPCT_DBPROLEN();
    maxsizedbposs = PMBatchProData::GetPMDataSizeUBTotal( maxdbposs );
    maxsizedbposs = SLC_MIN( maxsizedbposs, maxamountofleftover );
    maxsizedbposs = ALIGN_UP( maxsizedbposs, cszalnment );

    szdpdiag = GetSizeOfDPDiagScores( maxdbposs );
    szdpbottom = GetSizeOfDPBottomScores( maxdbposs );
    szdpmaxcoords = GetSizeOfDPMaxCoords( maxdbposs );
    szdpbtckdat = GetSizeOfDPBackTckData( nqyposs, maxdbposs );
    //{{phase-2 memory requirements
    maxdbposspass2 = (size_t)(maxdbposs * pass2memperc);
    maxdbposspass2 = ALIGN_DOWN( maxdbposspass2, CUL2CLINESIZE );
    maxndbprospass2 = maxdbposspass2 / CLOptions::GetDEV_EXPCT_DBPROLEN();
    szdp2diag = GetSizeOfDP2DiagScores( maxdbposspass2 );
    szdp2bottom = GetSizeOfDP2BottomScores( maxdbposspass2 );
    szdp2maxcoords = GetSizeOfDP2MaxCoords( maxdbposspass2 );
    szdp2btckdat = GetSizeOfDP2BackTckData( nqyposs, maxdbposspass2 );
    szdp2fwdmtx = GetSizeOfDP2FwdMtx( nqyposs, maxdbposspass2 );
    szss2data = GetSizeOfSS2Data( maxdbposspass2, maxndbprospass2 );
    szdp2alndata = GetSizeOfDP2AlnData( maxdbposspass2, maxndbprospass2 );
    szdp2alns = GetSizeOfDP2Alns( nqyposs, maxdbposspass2, maxndbprospass2 );
    sum = szdp2diag + szdp2bottom + szdp2maxcoords + szdp2btckdat + 
          szdp2fwdmtx + szss2data + szdp2alndata + szdp2alns;
    //}}
    szovlpos = szdpdiag + szdpbottom  + szdpmaxcoords + szdpbtckdat;
    if( szovlpos < sum )
        szovlpos = sum;

    if( maxdbposs < 1 || maxsizedbposs < 1 || 
        leftoversize <= maxsizedbposs + szovlpos)
    {
        char msgbuf[KBYTE];
        sprintf( msgbuf, "Insufficient amount of leftover device memory: %zu "
                    "(maxszdbpos= %zu, szdpdiag= %zu, szdpbottom= %zu, "
                    "szdpmaxcoords= %zu, szdpbtckdat= %zu; ndbpos= %zu; "
                    "szovlpos= %zu).", 
                 leftoversize, maxsizedbposs, szdpdiag, szdpbottom, 
                 szdpmaxcoords, szdpbtckdat, maxdbposs, szovlpos );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    //adjust n. of positions to take into account max size of profiles
    maxdbposs = (leftoversize - maxsizedbposs - szovlpos) / qsspace;
    maxdbposs = ALIGN_DOWN( maxdbposs, CUL2CLINESIZE );
    maxndbpros = maxdbposs / CLOptions::GetDEV_EXPCT_DBPROLEN();
    maxsizedbposs = PMBatchProData::GetPMDataSizeUBTotal( maxdbposs );
    maxsizedbposs = SLC_MIN( maxsizedbposs, maxamountofleftover );
    maxsizedbposs = ALIGN_UP( maxsizedbposs, cszalnment );

    szdpdiag = GetSizeOfDPDiagScores( maxdbposs );
    szdpbottom = GetSizeOfDPBottomScores( maxdbposs );
    szdpmaxcoords = GetSizeOfDPMaxCoords( maxdbposs );
    szdpbtckdat = GetSizeOfDPBackTckData( nqyposs, maxdbposs );
    //{{phase-2 memory requirements
    maxdbposspass2 = (size_t)(maxdbposs * pass2memperc);
    maxdbposspass2 = ALIGN_DOWN( maxdbposspass2, CUL2CLINESIZE );
    maxndbprospass2 = maxdbposspass2 / CLOptions::GetDEV_EXPCT_DBPROLEN();
    szdp2diag = GetSizeOfDP2DiagScores( maxdbposspass2 );
    szdp2bottom = GetSizeOfDP2BottomScores( maxdbposspass2 );
    szdp2maxcoords = GetSizeOfDP2MaxCoords( maxdbposspass2 );
    szdp2btckdat = GetSizeOfDP2BackTckData( nqyposs, maxdbposspass2 );
    szdp2fwdmtx = GetSizeOfDP2FwdMtx( nqyposs, maxdbposspass2 );
    szss2data = GetSizeOfSS2Data( maxdbposspass2, maxndbprospass2 );
    szdp2alndata = GetSizeOfDP2AlnData( maxdbposspass2, maxndbprospass2 );
    szdp2alns = GetSizeOfDP2Alns( nqyposs, maxdbposspass2, maxndbprospass2 );
    sum = szdp2diag + szdp2bottom + szdp2maxcoords + szdp2btckdat + 
          szdp2fwdmtx + szss2data + szdp2alndata + szdp2alns;
    //}}
    szovlpos = szdpdiag + szdpbottom  + szdpmaxcoords + szdpbtckdat;
    if( szovlpos < sum )
        szovlpos = sum;

    if( maxdbposs < 1 || maxsizedbposs < 1 )
    {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of leftover device memory: %zu "
                    "(maxszdbpos= %zu, ndbpos= %zu, szovlpos= %zu).", 
                 leftoversize, maxsizedbposs, maxdbposs, szovlpos );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    //adjust further to efficiently utilize memory
    size_t maxdbposs1;
    for( int i = 0; i < 4; i++ ) {
        MYMSGBEGl(5)
            char msgbuf[KBYTE];
            mystring strbuf = preamb + "Adjustment of total size for db profiles: ";
            sprintf(msgbuf, "it %d: %zu; maxdbpos= %zu; szdpdiag= %zu, szdpbottom= %zu, "
                            "szdpmaxcoords= %zu, szdpbtckdat= %zu; szovlpos= %zu; "
                            "leftoversize= %zu; denom= %zu", i, 
                    maxsizedbposs, maxdbposs, szdpdiag, szdpbottom, 
                    szdpmaxcoords, szdpbtckdat, szovlpos, leftoversize, qsspace );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),5);
        MYMSGENDl
        if( leftoversize <= maxsizedbposs + szovlpos)
            break;
        maxdbposs1 = (leftoversize - maxsizedbposs - szovlpos) / qsspace;
        if( maxdbposs1 == maxdbposs )
            break;
        maxdbposs = ALIGN_DOWN( maxdbposs1, CUL2CLINESIZE );
        maxndbpros = maxdbposs / CLOptions::GetDEV_EXPCT_DBPROLEN();
        maxsizedbposs = PMBatchProData::GetPMDataSizeUBTotal( maxdbposs1 );
        maxsizedbposs = SLC_MIN( maxsizedbposs, maxamountofleftover );
        maxsizedbposs = ALIGN_UP( maxsizedbposs, cszalnment );
        szdpdiag = GetSizeOfDPDiagScores( maxdbposs );
        szdpbottom = GetSizeOfDPBottomScores( maxdbposs );
        szdpmaxcoords = GetSizeOfDPMaxCoords( maxdbposs );
        szdpbtckdat = GetSizeOfDPBackTckData( nqyposs, maxdbposs );
        //{{phase-2 memory requirements
        maxdbposspass2 = (size_t)(maxdbposs * pass2memperc);
        maxdbposspass2 = ALIGN_DOWN( maxdbposspass2, CUL2CLINESIZE );
        maxndbprospass2 = maxdbposspass2 / CLOptions::GetDEV_EXPCT_DBPROLEN();
        szdp2diag = GetSizeOfDP2DiagScores( maxdbposspass2 );
        szdp2bottom = GetSizeOfDP2BottomScores( maxdbposspass2 );
        szdp2maxcoords = GetSizeOfDP2MaxCoords( maxdbposspass2 );
        szdp2btckdat = GetSizeOfDP2BackTckData( nqyposs, maxdbposspass2 );
        szdp2fwdmtx = GetSizeOfDP2FwdMtx( nqyposs, maxdbposspass2 );
        szss2data = GetSizeOfSS2Data( maxdbposspass2, maxndbprospass2 );
        szdp2alndata = GetSizeOfDP2AlnData( maxdbposspass2, maxndbprospass2 );
        szdp2alns = GetSizeOfDP2Alns( nqyposs, maxdbposspass2, maxndbprospass2 );
        sum = szdp2diag + szdp2bottom + szdp2maxcoords + szdp2btckdat + 
            szdp2fwdmtx + szss2data + szdp2alndata + szdp2alns;
        //}}
        szovlpos = szdpdiag + szdpbottom  + szdpmaxcoords + szdpbtckdat;
        if( szovlpos < sum )
            szovlpos = sum;
    }

// maxsizedbposs = PMBatchProData::GetPMDataSizeLBTotal( maxdbposs );//leaves some space unexploited 
// maxsizedbposs = ALIGN_UP( maxsizedbposs, cszalnment );

    MYMSGBEGl(3)
        char msgbuf[KBYTE];
        mystring strbuf = preamb + "Total size for db profiles ";
        sprintf(msgbuf, "(max %.1f of leftover): %zu; maxdbpos= %zu", 
                maxperc, maxsizedbposs, maxdbposs );
        strbuf += msgbuf;
        sprintf(msgbuf, " (dev alloc= %zu; dev leftover= %zu; "
                "szdpdiag= %zu, szdpbottom= %zu, "
                "szdpmaxcoords= %zu, szdpbtckdat= %zu; szovlpos= %zu)", 
                deviceallocsize_, leftoversize, szdpdiag, szdpbottom, 
                szdpmaxcoords, szdpbtckdat, szovlpos );
        strbuf += msgbuf;
        MYMSG(strbuf.c_str(),3);
    MYMSGENDl

    //size of score matrix
    const size_t szsmatrix = nqyposs * maxdbposs * sizeof(CUBSM_TYPE);

    SetCurrentMaxDbPos(maxdbposs);
    SetCurrentMaxNDbPros(maxndbpros);
    //{{phase 2
    SetCurrentMaxDbPosPass2( maxdbposspass2 );
    SetCurrentMaxNDbProsPass2( maxndbprospass2 );
    //}}

    sz_heapsections_[ddsEndOfDbChunk] = sz_heapsections_[ddsEndOfCached] + maxsizedbposs;
    sz_heapsections_[ddsEndOfOrgScores] = sz_heapsections_[ddsEndOfDbChunk] + szsmatrix;
    sz_heapsections_[ddsEndOfModScores] = sz_heapsections_[ddsEndOfOrgScores];
    if( GetModScoreMatrixInUse())
        sz_heapsections_[ddsEndOfModScores] += szsmatrix;
    //phase-1 sections
    sz_heapsections_[ddsEndOfDPDiagScores] = sz_heapsections_[ddsEndOfModScores] + szdpdiag;
    sz_heapsections_[ddsEndOfDPBottomScores] = sz_heapsections_[ddsEndOfDPDiagScores] + szdpbottom;
    sz_heapsections_[ddsEndOfDPMaxCoords] = sz_heapsections_[ddsEndOfDPBottomScores] + szdpmaxcoords;
    sz_heapsections_[ddsEndOfDPBackTckData] = sz_heapsections_[ddsEndOfDPMaxCoords] + szdpbtckdat;
    //phase-2 sections: reused previous occupations
    sz_heapsections_[ddsEndOfDP2DiagScores] = sz_heapsections_[ddsEndOfModScores] + szdp2diag;
    sz_heapsections_[ddsEndOfDP2BottomScores] = sz_heapsections_[ddsEndOfDP2DiagScores] + szdp2bottom;
    sz_heapsections_[ddsEndOfDP2MaxCoords] = sz_heapsections_[ddsEndOfDP2BottomScores] + szdp2maxcoords;
    sz_heapsections_[ddsEndOfDP2BackTckData] = sz_heapsections_[ddsEndOfDP2MaxCoords] + szdp2btckdat;
    sz_heapsections_[ddsEndOfDP2FwdMtx] = sz_heapsections_[ddsEndOfDP2BackTckData] + szdp2fwdmtx;
    sz_heapsections_[ddsEndOfSS2Data] = sz_heapsections_[ddsEndOfDP2FwdMtx] + szss2data;
    sz_heapsections_[ddsEndOfDP2AlnData] = sz_heapsections_[ddsEndOfSS2Data] + szdp2alndata;
    sz_heapsections_[ddsEndOfDP2Alns] = sz_heapsections_[ddsEndOfDP2AlnData] + szdp2alns;

    MsgAddressTable( preamb, 3 );

    return maxsizedbposs;
}

// -------------------------------------------------------------------------
// ComputeScoreMatrix: compute score matrices for part of query and database 
// profiles
//
void CuBatchProcessing::ProcessScoreMatrix(
    int qrysernr,
    std::unique_ptr<PMBatchProData> bdbC,
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend,
    char** bdbCpmbeg,
    char** bdbCpmend )
{
    MYMSG( "CuBatchProcessing::ComputeScoreMatrix", 6 );
    BatchProcessScoreMatrixDevice(
        qrysernr,
        std::move(bdbC),
        querypmbeg, querypmend,
        bdb1pmbeg, bdb1pmend,
        bdbCpmbeg, bdbCpmend );
}





// =========================================================================
// Helper code associated with device 
// 


//device addresses of vectors representing profile model data:
__constant__ char* dc_pm2dvfields_[fsg_sztotflds3];
//ANN weights for significance etimation:
__constant__ float dc_cuss_NN_weights_[sz_total_dc_cuss_NNwghts];


// -------------------------------------------------------------------------
// DestroyTextureObject: destroy texture object
inline
void CuBatchProcessing::DestroyTextureObject( cudaTextureObject_t& texObj )
{
    MYMSG( "CuBatchProcessing::DestroyTextureObject", 6 );
    if( texObj ) {
        MYCUDACHECK( cudaDestroyTextureObject( texObj ));
        MYCUDACHECKLAST;
        texObj = 0;
    }
}

// -------------------------------------------------------------------------
// FreeDevicePtr: free device pointer
inline
void CuBatchProcessing::FreeDevicePtr( char*& d_ptr )
{
    MYMSG( "CuBatchProcessing::FreeDevicePtr", 6 );
    if( d_ptr ) {
        MYCUDACHECK( cudaFree( d_ptr ));
        d_ptr = NULL;
        MYCUDACHECKLAST;
    }
}

// -------------------------------------------------------------------------
// CacheSSENNWeightsDevice: transfer NN weights for significance 
// estimation to device CMEM
// 
void CuBatchProcessing::CacheSSENNWeightsDevice()
{
    MYMSG( "CuBatchProcessing::CacheSSENNWeightsDevice", 4 );
    MYCUDACHECK( cudaMemcpyToSymbol( dc_cuss_NN_weights_, SSE_NN_WEIGHTS_, 
                        sz_total_dc_cuss_NNwghts*sizeof(float)));
    MYCUDACHECKLAST;
}

// -------------------------------------------------------------------------
// CacheSSSScoresDevice: cache a score table of secondary structure 
// predictions to a device
// 
void CuBatchProcessing::CacheSSSScoresDevice( const SSSScores& sssscores )
{
    MYMSG( "CuBatchProcessing::CacheSSSScoresDevice", 4 );
    const mystring preamb = "CuBatchProcessing::CacheSSSScoresDevice: ";

    SerializedScoresCtor<CUBSM_TYPE> ssss(
        sssscores.GetScores(), 
        sssscores.GetNoPrbLvs(), sssscores.GetNoTables(), sssscores.GetCardinality(), 
        sssscores.GetPrbLevels(), sssscores.GetLevels()
    );

    //size to align the section of SSS scores:
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = ssss.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[ddsEndOfPadding]) {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: %zu.", deviceallocsize_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofpad = d_heap_ + sz_heapsections_[ddsEndOfPadding];

    MYCUDACHECK( cudaMemcpy( d_eofpad, ssss.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    sz_heapsections_[ddsEndOfSSSTable] = sz_heapsections_[ddsEndOfPadding] + szscores;

    ssssattr_.ntotents_ = ssss.GetSizeAlloc() / sizeof(CUBSM_TYPE);
    ssssattr_.szalloc_ = ssss.GetSizeAlloc();
    ssssattr_.nplvs_ = ssss.GetNProbLevs();
    ssssattr_.nenos_ = ssss.GetNENOLevs();
    ssssattr_.card_ = ssss.GetCardinality();
    ssssattr_.ntbls_ = ssss.GetNTables();
    ssssattr_.nelems_ = ssss.GetNEntries();
}

// -------------------------------------------------------------------------
// CacheCVS2ScoresDevice: cache a map between two types of scores to a 
// device
// 
void CuBatchProcessing::CacheCVS2ScoresDevice( const CVS2Scores& cvs2scores )
{
    MYMSG( "CuBatchProcessing::CacheCVS2ScoresDevice", 4 );
    const mystring preamb = "CuBatchProcessing::CacheCVS2ScoresDevice: ";

    SerializedCVS2ScoresCtor<CUBSM_TYPE> cvs2s( cvs2scores );

    //size to align the section of SSS scores:
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = cvs2s.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[ddsEndOfSSSTable])
    {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: "
                "%zu (ssss+cvs2s, %zu).", deviceallocsize_,
                szscores + sz_heapsections_[ddsEndOfSSSTable]);
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofssssdata = d_heap_ + sz_heapsections_[ddsEndOfSSSTable];

    MYCUDACHECK( cudaMemcpy( d_eofssssdata, cvs2s.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    sz_heapsections_[ddsEndOfCVS2SMap] = sz_heapsections_[ddsEndOfSSSTable] + szscores;

    int card, shft;
    //boundary keys and values from the first map:
    CUBSM_TYPE key1first, value1first, key1last, value1last;
    cvs2s.Get1stCardShiftScale( &card, &shft, NULL, 
                                &key1first, &value1first, &key1last, &value1last );

    cvs2sattr_.ntotents_ = cvs2s.GetSizeAlloc() / sizeof(CUBSM_TYPE);
    cvs2sattr_.szalloc_ = cvs2s.GetSizeAlloc();
    cvs2sattr_.nenos_ = cvs2s.GetNENOLevs();
    cvs2sattr_.ntbls_ = cvs2s.GetNTables();
    cvs2sattr_.card_ = card;
    cvs2sattr_.shft_ = shft;
    cvs2sattr_.key1first_ = key1first;
    cvs2sattr_.value1first_ = value1first;
    cvs2sattr_.key1last_ = key1last;
    cvs2sattr_.value1last_ = value1last;
    cvs2sattr_.CVS_loKAPPA0_ = CVS.loKAPPA0;
    cvs2sattr_.CVS_PowerNU0_ = CVS.PowerNU0;
    cvs2sattr_.CVS_CTERM_ = CVS.CTERM;
}

// -------------------------------------------------------------------------
// CacheHDPScoresDevice: cache score tables of HDP cluster membership 
// predictions to a device
// 
void CuBatchProcessing::CacheHDPScoresDevice( const HDPscores& hdpscores )
{
    MYMSG( "CuBatchProcessing::CacheHDPScoresDevice", 4 );
    const mystring preamb = "CuBatchProcessing::CacheHDPScoresDevice: ";

    SerializedScoresCtor<CUBSM_TYPE> hdp1s(
        hdpscores.GetScores(),
        hdpscores.GetNoPrbLvs(), hdpscores.GetNoTables(), hdpscores.GetCardinality(),
        hdpscores.GetPrbLevels(), hdpscores.GetLevels()
    );

    //size to align the HDP scores in texture memory
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = hdp1s.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[ddsEndOfCVS2SMap])
    {
        char msgbuf[BUF_MAX];
//         sprintf( msgbuf, "Size of HDP scores exceeds allocated device memory size: %lu.", 
//                 deviceallocsize_ );
        sprintf( msgbuf, "Insufficient amount of allocated device memory: "
                "%zu (ssss+cvs2s+hdp1s, %zu).", deviceallocsize_,
                szscores + sz_heapsections_[ddsEndOfCVS2SMap]);
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofcvs2sdata = d_heap_ + sz_heapsections_[ddsEndOfCVS2SMap];

    MYCUDACHECK( cudaMemcpy( d_eofcvs2sdata, hdp1s.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    sz_heapsections_[ddsEndOfHDPscores] = sz_heapsections_[ddsEndOfCVS2SMap] + szscores;

    float eth1 = 0.0f, eth2 = 0.0f;
    hdp1s.GetEth12( &eth1, &eth2 );

    hdpsattr_.ntotents_ = hdp1s.GetSizeAlloc() / sizeof(CUBSM_TYPE);
    hdpsattr_.szalloc_ = hdp1s.GetSizeAlloc();
    hdpsattr_.nplvs_ = hdp1s.GetNProbLevs();
    hdpsattr_.nenos_ = hdp1s.GetNENOLevs();
    hdpsattr_.card_ = hdp1s.GetCardinality();
    hdpsattr_.ntbls_ = hdp1s.GetNTables();
    hdpsattr_.nelems_ = hdp1s.GetNEntries();
    hdpsattr_.eth1_ = eth1;
    hdpsattr_.eth2_ = eth2;

    if( (size_t)deviceProp_.maxTexture1DLinear_ < szscores / sizeof(CUBSM_TYPE)) {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Number of entries in the HDP scores buffer exceeds "
                "the maximum limit: %zu > %d.",
                szscores/sizeof(CUBSM_TYPE), deviceProp_.maxTexture1DLinear_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    DestroyTextureObject( hdp1sTexObj_ );
    //create texture object
    //first, specify texture resource descriptor
    struct cudaResourceDesc resDesc;
    memset( &resDesc, 0, sizeof(resDesc));
    //create channel descriptor
    //cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = d_eofcvs2sdata;
    resDesc.res.linear.sizeInBytes = szscores;
    resDesc.res.linear.desc = channelDesc;
    //next, specify texture descriptor (object parameters)
    struct cudaTextureDesc texDesc;
    memset( &texDesc, 0, sizeof(texDesc));
    //many fields are ignored when cudaResourceDesc::resType is cudaResourceTypeLinear
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 0;
    //finaly, create texture object
    MYCUDACHECK( cudaCreateTextureObject( &hdp1sTexObj_, &resDesc, &texDesc, NULL ));
    MYCUDACHECKLAST;
}

// -------------------------------------------------------------------------
// PackPMDataForDevice: pack profile model data into one buffer for 
// sending it to a device
// szpm2dvfields, an array of [fsg_sztotflds2+1] elements, sizes of each field 
// used to calculate addresses passed to a device;
// tmpbuf, temporary buffer for packing data;
// NOTE: the memory for this buffer is allocated in this function and the 
// control over it is passed to the calling function;
// szalloc, total size in bytes allocated for tmpbuf
//
void CuBatchProcessing::PackPMDataForDevice( 
    size_t szpm2dvfields[], char*& tmpbuf, size_t& szalloc,
    char** querypmbeg,//beginning address of queries
    char** querypmend,//terminal address of queries
    char** bdb1pmbeg,//beginning address of cached database profiles
    char** bdb1pmend,//terminal address of cached database profiles
    char** bdbCpmbeg,//beginning address of new profiles read from the database 
    char** bdbCpmend )//terminal address of new profiles read from the database 
{
    MYMSG( "CuBatchProcessing::PackPMDataForDevice", 4 );
    const mystring preamb = "CuBatchProcessing::PackPMDataForDevice: ";
    myruntime_error mre;

    const size_t szcalign = CUL2CLINESIZE;//alignment size wrt cache line size

    size_t i = 1, sz;
    size_t nflds = ( querypmbeg && querypmend )? fsg_sztotflds2: (size_t)pmv2DTotFlds;
    size_t dbfldsbeg = nflds - pmv2DTotFlds;//starting field index for db profiles

    memset( szpm2dvfields, 0, (nflds+1) * sizeof(size_t));

    if( querypmbeg && querypmend ) {
        for( i = 1; i <= pmv2DTotFlds; i++ ) {
#ifdef __DEBUG__
        if( querypmend[i-1] < querypmbeg[i-1] )
            throw MYRUNTIME_ERROR( preamb + "Invalid query pointers.");
#endif
            szpm2dvfields[i] = szpm2dvfields[i-1] + (size_t)(querypmend[i-1]-querypmbeg[i-1]);
            szpm2dvfields[i] = ALIGN_UP( szpm2dvfields[i], szcalign );
        }
    }

    for( ; i <= nflds; i++ ) {
        szpm2dvfields[i] = szpm2dvfields[i-1];
        if( bdb1pmend && bdb1pmbeg ) {
#ifdef __DEBUG__
            if( bdb1pmend[i-dbfldsbeg-1] < bdb1pmbeg[i-dbfldsbeg-1] )
                throw MYRUNTIME_ERROR( preamb + "Invalid database profile pointers.");
#endif
            szpm2dvfields[i] += 
                (size_t)(bdb1pmend[i-dbfldsbeg-1]-bdb1pmbeg[i-dbfldsbeg-1]);
        }
        if( bdbCpmend && bdbCpmbeg ) {
#ifdef __DEBUG__
            if( bdbCpmend[i-dbfldsbeg-1] < bdbCpmbeg[i-dbfldsbeg-1] )
                throw MYRUNTIME_ERROR( preamb + "Invalid additional database profile pointers.");
#endif
            szpm2dvfields[i] += 
                (size_t)(bdbCpmend[i-dbfldsbeg-1]-bdbCpmbeg[i-dbfldsbeg-1]);
        }
        //szpm2dvfields[i-1] has been aligned anyway
        szpm2dvfields[i] = ALIGN_UP( szpm2dvfields[i], szcalign );
    }

#ifdef __DEBUG__
    if( szpm2dvfields[pmv2DTotFlds] < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profile positions.");
#endif

    szalloc = szpm2dvfields[nflds];
    tmpbuf = NULL;

    if( szalloc < 1 )
        throw MYRUNTIME_ERROR( preamb + "No memory allocation.");

    tmpbuf = (char*)malloc(szalloc * sizeof(char));
    if( !tmpbuf )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

    try {
        //pack data into one buffer
        i = 0;
        if( querypmbeg && querypmend ) {
            for( ; i < pmv2DTotFlds; i++ ) {
                memcpy( tmpbuf+szpm2dvfields[i], querypmbeg[i], 
                    (size_t)(querypmend[i]-querypmbeg[i]));
            }
        }

        for( ; i < nflds; i++ ) {
            sz = 0;
            if( bdb1pmend && bdb1pmbeg ) {
                sz = (size_t)(bdb1pmend[i-dbfldsbeg]-bdb1pmbeg[i-dbfldsbeg]);
                memcpy( tmpbuf+szpm2dvfields[i], bdb1pmbeg[i-dbfldsbeg], sz );
            }
            if( bdbCpmend && bdbCpmbeg )
                memcpy( tmpbuf+szpm2dvfields[i]+sz, bdbCpmbeg[i-dbfldsbeg], 
                    (size_t)(bdbCpmend[i-dbfldsbeg]-bdbCpmbeg[i-dbfldsbeg])
                );
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// CopyPMDataToDevice: pack and copy profile model data to a device;
// querypmbeg, host pointer for query beginning addreses;
// querypmend, host pointer for query terminal addreses;
// bdb1pmbeg, host pointer for beginning addreses of cached db profiles;
// bdb1pmend, host pointer for terminal addreses of cached db profiles;
// bdbCpmbeg, host pointer for beginning addreses of read db profiles;
// bdbCpmend, host pointer for terminal addreses of read db profiles;
// dev_pckdpm, device pointer to (the end  of) packed data; (changes on exit);
// szmaxsize, maximum amount of memory allowed to be occupied for the data passed by arguments;
// cmndx, index of constant memory array where field addresses reside;
// NOTE: memory is allocated for device pointer dev_pckdpm 
//
void CuBatchProcessing::CopyPMDataToDevice(
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend,
    char** bdbCpmbeg,
    char** bdbCpmend,
    char*& dev_pckdpm,
    size_t szmaxsize,
    size_t cmndx )
{
    MYMSG( "CuBatchProcessing::CopyPMDataToDevice", 4 );
    const mystring preamb = "CuBatchProcessing::CopyPMDataToDevice: ";
    myruntime_error mre;

    //sizes of each field used to calculate addresses passed to a device
    size_t szpm2dvfields[fsg_sztotflds2+1];
    size_t szalloc = 0UL, i;
    //temporary buffer for packing data; this works much faster than 
    // conducting multiple transfers to a device
    char* tmpbuf = NULL;
    //device pointers to the fields of profile model data packed to one buffer
    char* d_pckdpmdatflds[fsg_sztotflds2];
    size_t nflds = ( querypmbeg && querypmend )? fsg_sztotflds2: (size_t)pmv2DTotFlds;

    try {
        PackPMDataForDevice( 
            szpm2dvfields, tmpbuf, szalloc, 
            querypmbeg, querypmend,
            bdb1pmbeg, bdb1pmend,
            bdbCpmbeg, bdbCpmend
        );

        if( szalloc < 1 ) {
            MYMSGBEGl(3)
                mystring strbuf = preamb + "No data transfered!";
                MYMSG(strbuf.c_str(),3);
            MYMSGENDl
            return;
        }

        if( szmaxsize < szalloc ) {
            char msgbuf[BUF_MAX];
            mystring strbuf = preamb + "Short of allocated device memory: ";
            sprintf( msgbuf, "%zuB requested vs %zuB allowed.", szalloc, szmaxsize );
            throw MYRUNTIME_ERROR( strbuf );
        }

//         FreeDevicePtr( dev_pckdpm );

//         MYCUDACHECK( cudaMalloc((void**)&dev_pckdpm, szalloc ));
//         MYCUDACHECKLAST;

        MYCUDACHECK( cudaMemcpy( dev_pckdpm, tmpbuf, szalloc, cudaMemcpyHostToDevice ));
        MYCUDACHECKLAST;

        for( i = 0; i < nflds; i++ )
            d_pckdpmdatflds[i] = dev_pckdpm + szpm2dvfields[i];

        dev_pckdpm += szalloc;

        //put device addresses pointing to locations of the same buffer to constant memory
        switch( cmndx ) {
            case ndx_dc_pm2dvfields_: 
                MYCUDACHECK( cudaMemcpyToSymbol( dc_pm2dvfields_, d_pckdpmdatflds, nflds*sizeof(char*)));
                break;
            case ndx_dc_newpm2dvfields_:
                MYCUDACHECK( cudaMemcpyToSymbol( dc_pm2dvfields_, d_pckdpmdatflds, nflds*sizeof(char*), 
                                        fsg_sztotflds2*sizeof(char*)/*offset*/));
                break;
            default: throw MYRUNTIME_ERROR( preamb + "Invalid constant memory array index.");
        }
        MYCUDACHECKLAST;

        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            mystring strbuf = preamb + "Data transfered:";
            if( querypmbeg && querypmend ) {
                sprintf(msgbuf, "  query size= %zu", szpm2dvfields[pmv2DTotFlds]);
                strbuf += msgbuf;
            }
            sprintf(msgbuf, "  total data size= %zu", szalloc );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( tmpbuf ) {
        free(tmpbuf);
        tmpbuf = NULL;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// CacheDataDevice: pack and transfer cached query and db profile model 
// data to a device
// 
inline
void CuBatchProcessing::CacheDataDevice(
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend )
{
    MYMSG( "CuBatchProcessing::CacheDataDevice", 4 );
    const mystring preamb = "CuBatchProcessing::CacheDataDevice: ";
    const size_t cszalnment = GetMemAlignment();

    char* d_eofcnstdata = d_heap_ + sz_heapsections_[ddsEndOfConstantData];
    char* d_eofcnstdataplus = d_eofcnstdata;

    if( deviceallocsize_ <= sz_heapsections_[ddsEndOfConstantData])
        throw MYRUNTIME_ERROR( preamb + "Insufficient device memory allocated to cache data.");

    //maxsize for data; TODO: to be changed
    size_t szmaxallowed = deviceallocsize_ - sz_heapsections_[ddsEndOfConstantData];

    CopyPMDataToDevice( 
        querypmbeg, querypmend, bdb1pmbeg, bdb1pmend,
        NULL/*bdbCpmbeg*/, NULL/*bdbCpmend*/,
        d_eofcnstdataplus, 
        szmaxallowed,
        ndx_dc_pm2dvfields_ );
    //
    size_t szdata = (size_t)(d_eofcnstdataplus-d_eofcnstdata);
    szdata = ALIGN_UP( szdata, cszalnment );
    sz_heapsections_[ddsEndOfCached] = sz_heapsections_[ddsEndOfConstantData] + szdata;
}

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
    std::unique_ptr<PMBatchProData> bdbC,
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend,
    char** bdbCpmbeg,
    char** bdbCpmend )
{
    MYMSG( "CuBatchProcessing::BatchProcessScoreMatrixDevice", 4 );
    const mystring preamb = "CuBatchProcessing::BatchProcessScoreMatrixDevice: ";
    myruntime_error mre;

#ifdef __DEBUG__
    if( !querypmbeg || !querypmend || !h_querypmbeg_)
        throw MYRUNTIME_ERROR( preamb + "Null query addresses.");
    if( d_heap_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null device heap.");
    if( cbsm_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch score matrix object.");
    if( cbdp_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch DP object.");
    if( cbss_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch SS object.");
    if( cbmapdp_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null batch MAP DP object.");
#endif

    unsigned int attrpassed[nTPassedPros] = {0,0,0};
    //attributes for the allowed number of passed profiles:
    unsigned int attrpassedallwd[nTPassedPros] = {0,0,0};

    size_t ndb1pros = 0;//total number of profiles in the first buffer
    size_t ndbCpros = 0;//total number of profiles in the complementary buffer
    size_t dbpro1len = 0;//the length of the first profile, which is the largest length

    size_t querprosOmtd = PMBatchProData::GetNoProsFromTo( h_querypmbeg_, querypmbeg );
    size_t ndb1prosOmtd = 0;//number of profiles omitted in the first buffer
    size_t ndbCprosOmtd = 0;//number of profiles omitted in the complementary buffer

//     uint querposoffset = (uint)PMBatchProData::GetNoProsFromTo( h_querypmbeg_, querypmbeg );
    size_t querposoffset = PMBatchProData::GetNoPositsFromTo( h_querypmbeg_, querypmbeg );
    size_t bdb1posoffset = 0;
    size_t bdbCposoffset = 0;

    const size_t nqyposs = PMBatchProData::GetNoPositsFromTo( querypmbeg, querypmend );
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
        if( !h_bdb1pmbeg_)
            throw MYRUNTIME_ERROR( preamb + "Null cached data address.");
#endif
        bdb1posoffset = PMBatchProData::GetNoPositsFromTo( h_bdb1pmbeg_, bdb1pmbeg );
        ndb1prosOmtd = PMBatchProData::GetNoProsFromTo( h_bdb1pmbeg_, bdb1pmbeg );
        //bdb1pmbeg..bdb1pmend come in front of bdbCpmbeg..bdbCpmend, 
        //the first profile length here is hence the largest
        dbpro1len = PMBatchProData::GetPMDataLen1At( bdb1pmbeg );
    }

    size_t ndbxposs = ndb1poss + ndbCposs;

    //padding along the x axis
    const size_t dbxpad = ALIGN_UP( ndbxposs, CUL2CLINESIZE ) - ndbxposs;
    SetCurrentDbxPadding((unsigned int)dbxpad);

#ifdef __DEBUG__
    size_t szsspace = nqyposs * ndbxposs;//search space

    if( ndb1pros + ndbCpros < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles of part of db.");
    if( szsspace < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid size of search space.");
    if( GetCurrentMaxDbPos() < ndb1poss + ndbCposs + dbxpad ) {
        char msgbuf[BUF_MAX];
        sprintf(msgbuf, "Invalid accumulated number of Db profile positions: %zu < %zu+%zu+%zu.",
            GetCurrentMaxDbPos(), ndb1poss, ndbCposs, dbxpad );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }
#endif

    cudaStream_t streamproc;
    MYCUDACHECK( cudaStreamCreate(&streamproc));

    try {
        if( bdbCpmbeg && bdbCpmend ) {
            //on exit of the call the address changes to point to the end of writtent data
            char* d_eofqueryandbdb1 = d_heap_ + sz_heapsections_[ddsEndOfCached];
            if( sz_heapsections_[ddsEndOfDbChunk] <= sz_heapsections_[ddsEndOfCached])
                throw MYRUNTIME_ERROR( preamb + "Unallocated memory for chunked db data.");
            CopyPMDataToDevice(
                NULL/*querypmbeg*/, NULL/*querypmend*/, 
                NULL/*bdb1pmbeg*/, NULL/*bdb1pmend*/,
                bdbCpmbeg, bdbCpmend,
                d_eofqueryandbdb1,
                sz_heapsections_[ddsEndOfDbChunk] - sz_heapsections_[ddsEndOfCached],//maxsize for data
                ndx_dc_newpm2dvfields_ );
        }

        cbsm_->ComputeScoreMatrixDevice(
            streamproc,
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfPadding])/*sssscores*/,
            ssssattr_,
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfSSSTable])/*cvs2scores*/,
            cvs2sattr_,
            hdp1sTexObj_,
            hdpsattr_,
            GetHDP1ScoresInUse(),
            nqyposs,
            ndb1poss,
            ndbCposs,
            dbxpad,
            querposoffset,
            bdb1posoffset,
            bdbCposoffset,
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDbChunk])/*outscores*/,
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfOrgScores])/*outmodscores*/
        );

        if( !GetMAPAlnInUse()) {
            //make sure results are not being copied before overwriting them:
            MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
            MYCUDACHECKLAST;
        }

        cbdp_->PerformCompleteDynProgDevice(
            streamproc,
            GetScoreThreshold(),
            !GetMAPAlnInUse(),
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
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDbChunk])/*scores [in]*/,
            (CUBDP_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfModScores])/*tmpdpdiagbuffers [out]*/,
            (CUBDP_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDPDiagScores])/*tmpdpbotbuffer [out]*/,
            (unsigned int*)(d_heap_ + sz_heapsections_[ddsEndOfDPBottomScores])/*maxcoordsbuf [out]*/,
            (char*)(d_heap_ + sz_heapsections_[ddsEndOfDPMaxCoords])/*btckdata [out]*/,
            //[out:]
            attrpassed
        );

        if( querylength_ && dblength_ && ndbsequences_ ) {
            cbsm_->ComputeLengthAdjustment( querylength_, dblength_, ndbsequences_ );
            ResetDbDetails();
        }

        MYCUDACHECK( cudaStreamSynchronize(streamproc));
        MYCUDACHECKLAST;

        //make sure results have been copied before overwriting them:
        MYCUDACHECK( cudaStreamSynchronize(streamcopyres_));
        MYCUDACHECKLAST;

        MYCUDACHECK( cudaMemcpyFromSymbol( attrpassed, d_gDPPassedPros, sizeof(d_gDPPassedPros)));
        MYCUDACHECKLAST;

        MYMSGBEGl(1)
            char msgbuf[KBYTE];
            sprintf(msgbuf,"%6c%u profiles passed to the realignment phase",' ',
                    attrpassed[dpppNPassedPros]);
            MYMSG(msgbuf,1);
        MYMSGENDl

        MYMSGBEGl(3)
            char msgbuf[KBYTE];
            mystring strbuf = preamb;
            sprintf(msgbuf,"%sPhase-1 results summary: ",NL);
            strbuf += msgbuf;
            sprintf(msgbuf,
                "# db_pros_passed_to_phase_2= %u (thld= %.1f), total positions= %u, max len= %u",
                attrpassed[dpppNPassedPros], GetScoreThreshold(), 
                attrpassed[dpppNPosits], attrpassed[dpppMaxProLen]);
            strbuf += msgbuf;
            //MYMSG(strbuf.c_str(),3);
            sprintf(msgbuf,
                "%s(phase-2 limit values: # profiles= %zu, # positions= %zu)",
                NL, GetCurrentMaxNDbProsPass2(), GetCurrentMaxDbPosPass2());
            //strbuf = preamb + msgbuf;
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

        for( int i = 0; i < nTPassedPros; i++)
            attrpassedallwd[i] = attrpassed[i];

        if( GetCurrentMaxNDbProsPass2() < attrpassed[dpppNPassedPros] ||
            GetCurrentMaxDbPosPass2() < attrpassed[dpppNPosits]) 
        {
            warning("RESULTS WILL NOT BE COMPLETE. "
                    "Increase the value of option --dev-pass2memp and rerun.");
            if( GetCurrentMaxNDbProsPass2() < attrpassed[dpppNPassedPros])
                attrpassedallwd[dpppNPassedPros] = (unsigned int)GetCurrentMaxNDbProsPass2();
            if( GetCurrentMaxDbPosPass2() < attrpassed[dpppNPosits])
                attrpassedallwd[dpppNPosits] = (unsigned int)GetCurrentMaxDbPosPass2();
        }

        cbss_->CalculateAlnStatisticsDevice(
            streamproc,
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
            (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDbChunk])/*scores [in]*/,
            (CUBDP_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfModScores])/*tmpdpdiagbuffers [in]*/,
            //(CUBDP_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDPDiagScores])/*tmpdpbotbuffer [in]*/,
            (float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2FwdMtx])/*tmpss2datbuffers [out]*/,
            (float*)(d_heap_ + sz_heapsections_[ddsEndOfSS2Data])/*dp2alndatbuffers [out]*/,
            attrpassedallwd/*[in]*/
        );

        //padding along the x axis for db profile data in the phase-2 processing
        const size_t dbxpad2 = ALIGN_UP( attrpassedallwd[dpppNPosits], CUL2CLINESIZE ) - 
                attrpassedallwd[dpppNPosits];
        SetCurrentDbxPaddingPass2((unsigned int)dbxpad2);
        //set db length for alignments, including padding
        size_t cbdbalnlen2 = 
            attrpassedallwd[dpppNPosits] + attrpassedallwd[dpppNPassedPros] * (nqyposs+1);
        cbdbalnlen2 = ALIGN_UP( cbdbalnlen2, CUL2CLINESIZE );
        SetCurrentDbAlnLengthWithPaddingPass2((unsigned int)cbdbalnlen2);

        if( GetMAPAlnInUse())
            cbmapdp_->PerformCompleteMAPDynProgDevice(
                streamproc,
                GetLogEThreshold(),
                GetModScoreMatrixInUse(),
                ndb1pros,
                ndbCpros,
                querprosOmtd,
                ndb1prosOmtd,
                ndbCprosOmtd,
                attrpassed[dpppMaxProLen],
                nqyposs,
                ndb1poss,
                ndbCposs,
                dbxpad,
                querposoffset,
                bdb1posoffset,
                bdbCposoffset,
                (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDbChunk])/*scores [in]*/,
                (CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfOrgScores])/*mod scores [in]*/,
                (float*)(d_heap_ + sz_heapsections_[ddsEndOfModScores])/*tmpdpdiagbuffers [tmp]*/,
                (float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2DiagScores])/*tmpdpbotbuffer [tmp]*/,
//(float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2BottomScores])/*tmpblockprobscales [tmp]*/,
//(float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2DiagBlkProbScales])/*tmpprobscales [tmp]*/,
                (unsigned int*)(d_heap_ + sz_heapsections_[ddsEndOfDP2BottomScores])/*maxcoordsbuf [tmp]*/,
                (char*)(d_heap_ + sz_heapsections_[ddsEndOfDP2MaxCoords])/*btckdata [tmp]*/,
                (float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2BackTckData])/*tmpfwdprobs [tmp]*/,
//(float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2FwdMtx])/*tmpbwdprobs [tmp]*/,
                (float*)(d_heap_ + sz_heapsections_[ddsEndOfDP2FwdMtx])/*tmpss2datbuffers [tmp]*/,
                (float*)(d_heap_ + sz_heapsections_[ddsEndOfSS2Data])/*dp2alndatbuffers [in]*/,
                (char*)(d_heap_ + sz_heapsections_[ddsEndOfDP2AlnData])/*output alignments [out]*/,
                attrpassedallwd,//[in]
                dbxpad2,
                cbdbalnlen2
            );

        MYCUDACHECK( cudaStreamSynchronize(streamproc));
        MYCUDACHECKLAST;

        TransferResultsDevice(
            qrysernr,
            nqyposs, qyeno,
            std::move(bdbC),
            querypmbeg, querypmend,
            bdb1pmbeg, bdb1pmend,
            bdbCpmbeg, bdbCpmend,
            ndb1pros,
            ndbCpros,
            querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
            querposoffset, bdb1posoffset, bdbCposoffset,
            attrpassedallwd[dpppNPosits],
            attrpassedallwd[dpppNPassedPros]);

// size_t szdata = (sz_heapsections_[ddsEndOfDPBackTckData]-sz_heapsections_[ddsEndOfDPMaxCoords])/10;
// char* tmpbuf = (char*)malloc(szdata * sizeof(char));
// if( !tmpbuf ) throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
// MYCUDACHECK( cudaMemcpy( tmpbuf, 
//     (char*)(d_heap_+sz_heapsections_[ddsEndOfDPMaxCoords]), szdata, cudaMemcpyDeviceToHost ));
// MYCUDACHECKLAST;
// tmpbuf[0]=tmpbuf[100]=23;
// free(tmpbuf);


if( nqyposs == /*26*/99999/*394*//*1074*//*967*//*1122*/ )
cbsm_->TESTPrintProProScores1(
    h_bdb1pmbeg_, h_bdb1pmend_,
    querypmbeg, querypmend,
    bdb1pmbeg, bdb1pmend, bdbCpmbeg, bdbCpmend, GetCurrentDbxPadding(),
    sz_heapsections_[ddsEndOfOrgScores] - sz_heapsections_[ddsEndOfDbChunk]
    ,(CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfDbChunk])/*outscores*/
    // ,(CUBSM_TYPE*)(d_heap_ + sz_heapsections_[ddsEndOfOrgScores])/*outmodscores*/
);

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
// bdbCpmbeg, beginning address of new profiles read from the database; 
// bdbCpmend, terminal address of new profiles read from the database;
// NOTE: query and cached db data addresses are relative to the beginning 
// addresses stored in host memebers h_querypmbeg_, h_querypmend_, 
// h_bdb1pmbeg_, and h_bdb1pmend_;
//
void CuBatchProcessing::TransferResultsDevice(
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
    unsigned int npros )
{
    MYMSG( "CuBatchProcessing::TransferResultsDevice", 3 );
    mystring preamb = "CuBatchProcessing::TransferResultsDevice: ";
    myruntime_error mre;

    size_t szresults = //size allocated for results
        sz_heapsections_[ddsEndOfDP2Alns] - 
        sz_heapsections_[ddsEndOfSS2Data];

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
                    h_results_, d_heap_ + sz_heapsections_[ddsEndOfSS2Data],
                    szresults, cudaMemcpyDeviceToHost, streamcopyres_));
                MYCUDACHECKLAST;
            }
            else {
                MYCUDACHECK( cudaMemcpy(
                    h_results_, d_heap_ + sz_heapsections_[ddsEndOfSS2Data],
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
            std::move(bdbC),
            querypmbeg, querypmend, 
            bdb1pmbeg, bdb1pmend,
            bdbCpmbeg, bdbCpmend,
            ndb1pros,
            ndbCpros,
            querprosOmtd, ndb1prosOmtd, ndbCprosOmtd,
            querposoffset, bdb1posoffset, bdbCposoffset,
            nposits,
            npros,
            h_results_,
            sz_heapsections_[ddsEndOfDP2AlnData]-sz_heapsections_[ddsEndOfSS2Data],
            sz_heapsections_[ddsEndOfDP2Alns]-sz_heapsections_[ddsEndOfDP2AlnData],
            GetCurrentDbAlnLengthWithPaddingPass2()
        );
        cbpfin_->Notify(CuBatchProcessingFinalizer::cubpthreadmsgFinalize);

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}
