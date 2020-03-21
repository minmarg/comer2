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

#include "extsp/psl.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/TCTXVECT.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cualn/Devices.h"
#include "libmycu/cupro/SerializedScores.cuh"
#include "libmycu/cupro/SerializedScoresCtor.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresCtor.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "CuBatchProcessing_com.cuh"
#include "CuDeviceMemory.cuh"

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
CuDeviceMemory::CuDeviceMemory(
        DeviceProperties dprop, 
        bool hdp1scoresinuse,
        bool mapalninuse,
        bool Xuninf,
        size_t deviceallocsize,
        int nareas )
:
    hdp1scoresinuse_( hdp1scoresinuse ),
    mapalninuse_( mapalninuse ),
    Xuninf_( Xuninf ),
    deviceallocsize_( deviceallocsize ),
    nareas_(nareas),
    curmaxdbpos_( 0UL ),
    curmaxndbpros_( 0UL ),
    curmaxdbpospass2_( 0UL ),
    curmaxndbprospass2_( 0UL ),
    h_querypmbeg_(NULL),
    h_querypmend_(NULL),
    h_bdb1pmbeg_(NULL),
    h_bdb1pmend_(NULL),
    //
    deviceProp_(dprop),
    hdp1sTexObj_(0),
    d_heap_(NULL),
    sz_heapsections_(NULL)
{
    MYMSG( "CuDeviceMemory::CuDeviceMemory", 4 );

    if( nareas_ < 1 )
        throw MYRUNTIME_ERROR("CuDeviceMemory::CuDeviceMemory: Invalid number of memory areas.");

    MYCUDACHECK( cudaSetDevice( deviceProp_.devid_ ));
    MYCUDACHECKLAST;

    sz_heapsections_ = (size_t(*)[nDevDataSectionsDiv2])
        malloc( nareas_ * nDevDataSectionsDiv2 * sizeof(size_t));

    if( sz_heapsections_ == NULL )
        throw MYRUNTIME_ERROR("CuDeviceMemory::CuDeviceMemory: Not enough memory.");

//     for( int a = 0; a < nareas_; a++ )
//         for( int i = 0; i < nDevDataSectionsDiv2; i++ )
//             sz_heapsections_[a][i] = 0UL;
    memset(sz_heapsections_, 0, nareas_ * nDevDataSectionsDiv2 * sizeof(size_t));

    if( deviceallocsize_ )
    {
        const size_t cszalnment = GetMemAlignment();

        MYCUDACHECK( cudaMalloc((void**)&d_heap_, deviceallocsize_ ));
        MYCUDACHECKLAST;

        size_t heapaligned = ALIGN_UP((size_t)d_heap_, cszalnment );
        for( int a = 0; a < nareas_; a++ )
            sz_heapsections_[a][ddsEndOfPadding] = (size_t)d_heap_ - heapaligned;
    }
}

// -------------------------------------------------------------------------
// default constructor
//
CuDeviceMemory::CuDeviceMemory()
{
    throw MYRUNTIME_ERROR("CuDeviceMemory::CuDeviceMemory: "
                "Default initialization is not allowed.");
}

// -------------------------------------------------------------------------
// destructor
//
CuDeviceMemory::~CuDeviceMemory()
{
    MYMSG( "CuDeviceMemory::~CuDeviceMemory", 4 );
    if(sz_heapsections_) {
        free(sz_heapsections_);
        sz_heapsections_ = NULL;
    }
    DestroyTextureObject( hdp1sTexObj_ );
    FreeDevicePtr( d_heap_ );
}

// -------------------------------------------------------------------------
// CacheCompleteData: transfer complete required data to device
// 
void CuDeviceMemory::CacheCompleteData(
    char** querypmbeg, char** querypmend,
    char** bdb1pmbeg, char** bdb1pmend)
{
    CacheSSENNWeights();
    CacheSSSScores( SSSSCORES );
    CacheCVS2Scores( CVS2SCORES );
    CacheHDPScores( HDPSCORES );
    CacheData(querypmbeg, querypmend,  bdb1pmbeg, bdb1pmend);
}





// -------------------------------------------------------------------------
// CacheSSENNWeights: cache NN weights for significance estimation
// 
void CuDeviceMemory::CacheSSENNWeights()
{
    MYMSG( "CuDeviceMemory::CacheSSENNWeights", 6 );
    if( MOptions::GetSSEMODEL() == 0 )
        return;
    CacheSSENNWeightsDevice();
}

// -------------------------------------------------------------------------
// CacheSSSScores: cache a score table of secondary structure predictions
// 
void CuDeviceMemory::CacheSSSScores( const SSSScores& sssscores )
{
    MYMSG( "CuDeviceMemory::CacheSSSScores", 6 );
    if( MOptions::GetSSSWGT() <= 0.0f )
        return;
    CacheSSSScoresDevice( sssscores );
}

// -------------------------------------------------------------------------
// CacheCVS2Scores: cache a map between cv scores and translated scores
// 
void CuDeviceMemory::CacheCVS2Scores( const CVS2Scores& cvs2scores )
{
    MYMSG( "CuDeviceMemory::CacheCVS2Scores", 6 );
    if( MOptions::GetCVSWGT() <= 0.0f )
        return;
    CacheCVS2ScoresDevice( cvs2scores );
}

// -------------------------------------------------------------------------
// CacheCVS2Scores: cache a map between cv scores and translated scores
// 
void CuDeviceMemory::CacheHDPScores( const HDPscores& hdpscores )
{
    MYMSG( "CuDeviceMemory::CacheHDPScores", 6 );
    if( !GetHDP1ScoresInUse() || MOptions::GetADJWGT() <= 0.0f )
        return;
    CacheHDPScoresDevice( hdpscores );
}


// -------------------------------------------------------------------------
// CacheData: cache query and db profile model data; host adresses are 
// saved for counting of positions later
// 
void CuDeviceMemory::CacheData(
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend )
{
    MYMSG( "CuDeviceMemory::CacheData", 6 );
    CacheDataDevice(
        h_querypmbeg_ = querypmbeg, 
        h_querypmend_ = querypmend, 
        h_bdb1pmbeg_ = bdb1pmbeg, 
        h_bdb1pmend_ = bdb1pmend );
}

// -------------------------------------------------------------------------
// TransferCPMData: transfer (prepare) a new chunk of db profile model data
// 
void CuDeviceMemory::TransferCPMData(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields )
{
    MYMSG( "CuDeviceMemory::TransferCPMData", 6 );
    TransferCPMDataDevice( bdbCpmbeg, bdbCpmend, szCpm2dvfields );
}





// =========================================================================
// CalcMaxDbDataChunkSize: delineate device memory sections given a query 
// length; the boundaries of device memory sections are calculated for each 
// memory area
//
size_t CuDeviceMemory::CalcMaxDbDataChunkSize( size_t nqyposs )
{
    MYMSG( "CuDeviceMemory::CalcMaxDbDataChunkSize", 4 );
    const mystring preamb = "CuDeviceMemory::CalcMaxDbDataChunkSize: ";
    char msgbuf[BUF_MAX];

    if( deviceallocsize_ <= sz_heapsections_[0][ddsEndOfCached]) {
        sprintf( msgbuf, "Insufficient amount of allocated device memory: %zu.", deviceallocsize_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    if( nareas_ < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of memory areas.");

    const size_t gapbtwas = 16UL*ONEM;//gap between areas
    const size_t leftoversize = deviceallocsize_ - sz_heapsections_[0][ddsEndOfCached];
    const size_t areasize = leftoversize / nareas_ - (nareas_-1) * gapbtwas;

    //TODO: currently each memory area allows for db chunk data and 
    // allocations of its size; however, percentage is not large wrt the 
    // total allocated memory;
    //this can be changed in future to abandon the reservations for all 
    // memory areas but the first;
    size_t chunkdatasize = CalcMaxDbDataChunkSizeHelper(nqyposs, areasize);

    MsgAddressTable( 0, preamb, 3 );

    //the first memory area has been configured;
    //derive memory section sizes for the other memory areas
    for(int a = 1; a < nareas_; a++) {
        sz_heapsections_[a][ddsEndOfDbChunk] = sz_heapsections_[0][ddsEndOfDbChunk];
        sz_heapsections_[a][ddsBegOfOrgScores] = sz_heapsections_[a-1][nDevDataSectionsDiv2-1]  +  gapbtwas;
        //sz_heapsections_[a][ddsBegOfOrgScores] = ALIGN_UP(sz_heapsections_[a][ddsBegOfOrgScores],4096);
        //
        for(int s = ddsBegOfOrgScores+1; s < nDevDataSectionsDiv2; s++) {
            sz_heapsections_[a][s] = sz_heapsections_[a][ddsBegOfOrgScores] + 
                (sz_heapsections_[0][s]-sz_heapsections_[0][ddsBegOfOrgScores]);
            //
            if( deviceallocsize_ < sz_heapsections_[a][s] &&
                (GetMAPAlnInUse() || s != ddsEndOfDP2FwdMtx))
            {
                sprintf( msgbuf, "Size of section %d of area %d exceeds "
                    "device memory allocations: %zu > %zu.",
                    s, a, sz_heapsections_[a][s], deviceallocsize_);
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }
        }
        MsgAddressTable( a, preamb, 3 );
        //TODO: once memory configuration has changed, notify 
        // processing children to check host memory synchronization necessity
        //CheckHostResultsSync();
    }

    return chunkdatasize;
}

// -------------------------------------------------------------------------
// CalcMaxDbDataChunkSizeHelper: recalculate and return a maximum allowed 
// number of positions and memory size for db profile model data (profiles) 
// so that the implied search space does not exceed allocated memory 
// (representing the maximum allowed limit);
// NOTE: results are saved for the first memory area;
// querypmbeg, relative beginning address of queries;
// querypmend, relative terminal address of queries;
// 
size_t CuDeviceMemory::CalcMaxDbDataChunkSizeHelper( size_t nqyposs, size_t leftoversize )
{
    MYMSG( "CuDeviceMemory::CalcMaxDbDataChunkSizeHelper", 4 );
    const mystring preamb = "CuDeviceMemory::CalcMaxDbDataChunkSizeHelper: ";
    const size_t cszalnment = GetMemAlignment();
    //
    if( nareas_ < 1 || sz_heapsections_ == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error.");
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
    size_t tmpszdp2maxcoords, szdp2maxcoords = 0UL;//maximum size for coordinates of max alignment scores in phase 2 
    size_t tmpszdp2btckdat, szdp2btckdat = 0UL;//maximum size for phase-2 DP backtracking data
    size_t tmpszdp2fwdmtx, szdp2fwdmtx = 0UL;//maximum size for phase-2 forward probability matrix
//size_t tmpszdp2bwdmtx, szdp2bwdmtx = 0UL;//maximum size for phase-2 backward probability matrix
    size_t tmpszss2data, szss2data = 0UL;//maximum size for phase-2 statistical significance data
    size_t tmpszdp2alndata, szdp2alndata = 0UL;//maximum size for phase-2 alignment data
    size_t tmpszdp2alns, szdp2alns = 0UL;//maximum size for phase-2 alignments themselves
    //
    size_t tmpszglbvars, szglbvars = 0UL;//maximum size for global variables
    //
    size_t tmpszovlpos, szovlpos = 0UL;//overall size of position-dependent data excluding maxsizedbposs
    size_t tmpdiff;

    size_t nmatrices = nDevMatrixSections;
    if( !GetModScoreMatrixInUse())
        nmatrices--;

    //seqch space divided by the n. of db profile positions
    //note: assuming matrices of type CUBSM_TYPE
    const size_t qsspace = nqyposs * nmatrices * sizeof(CUBSM_TYPE);
    //const size_t leftoversize = deviceallocsize_ - sz_heapsections_[0][ddsEndOfCached];
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
        //global variables
        &tmpszglbvars,
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
            szdp2maxcoords = tmpszdp2maxcoords; 
            szdp2btckdat = tmpszdp2btckdat;
            szdp2fwdmtx = tmpszdp2fwdmtx; 
//szdp2bwdmtx = tmpszdp2bwdmtx; 
            szss2data = tmpszss2data; 
            szdp2alndata = tmpszdp2alndata; 
            szdp2alns = tmpszdp2alns;
            szovlpos = tmpszovlpos;
            szglbvars = tmpszglbvars;
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
            &tmpszdp2maxcoords, &tmpszdp2btckdat,
            &tmpszdp2fwdmtx,// &tmpszdp2bwdmtx, 
            &tmpszss2data, &tmpszdp2alndata, &tmpszdp2alns,
            //global variables
            &tmpszglbvars,
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
        sprintf(msgbuf, " (dev alloc= %zu; dev leftover= %zu; maxqrylen= %zu%s "
                "szdpdiag= %zu, szdpbottom= %zu, "
                "szdpmaxcoords= %zu, szdpbtckdat= %zu; szovlpos= %zu)", 
                deviceallocsize_, leftoversize, nqyposs, NL, szdpdiag, szdpbottom, 
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

    sz_heapsections_[0][ddsEndOfDbChunk] = sz_heapsections_[0][ddsEndOfCached] + maxsizedbposs;
    sz_heapsections_[0][ddsBegOfOrgScores] = sz_heapsections_[0][ddsEndOfDbChunk];
    sz_heapsections_[0][ddsEndOfOrgScores] = sz_heapsections_[0][ddsEndOfDbChunk] + szsmatrix;
    sz_heapsections_[0][ddsEndOfModScores] = sz_heapsections_[0][ddsEndOfOrgScores];
    if( GetModScoreMatrixInUse())
        sz_heapsections_[0][ddsEndOfModScores] += szsmatrix;
    //phase-1 sections
    sz_heapsections_[0][ddsEndOfDPDiagScores] = sz_heapsections_[0][ddsEndOfModScores] + szdpdiag;
    sz_heapsections_[0][ddsEndOfDPBottomScores] = sz_heapsections_[0][ddsEndOfDPDiagScores] + szdpbottom;
    sz_heapsections_[0][ddsEndOfDPMaxCoords] = sz_heapsections_[0][ddsEndOfDPBottomScores] + szdpmaxcoords;
    sz_heapsections_[0][ddsEndOfDPBackTckData] = sz_heapsections_[0][ddsEndOfDPMaxCoords] + szdpbtckdat;
    //phase-2 sections: reused previous occupations
    sz_heapsections_[0][ddsEndOfDP2DiagScores] = sz_heapsections_[0][ddsEndOfModScores] + szdp2diag;
    sz_heapsections_[0][ddsEndOfDP2BottomScores] = sz_heapsections_[0][ddsEndOfDP2DiagScores] + szdp2bottom;
//sz_heapsections_[ddsEndOfDP2MaxCoords] = sz_heapsections_[ddsEndOfDP2PrbScales] + szdp2maxcoords;
    sz_heapsections_[0][ddsEndOfDP2MaxCoords] = sz_heapsections_[0][ddsEndOfDP2BottomScores] + szdp2maxcoords;
    sz_heapsections_[0][ddsEndOfDP2BackTckData] = sz_heapsections_[0][ddsEndOfDP2MaxCoords] + szdp2btckdat;
    sz_heapsections_[0][ddsEndOfDP2FwdMtx] = sz_heapsections_[0][ddsEndOfDP2BackTckData] + szdp2fwdmtx;
//sz_heapsections_[ddsEndOfDP2BwdMtx] = sz_heapsections_[ddsEndOfDP2FwdMtx] + szdp2bwdmtx;
//sz_heapsections_[ddsEndOfSS2Data] = sz_heapsections_[ddsEndOfDP2BwdMtx] + szss2data;
    //statistical significance and alignment data section:
    // the beginning depends on whether MAP alignment is inuse
    sz_heapsections_[0][ddsEndOfSS2Data] = 
        (GetMAPAlnInUse()? sz_heapsections_[0][ddsEndOfDP2FwdMtx]: sz_heapsections_[0][ddsEndOfDPBackTckData]) + 
        szss2data;
    sz_heapsections_[0][ddsEndOfDP2AlnData] = sz_heapsections_[0][ddsEndOfSS2Data] + szdp2alndata;
    sz_heapsections_[0][ddsEndOfDP2Alns] = sz_heapsections_[0][ddsEndOfDP2AlnData] + szdp2alns;
    //variables:
    sz_heapsections_[0][ddsBegOfGlobVars] = sz_heapsections_[0][ddsEndOfDP2Alns];
    sz_heapsections_[0][ddsEndOfGlobVars] = sz_heapsections_[0][ddsBegOfGlobVars] + szglbvars;

    return maxsizedbposs;
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
void CuDeviceMemory::DestroyTextureObject( cudaTextureObject_t& texObj )
{
    MYMSG( "CuDeviceMemory::DestroyTextureObject", 6 );
    if( texObj ) {
        MYCUDACHECK( cudaDestroyTextureObject( texObj ));
        MYCUDACHECKLAST;
        texObj = 0;
    }
}

// -------------------------------------------------------------------------
// FreeDevicePtr: free device pointer
inline
void CuDeviceMemory::FreeDevicePtr( char*& d_ptr )
{
    MYMSG( "CuDeviceMemory::FreeDevicePtr", 6 );
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
void CuDeviceMemory::CacheSSENNWeightsDevice()
{
    MYMSG( "CuDeviceMemory::CacheSSENNWeightsDevice", 4 );
    MYCUDACHECK( cudaMemcpyToSymbol( dc_cuss_NN_weights_, SSE_NN_WEIGHTS_, 
                        sz_total_dc_cuss_NNwghts*sizeof(float)));
    MYCUDACHECKLAST;
}

// -------------------------------------------------------------------------
// CacheSSSScoresDevice: cache a score table of secondary structure 
// predictions to a device
// 
void CuDeviceMemory::CacheSSSScoresDevice( const SSSScores& sssscores )
{
    MYMSG( "CuDeviceMemory::CacheSSSScoresDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CacheSSSScoresDevice: ";

    SerializedScoresCtor<CUBSM_TYPE> ssss(
        sssscores.GetScores(), 
        sssscores.GetNoPrbLvs(), sssscores.GetNoTables(), sssscores.GetCardinality(), 
        sssscores.GetPrbLevels(), sssscores.GetLevels()
    );

    //size to align the section of SSS scores:
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = ssss.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[0][ddsEndOfPadding]) {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: %zu.", deviceallocsize_ );
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofpad = d_heap_ + sz_heapsections_[0][ddsEndOfPadding];

    MYCUDACHECK( cudaMemcpy( d_eofpad, ssss.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    for(int a = 0; a < nareas_; a++)
        sz_heapsections_[a][ddsEndOfSSSTable] = sz_heapsections_[0][ddsEndOfPadding] + szscores;

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
void CuDeviceMemory::CacheCVS2ScoresDevice( const CVS2Scores& cvs2scores )
{
    MYMSG( "CuDeviceMemory::CacheCVS2ScoresDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CacheCVS2ScoresDevice: ";

    SerializedCVS2ScoresCtor<CUBSM_TYPE> cvs2s( cvs2scores );

    //size to align the section of SSS scores:
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = cvs2s.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[0][ddsEndOfSSSTable])
    {
        char msgbuf[BUF_MAX];
        sprintf( msgbuf, "Insufficient amount of allocated device memory: "
                "%zu (ssss+cvs2s, %zu).", deviceallocsize_,
                szscores + sz_heapsections_[0][ddsEndOfSSSTable]);
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofssssdata = d_heap_ + sz_heapsections_[0][ddsEndOfSSSTable];

    MYCUDACHECK( cudaMemcpy( d_eofssssdata, cvs2s.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    for(int a = 0; a < nareas_; a++)
        sz_heapsections_[a][ddsEndOfCVS2SMap] = sz_heapsections_[0][ddsEndOfSSSTable] + szscores;

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
void CuDeviceMemory::CacheHDPScoresDevice( const HDPscores& hdpscores )
{
    MYMSG( "CuDeviceMemory::CacheHDPScoresDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CacheHDPScoresDevice: ";

    SerializedScoresCtor<CUBSM_TYPE> hdp1s(
        hdpscores.GetScores(),
        hdpscores.GetNoPrbLvs(), hdpscores.GetNoTables(), hdpscores.GetCardinality(),
        hdpscores.GetPrbLevels(), hdpscores.GetLevels()
    );

    //size to align the HDP scores in texture memory
    const size_t cszalnment = GetMemAlignment();
    size_t szscores = hdp1s.GetSizeAlloc();

    if( deviceallocsize_ <= szscores + sz_heapsections_[0][ddsEndOfCVS2SMap])
    {
        char msgbuf[BUF_MAX];
//         sprintf( msgbuf, "Size of HDP scores exceeds allocated device memory size: %lu.", 
//                 deviceallocsize_ );
        sprintf( msgbuf, "Insufficient amount of allocated device memory: "
                "%zu (ssss+cvs2s+hdp1s, %zu).", deviceallocsize_,
                szscores + sz_heapsections_[0][ddsEndOfCVS2SMap]);
        throw MYRUNTIME_ERROR( preamb + msgbuf );
    }

    char* d_eofcvs2sdata = d_heap_ + sz_heapsections_[0][ddsEndOfCVS2SMap];

    MYCUDACHECK( cudaMemcpy( d_eofcvs2sdata, hdp1s.GetScores(), szscores, cudaMemcpyHostToDevice ));
    MYCUDACHECKLAST;

    szscores = ALIGN_UP( szscores, cszalnment );

    for(int a = 0; a < nareas_; a++)
        sz_heapsections_[a][ddsEndOfHDPscores] = sz_heapsections_[0][ddsEndOfCVS2SMap] + szscores;

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
void CuDeviceMemory::PackPMDataForDevice( 
    size_t szpm2dvfields[], char*& tmpbuf, size_t& szalloc,
    char** querypmbeg,//beginning address of queries
    char** querypmend,//terminal address of queries
    char** bdb1pmbeg,//beginning address of cached database profiles
    char** bdb1pmend,//terminal address of cached database profiles
    char** bdbCpmbeg,//beginning address of new profiles read from the database 
    char** bdbCpmend )//terminal address of new profiles read from the database 
{
    MYMSG( "CuDeviceMemory::PackPMDataForDevice", 4 );
    const mystring preamb = "CuDeviceMemory::PackPMDataForDevice: ";
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
        free(tmpbuf);
        tmpbuf = NULL;
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
//
void CuDeviceMemory::CopyPMDataToDevice(
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
    MYMSG( "CuDeviceMemory::CopyPMDataToDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CopyPMDataToDevice: ";
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

        if( szalloc < 1 || tmpbuf == NULL ) {
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
// CopyCPMDataToDevice: copy packed profile model data to a device;
// bdbCpmbeg, host pointer for beginning addreses of read db profiles;
// bdbCpmend, host pointer for terminal addreses of read db profiles;
// szCpm2dvfields, sizes of the fields occupied by bdbC data;
// dev_pckdpm, device pointer to (the end of) packed data;
// szmaxsize, maximum amount of memory allowed to be occupied for the data passed by arguments;
// cmndx, index of constant memory array where field addresses reside;
// NOTE: memory is allocated for device pointer dev_pckdpm 
//
void CuDeviceMemory::CopyCPMDataToDevice(
    char** bdbCpmbeg,
    char** bdbCpmend,
    size_t* szCpm2dvfields,
    char* dev_pckdpm,
    size_t szmaxsize,
    size_t /*cmndx*/ )
{
    MYMSG( "CuDeviceMemory::CopyCPMDataToDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CopyCPMDataToDevice: ";
    myruntime_error mre;

    if( bdbCpmbeg == NULL || bdbCpmend == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null data structures.");
    if( szCpm2dvfields == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null data size structure.");

    //sizes of each field used to calculate addresses passed to a device
    size_t szalloc = szCpm2dvfields[pmv2DTotFlds], i;
    //device pointers to the fields of profile model data packed to one buffer
    char* d_pckdpmdatflds[fsg_sztotflds2];
    size_t nflds = (size_t)pmv2DTotFlds;

    try {
        if( szalloc < 1 ) {
            warning((preamb + "No data transfered!").c_str());
            return;
        }

        if( szmaxsize < szalloc ) {
            char msgbuf[BUF_MAX];
            mystring strbuf = preamb + "Short of allocated device memory: ";
            sprintf( msgbuf, "%zuB requested vs %zuB allowed.", szalloc, szmaxsize );
            throw MYRUNTIME_ERROR(strbuf + msgbuf);
        }

        //synchronous copy; have to wait for data to be present on device
        MYCUDACHECK( cudaMemcpy( dev_pckdpm, bdbCpmbeg[0], szalloc, cudaMemcpyHostToDevice ));
        MYCUDACHECKLAST;

        for( i = 0; i < nflds; i++ )
            d_pckdpmdatflds[i] = dev_pckdpm + szCpm2dvfields[i];

        //put device addresses pointing to locations of the same buffer to constant memory
        MYCUDACHECK( cudaMemcpyToSymbol( dc_pm2dvfields_, d_pckdpmdatflds, nflds*sizeof(char*), 
                            fsg_sztotflds2*sizeof(char*)/*offset*/));
        MYCUDACHECKLAST;

        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            mystring strbuf = preamb + "Data transfered:";
            sprintf(msgbuf, "  total data size= %zu", szalloc );
            strbuf += msgbuf;
            MYMSG(strbuf.c_str(),3);
        MYMSGENDl

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// CacheDataDevice: pack and transfer cached query and db profile model 
// data to a device
// 
inline
void CuDeviceMemory::CacheDataDevice(
    char** querypmbeg,
    char** querypmend,
    char** bdb1pmbeg,
    char** bdb1pmend )
{
    MYMSG( "CuDeviceMemory::CacheDataDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CacheDataDevice: ";
    const size_t cszalnment = GetMemAlignment();

    char* d_eofcnstdata = d_heap_ + sz_heapsections_[0][ddsEndOfConstantData];
    char* d_eofcnstdataplus = d_eofcnstdata;

    if( deviceallocsize_ <= sz_heapsections_[0][ddsEndOfConstantData])
        throw MYRUNTIME_ERROR( preamb + "Insufficient device memory allocated to cache data.");

    //maxsize for data; TODO: to be changed
    size_t szmaxallowed = deviceallocsize_ - sz_heapsections_[0][ddsEndOfConstantData];

    CopyPMDataToDevice( 
        querypmbeg, querypmend, bdb1pmbeg, bdb1pmend,
        NULL/*bdbCpmbeg*/, NULL/*bdbCpmend*/,
        d_eofcnstdataplus, 
        szmaxallowed,
        ndx_dc_pm2dvfields_ );
    //
    size_t szdata = (size_t)(d_eofcnstdataplus-d_eofcnstdata);
    szdata = ALIGN_UP( szdata, cszalnment );
    for(int a = 0; a < nareas_; a++)
        sz_heapsections_[a][ddsEndOfCached] = sz_heapsections_[0][ddsEndOfConstantData] + szdata;
}

// -------------------------------------------------------------------------
// TransferCPMDataDevice: transfer a new chunk of db profile model data to a
// device
// 
inline
void CuDeviceMemory::TransferCPMDataDevice(
        char** bdbCpmbeg,
        char** bdbCpmend,
        size_t* szCpm2dvfields )
{
    MYMSG( "CuDeviceMemory::CacheCPMDataDevice", 4 );
    const mystring preamb = "CuDeviceMemory::CacheCPMDataDevice: ";

    if( sz_heapsections_[0][ddsEndOfDbChunk] <= sz_heapsections_[0][ddsEndOfCached])
        throw MYRUNTIME_ERROR( preamb + "Unallocated memory for chunked db data.");

    //eof cached data and beginning of variable data
    char* d_eofcacheddata = d_heap_ + sz_heapsections_[0][ddsEndOfCached];

    //maxsize for data
    size_t szmaxallowed = sz_heapsections_[0][ddsEndOfDbChunk] - sz_heapsections_[0][ddsEndOfCached];

    CopyCPMDataToDevice(
        bdbCpmbeg, bdbCpmend, szCpm2dvfields,
        d_eofcacheddata,
        szmaxallowed,
        ndx_dc_newpm2dvfields_ );
}
