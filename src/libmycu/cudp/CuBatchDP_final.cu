/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/mymutex.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "CuBatchDP_com.h"
#include "CuBatchDP_init_btck.cuh"
#include "CuBatchDP_final.cuh"

// #define CUDP_FINAL_TESTPRINT 1682 //1286 //-1 //0//3024//4899

// =========================================================================

//attributes of passed profiles
__device__ unsigned int d_gDPPassedPros[nTPassedPros];
//mutex for updating the attributes of passed profiles
__device__ volatile int d_dpfinal_mutex = 0;

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for finalizing dynamic programming;
// NOTE: thread block is supposed to consist of one warp;
// blkdiagnum, block diagonal serial number;
// lastydiagnum, last block diagonal serial number along y axis 
// (starting at x=-CUDP_2DCACHE_DIM);
// scorethld, score threshold for alignments to pass to the next phase;
// ndb1pros, number of profiles in the first profile data buffer db1;
// ndbCpros, number of profiles in the second profile data buffer dbC;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// dbxpad, number of padded positions for memory alignment;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
//

// -------------------------------------------------------------------------
// FinalizeDPPro: device code for finalizing dynamic programming in parallel
// over profiles;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers of calculated diagonal scores;
// tmpdpbotbuffer, this buffer is used to store now profile-specific 
//  information once the phase-1 DP has completed;
// maxscoordsbuf, coordinates of maximum alignment scores;
// btckdata, backtracking information data;
// 
__global__ void FinalizeDPPro(
    CUBSM_TYPE scorethld,
    uint ndb1pros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
//     CUBDP_TYPE* __restrict__ tmpdpbotbuffer,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata )
{
    /*__shared__ */LNTYPE dbprodstCache;//distance in positions to db profile blockIdx.x
    /*__shared__ */INTYPE dbprolenCache;//length of profile blockIdx.x
    //NOTE: keep SMEM below 2KB to consistently ensure high occupancy;
    //NOTE: assumed CUDP_CORR_NSCORES < CUDP_2DCACHE_DIM_D; 
    //reuse maxscCache
    //__shared__ CUBSM_TYPE
    //        scoreCache[CUDP_2DCACHE_DIM_D+CUDP_CORR_NSCORES];
    //cache of maximum scores from the last processed diagonals
    // (add 1 to dimensions to avoid bank conflicts):
    __shared__ CUBDP_TYPE maxscCache[TIMES2(CUDP_2DCACHE_DIM_D)+1];
    __shared__ int maxscndxCache[TIMES2(CUDP_2DCACHE_DIM_D)];//indices
    __shared__ CUBDP_TYPE outputCache[nTDP2OutputPhase1];
    //
    CUBDP_TYPE maxsc;
    int ndx;
    //
    // blockIdx.x is the profile serial number;
    uint dbfldsndx;
    uint pronr;
    //NOTE: protection against overflow ensured on the host side
    if( blockIdx.x < ndb1pros ) { pronr = blockIdx.x + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = blockIdx.x - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( threadIdx.x == 0 ) {
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[pronr];
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[pronr];
    }

    //{{use registers efficiently
    dbprodstCache = __shfl_sync(0xffffffff, dbprodstCache, 0);
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}
    //...or use shared memory; bank conflicts arise when accessed
//     __syncthreads();


    const int dblen = ndb1poss + ndbCposs + dbxpad;
    //dbpos is the beginning x position in the score matrix;
    const int dbpos = (blockIdx.x < ndb1pros)? dbprodstCache - bdb1posoffset: dbprodstCache + ndb1poss;
    //db profle position calculated from the end, including the offset determined by tid:
    int dbposoff = dbpos + dbprolenCache-1-threadIdx.x;

    int qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
    maxscCache[threadIdx.x] = CUBDP_Q(0);
    maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D] = CUBDP_Q(0);
    //read maximum alignment scores at the end
    if( threadIdx.x < dbprolenCache)
        maxscCache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbposoff];
    if( threadIdx.x+CUDP_2DCACHE_DIM_D < dbprolenCache )
        maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D] = 
            tmpdpdiagbuffers[qpos + dbposoff-CUDP_2DCACHE_DIM_D];

    __syncthreads();

    //find the maximum score and its index (by index propagation)
    dpfinmaxndxinit( maxscCache, maxscndxCache, threadIdx.x, threadIdx.x+CUDP_2DCACHE_DIM_D );
    __syncthreads();
    for(int i = CUDP_2DCACHE_DIM_D >> 1; i >= 1; i >>= 1 ) {
        dpfinmaxndx( maxscCache, maxscndxCache, threadIdx.x, threadIdx.x+i );
        __syncthreads();
    }

    maxsc = maxscCache[0];
    ndx = maxscndxCache[0];

    if( maxsc <= scorethld ) {
        //useless alignment score; 
        //write the result, and all the block's threads exit
#if !defined(UNSRT_PASSED_PROFILES)
        if( threadIdx.x == 0 ) {
            //tmpdpbotbuffer[blockIdx.x*nTDP2OutputPhase1] = maxsc;//CUBDP_Q(0)
            tmpdpdiagbuffers[blockIdx.x*nTDP2OutputPhase1] = maxsc;//CUBDP_Q(0)
        }
#endif
        return;
    }

    int xx, x, y;
    int psts = 0;//number of positive substitution scores
    char btck = dpbtckDIAG;
    //read the coordinates of the maximum score
    if( threadIdx.x == 0 )
        xx = maxscoordsbuf[dbposoff-ndx];
    //coordinates xx, x, and y are used only by tid 0
    //xx = __shfl_sync(0xffffffff, xx, 0);
    y = GetCoordY(xx);
    x = GetCoordX(xx);

#ifdef CUDP_FINAL_TESTPRINT
    if((int)pronr==CUDP_FINAL_TESTPRINT){
        printf(" Fin: bid= %u tid= %u: len= %d addr= %u maxsc= %.6f ndx= %d y= %d x= %d\n",
            blockIdx.x,threadIdx.x,
            dbprolenCache,dbprodstCache,maxsc,ndx,y,x
        );
    }
#endif

    CUBDP_TYPE corrsc, sumcorrsc = CUBDP_Q(0);

    if( threadIdx.x < CUDP_CORR_NSCORES )
        //initialize in parallel
        maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D] = CUBDP_Q(0);
    __syncthreads();

    //backtrace over the alignment to calculate correlation score;
    while( btck != dpbtckSTOP ) {
        //first, cache a window of scores;
        if( threadIdx.x < CUDP_CORR_NSCORES )
            //copy the most recent scores to the end (physically beginning)
            // in parallel
            maxscCache[threadIdx.x] = maxscCache[threadIdx.x+CUDP_2DCACHE_DIM_D];
        if( threadIdx.x == 0 ) {
            //only thread 0 caches scores
            for( ndx = CUDP_CORR_NSCORES; ndx < CUDP_2DCACHE_DIM_D+CUDP_CORR_NSCORES; ) {
                if( x < 0 || y < 0 ) {
                    btck = dpbtckSTOP;
                    break;
                }
                qpos = y * dblen + dbpos + x;
                btck = btckdata[qpos];
                if( btck == dpbtckSTOP )
                    break;
                if( btck == dpbtckUP ) { y--; continue; }
                else if( btck == dpbtckLEFT ) { x--; continue; }
                //(btck == dpbtckDIAG)
                maxscCache[ndx] = scores[qpos];
                if( maxscCache[ndx] > 0.0f )
                    psts++;
                x--; y--;
                ndx++;
            }
        }
        __syncwarp();
        btck = __shfl_sync(0xffffffff, btck, 0);
        ndx = __shfl_sync(0xffffffff, ndx, 0);
        //next, parallel calculation and reduction
        if( ndx <= threadIdx.x+CUDP_CORR_NSCORES )
            corrsc = CUBDP_Q(0);
        else
            corrsc = maxscCache[threadIdx.x+CUDP_CORR_NSCORES] * 
                    ( maxscCache[threadIdx.x] + maxscCache[threadIdx.x+1] +
                    maxscCache[threadIdx.x+2] + maxscCache[threadIdx.x+3] );
        corrsc += __shfl_xor_sync(0xffffffff, corrsc, 16);
        corrsc += __shfl_xor_sync(0xffffffff, corrsc, 8);
        corrsc += __shfl_xor_sync(0xffffffff, corrsc, 4);
        corrsc += __shfl_xor_sync(0xffffffff, corrsc, 2);
        corrsc += __shfl_xor_sync(0xffffffff, corrsc, 1);
        sumcorrsc += corrsc;
#ifdef CUDP_FINAL_TESTPRINT
        if((int)pronr==CUDP_FINAL_TESTPRINT){
            printf(" Fin loop: sumcorrsc= %.6f y= %d x= %d\n", sumcorrsc,y,x );
        }
#endif
    }

    int npr = -1;
    uint dst;

    //write the maximum alignment score
    if( threadIdx.x == 0 ) {
        sumcorrsc *= 0.5f;
        sumcorrsc = sumcorrsc<0? -sqrtf(-sumcorrsc): sqrtf(sumcorrsc);
#ifdef CUDP_FINAL_TESTPRINT
        if((int)pronr==CUDP_FINAL_TESTPRINT){
            printf(" Fin0: maxsc= %.6f sumcorrsc= %.6f added= %.6f y= %d x= %d psts= %d\n",
                maxsc,sumcorrsc,maxsc+sumcorrsc,y,x,psts );
        }
#endif

        maxsc += sumcorrsc;

        if( scorethld < maxsc ) {
#ifdef UNSRT_PASSED_PROFILES
            //get profile serial number and address for this profile;
            //NOTE: perform the operations locked by a mutex so that
            // other thread blocks do not intervene in between!
            LOCK( &d_dpfinal_mutex );
                npr = atomicAdd( d_gDPPassedPros+dpppNPassedPros, 1 );
                dst = atomicAdd( d_gDPPassedPros+dpppNPosits, dbprolenCache );
            UNLOCK( &d_dpfinal_mutex );
            atomicMax( d_gDPPassedPros+dpppMaxProLen, dbprolenCache );
#else
            npr = blockIdx.x;
#endif
            outputCache[dp2oadScore] = maxsc;//score
            *(uint*)(outputCache+dp2oadOrgProNo) = blockIdx.x;//original profile number
            *(uint*)(outputCache+dp2oadProNewDst) = dst;//profile new distance from the beginning
            outputCache[dp2oadPstvs] = psts;//# positive scores
            //alignment beginning coordinates: beg state is always a match, add 1
            *(uint*)(outputCache+dp2oadBegCoords) = CombineCoords(x+1,y+1);
            *(uint*)(outputCache+dp2oadEndCoords) = xx;//end coordinates of the alignment
        }
#if !defined(UNSRT_PASSED_PROFILES)
        else {
            //tmpdpbotbuffer[nTDP2OutputPhase1*blockIdx.x] = maxsc;//CUBDP_Q(0)
            tmpdpdiagbuffers[nTDP2OutputPhase1*blockIdx.x] = maxsc;//CUBDP_Q(0)
        }
#endif
    }

    //NOTE:this is valid only if the block dimension matches the warp size!
    //otherwise, use SMEM and __syncthreads()
    __syncwarp();
    npr = __shfl_sync(0xffffffff, npr, 0);

    if( 0 <= npr && threadIdx.x < nTDP2OutputPhase1 )
        //tmpdpbotbuffer[nTDP2OutputPhase1*npr+threadIdx.x] = outputCache[threadIdx.x];
        tmpdpdiagbuffers[nTDP2OutputPhase1*npr+threadIdx.x] = outputCache[threadIdx.x];
}

// -------------------------------------------------------------------------
// FinalizeDP: device code for finalizing dynamic programming;
// NOTE: This kernel is the parent kernel to correspond to one thread!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers of calculated diagonal scores used to 
// store now profile-specific information;
// maxscoordsbuf, coordinates of maximum alignment scores;
// btckdata, backtracking information data;
//
__global__ void FinalizeDP(
    CUBSM_TYPE scorethld,
    uint ndb1pros, uint ndbCpros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ /*tmpdpbotbuffer*/,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata )
{
    if( blockIdx.x || threadIdx.x )
        return;

    uint npros = ndb1pros+ndbCpros;

#ifdef UNSRT_PASSED_PROFILES
    d_gDPPassedPros[dpppNPassedPros] = 0;
    d_gDPPassedPros[dpppNPosits] = 0;
    d_gDPPassedPros[dpppMaxProLen] = 0;
#endif

    UNLOCK( &d_dpfinal_mutex );

    FinalizeDPPro<<<npros,CUDP_2DCACHE_DIM_D>>>( 
        scorethld,
        ndb1pros,
        ndb1prosOmtd, ndbCprosOmtd,
        ndb1poss, ndbCposs, dbxpad,
        bdb1posoffset, bdbCposoffset,
        scores,
        tmpdpdiagbuffers, 
//         tmpdpbotbuffer,
        maxscoordsbuf,
        btckdata
    );

#if !defined(UNSRT_PASSED_PROFILES)
    MYCUDACHECK( cudaDeviceSynchronize());
    MYCUDACHECKLAST;

    //compact data of profiles that require further processing, 
    // once DP has been finalized
    const int dblen = ndb1poss + ndbCposs + dbxpad;
    uint npr = 0;//current write position (number of profiles)
    uint totposs = 0;//total number of positions

    //NOTE: length is not included due to the output buffer
    // dimensions constraints
    for( uint i = 0; i < npros; i++ ) {
        //CUBDP_TYPE alnsc = tmpdpbotbuffer[nTDP2OutputPhase1*i];
        CUBDP_TYPE alnsc = tmpdpdiagbuffers[nTDP2OutputPhase1*i];
        uint len = 0;
        if( alnsc <= scorethld )
            continue;
        ///len = (uint)tmpdpbotbuffer[nTDP2OutputPhase1*i+];
        ///len = (uint)tmpdpdiagbuffers[nTDP2OutputPhase1*i+];
        totposs += len;
        if( npr == 0 )
            d_gDPPassedPros[dpppMaxProLen] = len;
        if( i != npr ) {
            tmpdpdiagbuffers[nTDP2OutputPhase1*npr+dp2oadScore] = alnsc;
            *(uint)(tmpdpdiagbuffers + nTDP2OutputPhase1*npr+dp2oadOrgProNo) = i;//profile number
            ///tmpdpbotbuffer[nTDP2OutputPhase1*npr+] = len;//profile length
            ///tmpdpdiagbuffers[nTDP2OutputPhase1*npr+] = len;//profile length
        }
        npr++;
    }

    d_gDPPassedPros[dpppNPassedPros] = npr;
    d_gDPPassedPros[dpppNPosits] = totposs;
#endif

#ifdef CUDP_FINAL_TESTPRINT
# ifdef UNSRT_PASSED_PROFILES
    const int dblen = ndb1poss + ndbCposs + dbxpad;
# endif
    for( uint i = 0; i < d_gDPPassedPros[dpppNPassedPros]; i++ ) {
        printf(" FinalizeDP %u SC= %.6f SNR= %u LEN= %u TOTPOSS= %u MAXLEN= %u\n",
            i, tmpdpdiagbuffers/*tmpdpbotbuffer*/[nTDP2OutputPhase1*i],
            *(uint*)(tmpdpdiagbuffers/*tmpdpbotbuffer*/ + nTDP2OutputPhase1*i+dp2oadOrgProNo),
            0,//(uint)tmpdpdiagbuffers/*tmpdpbotbuffer*/[nTDP2OutputPhase1*i+], 
            d_gDPPassedPros[dpppNPosits], d_gDPPassedPros[dpppMaxProLen]
        );
    }
#endif
}
