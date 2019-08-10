/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init.cuh"
#include "libmycu/cudp/CuBatchDP_init_btck.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "CuBatchMAPDP_fwd.cuh"
#include "CuBatchMAPDP_bwd.cuh"
#include "CuBatchMAPDP_map_btck.cuh"

// #define CUMAPDP_MAP_BTCK_TESTPRINT 453 //0 //1286 //0 //3024//4899

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming with backtracking 
// information;
// NOTE: Version for thread block of one warp!
// alnextreg, alignment extension regularizer;
// alnextgapreg, regularizer for gaps during alignment extension;
// blkdiagnum, block diagonal serial number;
// lastydiagnum, last block diagonal serial number along y axis 
// (starting at x=-CUMAPDP_2DCACHE_DIM);
// ndb1pros, number of profiles in the first profile data buffer db1;
// querprosOmtd, number of query profiles up to this query profile;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// nqyposs, number of query positions to process;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
// dblen2, total number of positions of profiles passed to phase 2 PLUS 
// padding for this data;
// dbxpad2, number of padded positions for memory alignment;
//

// DP processing layout:
// +---------===-+--------+--+-
// |  /  /  /  / |  /  /  |  |
// | /  /  /  /  | /  /  /| /|
// |/  /  /  /  /|/  /  / |/ |
// +---======----+---=====+--+-
// |  /  /  /  / |  /  /  |  |
// +=====--------+=====---+--+-
// (double line indicates current parallel processing)

// -------------------------------------------------------------------------
// ExecMAPDP_Bwd_Unroll32x: device code for executing dynamic programming to 
// calculate backward probabilities using shared memory and 32-fold 
// unrolling along the diagonal of dimension CUMAPDP_2DCACHE_DIM;
// NOTE: memory pointers should be aligned!
// tmpfwdprobs, forward/backward probabilities;
// dp2alndatbuffers, profile-specific statistics from (previous) phase 1;
// tmpdpdiagbuffers, temporary buffers for last calculated diagonal scores;
// tmpdpbotbuffer, temporary buffers for last calculated bottom scores;
// tmpss2datbuffers, temporary buffer of the total probability for each 
// pair of profiles;
// maxscoordsbuf, coordinates of maximum alignment scores;
// btckdata, backtracking information data;
//
__global__ void ExecMAPDP_MAP_Btck_Unroll32x(
    const float alnextreg,
    float alnextgapreg,
    uint blkdiagnum,
    uint lastydiagnum,
    float logevthld,
    uint ndb1pros,
    uint /*querprosOmtd*/, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, //uint ndb1poss, uint ndbCposs, uint dbxpad,
    //uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    uint dblen2,
    uint dbxpad2,
    const float* __restrict__ tmpfwdprobs,
    const float* __restrict__ dp2alndatbuffers,
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpss2datbuffers,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata )
{
    /*__shared__ */INTYPE dbprolenCache;//length of profile
    //cache of scores: use padding of +1 to eliminate bank conflicts when accessing the cache
    __shared__ float
            scoreCache[CUMAPDP_2DCACHE_DIM_D][CUMAPDP_2DCACHE_DIM_X+CUDPIDS];
    __shared__ float
            diag1Cache[CUMAPDP_2DCACHE_DIM_D+1],//cache for state scores of the 1st diagonal
            diag2Cache[CUMAPDP_2DCACHE_DIM_D+1],//cache for state scores of the last diagonal
            bottmCache[CUMAPDP_2DCACHE_DIM_X];//cache for state scores of the top of the diagonals
    float maxscCache;
    //SECTION for backtracking information
    uint maxscCoords;//coordinates of the maximum alignment score maxscCache
    __shared__ char btckCache[CUMAPDP_2DCACHE_DIM_D][CUMAPDP_2DCACHE_DIM_X+CUDPIDS];//backtracking cache
    //
    //log total forward probability
    float logtotfwdprob;
    //using SMEM (even with bank conflicts) works faster than using (with constant indices) and 
    // synchronizing registers here (small number of accesses):
    __shared__ uint proattr[2];//original profile NUMBER and new DISTANCE from the beginning
    __shared__ float evals[2];//log of e-value and pair e-value
//     float eval;//log e-value
    uint dbfldsndx;
    uint pronr;


    // blockIdx.y is the profile serial number in phase 2;
    //{{get the original profile number and e-value too
    if( /*threadIdx.y == 0 && */threadIdx.x < 2 ) {
        proattr[threadIdx.x] = 
            *(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.y+dp2oadOrgProNo+threadIdx.x);
        evals[threadIdx.x] = dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.y+dp2oadEvalue+threadIdx.x];
//         if( threadIdx.x == 0 )
//            eval = dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.y+dp2oadEvalue];
    }
//     eval = __shfl_sync(0xffffffff, eval, 0);
    //}}

    __syncthreads();

    //all threads exit if calculated e-value exceeds the threshold
    if( logevthld < evals[0]/*eval*/)
        return;


    //NOTE: protection against overflow ensured on the host side
    if( proattr[0] < ndb1pros ) { pronr = proattr[0] + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = proattr[0] - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( /*threadIdx.y == 0 && */threadIdx.x == 0 ) {
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[pronr];
    }


    // blockIdx.x is block serial number s within diagonal blkdiagnum;
    // (x,y) is the bottom-left corner (x,y) coordinates for profile proattr[0];
    int x, y;
    if( blkdiagnum <= lastydiagnum) {
        //x=-!(d%2)w+2ws; y=dw/2+w-sw -1 (-1, zero-based indices); [when w==b]
        //(b, block's length; w, block's width)
        x = (2*blockIdx.x - (!(blkdiagnum & 1))) * CUMAPDP_2DCACHE_DIM_D;
        y = ((blkdiagnum>>1) + 1 - blockIdx.x) * CUMAPDP_2DCACHE_DIM_D - 1;
    } else {
        //x=-w+(d-d_l)w+2ws; y=dw/2+w-sw -1; [when w==b]
        x = (2*blockIdx.x + (blkdiagnum-lastydiagnum-1)) * CUMAPDP_2DCACHE_DIM_D;
        y = ((lastydiagnum>>1) + 1 - blockIdx.x) * CUMAPDP_2DCACHE_DIM_D - 1;
    }


    //{{use registers efficiently
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}


    if(dblen2 < proattr[1]+dbprolenCache+dbxpad2)
        return;


    //number of iterations for this block to perform;
    //uint makes the compiler to reduce #registers
    uint ilim = GetMaqxNoIterations( x, y, nqyposs, dbprolenCache, CUMAPDP_2DCACHE_DIM_X );

    if( y < 0 || (int)nqyposs <= (int)(y+1 - CUMAPDP_2DCACHE_DIM_D) || 
        (int)dbprolenCache <= x /*+ CUMAPDP_2DCACHE_DIM_DpX */ ||
        ilim < 1 )
        //block does not participate: out of profile boundaries
        return;

    //{{NOTE: do not align the distant corners;
    //get the edge length of a triangular area to be ignored:
    int crnedgelen = GetCornerMaxlen(nqyposs, dbprolenCache);
    //}}



    int qpos = y - threadIdx.x;//going upwards



    //x is now the position this thread will process
    x += threadIdx.x;

    //dbpos2 is the x position of the diagonal block's bottom-left 
    // corner in the phase-2 diagonal and bottom structures plus the offset 
    // determined by thread id:
    //(NOTE: int: negative x may complicate an overflow check)
    int dbpos2 = x + proattr[1];



    //get log of the total forward probability
    logtotfwdprob = MAPDP_GetLogTotalFwdProbability( 
        blkdiagnum == 0,
        blockIdx.y/*pronr2*/,
        tmpss2datbuffers );



    //cache forward/backward probabilities over the oblique diagonal block
    #pragma unroll 4
    for( int i = 0; i < CUMAPDP_2DCACHE_DIM_D; i++ ) {
        //going upwards
        scoreCache[i][threadIdx.x] = CUBDP_Q(0);
        if( CellYXinValidArea(nqyposs, dbprolenCache, crnedgelen, y-i, x+i )) {
            //starting position of line i of the oblq. diagonal block in the 
            // forward probability matrix: 
            qpos = (y-i) * dblen2 + i;
            scoreCache[i][threadIdx.x] = tmpfwdprobs[qpos + dbpos2];
        }
    }


    //cache two diagonals from the previous (along the x axis) diagonal block;
    //the structure of tmpdpdiagbuffers is position-specific (1D, along x-axis)
    if( 0 <= x-1 && x < (int)dbprolenCache ) {
        qpos = dpdssDiag2 * nTDPDiagScoreSubsections * dblen2;
        diag1Cache[threadIdx.x] = tmpdpdiagbuffers[dbpos2-1];
        diag2Cache[threadIdx.x+1] = tmpdpdiagbuffers[qpos + dbpos2-1];
        //cache the buffer of maximum probabilities for a query position
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen2;
        maxscCache = tmpdpdiagbuffers[qpos + dbpos2 - 1];
        maxscCoords = maxscoordsbuf[dbpos2-1];
    }
    else {
        diag1Cache[threadIdx.x] = CUBDP_Q(0);
        diag2Cache[threadIdx.x+1] = CUBDP_Q(0);
        maxscCache = CUBDP_Q(0);
    }


    //cache the bottom of the above diagonal blocks;
    //the structure of tmpdpbotbuffer is position-specific (1D, along x-axis)
    if( CUMAPDP_2DCACHE_DIM_D <= y && 
        0 <= x+CUMAPDP_2DCACHE_DIM_D-1 && x+CUMAPDP_2DCACHE_DIM_D-1 < (int)dbprolenCache )
        bottmCache[threadIdx.x] = tmpdpbotbuffer[dbpos2+CUMAPDP_2DCACHE_DIM_D-1];
    else
        bottmCache[threadIdx.x] = CUBDP_Q(0);


//     //dynamic gap extension multiplier
//     float dmul;
//     MAPDPGetDynGapExtParams( evals[1], &dmul );
//     alnextgapreg *= dmul;

    __syncthreads();


    float *pdiag1 = diag1Cache;
    float *pdiag2 = diag2Cache;

    //start calculations for this position with 32x unrolling
    //NOTE: sync inside: do not branch;
    for( int i = 0; i < ilim/*CUMAPDP_2DCACHE_DIM_X*/; i++ ) {
        float val1;
        int btck = dpbtckSTOP;
        //
        if( threadIdx.x == blockDim.x-1 ) {
            pdiag1[threadIdx.x+1] = bottmCache[i];
        }

        //match state (diagonal direction)
        val1 = pdiag2[threadIdx.x+1] + 
            MAPDPexpf(scoreCache[threadIdx.x][i] - logtotfwdprob) -
            alnextreg;
        if( val1 > 0.0f ) btck = dpbtckDIAG; else val1 = 0.0f;

        __syncthreads();

        //NOTE: coord. overflow will be managed below while writing data
        dpmaxandcoords( maxscCache, val1, maxscCoords, x+i, y-threadIdx.x );

        //delete/insert state (up direction)
        myhdmaxassgn( val1, pdiag1[threadIdx.x+1] - alnextgapreg, btck, (int)dpbtckUP );

        //delete/insert state (left direction)
        myhdmaxassgn( val1, pdiag1[threadIdx.x] - alnextgapreg, btck, (int)dpbtckLEFT );

        pdiag2[threadIdx.x] = val1;//WRITE
        btckCache[threadIdx.x][i] = btck;

        if( threadIdx.x == 0 ) {
            //WRITE: this position is not used by other threads in the current iteration
            bottmCache[i] = pdiag2[threadIdx.x];
        }

#ifdef CUMAPDP_MAP_BTCK_TESTPRINT
        if(pronr==CUMAPDP_MAP_BTCK_TESTPRINT){
            __syncthreads();
            printf(" d=%u(%u) s=%u i%02d/%u (t%02u): len= %d (yx: %d,%d) "
                    "MM= %.6f  MAX= %.6f COORD= %x BTCK= %d  logs: F+B= %.4f FE= %.4f\n",
                    blkdiagnum,lastydiagnum,blockIdx.x,i,ilim,threadIdx.x,
                    dbprolenCache, y-threadIdx.x,x+i,
                    pdiag2[threadIdx.x],
                    maxscCache, maxscCoords, btck,
                    scoreCache[threadIdx.x][i], logtotfwdprob
            );
            for(size_t _k=0;_k<1000000000UL;_k++)clock();
            for(size_t _k=0;_k<1000000000UL;_k++)clock();
        }
#endif

        myhdswap( pdiag1, pdiag2 );

        __syncthreads();
    }


#ifdef CUMAPDP_MAP_BTCK_TESTPRINT
        if(pronr==CUMAPDP_MAP_BTCK_TESTPRINT)
            printf(" >>> d=%u(%u) s=%u (t%02u): pronr= %u MAX= %.6f COORD= %x wrt= %d xpos= %d\n", 
                    blkdiagnum,lastydiagnum,blockIdx.x,threadIdx.x, pronr,
                    maxscCache, maxscCoords, x+CUMAPDP_2DCACHE_DIM_X-1<(int)dbprolenCache, 
                    x+CUMAPDP_2DCACHE_DIM_X-1);
#endif


    //write the result of calculations for the upcoming blocks;
    //write two diagonals;
    if( 0 <= x+CUMAPDP_2DCACHE_DIM_X-1 && x+CUMAPDP_2DCACHE_DIM_X-1 < (int)dbprolenCache ) {
        qpos = dpdssDiag2 * nTDPDiagScoreSubsections * dblen2;
        tmpdpdiagbuffers[dbpos2+CUMAPDP_2DCACHE_DIM_X-1] = pdiag1[threadIdx.x];
        tmpdpdiagbuffers[qpos + dbpos2+CUMAPDP_2DCACHE_DIM_X-1] = pdiag2[threadIdx.x+1];
    }

    int ndx = CUMAPDP_2DCACHE_DIM_X-1;
    int xtl = x-threadIdx.x+blockDim.x-1;//top-left x coordinate of the diagonal block
    if( xtl >= dbprolenCache )
        ndx = 0;
    else if( xtl+ndx >= dbprolenCache )
        ndx = dbprolenCache - xtl - 1;

    //write the buffer of maximum probabilities
    if( 0 <= x+ndx && x+ndx < (int)dbprolenCache ) {
        //NOTE: make an additional check as zero costs (alnextreg) may 
        // induce max score propagation out of boundaries
        //if((uint)dbprolenCache <= GetCoordX(maxscCoords) || 
        //    nqyposs <= GetCoordY(maxscCoords)) {
        if(!CellYXinValidArea(nqyposs, dbprolenCache, crnedgelen, 
                GetCoordY(maxscCoords), GetCoordX(maxscCoords))) {
            maxscCache = 0.0f;
            maxscCoords = 0;
        }
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen2;
        //tmpdpdiagbuffers[qpos + dbpos2+ndx] = maxscCache[threadIdx.x];
        if( CUMAPDP_2DCACHE_DIM_D <= y ) {
            if( tmpdpdiagbuffers[qpos + dbpos2+ndx] < maxscCache ) {
                tmpdpdiagbuffers[qpos + dbpos2+ndx] = maxscCache;
                maxscoordsbuf[dbpos2+ndx] = maxscCoords;
            }
        } else {
            tmpdpdiagbuffers[qpos + dbpos2+ndx] = maxscCache;
            maxscoordsbuf[dbpos2+ndx] = maxscCoords;
        }
    }

    //write the bottom of the diagonal block;
    if( 0 <= x && x < dbprolenCache )
        tmpdpbotbuffer[dbpos2] = bottmCache[threadIdx.x];

    //write backtracking information with unrolling
    #pragma unroll 8
    for( int i = 0; i < CUMAPDP_2DCACHE_DIM_D; i++ ) {
        //going upwards
        if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x+i && x+i < (int)dbprolenCache ) {
            //starting position of line i of the oblq. diagonal block in the matrix:
            qpos = (y-i) * dblen2 + i;
            btckdata[qpos + dbpos2] = btckCache[i][threadIdx.x];
        }
    }
}
