/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"
#include "CuBatchDP_init.cuh"
#include "CuBatchDP_init_corr.cuh"

// #define CUDP_CORR_INIT_TESTPRINT

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming;
// NOTE: Version for thread block of one warp!
// NOTE: This version also calculates correlated match scores along the 
// alignments inline!
// blkdiagnum, block diagonal serial number;
// lastydiagnum, last block diagonal serial number along y axis 
// (starting at x=-CUDP_2DCACHE_DIM);
// ndb1pros, number of profiles in the first profile data buffer db1;
// querprosOmtd, number of query profiles up to this query profile;
// ndb1prosOmtd, number of profiles missed up to the first one in db1;
// ndbCprosOmtd, number of profiles missed up to the first one in dbC;
// nqyposs, number of query positions to process;
// ndb1poss, number of cached db profile positions to process;
// ndbCposs, number of new db profile positions to process;
// dbxpad, number of padded positions for memory alignment;
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
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
// ExecDP_Corr_Unroll32x: device code for executing dynamic programming 
// using shared memory and 32-fold unrolling along the diagonal of dimension 
// CUDP_2DCACHE_DIM;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers for last calculated diagonal scores;
// tmpdpbotbuffer, temporary buffers for last calculated bottom scores;
// 
__global__ void ExecDP_Corr_Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer )
{
    /*__shared__ */LNTYPE dbprodstCache;//distance in positions to db profile blockIdx.y
    /*__shared__ */INTYPE dbprolenCache;//length of profile blockIdx.y
    __shared__ FPTYPE 
            qrtrpCache[nTDPUsedProTrn][CUDP_2DCACHE_DIM_D],//cache for query transition probabilities
            dbtrpCache[nTDPUsedProTrn][CUDP_2DCACHE_DIM_DpX];//cache for db transition probabilities
    //cache of scores: use padding of +1 to eliminate bank conflicts when accessing the cache
    __shared__ CUBSM_TYPE
            scoreCache[CUDP_2DCACHE_DIM_D][CUDP_2DCACHE_DIM_X+CUDPIDS];
    __shared__ CUBDP_TYPE
            diag1Cache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_D+1],//cache for state scores of the 1st diagonal
            diag2Cache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_D+1],//cache for state scores of the last diagonal
            bottmCache[nTDPDiagScoreSubsections][CUDP_2DCACHE_DIM_X];//cache for state scores of the bottom of the diagonals
    //__shared__ CUBDP_TYPE
    //        maxscCache[CUDP_2DCACHE_DIM_D];//cache for maximum scores of the last processed diagonal
    CUBDP_TYPE maxscCache;
    //SECTION of cache for correlated scores:
    __shared__ CUBDP_TYPE
            mscodg1Cache[CUDP_CORR_NSCORES][CUDP_2DCACHE_DIM_D],//cache of match scores of the 1st diagonal
            mscodg2Cache[CUDP_CORR_NSCORES][CUDP_2DCACHE_DIM_D],//cache of match scores of the last diagonal
            mscobotCache[CUDP_CORR_NSCORES][CUDP_2DCACHE_DIM_X];//cache of match scores of the block's bottom
    __shared__ CUBDP_TYPE
            //cache of correlation scores of the 1st diagonal
            corrsdg1Cache[CUDP_2DCACHE_DIM_D],
            //cache of correlation scores of the last diagonal
            corrsdg2Cache[CUDP_2DCACHE_DIM_D];
            //cache of correlation scores of the block's bottom
    __shared__ CUBDP_TYPE
            corrsbotCache[CUDP_2DCACHE_DIM_X];
    //use registers instead of SMEM and increase occupancy; needs writing to other threads
    //CUBDP_TYPE corrsbotCache;
    CUBDP_TYPE maxcorrscCache;//correlation score corresponding to the max DP aln score 
    //
    // blockIdx.y is the profile serial number;
    uint dbfldsndx;
    uint pronr;
    //NOTE: protection against overflow ensured on the host side
    if( blockIdx.y < ndb1pros ) { pronr = blockIdx.y + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = blockIdx.y - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( /*threadIdx.y == 0 && */threadIdx.x == 0 ) {
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[pronr];
        dbprolenCache = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DLen]))[pronr];
    }

    // blockIdx.x is block serial number s within diagonal blkdiagnum;
    // (x,y) is the bottom-left corner (x,y) coordinates for profile blockIdx.y
    int x, y;
    if( blkdiagnum <= lastydiagnum) {
        //x=-!(d%2)w+2ws; y=dw/2+w-sw -1 (-1, zero-based indices); [when w==b]
        //(b, block's length; w, block's width)
#ifdef CUDP_2DCACHE_DIM_DequalsX
        x = (2*blockIdx.x - (!(blkdiagnum & 1))) * CUDP_2DCACHE_DIM_D;
        y = ((blkdiagnum>>1) + 1 - blockIdx.x) * CUDP_2DCACHE_DIM_D - 1;
#else
        //x=-w+(d*b)%(w+b)+(w+b)s; y=(d*b/(w+b))w+w-sw -1 (-1, zero-based indices);
//         x = CUDP_2DCACHE_DIM_DpX * blockIdx.x -
//             CUDP_2DCACHE_DIM_D + (blkdiagnum * CUDP_2DCACHE_DIM_X) % CUDP_2DCACHE_DIM_DpX;
//         y = ((blkdiagnum * CUDP_2DCACHE_DIM_X) / CUDP_2DCACHE_DIM_DpX + 1 - blockIdx.x) * 
//             CUDP_2DCACHE_DIM_D - 1;
        //the below is more universal 
        x = (blkdiagnum * CUDP_2DCACHE_DIM_X) % CUDP_2DCACHE_DIM_DpX - CUDP_2DCACHE_DIM_D;
        if( x > 0 ) {
            x -= CUDP_2DCACHE_DIM_DpX;
            y = CUDP_2DCACHE_DIM_D;
        } else
            y = 0;
        x += CUDP_2DCACHE_DIM_DpX * blockIdx.x;
        y += ((blkdiagnum * CUDP_2DCACHE_DIM_X) / CUDP_2DCACHE_DIM_DpX + 1 - blockIdx.x) * 
            CUDP_2DCACHE_DIM_D - 1;
#endif
    } else {
#ifdef CUDP_2DCACHE_DIM_DequalsX
        //x=-w+(d-d_l)w+2ws; y=dw/2+w-sw -1; [when w==b]
        x = (2*blockIdx.x + (blkdiagnum-lastydiagnum-1)) * CUDP_2DCACHE_DIM_D;
        y = ((lastydiagnum>>1) + 1 - blockIdx.x) * CUDP_2DCACHE_DIM_D - 1;
#else
        //x=-w+(d-d_l)w+2ws; y=(d*b/(w+b))w+w-sw -1;
//         x = CUDP_2DCACHE_DIM_DpX * blockIdx.x -
//             CUDP_2DCACHE_DIM_D + (blkdiagnum-lastydiagnum) * CUDP_2DCACHE_DIM_X;
//         y = ((lastydiagnum * CUDP_2DCACHE_DIM_X) / CUDP_2DCACHE_DIM_DpX + 1 - blockIdx.x) * 
//             CUDP_2DCACHE_DIM_D - 1;
        //the below is more universal 
        x = (lastydiagnum * CUDP_2DCACHE_DIM_X) % CUDP_2DCACHE_DIM_DpX - CUDP_2DCACHE_DIM_D;
        if( x > 0 ) {
            x -= CUDP_2DCACHE_DIM_DpX;
            y = CUDP_2DCACHE_DIM_D;
        } else
            y = 0;
        x += CUDP_2DCACHE_DIM_DpX * blockIdx.x + (blkdiagnum-lastydiagnum) * CUDP_2DCACHE_DIM_X;
        y += ((lastydiagnum * CUDP_2DCACHE_DIM_X) / CUDP_2DCACHE_DIM_DpX + 1 - blockIdx.x) * 
            CUDP_2DCACHE_DIM_D - 1;
#endif
    }

    //{{use registers efficiently
    dbprodstCache = __shfl_sync(0xffffffff, dbprodstCache, 0);
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}
    //...or use shared memory; bank conflicts arise when accessed
//     __syncthreads();


    //number of iterations for this block to perform;
    //uint makes the compiler to reduce #registers
    uint ilim = GetMaqxNoIterations( x, y, nqyposs, dbprolenCache );

    if( y < 0 || (int)nqyposs <= (int)(y+1 - blockDim.x) || 
        (int)dbprolenCache <= x /*+ CUDP_2DCACHE_DIM_DpX */ ||
        ilim < 1 )
        //block does not participate: out of profile boundaries
        return;

    int qpos = y - threadIdx.x;//going upwards

    if( 0 <= qpos && qpos < (int)nqyposs )
        // (add querprosOmtd vectors of beginning transitions):
        DPLocCacheTrnProbs( qrtrpCache, 0, querposoffset + qpos + querprosOmtd );
    else
        DPLocInitTrnProbs( qrtrpCache );

    //x is now the position this thread will process
    x += threadIdx.x;

    //db profile position corresponding to the diagonal block's bottom-left 
    // corner in the buffers dc_pm2dvfields_:
    int dbpos = x + dbprodstCache;//going right
    int dblen = ndb1poss + ndbCposs + dbxpad;

    if( 0 <= x && x < (int)dbprolenCache )
        // (add pronr vectors of beginning transitions):
        DPLocCacheTrnProbs( dbtrpCache, 0, dbfldsndx, dbpos + pronr );
    else
        DPLocInitTrnProbs( dbtrpCache, 0 );
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( 0 <= (int)(x+blockDim.x) && (int)(x+blockDim.x) < (int)dbprolenCache )
            DPLocCacheTrnProbs( dbtrpCache, blockDim.x, dbfldsndx, dbpos + pronr + blockDim.x );
        else
            DPLocInitTrnProbs( dbtrpCache, blockDim.x );
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif

    //dbpos is now the x position of the diagonal block's bottom-left 
    // corner in the score matrix plus the offset determined by thread id:
    dbpos = (blockIdx.y < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;

    //cache scores over the oblique diagonal block
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        #pragma unroll 4
        for( int i = 0; i < CUDP_2DCACHE_DIM_D; i++ ) {
            //going upwards
            scoreCache[i][threadIdx.x] = CUBDP_Q(0);
            if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x+i && x+i < (int)dbprolenCache ) {
                //starting position of line i of the oblq. diagonal block in the 
                // score matrix: 
                qpos = (y-i) * dblen + i;
                scoreCache[i][threadIdx.x] = scores[qpos + dbpos];
            }
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif

    //cache two diagonals from the previous (along the x axis) diagonal block;
    //the structure of tmpdpdiagbuffers is position-specific (1D, along x-axis)
    if( 0 <= x-1 && x < (int)dbprolenCache ) {
        DPLocCacheBuffer<CUDP_2DCACHE_DIM_D>( diag1Cache, tmpdpdiagbuffers, dbpos-1, 0, dblen );
        DPLocCacheBuffer<CUDP_2DCACHE_DIM_D>( diag2Cache, tmpdpdiagbuffers, dbpos-1, 
                          dpdssDiag2 * nTDPDiagScoreSubsections * dblen, 
                          dblen, 
                          1/*shift*/);
        //cache the buffer of maximum scores here
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
        //maxscCache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbpos - 1];
        maxscCache = tmpdpdiagbuffers[qpos + dbpos - 1];
        //
        //buffers for calculating CORRELATION SCORES
        //NOTE: indices of the SMEM buffers are one less (to save SMEM)
        qpos += dblen;
        maxcorrscCache = tmpdpdiagbuffers[qpos + dbpos - 1];
        qpos += dblen;
        DPLocCacheCorrBuffer( mscodg1Cache, tmpdpdiagbuffers, dbpos-1, qpos, dblen );
        qpos += CUDP_CORR_NSCORES * dblen;
        DPLocCacheCorrBuffer( mscodg2Cache, tmpdpdiagbuffers, dbpos-1, qpos, dblen );
        qpos += CUDP_CORR_NSCORES * dblen;
        corrsdg1Cache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbpos - 1];
        qpos += dblen;
        corrsdg2Cache[threadIdx.x] = tmpdpdiagbuffers[qpos + dbpos - 1];
    }
    else {
        DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag1Cache);
        DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag2Cache, 1/*shift*/);
        //maxscCache[threadIdx.x] = CUBDP_Q(0);
        maxscCache = CUBDP_Q(0);
        //buffers for calculating CORRELATION SCORES
        DPLocInitCorrCache(mscodg1Cache);
        DPLocInitCorrCache(mscodg2Cache);
        corrsdg1Cache[threadIdx.x] = CUBDP_Q(0);
        corrsdg2Cache[threadIdx.x] = CUBDP_Q(0);
        maxcorrscCache = CUBDP_Q(0);
    }

    //cache the bottom of the above diagonal blocks;
    //the structure of tmpdpbotbuffer is position-specific (1D, along x-axis)
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( CUDP_2DCACHE_DIM_D <= y && 
            0 <= x+CUDP_2DCACHE_DIM_D-1 && x+CUDP_2DCACHE_DIM_D-1 < (int)dbprolenCache ) {
            DPLocCacheBuffer<CUDP_2DCACHE_DIM_X>( 
                bottmCache, tmpdpbotbuffer, dbpos+CUDP_2DCACHE_DIM_D-1, 0, 
                            dblen );
            //buffers for calculating CORRELATION SCORES
            qpos = nTDPDiagScoreSubsections * dblen;
            DPLocCacheCorrBuffer( mscobotCache, tmpdpbotbuffer, dbpos+CUDP_2DCACHE_DIM_D-1, 
                              qpos, dblen );
            qpos += CUDP_CORR_NSCORES * dblen;
            corrsbotCache[threadIdx.x] = tmpdpbotbuffer[qpos + dbpos+CUDP_2DCACHE_DIM_D-1];
            //corrsbotCache = tmpdpbotbuffer[qpos + dbpos+CUDP_2DCACHE_DIM_D-1];
        }
        else {
            DPLocInitCache<CUDP_2DCACHE_DIM_X>(bottmCache);
            //buffers for calculating CORRELATION SCORES
            DPLocInitCorrCache(mscobotCache);
            corrsbotCache[threadIdx.x] = CUBDP_Q(0);
            //corrsbotCache = CUBDP_Q(0);
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif


    __syncthreads();


    CUBDP_TYPE (*pdiag1)[CUDP_2DCACHE_DIM_D+1] = diag1Cache;
    CUBDP_TYPE (*pdiag2)[CUDP_2DCACHE_DIM_D+1] = diag2Cache;

    CUBDP_TYPE (*mscdg1)[CUDP_2DCACHE_DIM_D] = mscodg1Cache;
    CUBDP_TYPE (*mscdg2)[CUDP_2DCACHE_DIM_D] = mscodg2Cache;

    CUBDP_TYPE* cosdg1 = corrsdg1Cache;
    CUBDP_TYPE* cosdg2 = corrsdg2Cache;


    //start calculations for this position with 32x unrolling
    //NOTE: sync inside: do not branch;
    for( int i = 0; i < ilim/*CUDP_2DCACHE_DIM_X*/; i++ ) {
        CUBDP_TYPE val1, val2, corr;
        if( threadIdx.x == blockDim.x-1 ) {
            pdiag1[dpdsssStateMM][threadIdx.x+1] = bottmCache[dpdsssStateMM][i];
            pdiag1[dpdsssStateMI][threadIdx.x+1] = bottmCache[dpdsssStateMI][i];
            pdiag1[dpdsssStateIM][threadIdx.x+1] = bottmCache[dpdsssStateIM][i];
            pdiag1[dpdsssStateDN][threadIdx.x+1] = bottmCache[dpdsssStateDN][i];
            pdiag1[dpdsssStateND][threadIdx.x+1] = bottmCache[dpdsssStateND][i];
            //
            //indices of diagonal buffers (mscdg1 and cosdg1) are one less
            mscdg1[0][threadIdx.x] = mscobotCache[0][i];
            mscdg1[1][threadIdx.x] = mscobotCache[1][i];
            mscdg1[2][threadIdx.x] = mscobotCache[2][i];
            mscdg1[3][threadIdx.x] = mscobotCache[3][i];
            cosdg1[threadIdx.x] = corrsbotCache[i];
        }

        //MM state update (diagonal direction)
        val1 = pdiag2[dpdsssStateMI][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIMp][threadIdx.x+i];
        val2 = pdiag2[dpdsssStateIM][threadIdx.x+1] + 
            qrtrpCache[dptrIMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateDN][threadIdx.x+1] + 
            qrtrpCache[dptrDMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateND][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrDMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        val2 = pdiag2[dpdsssStateMM][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        if( val1 == val2 )
            //mark that the MM state has a maximum value
            //(reuse val2)
            val2 = CUDP_FLOAT_DELM;
        //accessing scores is without bank conflicts because of padding
        val1 = myhdmax/*fmaxf*/( val1 + scoreCache[threadIdx.x][i], CUBDP_Q(0));
        __syncthreads();//last point pdiag2 is used for reading in this iteration


        //{{calculate CORRELATION scores
        if( val1 == CUBDP_Q(0)) {
            if( threadIdx.x ) {
                DPLocInitCorrCache(mscdg2, -1/*shift*/);//WRITE
                cosdg2[threadIdx.x-1] = CUBDP_Q(0);//WRITE
            } else {
                mscobotCache[0][i] = mscobotCache[1][i] =
                mscobotCache[2][i] = mscobotCache[3][i] = CUBDP_Q(0);//WRITE
                corrsbotCache[i] = CUBDP_Q(0);//WRITE
            }
        }
        else if( val2 == CUDP_FLOAT_DELM ) {
            //MM state has a maximum score;
            //update match score buffers and calculate correlation score;
            //NOTE: match and correlation score diagonals represent vectors with 
            // indices one greater!
            //cosdg2 will be bnecome cosdg1 in the next iteration
            corr = cosdg2[threadIdx.x] + scoreCache[threadIdx.x][i] *
                    ( mscdg2[0][threadIdx.x] + mscdg2[1][threadIdx.x] +
                      mscdg2[2][threadIdx.x] + mscdg2[3][threadIdx.x] );//WRITE
        }
        //sync in case of update of correlation scores;
        //avoid using divergent code here
        __syncthreads();
        if( val1 && val2 == CUDP_FLOAT_DELM ) {
            if( threadIdx.x ) {
                cosdg2[threadIdx.x-1] = corr;//WRITE
                mscdg2[0][threadIdx.x-1] = mscdg2[1][threadIdx.x];//WRITE
                mscdg2[1][threadIdx.x-1] = mscdg2[2][threadIdx.x];//WRITE
                mscdg2[2][threadIdx.x-1] = mscdg2[3][threadIdx.x];//WRITE
                mscdg2[3][threadIdx.x-1] = scoreCache[threadIdx.x][i];//WRITE
            } else {
                corrsbotCache[i] = corr;//WRITE
                mscobotCache[0][i] = mscdg2[1][threadIdx.x];//WRITE
                mscobotCache[1][i] = mscdg2[2][threadIdx.x];//WRITE
                mscobotCache[2][i] = mscdg2[3][threadIdx.x];//WRITE
                mscobotCache[3][i] = scoreCache[threadIdx.x][i];//WRITE
            }
        }
        //}}


        //__syncthreads();//last point pdiag2 is used for reading in this iteration
        pdiag2[dpdsssStateMM][threadIdx.x] = val1;//WRITE
        //maxscCache[threadIdx.x] = myhdmax/*fmaxf*/( maxscCache[threadIdx.x], val1 );
        //maxscCache = myhdmax/*fmaxf*/( maxscCache, val1 );
        myhdmaxassgn(maxscCache, val1, maxcorrscCache, corr);


        //MI state update (up direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMIc][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateMI][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIIc][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateMI][threadIdx.x] = val1;//WRITE

        //IM state update (left direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateIM][threadIdx.x] + 
            qrtrpCache[dptrIIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateIM][threadIdx.x] = val1;//WRITE

        //DN state update (up)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + qrtrpCache[dptrMDp][threadIdx.x];
        val2 = pdiag1[dpdsssStateDN][threadIdx.x+1] + qrtrpCache[dptrDDp][threadIdx.x];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateDN][threadIdx.x] = val1;//WRITE

        //ND state update (left)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + dbtrpCache[dptrMDp][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateND][threadIdx.x] + dbtrpCache[dptrDDp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateND][threadIdx.x] = val1;//WRITE

        if( threadIdx.x == 0 ) {
            //WRITE
            //this position is not used by other threads in the current iteration
            bottmCache[dpdsssStateMM][i] = pdiag2[dpdsssStateMM][threadIdx.x];
            bottmCache[dpdsssStateMI][i] = pdiag2[dpdsssStateMI][threadIdx.x];
            bottmCache[dpdsssStateIM][i] = pdiag2[dpdsssStateIM][threadIdx.x];
            bottmCache[dpdsssStateDN][i] = pdiag2[dpdsssStateDN][threadIdx.x];
            bottmCache[dpdsssStateND][i] = pdiag2[dpdsssStateND][threadIdx.x];
        }

#ifdef CUDP_CORR_INIT_TESTPRINT
        if(pronr==3024/*4899*/){
            printf(" d=%u(%u) s=%u i%02d/%u (t%02u): len= %d addr= %u SC= %.4f (yx: %d,%d) "
                    "MM= %.6f MI= %.6f IM= %.6f DN= %.6f ND= %.6f  MAX= %.6f\n"
                    "  >qMM= %.4f qMIc= %.4f qMD= %.4f qIM= %.4f qIIc= %.4f qDM= %.4f qDD= %.4f\n"
                    "  >dMM= %.4f dMIc= %.4f dMD= %.4f dIM= %.4f dIIc= %.4f dDM= %.4f dDD= %.4f\n",
                    blkdiagnum,lastydiagnum,blockIdx.x,i,ilim,threadIdx.x,
                    dbprolenCache,dbprodstCache,scoreCache[threadIdx.x][i], y-threadIdx.x,x+i,
                    pdiag2[dpdsssStateMM][threadIdx.x],
                    pdiag2[dpdsssStateMI][threadIdx.x],
                    pdiag2[dpdsssStateIM][threadIdx.x],
                    pdiag2[dpdsssStateDN][threadIdx.x],
                    pdiag2[dpdsssStateND][threadIdx.x],
                    //maxscCache[threadIdx.x],
                    maxscCache,
                    qrtrpCache[dptrMMp][threadIdx.x],qrtrpCache[dptrMIc][threadIdx.x],
                    qrtrpCache[dptrMDp][threadIdx.x],qrtrpCache[dptrIMp][threadIdx.x],
                    qrtrpCache[dptrIIc][threadIdx.x],qrtrpCache[dptrDMp][threadIdx.x],
                    qrtrpCache[dptrDDp][threadIdx.x],
                    dbtrpCache[dptrMMp][threadIdx.x+i],dbtrpCache[dptrMIc][threadIdx.x+i],
                    dbtrpCache[dptrMDp][threadIdx.x+i],dbtrpCache[dptrIMp][threadIdx.x+i],
                    dbtrpCache[dptrIIc][threadIdx.x+i],dbtrpCache[dptrDMp][threadIdx.x+i],
                    dbtrpCache[dptrDDp][threadIdx.x+i]
            );
            for(size_t _k=0;_k<1000000000UL;_k++)clock();
            for(size_t _k=0;_k<10000000UL;_k++)clock();
        }
#endif

        myhdswap( pdiag1, pdiag2 );
        //swap data for correlation scores
        myhdswap( mscdg1, mscdg2 );
        myhdswap( cosdg1, cosdg2 );

        __syncthreads();
    }


//     //XOR mode warp sync for butterfly reduction;
//     //NOTE: too expensive to perform for each block; 
//     // this synchronization is only needed for the last block;
//     // thus, skip this step altogether;
//     for(int i = 16/*warpSize>>1*/; 1 <= i; i /= 2 )
//         maxscCache = myhdmax( maxscCache, __shfl_xor_sync(0xffffffff, maxscCache, i/*, 32*//*warpsize*/));

#ifdef CUDP_CORR_INIT_TESTPRINT
        if(pronr==3024/*4899*/)
            printf(" >>> d=%u(%u) s=%u (t%02u): MAX= %.6f\n", 
                    blkdiagnum,lastydiagnum,blockIdx.x,threadIdx.x, maxscCache );
#endif


    //write the result of calculations for following blocks;
    //write two diagonals;
    if( 0 <= x+CUDP_2DCACHE_DIM_X-1 && x+CUDP_2DCACHE_DIM_X-1 < (int)dbprolenCache ) {
        DPLocWriteBuffer<CUDP_2DCACHE_DIM_D>( diag1Cache, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, 0, 
                          dblen );
        DPLocWriteBuffer<CUDP_2DCACHE_DIM_D>( diag2Cache, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, 
                          dpdssDiag2 * nTDPDiagScoreSubsections * dblen, 
                          dblen,
                          1/*shift*/);
        //write the buffer of maximum scores
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
        //tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxscCache[threadIdx.x];
        if( CUDP_2DCACHE_DIM_D <= y ) {
            if( tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] < maxscCache ) {
                tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxscCache;
                qpos += dblen;
                tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxcorrscCache;
            }
        } else {
            tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxscCache;
            qpos += dblen;
            tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = maxcorrscCache;
        }
        //write buffers related to CORRELATION SCORES
        qpos += dblen;
        DPLocWriteCorrBuffer( mscodg1Cache, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, qpos, dblen );
        qpos += CUDP_CORR_NSCORES * dblen;
        DPLocWriteCorrBuffer( mscodg2Cache, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, qpos, dblen );
        qpos += CUDP_CORR_NSCORES * dblen;
        tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = corrsdg1Cache[threadIdx.x];
        qpos += dblen;
        tmpdpdiagbuffers[qpos + dbpos+CUDP_2DCACHE_DIM_X-1] = corrsdg2Cache[threadIdx.x];
    }

    //write the bottom of the diagonal block;
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( 0 <= x && x < (int)dbprolenCache ) {
            DPLocWriteBuffer<CUDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer, dbpos, 0, dblen );
            //write buffers related to CORRELATION SCORES
            qpos = nTDPDiagScoreSubsections * dblen;
            DPLocWriteCorrBuffer( mscobotCache, tmpdpbotbuffer, dbpos, qpos, dblen );
            qpos += CUDP_CORR_NSCORES * dblen;
            tmpdpbotbuffer[qpos + dbpos] = corrsbotCache[threadIdx.x];
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif
}
