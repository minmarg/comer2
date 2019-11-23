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
#include "CuBatchDP_init_btck.cuh"

// #define CUDP_INIT_BTCK_TESTPRINT 0 //3024//4899

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming with backtracking 
// information;
// NOTE: Version for thread block of one warp!
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
// ExecDP_Btck_Unroll32x: device code for executing dynamic programming with 
// backtracking information using shared memory and 32-fold unrolling 
// along the diagonal of dimension CUDP_2DCACHE_DIM;
// NOTE: memory pointers should be aligned!
// scores, calculated scores used as input;
// tmpdpdiagbuffers, temporary buffers for last calculated diagonal scores;
// tmpdpbotbuffer, temporary buffers for last calculated bottom scores;
// maxscoordsbuf, coordinates of maximum alignment scores;
// btckdata, backtracking information data;
// 
__global__ void ExecDP_Btck_Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata )
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
    //SECTION for backtracking information
    uint maxscCoords;//coordinates of the maximum alignment score maxscCache
    __shared__ char btckCache[CUDP_2DCACHE_DIM_D][CUDP_2DCACHE_DIM_X+CUDPIDS];//backtracking cache
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

    //{{NOTE: do not align the distant corners;
    //get the edge length of a triangular area to be ignored:
    int crnedgelen = GetCornerMaxlen(nqyposs, dbprolenCache);
    //}}



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
            scoreCache[i][threadIdx.x] = CUBDP_Q(-9);
            if( CellYXinValidArea(nqyposs, dbprolenCache, crnedgelen, y-i, x+i )) {
            //if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x+i && x+i < (int)dbprolenCache ) {
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
        maxscCoords = maxscoordsbuf[dbpos-1];
    }
    else {
        DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag1Cache);
        DPLocInitCache<CUDP_2DCACHE_DIM_D>(diag2Cache, 1/*shift*/);
        //maxscCache[threadIdx.x] = CUBDP_Q(0);
        maxscCache = CUBDP_Q(0);
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
        }
        else {
            DPLocInitCache<CUDP_2DCACHE_DIM_X>(bottmCache);
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif


    __syncthreads();

    CUBDP_TYPE (*pdiag1)[CUDP_2DCACHE_DIM_D+1] = diag1Cache;
    CUBDP_TYPE (*pdiag2)[CUDP_2DCACHE_DIM_D+1] = diag2Cache;

    //start calculations for this position with 32x unrolling
    //NOTE: sync inside: do not branch;
    for( int i = 0; i < ilim/*CUDP_2DCACHE_DIM_X*/; i++ ) {
        CUBDP_TYPE val1, val2, max;
        int btck = dpbtckSTOP;
        //
        if( threadIdx.x == blockDim.x-1 ) {
            pdiag1[dpdsssStateMM][threadIdx.x+1] = bottmCache[dpdsssStateMM][i];
            pdiag1[dpdsssStateMI][threadIdx.x+1] = bottmCache[dpdsssStateMI][i];
            pdiag1[dpdsssStateIM][threadIdx.x+1] = bottmCache[dpdsssStateIM][i];
            pdiag1[dpdsssStateDN][threadIdx.x+1] = bottmCache[dpdsssStateDN][i];
            pdiag1[dpdsssStateND][threadIdx.x+1] = bottmCache[dpdsssStateND][i];
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
        //accessing scores is without bank conflicts because of padding
        val1 = myhdmax/*fmaxf*/( val1 + scoreCache[threadIdx.x][i], CUBDP_Q(0));
        max = val1;
        if(max) btck = dpbtckDIAG;

        __syncthreads();//last point pdiag2 is used for reading in this iteration
        pdiag2[dpdsssStateMM][threadIdx.x] = val1;//WRITE
        //maxscCache[threadIdx.x] = myhdmax/*fmaxf*/( maxscCache[threadIdx.x], val1 );
        //maxscCache = myhdmax/*fmaxf*/( maxscCache, val1 );
        //NOTE: coord. overflow will be managed below while writing data
        dpmaxandcoords( maxscCache, val1, maxscCoords, x+i, y-threadIdx.x );

        //MI state update (up direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMIc][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateMI][threadIdx.x+1] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIIc][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateMI][threadIdx.x] = val1;//WRITE
        myhdmaxassgn( max, val1, btck, (int)dpbtckUP );

        //IM state update (left direction)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateIM][threadIdx.x] + 
            qrtrpCache[dptrIIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateIM][threadIdx.x] = val1;//WRITE
        myhdmaxassgn( max, val1, btck, (int)dpbtckLEFT );

        //DN state update (up)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x+1] + qrtrpCache[dptrMDp][threadIdx.x];
        val2 = pdiag1[dpdsssStateDN][threadIdx.x+1] + qrtrpCache[dptrDDp][threadIdx.x];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateDN][threadIdx.x] = val1;//WRITE
        myhdmaxassgn( max, val1, btck, (int)dpbtckUP );

        //ND state update (left)
        val1 = pdiag1[dpdsssStateMM][threadIdx.x] + dbtrpCache[dptrMDp][threadIdx.x+i];
        val2 = pdiag1[dpdsssStateND][threadIdx.x] + dbtrpCache[dptrDDp][threadIdx.x+i];
        val1 = myhdmax/*fmaxf*/( val1, val2 );
        pdiag2[dpdsssStateND][threadIdx.x] = val1;//WRITE
        myhdmaxassgn( max, val1, btck, (int)dpbtckLEFT );

        btckCache[threadIdx.x][i] = btck;

        if( threadIdx.x == 0 ) {
            //WRITE
            //this position is not used by other threads in the current iteration
            bottmCache[dpdsssStateMM][i] = pdiag2[dpdsssStateMM][threadIdx.x];
            bottmCache[dpdsssStateMI][i] = pdiag2[dpdsssStateMI][threadIdx.x];
            bottmCache[dpdsssStateIM][i] = pdiag2[dpdsssStateIM][threadIdx.x];
            bottmCache[dpdsssStateDN][i] = pdiag2[dpdsssStateDN][threadIdx.x];
            bottmCache[dpdsssStateND][i] = pdiag2[dpdsssStateND][threadIdx.x];
        }

#ifdef CUDP_INIT_BTCK_TESTPRINT
        if(pronr==CUDP_INIT_BTCK_TESTPRINT){
            printf(" d=%u(%u) s=%u i%02d/%u (t%02u): len= %d addr= %u SC= %.4f (yx: %d,%d) "
                    "MM= %.6f MI= %.6f IM= %.6f DN= %.6f ND= %.6f  MAX= %.6f COORD= %x\n"// BTCK= %d\n"
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
                    maxscCache, maxscCoords,// btck,
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
            for(size_t _k=0;_k<1000000000UL;_k++)clock();
        }
#endif

        myhdswap( pdiag1, pdiag2 );

        __syncthreads();
    }


//     //XOR mode warp sync for butterfly reduction;
//     //NOTE: too expensive to perform for each block; 
//     // this synchronization is only needed for the last block;
//     // thus, skip this step altogether;
//     for(int i = 16/*warpSize>>1*/; 1 <= i; i >>= 1 )
//         maxscCache = myhdmax( maxscCache, __shfl_xor_sync(0xffffffff, maxscCache, i/*, 32*//*warpsize*/));

#ifdef CUDP_INIT_BTCK_TESTPRINT
        if(pronr==CUDP_INIT_BTCK_TESTPRINT)
            printf(" >>> d=%u(%u) s=%u (t%02u): len= %d MAX= %.6f COORD= %x (yx %u %u) wrt= %d xpos= %d\n", 
                    blkdiagnum,lastydiagnum,blockIdx.x,threadIdx.x, dbprolenCache,
                    maxscCache, maxscCoords, GetCoordY(maxscCoords),GetCoordX(maxscCoords),
                    x+CUDP_2DCACHE_DIM_X-1<(int)dbprolenCache, x+CUDP_2DCACHE_DIM_X-1);
#endif


    //write the result of calculations for following blocks;
    //write two diagonals;
    if( 0 <= x+CUDP_2DCACHE_DIM_X-1 && x+CUDP_2DCACHE_DIM_X-1 < (int)dbprolenCache ) {
        DPLocWriteBuffer<CUDP_2DCACHE_DIM_D>( pdiag1/*diag1Cache*/, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, 0, 
                          dblen );
        DPLocWriteBuffer<CUDP_2DCACHE_DIM_D>( pdiag2/*diag2Cache*/, tmpdpdiagbuffers, dbpos+CUDP_2DCACHE_DIM_X-1, 
                          dpdssDiag2 * nTDPDiagScoreSubsections * dblen, 
                          dblen,
                          1/*shift*/);
    }

    int ndx = CUDP_2DCACHE_DIM_X-1;
    int xtl = x-threadIdx.x+blockDim.x-1;//top-left x coordinate of the diagonal block
    if( xtl >= dbprolenCache )
        ndx = 0;
    else if( xtl+ndx >= dbprolenCache )
        ndx = dbprolenCache - xtl - 1;

    //write the buffer of maximum scores
    if( 0 <= x+ndx && x+ndx < (int)dbprolenCache ) {
//         //{{NOTE:no need to use the check as the scores in the triangle areas are highly negative.
//         //if((uint)dbprolenCache <= GetCoordX(maxscCoords) || 
//         //    nqyposs <= GetCoordY(maxscCoords)) {
//         if(!CellYXinValidArea(nqyposs, dbprolenCache, crnedgelen, 
//                 GetCoordY(maxscCoords), GetCoordX(maxscCoords))) {
//             maxscCache = CUBDP_Q(0);
//             maxscCoords = 0;
//         }
//         //}}
        qpos = dpdssDiagM * nTDPDiagScoreSubsections * dblen;
        //tmpdpdiagbuffers[qpos + dbpos+ndx] = maxscCache[threadIdx.x];
        if( CUDP_2DCACHE_DIM_D <= y ) {
            if( tmpdpdiagbuffers[qpos + dbpos+ndx] < maxscCache ) {
                tmpdpdiagbuffers[qpos + dbpos+ndx] = maxscCache;
                maxscoordsbuf[dbpos+ndx] = maxscCoords;
            }
        } else {
            tmpdpdiagbuffers[qpos + dbpos+ndx] = maxscCache;
            maxscoordsbuf[dbpos+ndx] = maxscCoords;
        }
    }

    //write the bottom of the diagonal block;
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        if( 0 <= x && x < (int)dbprolenCache )
            DPLocWriteBuffer<CUDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer, dbpos, 0, dblen );
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif

    //write backtracking information
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    if( threadIdx.x < CUDP_2DCACHE_DIM_X ) {
#endif
        //unroll since #registers is not a problem over all compute capabilities
        #pragma unroll 8
        for( int i = 0; i < CUDP_2DCACHE_DIM_D; i++ ) {
            //going upwards
            if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x+i && x+i < (int)dbprolenCache ) {
                //starting position of line i of the oblq. diagonal block in the matrix:
                qpos = (y-i) * dblen + i;
                btckdata[qpos + dbpos] = btckCache[i][threadIdx.x];
            }
        }
#if !defined(CUDP_2DCACHE_DIM_DequalsX)
    }
#endif
}
