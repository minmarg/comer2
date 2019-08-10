/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "libmycu/cudp/CuBatchDP_init.cuh"
#include "libmycu/cuss/CuBatchSS_com.h"
#include "CuBatchMAPDP_fwd.cuh"
#include "CuBatchMAPDP_bwd.cuh"

// #define CUMAPDP_BWD_TESTPRINT 1682 //27 //0 //1286 //0 //3024//4899

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for executing dynamic programming with backtracking 
// information;
// NOTE: Version for thread block of one warp!
// blkdiagnum, block diagonal serial number;
// lastydiagnum, last block diagonal serial number along y axis 
// (starting at x=-CUMAPDP_2DCACHE_DIM);
// logevthld, log e-value threshold;
// modscoresinuse, indication of using the modular scores;
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
// dblen2, total number of db profile positions + padding in phase 2;
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
// scores, calculated scores used as input;
// modscores, input modular scores;
// dp2alndatbuffers, profile-specific statistics from (previous) phase 1;
// tmpdpdiagbuffers, temporary buffers for last calculated diagonal scores;
// tmpdpbotbuffer, temporary buffers for last calculated bottom scores;
// tmpfwdprobs, forward probabilities used to calculate partition function 
// (total probability) while executing the backward calculations;
// // tmpbwdprobs, output backward probabilities to be calculated;
// tmpss2datbuffers, temporary buffer for updating the total probability for 
// each pair of profiles;
// NOTE: 13824B SMEM is a limit of reduced occupancy
// 
__global__ void ExecMAPDP_Bwd_Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    float logevthld,
    bool modscoresinuse,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint /*bdbCposoffset*/,
    uint dblen2,
    uint dbxpad2,
    const CUBSM_TYPE* __restrict__ scores,
    const CUBSM_TYPE* __restrict__ modscores,
    const float* __restrict__ dp2alndatbuffers,
    float* __restrict__ tmpdpdiagbuffers,
    float* __restrict__ tmpdpbotbuffer,
    float* __restrict__ tmpfwdprobs,
//float* __restrict__ tmpbwdprobs,
    float* __restrict__ tmpss2datbuffers )
{
    /*__shared__ */LNTYPE dbprodstCache;//distance in positions to db profile
    /*__shared__ */INTYPE dbprolenCache;//length of profile
    __shared__ FPTYPE 
            qrtrpCache[nTDPUsedProTrn][CUMAPDP_2DCACHE_DIM_D],//cache for query transition probabilities
            dbtrpCache[nTDPUsedProTrn][CUMAPDP_2DCACHE_DIM_DpX];//cache for db transition probabilities
    //cache of scores: use padding of +1 to eliminate bank conflicts when accessing the cache
    __shared__ float
            scoreCache[CUMAPDP_2DCACHE_DIM_D][CUMAPDP_2DCACHE_DIM_X+CUDPIDS];
    __shared__ float
            diag1Cache[nTDPDiagScoreSubsections][CUMAPDP_2DCACHE_DIM_D+1],//cache for state scores of the 1st diagonal
            diag2Cache[nTDPDiagScoreSubsections][CUMAPDP_2DCACHE_DIM_D+1],//cache for state scores of the last diagonal
            bottmCache[nTDPDiagScoreSubsections][CUMAPDP_2DCACHE_DIM_X];//cache for state scores of the top of the diagonals
    //
    //maximum log forward probability
    float maxlogprob;
    //using SMEM (even with bank conflicts) works faster than using (with constant indices) and 
    // synchronizing registers here (small number of accesses):
    __shared__ uint proattr[2];//original profile NUMBER and new DISTANCE from the beginning
    __shared__ float evals[2];//log of e-value and pair e-value
    uint dbfldsndx;
    uint pronr;


    // blockIdx.y is the profile serial number in phase 2;
    //{{get the original profile number and e-values too
    if( /*threadIdx.y == 0 && */threadIdx.x < 2 ) {
        evals[threadIdx.x] = dp2alndatbuffers[nTDP2OutputAlnData*blockIdx.y+dp2oadEvalue+threadIdx.x];
        proattr[threadIdx.x] = 
            *(uint*)(dp2alndatbuffers + nTDP2OutputAlnData*blockIdx.y+dp2oadOrgProNo+threadIdx.x);
    }
//     proattr[0] = __shfl_sync(0xffffffff, proattr[0], 0);
//     proattr[1] = __shfl_sync(0xffffffff, proattr[1], 1);
//     evals[0] = __shfl_sync(0xffffffff, evals[0], 0);
//     evals[1] = __shfl_sync(0xffffffff, evals[1], 1);
    //}}

    __syncthreads();

    //all threads exit if calculated e-value exceeds the threshold
    if( logevthld < evals[0])
        return;


    //NOTE: protection against overflow ensured on the host side
    if( proattr[0] < ndb1pros ) { pronr = proattr[0] + ndb1prosOmtd;
                dbfldsndx = pmv2DTotFlds;
    } else {    pronr = proattr[0] - ndb1pros + ndbCprosOmtd;//jump to section ndbCposs
                //bdb1posoffset = bdbCposoffset;
                dbfldsndx = TIMES2(pmv2DTotFlds);
    }

    if( /*threadIdx.y == 0 && */threadIdx.x == 0 ) {
        dbprodstCache = ((LNTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DDist]))[pronr];
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
    dbprodstCache = __shfl_sync(0xffffffff, dbprodstCache, 0);
    dbprolenCache = __shfl_sync(0xffffffff, dbprolenCache, 0);
    //}}
    //...or use shared memory; bank conflicts arise when accessed
//     __syncthreads();


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

    int qpos = y - threadIdx.x;//going upwards



    if( 0 <= qpos && qpos < (int)nqyposs )
        // (add querprosOmtd vectors of beginning transitions):
        MAPDP_Bwd_LocCacheTrnProbs( qrtrpCache, 0, querposoffset + qpos + querprosOmtd );
    else
        MAPDPLocInitTrnProbs( qrtrpCache );

    //x is now the position this thread will process
    x += threadIdx.x;

    //db profile position corresponding to the diagonal block's bottom-left 
    // corner in the buffers dc_pm2dvfields_:
    int dbpos = x + dbprodstCache;//going right
    int dblen = ndb1poss + ndbCposs + dbxpad;

    if( 0 <= x && x < (int)dbprolenCache )
        // (add pronr vectors of beginning transitions):
        MAPDP_Bwd_LocCacheTrnProbs( dbtrpCache, 0, dbfldsndx, dbpos + pronr );
    else
        MAPDPLocInitTrnProbs( dbtrpCache, 0 );

    if( 0 <= (int)(x+blockDim.x) && (int)(x+blockDim.x) < (int)dbprolenCache )
        MAPDP_Bwd_LocCacheTrnProbs( dbtrpCache, blockDim.x, dbfldsndx, dbpos + pronr + blockDim.x );
    else
        MAPDPLocInitTrnProbs( dbtrpCache, blockDim.x );



    //dbpos is now the x position of the diagonal block's bottom-left 
    // corner in the score matrix plus the offset determined by thread id:
    dbpos = (proattr[0] < ndb1pros)? dbpos - bdb1posoffset: dbpos + ndb1poss;
    //dbpos2 is the x position of the diagonal block's bottom-left 
    // corner in the phase-2 diagonal and bottom structures plus the offset 
    // determined by thread id:
    //(NOTE: use int because of an overflow check required on every chenge of x)
    int dbpos2 = x + proattr[1];



    //flag of the first block diagonal being processed (in reversed order):
    bool firstprocdiagonal = (
        blkdiagnum+1 == 
        GetNoBlockDiagonalsProfile(nqyposs, dbprolenCache, CUMAPDP_2DCACHE_DIM_X, CUMAPDP_2DCACHE_DIM_D)
    );
    //reduce or read the max log forward probability
    maxlogprob = MAPDP_GetMaxLogFwdProbability(
        firstprocdiagonal,
        blockIdx.y/*pronr2*/,
        dbprolenCache,
        dblen2,
        proattr[1],
        tmpdpdiagbuffers,
        tmpss2datbuffers );
    //update sum of exponentials for computing log of the total probability
    MAPDP_UpdateForTotalProbability(
        blockIdx.y,
        nqyposs,
        dbprolenCache,
        dblen2,
        x, y,
        dbpos2,
        maxlogprob,
        tmpfwdprobs,
        tmpss2datbuffers );



    //dynamic score offset and multiplier
    float dyno, dmul;
    MAPDPGetDynScoreParams( evals[1], &dyno, &dmul );

    //cache scores over the oblique diagonal block
    #pragma unroll 4
    for( int i = 0; i < CUMAPDP_2DCACHE_DIM_D; i++ ) {
        //going upwards
        scoreCache[i][threadIdx.x] = MAPDPLocGetLogProb0();
        //NOTE: cache a score block SHIFTED right and down
        if( 0 <= y+1-i && y+1-i < (int)nqyposs && 0 <= x+1+i && x+1+i < (int)dbprolenCache ) {
            //starting position of line i of the oblq. diagonal block in the 
            // score matrix: 
            qpos = (y+1-i) * dblen + i+1;
            scoreCache[i][threadIdx.x] = (
                (dmul && modscoresinuse) ? 
                    scores[qpos+dbpos] - dyno + modscores[qpos+dbpos] * dmul
                :   scores[qpos+dbpos] - dyno
            );
        }
    }



    //cache two diagonals from the previous (i.e., farther along the x axis) diagonal block;
    //the structure of tmpdpdiagbuffers is position-specific (1D, along x-axis)
    if( !firstprocdiagonal &&
        0 <= x+CUMAPDP_2DCACHE_DIM_X && x+CUMAPDP_2DCACHE_DIM_X < (int)dbprolenCache ) {
        DPLocCacheBuffer<CUMAPDP_2DCACHE_DIM_D>( diag1Cache, tmpdpdiagbuffers,
                          dbpos2+CUMAPDP_2DCACHE_DIM_X, 0, dblen2, 1/*shift*/);
        DPLocCacheBuffer<CUMAPDP_2DCACHE_DIM_D>( diag2Cache, tmpdpdiagbuffers,
                          dbpos2+CUMAPDP_2DCACHE_DIM_X,
                          dpdssDiag2 * nTDPDiagScoreSubsections * dblen2,
                          dblen2,
                          1/*shift*/);
    }
    else {
        DPLocInitCache<CUMAPDP_2DCACHE_DIM_D>(diag1Cache, 1/*shift*/, MAPDPLocGetLogProb0());
        DPLocInitCache<CUMAPDP_2DCACHE_DIM_D>(diag2Cache, 1/*shift*/, MAPDPLocGetLogProb0());
    }


    //cache the top of the below diagonal blocks;
    //the structure of tmpdpbotbuffer is position-specific (1D, along x-axis)
    if( threadIdx.x == 0 ) {
        //NOTE: cache one last top value into the diag2Cache buffer!
        if( !firstprocdiagonal &&
            y+1 < (int)nqyposs && 
            0 <= x+CUMAPDP_2DCACHE_DIM_X && x+CUMAPDP_2DCACHE_DIM_X < (int)dbprolenCache )
            DPLocCacheBuffer<CUMAPDP_2DCACHE_DIM_D>( diag2Cache, tmpdpbotbuffer,
                            dbpos2+CUMAPDP_2DCACHE_DIM_X,
                            0,
                            dblen2, 
                            0/*shift*/);
        else
            DPLocInitCache<CUMAPDP_2DCACHE_DIM_D>(diag2Cache, 0/*shift*/, MAPDPLocGetLogProb0());
    }
    if( y+1 < (int)nqyposs && 0 <= x && x < (int)dbprolenCache ) {
        DPLocCacheBuffer<CUMAPDP_2DCACHE_DIM_X>( bottmCache, tmpdpbotbuffer,
                        dbpos2, 0, dblen2 );
    }
    else {
        DPLocInitCache<CUMAPDP_2DCACHE_DIM_X>(bottmCache, MAPDPLocGetLogProb0());
    }


    __syncthreads();


    float (*pdiag1)[CUMAPDP_2DCACHE_DIM_D+1] = diag1Cache;
    float (*pdiag2)[CUMAPDP_2DCACHE_DIM_D+1] = diag2Cache;

    //start calculations for this position with 32x unrolling
    //NOTE: sync inside: do not branch;
    for( int i = ilim-1; i >= 0; i-- ) {
        float val1, val2, val3, val4, val5, max, sco;
        //
        if( threadIdx.x == 0 ) {
            //indices of pdiag1 and pdiag2 are SHIFTED down!
            pdiag1[dpdsssStateMM][threadIdx.x] = bottmCache[dpdsssStateMM][i];
            pdiag1[dpdsssStateMI][threadIdx.x] = bottmCache[dpdsssStateMI][i];
            pdiag1[dpdsssStateIM][threadIdx.x] = bottmCache[dpdsssStateIM][i];
            pdiag1[dpdsssStateDN][threadIdx.x] = bottmCache[dpdsssStateDN][i];
            pdiag1[dpdsssStateND][threadIdx.x] = bottmCache[dpdsssStateND][i];
        }

        __syncthreads();

        //MM state update (diagonal direction)
        max = 0.0f;//log 1
        //padding causes conflict-free access of scores whose 
        // true indices are SHIFTED right and down;
        //indices of pdiag1 and pdiag2 are SHIFTED down;
        //store score value in register;
        //transition probabilities represent the current position
        sco = scoreCache[threadIdx.x][i];
        //MI state (up direction)
        val1 = pdiag1[dpdsssStateMI][threadIdx.x] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMIc][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( max, val1 );
        //IM state (left direction)
        val2 = pdiag1[dpdsssStateIM][threadIdx.x+1] + 
            qrtrpCache[dptrMIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( max, val2 );
        //DN state (up)
        val3 = pdiag1[dpdsssStateDN][threadIdx.x] + qrtrpCache[dptrMDp][threadIdx.x];
        max = myhdmax/*fmaxf*/( max, val3 );
        //ND state (left)
        val4 = pdiag1[dpdsssStateND][threadIdx.x+1] + dbtrpCache[dptrMDp][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( max, val4 );
        //MM state (diagonal)
        val5 = pdiag2[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i] +
            sco;
        max = myhdmax/*fmaxf*/( max, val5 );
        val1 = MAPDPexpf(val1-max) + MAPDPexpf(val2-max) + MAPDPexpf(val3-max)+
               MAPDPexpf(val4-max) + MAPDPexpf(val5-max) + MAPDPexpf(-max)/*1*/;
        val1 = MAPDPlogf(val1) + max;


        //save log backward probabilities in the score cache:
        // the diagonal of the current iteration will not be used in the loop any more
        scoreCache[threadIdx.x][i] = val1;//WRITE


        //MI state update (up direction)
        //MM state (diagonal direction)
        val2 = pdiag2[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIMp][threadIdx.x+i] +
            sco;
        //MI state (up direction)
        val3 = pdiag1[dpdsssStateMI][threadIdx.x] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrIIc][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( val2, val3 );
        val2 = MAPDPexpf(val2-max) + MAPDPexpf(val3-max);
        val2 = MAPDPlogf(val2) + max;


        //IM state update (left direction)
        //MM state (diagonal direction)
        val3 = pdiag2[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrIMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i] +
            sco;
        //IM state (left direction)
        val4 = pdiag1[dpdsssStateIM][threadIdx.x+1] + 
            qrtrpCache[dptrIIc][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( val3, val4 );
        val3 = MAPDPexpf(val3-max) + MAPDPexpf(val4-max);
        val3 = MAPDPlogf(val3) + max;


        //DN state update (up)
        //MM state (diagonal direction)
        val4 = pdiag2[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrDMp][threadIdx.x] + dbtrpCache[dptrMMp][threadIdx.x+i] +
            sco;
        //DN state (up)
        val5 = pdiag1[dpdsssStateDN][threadIdx.x] + qrtrpCache[dptrDDp][threadIdx.x];
        max = myhdmax/*fmaxf*/( val4, val5 );
        val4 = MAPDPexpf(val4-max) + MAPDPexpf(val5-max);
        val4 = MAPDPlogf(val4) + max;


        //ND state update (left)
        //MM state (diagonal direction)
        val5 = pdiag2[dpdsssStateMM][threadIdx.x] + 
            qrtrpCache[dptrMMp][threadIdx.x] + dbtrpCache[dptrDMp][threadIdx.x+i] +
            sco;
        //ND state (left)
        sco = pdiag1[dpdsssStateND][threadIdx.x+1] + dbtrpCache[dptrDDp][threadIdx.x+i];
        max = myhdmax/*fmaxf*/( val5, sco );
        val5 = MAPDPexpf(val5-max) + MAPDPexpf(sco-max);
        val5 = MAPDPlogf(val5) + max;


        //WRITE...
        pdiag2[dpdsssStateMM][threadIdx.x+1] = val1;//WRITE
        pdiag2[dpdsssStateMI][threadIdx.x+1] = val2;//WRITE
        pdiag2[dpdsssStateIM][threadIdx.x+1] = val3;//WRITE
        pdiag2[dpdsssStateDN][threadIdx.x+1] = val4;//WRITE
        pdiag2[dpdsssStateND][threadIdx.x+1] = val5;//WRITE

        if( threadIdx.x == blockDim.x-1 ) {
            //WRITE the top row of the diagonal block
            //this position is not used by other threads in the current iteration
            bottmCache[dpdsssStateMM][i] = val1;
            bottmCache[dpdsssStateMI][i] = val2;
            bottmCache[dpdsssStateIM][i] = val3;
            bottmCache[dpdsssStateDN][i] = val4;
            bottmCache[dpdsssStateND][i] = val5;
        }

        //SYNCHRONIZATION on the NEXT iteration!

#ifdef CUMAPDP_BWD_TESTPRINT
        if(pronr==CUMAPDP_BWD_TESTPRINT){
            printf(" d=%u(%u) s=%u i%02d/%u (t%02u): len= %d addr= %u SC= %.4f (yx: %d,%d) "
                    "MM= %.6f MI= %.6f IM= %.6f DN= %.6f ND= %.6f\n"
                    "  >qMM= %.4f qMIc= %.4f qMD= %.4f qIM= %.4f qIIc= %.4f qDM= %.4f qDD= %.4f\n"
                    "  >dMM= %.4f dMIc= %.4f dMD= %.4f dIM= %.4f dIIc= %.4f dDM= %.4f dDD= %.4f\n",
                    blkdiagnum,lastydiagnum,blockIdx.x,i,ilim,threadIdx.x,
                    dbprolenCache,dbprodstCache,scoreCache[threadIdx.x][i], y-threadIdx.x,x+i,
                    val1,//pdiag2[dpdsssStateMM][threadIdx.x+1],
                    val2,//pdiag2[dpdsssStateMI][threadIdx.x+1],
                    val3,//pdiag2[dpdsssStateIM][threadIdx.x+1],
                    val4,//pdiag2[dpdsssStateDN][threadIdx.x+1],
                    val5,//pdiag2[dpdsssStateND][threadIdx.x+1],
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


#ifdef CUMAPDP_BWD_TESTPRINT
        if(pronr==CUMAPDP_BWD_TESTPRINT) {
//             if(threadIdx.x==0)
//                 for(float x=-23.0f;x<=1.0f;x+=0.2) 
//                     printf(" ><> exp(%.1f)= %.6f %.6f T=%.6f   log(%.1f)= %.6f T=%.6f\n",
//                         x,myfastapexpf(x),myNJfastapexpf(x),__expf(x),
//                         x+23.1f,myfastaplogf(x+23.1f),__logf(x+23.1f));
            printf(" >>> d=%u(%u) s=%u (t%02u): MAXLP= %.6f wrt= %d MAXLP= %.4f SUMEXP= %.6f  "
                   "Log_E-value= %.6f Log_Pair_E-value= %.6f\n", 
                    blkdiagnum,lastydiagnum,blockIdx.x,threadIdx.x,
                    maxlogprob, x+CUMAPDP_2DCACHE_DIM_X-1<(int)dbprolenCache,
                    tmpss2datbuffers[blockIdx.y*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES+1],
                    tmpss2datbuffers[blockIdx.y*CUSS_ALIGNED_N_DIFF_TOTAL_SCORES],
                    evals[0], evals[1]
            );
        }
#endif



    //write the result of calculations for upcoming blocks;
    //write two diagonals;
    if( 0 <= x && x < (int)dbprolenCache ) {
        DPLocWriteBuffer<CUMAPDP_2DCACHE_DIM_D>( pdiag1, tmpdpdiagbuffers, 
                          dbpos2, 0, 
                          dblen2, 1/*shift*/);
        DPLocWriteBuffer<CUMAPDP_2DCACHE_DIM_D>( pdiag2, tmpdpdiagbuffers, 
                          dbpos2, 
                          dpdssDiag2 * nTDPDiagScoreSubsections * dblen2, 
                          dblen2,
                          1/*shift*/);
    }


    //write the top of the diagonal block;
    if( 0 <= x+CUMAPDP_2DCACHE_DIM_D-1 && x+CUMAPDP_2DCACHE_DIM_D-1 < (int)dbprolenCache )
        DPLocWriteBuffer<CUMAPDP_2DCACHE_DIM_X>( 
            bottmCache, tmpdpbotbuffer, dbpos2+CUMAPDP_2DCACHE_DIM_D-1, 0, 
                        dblen2 );


    //write log backward probabilities written in scoreCache for SMEM reuse;
    //unroll since #registers is not a problem over all compute capabilities;
    #pragma unroll 8
    for( int i = 0; i < CUMAPDP_2DCACHE_DIM_D; i++ ) {
        //going upwards
        if( 0 <= y-i && y-i < (int)nqyposs && 0 <= x+i && x+i < (int)dbprolenCache ) {
            //starting position of line i of the oblq. diagonal block in the matrix:
            qpos = (y-i) * dblen2 + i;
            //tmpbwdprobs[qpos + dbpos2] = scoreCache[i][threadIdx.x];
            //NOTE: reuse tmpfwdprobs memory space for backward probabilities:
            // as forward probabilities will not be used for this diagonal block 
            // (they've been read above), add forward backward probabilities in place
            //(NOTE: dbpos2+i should be used here if dbpos2 were uint)
            atomicAdd( tmpfwdprobs + qpos + dbpos2, scoreCache[i][threadIdx.x]);
        }
    }
}
