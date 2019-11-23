/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/SerializedScoresSM.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_ssss.cuh"
#include "CuBatchScoreMatrix_cvs2s.cuh"
#include "CuBatchScoreMatrix_ssss_cvs2s.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// sssscores, serialized secondary structure scores
// ssattr, attributes of serialized SS scores
// ssswgt, weight of scoring secondary structure state predictions
// cvs2scores, serialized context vector scores
// cvattr, attributes of serialized scores
// cvswgt, weight of scoring CVS2S scores
// nqyposs, number of query positions to process
// ndb1poss, number of cached db profile positions to process
// ndbCposs, number of new db profile positions to process
// dbxpad, number of padded positions for memory alignment
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
//

// -------------------------------------------------------------------------
// CalcSM_SSSS_CVS2S_SMEMUnroll2x_2: device code for calculating pairwise 
// SS and context vector scores at once;
// NOTE: the kernel is valid only for CVS.AVGLEN == 1!
// matrix using shared memory and twofold unrolling along the x axis;
// each block processes two blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses more shared memory but allows for greater occupancy;
// although bank conflicts cannot be avoided (due to random acces to the 
// SMEM), using SMEM garantees much greater efficiency in comparison to 
// other types of memory;
// NOTE: output pointers should be aligned!
// NOTE: results add to outscores and rewrite outmodscores;
// 
template <typename Func>
__global__ void CalcSM_SSSS_CVS2S_SMEMUnroll2x_2(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr ssattr,
    float ssswgt,
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr cvattr,
    float cvswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    Func roundfunc,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores )
{
    extern __shared__ CUBSM_TYPE sssandcvsCache[];//cached serialized scores, both ssss and cvs2s
    //
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
    __shared__ FPTYPE 
            qrenoCache,//cache for query ENO
            dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
#endif
    __shared__ FPTYPE 
            qrposCache[pmv2DNoSSSps][SMINIT_2DCACHE_DIM],//cache for a tile of SSS probabilities at query positions
            dbposCache[TIMES2(pmv2DNoSSSps)][SMINIT_2DCACHE_DIM];//cache for a tile of SSS probs at db profile positions
    //
    __shared__ FPTYPE 
            qrcvpriorCache[SMINIT_2DCACHE_DIM],//cache for CV prior probabilities at query positions
            dbcvpriorCache[2][SMINIT_2DCACHE_DIM];//cache for CV prior probs at db profile positions
    __shared__ FPTYPE 
            qrcvnm2Cache[SMINIT_2DCACHE_DIM],//cache for CV squared norms calculated at query positions
            dbcvnm2Cache[2][SMINIT_2DCACHE_DIM];//cache for CV squared norms calculated at db profile positions
    __shared__ FPTYPE 
            qrcvCache[pmv2DNoCVEls][SMINIT_2DCACHE_DIM],//cache for a tile of CVs at query positions
            dbcvCache[TIMES2(pmv2DNoCVEls)][SMINIT_2DCACHE_DIM];//cache for a tile of CVs at db profile positions
    //
    uint blockbeg_y = blockIdx.y * blockDim.y;
    uint row = blockbeg_y + threadIdx.y;
    const uint col = blockIdx.x * blockDim.x * 2 + threadIdx.x;//logical column
    const uint col2 = col + blockDim.x;//logical column
    //physical indices:
    uint db1pos;
    uint db1pos2;
    uint dbfldsndx;
    uint dbfldsndx2;
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    if( col2 < ndb1poss ) { db1pos2 = col2 + bdb1posoffset;
                            dbfldsndx2 = pmv2DTotFlds;
    } else {                db1pos2 = col2 - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx2 = TIMES2(pmv2DTotFlds);
    }
    //
    uint dbpronr;
    //
    //{{CACHE SERIALIZED SCORES using coalescent access (reuse registers)
    //NOTE: IMPORTANT: total number of entries (attr.ntotents_) is 
    // assumed to be less than the block size;
    //if this is not the case, uncomment the for loop for n-fold caching!
    dbpronr = threadIdx.y * blockDim.x + threadIdx.x;
    if( dbpronr < ssattr.ntotents_ )
        sssandcvsCache[dbpronr] = sssscores[dbpronr];
    else if( dbpronr < ssattr.ntotents_ + cvattr.ntotents_ )
        sssandcvsCache[dbpronr] = cvs2scores[dbpronr-ssattr.ntotents_];
    //for( dbpronr += blockDim.x; dbpronr < ssattr.ntotents_; dbpronr += blockDim.x )
    //    sssandcvsCache[dbpronr] = sssscores[dbpronr];
    //for( ; dbpronr < ssattr.ntotents_ + cvattr.ntotents_; dbpronr += blockDim.x )
    //    sssandcvsCache[dbpronr] = cvs2scores[dbpronr-ssattr.ntotents_];
    //}}
    //
    //for comments, see the CalcSMInit... kernels
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            qrcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVentrs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            if( threadIdx.y < pmv2DNoSSSps ) {
                qrposCache[threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DSSsprbs+threadIdx.y]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
                if( threadIdx.y < 1 ) {
                    //read query log priors and squared norms
                    qrcvpriorCache[threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVprior]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
                    qrcvnm2Cache[threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVnorm2]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    //read query ENO
                    //NOTE: valid when processing one query at a time
                    if( threadIdx.x < 1 ) {
                        uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                        //read only one element per block (blockDim.y x blockDim.x)
                        qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr];
                    }
#endif
                }
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVentrs+threadIdx.y]))[db1pos];
            if( threadIdx.y < pmv2DNoSSSps ) {
                dbposCache[threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSsprbs+threadIdx.y]))[db1pos];
                if( threadIdx.y < 1 ) {
                    dbcvpriorCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVprior/*+threadIdx.y*/]))[db1pos];
                    dbcvnm2Cache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
                    dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                        ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
                }
            }
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[pmv2DNoCVEls+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVentrs+threadIdx.y]))[db1pos2];
            if( threadIdx.y < pmv2DNoSSSps ) {
                dbposCache[pmv2DNoSSSps+threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSsprbs+threadIdx.y]))[db1pos2];
                if( threadIdx.y < 1 ) {
                    dbcvpriorCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVprior/*+threadIdx.y*/]))[db1pos2];
                    dbcvnm2Cache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos2];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
                    dbenoCache[1/*+threadIdx.y*/][threadIdx.x] = 
                        ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
                }
            }
        }
    }
    //the warp reads data written by other warps, sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    //reuse registers
    dbpronr = (col2 < (ndb1poss + ndbCposs));
    //
    float score1, score2;
    float s1, s2;

    CalcSSSScoreSMEMUnroll2x( 
        ssattr,
        ssswgt,
        sssandcvsCache,
#if !defined(SSSS_SCORESP1E1)
        qrenoCache,
        dbenoCache,
#endif
        qrposCache,
        dbposCache,
        dbpronr/*col2 < (ndb1poss + ndbCposs)*/,
        score1, score2
    );

    CalcCVS2ScoreSMEMUnroll2x( 
        cvattr,
        cvswgt,
        sssandcvsCache + ssattr.ntotents_,
#if !defined(CVS2S_SCORESE1S1Step2)
        qrenoCache,
        dbenoCache,
#endif
        qrcvpriorCache,
        dbcvpriorCache,
        qrcvnm2Cache,
        dbcvnm2Cache,
        qrcvCache,
        dbcvCache,
        roundfunc,
        dbpronr/*col2 < (ndb1poss + ndbCposs)*/,
        s1, s2
    );

    score1 += s1;

    row = row * (ndb1poss + ndbCposs + dbxpad);

    //perform coalescent write of scores
    //atomicAdd is faster than += when we need coalescent write performed once
    atomicAdd( &outscores[row + col], score1 );
    outmodscores[row + col] = score1 - CONSTCVSSHIFT * cvswgt;

    if( dbpronr/*col2 < (ndb1poss + ndbCposs)*/) {
        score2 += s2;
        atomicAdd( &outscores[row + col2], score2 );
        outmodscores[row + col2] = score2 - CONSTCVSSHIFT * cvswgt;
    }
}



// -------------------------------------------------------------------------
// CalcSM_SSSS_CVS2S_SMEMUnroll2x_1: same as CalcSM_SSSS_CVS2S_SMEMUnroll2x_2
// except that this version rearranges the order of the final 
// operations for calculations and write
template <typename Func>
__global__ void CalcSM_SSSS_CVS2S_SMEMUnroll2x_1(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr ssattr,
    float ssswgt,
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr cvattr,
    float cvswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    Func roundfunc,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores )
{
    extern __shared__ CUBSM_TYPE sssandcvsCache[];//cached serialized scores, both ssss and cvs2s
    //
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
    __shared__ FPTYPE 
            qrenoCache,//cache for query ENO
            dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
#endif
    __shared__ FPTYPE 
            qrposCache[pmv2DNoSSSps][SMINIT_2DCACHE_DIM],//cache for a tile of SSS probabilities at query positions
            dbposCache[TIMES2(pmv2DNoSSSps)][SMINIT_2DCACHE_DIM];//cache for a tile of SSS probs at db profile positions
    //
    __shared__ FPTYPE 
            qrcvpriorCache[SMINIT_2DCACHE_DIM],//cache for CV prior probabilities at query positions
            dbcvpriorCache[2][SMINIT_2DCACHE_DIM];//cache for CV prior probs at db profile positions
    __shared__ FPTYPE 
            qrcvnm2Cache[SMINIT_2DCACHE_DIM],//cache for CV squared norms calculated at query positions
            dbcvnm2Cache[2][SMINIT_2DCACHE_DIM];//cache for CV squared norms calculated at db profile positions
    __shared__ FPTYPE 
            qrcvCache[pmv2DNoCVEls][SMINIT_2DCACHE_DIM],//cache for a tile of CVs at query positions
            dbcvCache[TIMES2(pmv2DNoCVEls)][SMINIT_2DCACHE_DIM];//cache for a tile of CVs at db profile positions
    //
    uint blockbeg_y = blockIdx.y * blockDim.y;
    uint row = blockbeg_y + threadIdx.y;
    const uint col = blockIdx.x * blockDim.x * 2 + threadIdx.x;//logical column
    const uint col2 = col + blockDim.x;//logical column
    //physical indices:
    uint db1pos;
    uint db1pos2;
    uint dbfldsndx;
    uint dbfldsndx2;
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    if( col2 < ndb1poss ) { db1pos2 = col2 + bdb1posoffset;
                            dbfldsndx2 = pmv2DTotFlds;
    } else {                db1pos2 = col2 - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx2 = TIMES2(pmv2DTotFlds);
    }
    //
    uint dbpronr;
    //
    //{{CACHE SERIALIZED SCORES using coalescent access (reuse registers)
    //NOTE: IMPORTANT: total number of entries (attr.ntotents_) is 
    // assumed to be less than the block size;
    //if this is not the case, uncomment the for loop for n-fold caching!
    dbpronr = threadIdx.y * blockDim.x + threadIdx.x;
    if( dbpronr < ssattr.ntotents_ )
        sssandcvsCache[dbpronr] = sssscores[dbpronr];
    else if( dbpronr < ssattr.ntotents_ + cvattr.ntotents_ )
        sssandcvsCache[dbpronr] = cvs2scores[dbpronr-ssattr.ntotents_];
    //for( dbpronr += blockDim.x; dbpronr < ssattr.ntotents_; dbpronr += blockDim.x )
    //    sssandcvsCache[dbpronr] = sssscores[dbpronr];
    //for( ; dbpronr < ssattr.ntotents_ + cvattr.ntotents_; dbpronr += blockDim.x )
    //    sssandcvsCache[dbpronr] = cvs2scores[dbpronr-ssattr.ntotents_];
    //}}
    //
    //for comments, see the CalcSMInit... kernels
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            qrcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVentrs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            if( threadIdx.y < pmv2DNoSSSps ) {
                qrposCache[threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DSSsprbs+threadIdx.y]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
                if( threadIdx.y < 1 ) {
                    //read query log priors and squared norms
                    qrcvpriorCache[threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVprior]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
                    qrcvnm2Cache[threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVnorm2]))[
                        blockbeg_y + querposoffset + threadIdx.x
                    ];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    //read query ENO
                    //NOTE: valid when processing one query at a time
                    if( threadIdx.x < 1 ) {
                        uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                        //read only one element per block (blockDim.y x blockDim.x)
                        qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr];
                    }
#endif
                }
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVentrs+threadIdx.y]))[db1pos];
            if( threadIdx.y < pmv2DNoSSSps ) {
                dbposCache[threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSsprbs+threadIdx.y]))[db1pos];
                if( threadIdx.y < 1 ) {
                    dbcvpriorCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVprior/*+threadIdx.y*/]))[db1pos];
                    dbcvnm2Cache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
                    dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                        ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
                }
            }
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        //NOTE: assuming pmv2DNoSSSps < pmv2DNoCVEls
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[pmv2DNoCVEls+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVentrs+threadIdx.y]))[db1pos2];
            if( threadIdx.y < pmv2DNoSSSps ) {
                dbposCache[pmv2DNoSSSps+threadIdx.y][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSsprbs+threadIdx.y]))[db1pos2];
                if( threadIdx.y < 1 ) {
                    dbcvpriorCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVprior/*+threadIdx.y*/]))[db1pos2];
                    dbcvnm2Cache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos2];
#if !defined(CVS2S_SCORESE1S1Step2) || !defined(SSSS_SCORESP1E1)
                    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
                    dbenoCache[1/*+threadIdx.y*/][threadIdx.x] = 
                        ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
                }
            }
        }
    }
    //the warp reads data written by other warps, sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;

    float score, s;

    CalcSSSScoreSMEMUnroll2x_1<0>( 
        ssswgt,
        sssandcvsCache,
        qrposCache,
        dbposCache,
        score
    );

    CalcCVS2ScoreSMEMUnroll2x_1<Func,0,0>( 
        cvattr,
        cvswgt,
        sssandcvsCache + ssattr.ntotents_,
        qrcvpriorCache,
        dbcvpriorCache,
        qrcvnm2Cache,
        dbcvnm2Cache,
        qrcvCache,
        dbcvCache,
        roundfunc,
        s
    );

    score += s;

    row = row * (ndb1poss + ndbCposs + dbxpad);

    //perform coalescent write of scores
    //atomicAdd is faster than += when we need coalescent write performed once
    atomicAdd( &outscores[row + col], score );
    outmodscores[row + col] = score - CONSTCVSSHIFT * cvswgt;


    if( ndb1poss + ndbCposs <= col2 )
        return;


    CalcSSSScoreSMEMUnroll2x_1<pmv2DNoSSSps>( 
        ssswgt,
        sssandcvsCache,
        qrposCache,
        dbposCache,
        score
    );

    CalcCVS2ScoreSMEMUnroll2x_1<Func,pmv2DNoCVEls,1>( 
        cvattr,
        cvswgt,
        sssandcvsCache + ssattr.ntotents_,
        qrcvpriorCache,
        dbcvpriorCache,
        qrcvnm2Cache,
        dbcvnm2Cache,
        qrcvCache,
        dbcvCache,
        roundfunc,
        s
    );

    score += s;

    atomicAdd( &outscores[row + col2], score );
    outmodscores[row + col2] = score - CONSTCVSSHIFT * cvswgt;
}
