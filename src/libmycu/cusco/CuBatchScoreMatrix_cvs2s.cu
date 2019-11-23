/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/SerializedCVS2Scores.cuh"
#include "libmycu/cupro/SerializedCVS2ScoresAttr.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_cvs2s.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// cvs2scores, serialized scores
// attr, attributes of serialized scores
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

// #define TEST_CalcSM_CVS2S_SMEMUnroll2x
// -------------------------------------------------------------------------
// CalcSM_CVS2S_SMEMUnroll2x: device code for calculating pairwise context 
// vector score;
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
// NOTE: results add to both outscores and outmodscores;
//
template <typename Func>
__global__ void CalcSM_CVS2S_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ cvs2scores,
    SerializedCVS2ScoresAttr attr,
    float cvswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    Func roundfunc,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores )
{
    extern __shared__ CUBSM_TYPE cvs2sCache[];//cached serialized scores
    //
#if !defined(CVS2S_SCORESE1S1Step2)
    __shared__ FPTYPE 
            qrenoCache,//cache for query ENO
            dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
#endif
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
    if( dbpronr < attr.ntotents_ )
        cvs2sCache[dbpronr] = cvs2scores[dbpronr];
    //for( dbpronr += blockDim.x; dbpronr < attr.ntotents_; dbpronr += blockDim.x )
    //    cvs2sCache[dbpronr] = cvs2scores[dbpronr];
    //}}
    //
    //for comments, see the CalcSMInit... kernels
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        if( threadIdx.y < pmv2DNoCVEls ) {
            qrcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVentrs+threadIdx.y]))[
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
#if !defined(CVS2S_SCORESE1S1Step2)
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
    //
    if( col < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVentrs+threadIdx.y]))[db1pos];
            if( threadIdx.y < 1 ) {
                dbcvpriorCache[0/*+threadIdx.y*/][threadIdx.x] = 
                   ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVprior/*+threadIdx.y*/]))[db1pos];
                dbcvnm2Cache[0/*+threadIdx.y*/][threadIdx.x] = 
                   ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos];
#if !defined(CVS2S_SCORESE1S1Step2)
                dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
                dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
            }
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < pmv2DNoCVEls ) {
            dbcvCache[pmv2DNoCVEls+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVentrs+threadIdx.y]))[db1pos2];
            if( threadIdx.y < 1 ) {
                dbcvpriorCache[1/*+threadIdx.y*/][threadIdx.x] = 
                   ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVprior/*+threadIdx.y*/]))[db1pos2];
                dbcvnm2Cache[1/*+threadIdx.y*/][threadIdx.x] = 
                   ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVnorm2/*+threadIdx.y*/]))[db1pos2];
#if !defined(CVS2S_SCORESE1S1Step2)
                dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
                dbenoCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
            }
        }
    }
    //the warp reads data written by other warps, sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    float score1 = 0.0f, score2 = 0.0f;
    FPTYPE f1;
#if !defined(CVS2S_SCORESE1S1Step2)
    FPTYPE e1 = qrenoCache;
    SerializedCVS2Scores<CUBSM_TYPE> cvs2s(
        cvs2sCache, attr.szalloc_,
        attr.nenos_, attr.ntbls_
    );
#endif
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoCVEls; i++ ) {
        f1 = qrcvCache[i][threadIdx.y];
        score1 += f1 * dbcvCache[i][threadIdx.x];
        if( col2 < (ndb1poss + ndbCposs))
            score2 += f1 * dbcvCache[pmv2DNoCVEls+i][threadIdx.x];
    }

    score1 += attr.CVS_loKAPPA0_;//add off-diagonal entry of C_2 to the dot product
    //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
    f1 = qrcvnm2Cache[threadIdx.y];
    score1 = (f1 + 1.0f + attr.CVS_loKAPPA0_) *
             (dbcvnm2Cache[0][threadIdx.x] + 1.0f + attr.CVS_loKAPPA0_) - 
             score1 * score1;
    MYASSERT( score1 > 0.0f, "Non-positive determinant 1.");
    //log odds score of normal vectors
    //NOTE: test (wrt sens and aq) and REMOVE typecast to int: CHANGED!
    score1 = attr.CVS_CTERM_ + attr.CVS_PowerNU0_ * __logf(score1) - 
             roundfunc(qrcvpriorCache[threadIdx.y]) - roundfunc(dbcvpriorCache[0][threadIdx.x]);
    //translate score
    score1 = 
#ifdef CVS2S_SCORESE1S1Step2
        //SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2( cvs2sCache, score1, attr.card_, attr.shft_ );
        SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2Boundary( 
            cvs2sCache, score1, attr.card_, attr.shft_,
            attr.key1first_, attr.value1first_, attr.key1last_, attr.value1last_);
#else
        //cvs2s.GetScore( score1, e1, dbenoCache[0][threadIdx.x] );
        cvs2s.GetScoreE1S1Step2( score1, attr.card_, attr.shft_, 0.f, 0.f );
#endif

    if( col2 < (ndb1poss + ndbCposs)) {
        score2 += attr.CVS_loKAPPA0_;//add off-diagonal entry of C_2 to the dot product
        //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
        score2 = (f1 + 1.0f + attr.CVS_loKAPPA0_) *
                (dbcvnm2Cache[1][threadIdx.x] + 1.0f + attr.CVS_loKAPPA0_) - 
                score2 * score2;
        MYASSERT( score2 > 0.0f, "Non-positive determinant 2.");
        //log odds score of normal vectors
        //NOTE: test (wrt sens and aq) and REMOVE typecast to int: CHANGED!
        score2 = attr.CVS_CTERM_ + attr.CVS_PowerNU0_ * __logf(score2) - 
                roundfunc(qrcvpriorCache[threadIdx.y]) - roundfunc(dbcvpriorCache[1][threadIdx.x]);
        //translate score
        score2 = 
#ifdef CVS2S_SCORESE1S1Step2
            //SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2( cvs2sCache, score2, attr.card_, attr.shft_ );
            SerializedCVS2Scores<CUBSM_TYPE>::GetScoreE1S1Step2Boundary(
                cvs2sCache, score2, attr.card_, attr.shft_,
                attr.key1first_, attr.value1first_, attr.key1last_, attr.value1last_);
#else
            //cvs2s.GetScore( score2, e1, dbenoCache[1][threadIdx.x] );
            cvs2s.GetScoreE1S1Step2( score2, attr.card_, attr.shft_, 0.f, 0.f );
#endif
    }

#ifdef TEST_CalcSM_CVS2S_SMEMUnroll2x
    size_t querypos = row + querposoffset;
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
# if !defined(CVS2S_SCORESE1S1Step2)
    uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    MYASSERT( qrenoCache == ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr], "Inconsistency.");
    MYASSERT( dbenoCache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO]))[dbpronr], "Inconsistency.");
# endif
    MYASSERT( qrcvpriorCache[threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVprior]))[querypos], "Inconsistency.");
    MYASSERT( qrcvnm2Cache[threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVnorm2]))[querypos], "Inconsistency.");
    MYASSERT( dbcvpriorCache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVprior]))[db1pos], "Inconsistency.");
    MYASSERT( dbcvnm2Cache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVnorm2]))[db1pos], "Inconsistency.");
    for( int i = 0; i < pmv2DNoCVEls; i++ ) {
        MYASSERT( qrcvCache[i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DCVentrs+i]))[querypos], "Inconsistency.");
        MYASSERT( dbcvCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DCVentrs+i]))[db1pos], "Inconsistency.");
    }
    if( col2 < (ndb1poss + ndbCposs)) {
# if !defined(CVS2S_SCORESE1S1Step2)
        size_t dbpronr2 = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        MYASSERT( dbenoCache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO]))[dbpronr2], "Inconsistency.");
# endif
        MYASSERT( dbcvpriorCache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVprior]))[db1pos2], "Inconsistency.");
        MYASSERT( dbcvnm2Cache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVnorm2]))[db1pos2], "Inconsistency.");
        for( int i = 0; i < pmv2DNoCVEls; i++ ) {
            MYASSERT( dbcvCache[pmv2DNoCVEls+i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DCVentrs+i]))[db1pos2], "Inconsistency.");
        }
        atomicAdd(&d_testvar,1UL);
    }
    atomicAdd(&d_testvar,1UL);
    if( querposoffset < 1)
    if((size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset) <= d_testvar )
        printf(" Counter %zu col2= %u qcnt= %u(%lu) dbcnt= %u(%lu)   ss= %lu\n", 
               d_testvar, col2, nqyposs, querposoffset, ndb1poss + ndbCposs, bdb1posoffset + bdbCposoffset,
               (size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset)
        );
#endif

    row = row * (ndb1poss + ndbCposs + dbxpad);

    //perform coalescent write of scores
    if( score1 ) {
        //atomicAdd is faster than += when we need coalescent write performed once
        atomicAdd( &outscores[row + col], score1 * cvswgt );
        atomicAdd( &outmodscores[row + col], (score1 - CONSTCVSSHIFT) * cvswgt );
    }

    if( col2 < (ndb1poss + ndbCposs) && score2 ) {
        atomicAdd( &outscores[row + col2], score2 * cvswgt );
        atomicAdd( &outmodscores[row + col2], (score2 - CONSTCVSSHIFT) * cvswgt );
    }
}
