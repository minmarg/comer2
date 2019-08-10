/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/SerializedScoresTM.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_hdp1s.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// hdp1sTexo, serialized scores represented by a texture object
// attr, attributes of serialized scores
// hdp1swgt, weight of scoring HDP cluster membership match
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

// #define TEST_CalcSM_HDP1S_SMEMUnroll2x
// -------------------------------------------------------------------------
// CalcSM_HDP1S_SMEMUnroll2x: device code for calculating HDP scores 
// using shared memory and twofold unrolling along the x axis;
// each block processes two blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses texture memory (for big data), which due to random acces 
// provides greater efficiency in comparison to other types of memory;
// NOTE: output pointers should be aligned!
// NOTE: results add to outscores;
// 
__global__ void CalcSM_HDP1S_SMEMUnroll2x(
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr attr,
    float hdp1swgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* /*outmodscores*/ )
{
    __shared__ FPTYPE 
            qrenoCache,//cache for query ENO
            dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
    __shared__ FPTYPE 
            qrhdp1prbCache[SMINIT_2DCACHE_DIM],//cache for HDP probabilities at query positions
            dbhdp1prbCache[2][SMINIT_2DCACHE_DIM];//cache for HDP probs at db profile positions
    __shared__ INTYPE 
            qrhdp1ndxCache[SMINIT_2DCACHE_DIM],//cache for HDP cluster indices at query positions
            dbhdp1ndxCache[2][SMINIT_2DCACHE_DIM];//cache for HDP cluster indices at db profile positions
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
    uint qpronr, dbpronr;
    //
    //for comments, see the CalcSMInit... kernels
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        if( threadIdx.y < 1 ) {
            qrhdp1prbCache[threadIdx.x] = 
            ((FPTYPE*)(dc_pm2dvfields_[pmv2DHDP1prb]))[
                blockbeg_y + querposoffset + threadIdx.x
            ];
            qrhdp1ndxCache[threadIdx.x] = 
            ((INTYPE*)(dc_pm2dvfields_[pmv2DHDP1ind]))[
                blockbeg_y + querposoffset + threadIdx.x
            ];
            //read query ENO
            //NOTE: valid when processing one query at a time
            if( threadIdx.x < 1 ) {
                qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                //read only one element per block (blockDim.y x blockDim.x)
                qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr];
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < 1 ) {
            dbhdp1prbCache[0/*+threadIdx.y*/][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1prb/*+threadIdx.y*/]))[db1pos];
            dbhdp1ndxCache[0/*+threadIdx.y*/][threadIdx.x] = 
                ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1ind/*+threadIdx.y*/]))[db1pos];
            dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
            dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[dbpronr];
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < 1 ) {
            dbhdp1prbCache[1/*+threadIdx.y*/][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1prb/*+threadIdx.y*/]))[db1pos2];
            dbhdp1ndxCache[1/*+threadIdx.y*/][threadIdx.x] = 
                ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1ind/*+threadIdx.y*/]))[db1pos2];
            dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
            dbenoCache[1/*+threadIdx.y*/][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO/*+threadIdx.y*/]))[dbpronr];
        }
    }
    //the warp reads data written by other warps, sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    float score1 = 0.0f, score2 = 0.0f;
    FPTYPE e1 = qrenoCache;
    INTYPE n1 = qrhdp1ndxCache[threadIdx.y];
    FPTYPE f1 = qrhdp1prbCache[threadIdx.y];
    //
    //using registers is a bit more effective although access to 
    //SMEM without bank conflicts (below) is almost as effective
    FPTYPE e21 = dbenoCache[0][threadIdx.x];
    FPTYPE e22 = dbenoCache[1][threadIdx.x];
    INTYPE n21 = dbhdp1ndxCache[0][threadIdx.x];
    INTYPE n22 = dbhdp1ndxCache[1][threadIdx.x];
    FPTYPE f21 = dbhdp1prbCache[0][threadIdx.x];
    FPTYPE f22 = dbhdp1prbCache[1][threadIdx.x];
    //
    qpronr = col2 < (ndb1poss + ndbCposs);
    //
#ifdef HDP1S_SCORESP1E2
    if( 0 <= n1 ) {
        if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
            score1 = SerializedScoresTM<CUBSM_TYPE>::GetScoreP1E2( 
                hdp1sTexo, attr.eth1_, attr.eth2_, attr.nelems_,
                n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
                e1, e21//dbenoCache[0][threadIdx.x]
            );
        if( qpronr && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
            score2 = SerializedScoresTM<CUBSM_TYPE>::GetScoreP1E2( 
                hdp1sTexo, attr.eth1_, attr.eth2_, attr.nelems_,
                n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
                e1, e22//dbenoCache[1][threadIdx.x]
            );
    }
#else
    if( 0 <= n1 ) {
        SerializedScoresTM<CUBSM_TYPE> hdp1s(
            hdp1sTexo, attr.szalloc_,
            attr.nplvs_, attr.nenos_, attr.card_,
            attr.ntbls_, attr.nelems_
        );
//         if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
//             score1 = hdp1s.GetScore( 
//                 n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
//                 f1, f21,//dbhdp1prbCache[0][threadIdx.x], 
//                 e1, e21//dbenoCache[0][threadIdx.x]
//             );
//         if( qpronr && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
//             score2 = hdp1s.GetScore( 
//                 n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
//                 f1, f22,//dbhdp1prbCache[1][threadIdx.x], 
//                 e1, e22//dbenoCache[1][threadIdx.x]
//             );
        if( 0 <= n21 )//dbhdp1ndxCache[0][threadIdx.x])
            score1 = hdp1s.GetScoreP1( 
                n1, n21,//dbhdp1ndxCache[0][threadIdx.x], 
                0.0f, 0.0f, 
                e1, e21//dbenoCache[0][threadIdx.x]
            );
        if( qpronr && 0 <= n22 )//dbhdp1ndxCache[1][threadIdx.x])
            score2 = hdp1s.GetScoreP1( 
                n1, n22,//dbhdp1ndxCache[1][threadIdx.x], 
                0.0f, 0.0f, 
                e1, e22//dbenoCache[1][threadIdx.x]
            );
    }
#endif
    //

#ifdef TEST_CalcSM_HDP1S_SMEMUnroll2x
    size_t querypos = row + querposoffset;
    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    MYASSERT( qrenoCache == ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr], "Inconsistency.");
    MYASSERT( dbenoCache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO]))[dbpronr], "Inconsistency.");
    MYASSERT( qrhdp1prbCache[threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DHDP1prb]))[querypos], "Inconsistency.");
    MYASSERT( qrhdp1ndxCache[threadIdx.y] == ((INTYPE*)(dc_pm2dvfields_[pmv2DHDP1ind]))[querypos], "Inconsistency.");
    MYASSERT( dbhdp1prbCache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1prb]))[db1pos], "Inconsistency.");
    MYASSERT( dbhdp1ndxCache[0][threadIdx.x] == ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1ind]))[db1pos], "Inconsistency.");
    if( col2 < (ndb1poss + ndbCposs)) {
        size_t dbpronr2 = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        MYASSERT( dbenoCache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO]))[dbpronr2], "Inconsistency.");
        MYASSERT( dbhdp1prbCache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1prb]))[db1pos2], "Inconsistency.");
        MYASSERT( dbhdp1ndxCache[1][threadIdx.x] == ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1ind]))[db1pos2], "Inconsistency.");
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
    if( -99.0f < score1 && score1 ) {
        score1 *= GetFactorWrtENOAndHDPPrb(
                        e1, e21,//dbenoCache[0][threadIdx.x],
                        f1, f21);//dbhdp1prbCache[0][threadIdx.x]);
        score1 *= hdp1swgt;
        //atomicAdd is faster than += when we need coalescent write performed once
        atomicAdd( &outscores[row + col], score1 );
    }

    if( qpronr/*col2 < (ndb1poss + ndbCposs)*/ && -99.0f < score2 && score2 ) {
        score2 *= GetFactorWrtENOAndHDPPrb(
                        e1, e22,//dbenoCache[1][threadIdx.x],
                        f1, f22);//dbhdp1prbCache[1][threadIdx.x]);
        score2 *= hdp1swgt;
        atomicAdd( &outscores[row + col2], score2 );
    }
}

