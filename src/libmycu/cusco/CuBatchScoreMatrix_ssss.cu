/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/SerializedScoresSM.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_ssss.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// sssscores, serialized scores
// attr, attributes of serialized scores
// ssswgt, weight of scoring secondary structure state predictions
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

// #define TEST_CalcSM_SSSS_SMEMUnroll2x
// -------------------------------------------------------------------------
// CalcSM_SSSS_SMEMUnroll2x: device code for calculating SSS score 
// matrix using shared memory and twofold unrolling along the x axis;
// each block processes two blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses more shared memory but allows for greater occupancy;
// although bank conflicts cannot be avoided (due to random acces to the 
// SMEM), using SMEM garantees much greater efficiency in comparison to 
// other types of memory;
// reduced load bank conflicts!
// NOTE: output pointers should be aligned!
// NOTE: results OVERWRITE outmodscores while they are added to at the 
// address of outscores;
// 
__global__ void CalcSM_SSSS_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ sssscores,
    SerializedScoresAttr attr,
    float ssswgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores )
{
    extern __shared__ CUBSM_TYPE ssssCache[];//cached serialized scores
    //
    //__shared__ CHTYPE qrsssCache[SMINIT_2DCACHE_DIM];//cache for a tile of SS states at query positions
    //__shared__ CHTYPE dbsssCache[2][SMINIT_2DCACHE_DIM];//cache for a tile of SS states at db profile positions
#if !defined(SSSS_SCORESP1E1)
    __shared__ FPTYPE 
            qrenoCache,//cache for query ENO
            dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
#endif
    __shared__ FPTYPE 
            qrposCache[pmv2DNoSSSps][SMINIT_2DCACHE_DIM],//cache for a tile of query positions
            dbposCache[TIMES2(pmv2DNoSSSps)][SMINIT_2DCACHE_DIM];//cache for a tile of db profile positions
    //
    uint blockbeg_y = blockIdx.y * blockDim.y;
    //uint blockbeg_x = blockIdx.x * blockDim.x * 2;
    uint row = blockbeg_y + threadIdx.y;
    //uint col = blockbeg_x + threadIdx.x;
    const uint col = blockIdx.x * blockDim.x * 2 + threadIdx.x;//logical column
    const uint col2 = col + blockDim.x;//logical column
    //physical indices:
    uint db1pos;// = col < ndb1poss ? col + bdb1posoffset: (col - ndb1poss + bdbCposoffset)/*jump to section ndbCposs*/;
    uint db1pos2;// = col2 < ndb1poss ? col2 + bdb1posoffset: (col2 - ndb1poss + bdbCposoffset)/*section ndbCposs*/;
    uint dbfldsndx;// = col < ndb1poss ? pmv2DTotFlds: TIMES2(pmv2DTotFlds);//db index for dc_pm2dvfields_
    uint dbfldsndx2;// = col2 < ndb1poss ? pmv2DTotFlds: TIMES2(pmv2DTotFlds);//dbC index for dc_pm2dvfields_
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
    //{{CACHE SERIALIZED SCORES using coalescent access (reuse registers)
    //NOTE: IMPORTANT: total number of entries (attr.ntotents_) is 
    // assumed to be less than the block size;
    //if this is not the case, uncomment the for loop for n-fold caching!
    dbpronr = threadIdx.y * blockDim.x + threadIdx.x;
    if( dbpronr < attr.ntotents_ )
        ssssCache[dbpronr] = sssscores[dbpronr];
    //for( dbpronr += blockDim.x; dbpronr < attr.ntotents_; dbpronr += blockDim.x )
    //    ssssCache[dbpronr] = sssscores[dbpronr];
    //}}
    //
    //for comments, see the CalcSMInit... kernels
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        if( threadIdx.y < pmv2DNoSSSps ) {
            qrposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DSSsprbs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            if( threadIdx.y < 1 ) {
                //read query SS states
                //qrsssCache[threadIdx.x] = 
                //    ((CHTYPE*)(dc_pm2dvfields_[pmv2DSSstate]))[
                //        blockbeg_y + querposoffset + threadIdx.x
                //    ];
#if !defined(SSSS_SCORESP1E1)
                //read query ENO
                //NOTE: valid when processing one query at a time
                if( threadIdx.x < 1 ) {
                    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                    //read only one element per block (blockDim.y x blockDim.x)
                    qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr];
                }
#endif
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
        if( threadIdx.y < pmv2DNoSSSps ) {
            dbposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSsprbs+threadIdx.y]))[db1pos];
            if( threadIdx.y < 1 ) {
                //dbsssCache[0/*+threadIdx.y*/][threadIdx.x] = 
                //    ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSstate/*+threadIdx.y*/]))[db1pos];
#if !defined(SSSS_SCORESP1E1)
                dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[dbpronr];
#endif
            }
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        if( threadIdx.y < pmv2DNoSSSps ) {
            dbposCache[pmv2DNoSSSps+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSsprbs+threadIdx.y]))[db1pos2];
            if( threadIdx.y < 1 ) {
                //dbsssCache[1/*+threadIdx.y*/][threadIdx.x] = 
                //    ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSstate/*+threadIdx.y*/]))[db1pos2];
#if !defined(SSSS_SCORESP1E1)
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
    FPTYPE f1, f2;
    //CHTYPE s1 = qrsssCache[threadIdx.y];
#if !defined(SSSS_SCORESP1E1)
    FPTYPE e1 = qrenoCache;
#endif
    //
    //
#if !defined(SSSS_SCORESP1E1)
    SerializedScoresSM<CUBSM_TYPE> ssss(
        ssssCache, attr.szalloc_,
        attr.nplvs_, attr.nenos_, attr.card_,
        attr.ntbls_, attr.nelems_
    );
#endif
    //
    //unrolling behaviour is default, and here it gives substantial speedup
    #pragma unroll
    for( int i = 0; i < pmv2DNoSSSps; i++ ) {
        //exchanged access for qrposCache
        f1 = qrposCache[i][threadIdx.y];
        qpronr = GetSSSSTableIndex( i/*s1*/, f1 );//reuse registers
        #pragma unroll
        for( int j = 0; j < pmv2DNoSSSps; j++ ) {
            f2 = dbposCache[j][threadIdx.x];
            dbpronr = GetSSSSTableIndex( j/*dbsssCache[0][threadIdx.x]*/, f2 );//reuse registers
            score1 += f1 * f2 * 
#ifdef SSSS_SCORESP1E1
                SerializedScoresSM<CUBSM_TYPE>::GetScoreP1E1( ssssCache, qpronr, dbpronr );
#else
                //ssss.GetScoreP1( qpronr, dbpronr, 0.f, 0.f, e1, dbenoCache[0][threadIdx.x] );
                ssss.GetScoreP1E1( qpronr, dbpronr, 0.f, 0.f, 0.f, 0.f );
#endif
            //
            if( col2 < (ndb1poss + ndbCposs)) {
                f2 = dbposCache[pmv2DNoSSSps+j][threadIdx.x];
                dbpronr = GetSSSSTableIndex( j/*dbsssCache[1][threadIdx.x]*/, f2 );//reuse registers
                score2 += f1 * f2 * 
#ifdef SSSS_SCORESP1E1
                    SerializedScoresSM<CUBSM_TYPE>::GetScoreP1E1( ssssCache, qpronr, dbpronr );
#else
                    //ssss.GetScoreP1( qpronr, dbpronr, 0.f, 0.f, e1, dbenoCache[1][threadIdx.x] );
                    ssss.GetScoreP1E1( qpronr, dbpronr, 0.f, 0.f, 0.f, 0.f );
#endif
            }
        }
    }

#ifdef TEST_CalcSM_SSSS_SMEMUnroll2x
    size_t querypos = row + querposoffset;
    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
# if !defined(SSSS_SCORESP1E1)
    MYASSERT( qrenoCache == ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr], "Inconsistency.");
    MYASSERT( dbenoCache[0][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO]))[dbpronr], "Inconsistency.");
# endif
    //MYASSERT( qrsssCache[threadIdx.y] == ((CHTYPE*)(dc_pm2dvfields_[pmv2DSSstate]))[querypos], "Inconsistency.");
    //MYASSERT( dbsssCache[0][threadIdx.x] == ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSstate]))[db1pos], "Inconsistency.");
    for( int i = 0; i < pmv2DNoSSSps; i++ ) {
        MYASSERT( qrposCache[i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DSSsprbs+i]))[querypos], "Inconsistency.");
        MYASSERT( dbposCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DSSsprbs+i]))[db1pos], "Inconsistency.");
    }
    if( col2 < (ndb1poss + ndbCposs)) {
# if !defined(SSSS_SCORESP1E1)
        size_t dbpronr2 = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        MYASSERT( dbenoCache[1][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO]))[dbpronr2], "Inconsistency.");
# endif
        //MYASSERT( dbsssCache[1][threadIdx.x] == ((CHTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSstate]))[db1pos2], "Inconsistency.");
        for( int i = 0; i < pmv2DNoSSSps; i++ ) {
            MYASSERT( dbposCache[pmv2DNoSSSps+i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DSSsprbs+i]))[db1pos2], "Inconsistency.");
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

    f1 = qrposCache[pmv2DNoSSSps-1][threadIdx.y];

    row = row * (ndb1poss + ndbCposs + dbxpad);
    
    //perform coalescent write of scores
    if( -99.f < score1 ) {
        //{{NOTE: test (wrt sens and aq) and REMOVE the block below; CHECKED!
        f2 = dbposCache[pmv2DNoSSSps-1][threadIdx.x];
        score1 += (f1 + 0.1f) * (f2 + 0.1f) * 0.1111111f;//(/9)
        //}}
        score1 -= CONSTSSSSSHIFT;
        score1 *= ssswgt;
        //atomicAdd is faster than += when we need coalescent write performed once
        atomicAdd( &outscores[row + col], score1 );
        outmodscores[row + col] = score1;
    }

    if( col2 < (ndb1poss + ndbCposs) && -99.f < score2 ) {
        //{{NOTE: test (wrt sens and aq) and REMOVE the block below; CHECKED!
        f2 = dbposCache[pmv2DNoSSSps+pmv2DNoSSSps-1][threadIdx.x];
        score2 += (f1 + 0.1f) * (f2 + 0.1f) * 0.1111111f;//(/9)
        //}}
        score2 -= CONSTSSSSSHIFT;
        score2 *= ssswgt;
        atomicAdd( &outscores[row + col2], score2 );
        outmodscores[row + col2] = score2;
    }
}

