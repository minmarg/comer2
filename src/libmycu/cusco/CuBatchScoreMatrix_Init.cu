/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_Init.cuh"

__device__ unsigned long long d_testvar = 0;

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
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

// #define TEST_CalcSMInitSMEMUnroll2x
// -------------------------------------------------------------------------
// CalcSMInitSMEMUnroll2x: device code for calculating initial score 
// matrix using shared memory and twofold unrolling along the x axis;
// each block processes two blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses more shared memory but allows for greater occupancy and 
// reduced load bank conflicts!
// NOTE: output pointers should be aligned!
// NOTE: results are not written to outmodscores;
// 
__global__ void CalcSM_Init_SMEMUnroll2x(
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* /*outmodscores*/ )
{
    __shared__ FPTYPE qrproCache[pmv2DNoElems],//cache for query background probabilities
            dbproCache[TIMES2(pmv2DNoElems)][SMINIT_2DCACHE_DIM];//cache for background probabilities of db profiles
    __shared__ FPTYPE 
            qrposCache[pmv2DNoElems][SMINIT_2DCACHE_DIM],//cache for a tile of query positions
            dbposCache[TIMES2(pmv2DNoElems)][SMINIT_2DCACHE_DIM];//cache for a tile of db profile positions
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
    //for comments, see CalcSMInitSMEM
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        if( threadIdx.y < pmv2DNoElems ) {
            qrposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            //read query background probabilities
            //NOTE: valid when processing one query at a time
            if( threadIdx.x < 1 ) {
                qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                //read only one element and once per block line (thread.y)
                qrproCache[threadIdx.y] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+threadIdx.y]))[qpronr];
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+threadIdx.y]))[db1pos];
            dbproCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+threadIdx.y]))[dbpronr];
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[pmv2DNoElems+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DTrgFrqs+threadIdx.y]))[db1pos2];
            dbproCache[pmv2DNoElems+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DBkgPrbs+threadIdx.y]))[dbpronr];
        }
    }
    //the warp reads data written by other warps, sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    float score1 = 0.0f, score2 = 0.0f;
    float f1, f2;
    //
    //
    #pragma unroll
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        //exchanged access for qrposCache
        f1 = qrposCache[i][threadIdx.y];
        f2 = qrproCache[i];
        score1 += __fdividef( 
            ( f1 * dbposCache[i][threadIdx.x]), 
            ( 0.5f * (f2 + dbproCache[i][threadIdx.x]))
        );
        if( col2 < (ndb1poss + ndbCposs)) {
            score2 += __fdividef( 
                ( f1 * dbposCache[pmv2DNoElems+i][threadIdx.x]), 
                ( 0.5f * (f2 + dbproCache[pmv2DNoElems+i][threadIdx.x]))
            );
        }
    }

#ifdef TEST_CalcSMInitSMEMUnroll2x
    size_t querypos = row + querposoffset;
    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        MYASSERT( qrproCache[i] == ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr], "Inconsistency.");
        MYASSERT( qrposCache[i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos], "Inconsistency.");
        MYASSERT( dbproCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr], "Inconsistency.");
        MYASSERT( dbposCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+i]))[db1pos], "Inconsistency.");
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        size_t dbpronr2 = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        for( uint i = 0; i < pmv2DNoElems; i++ ) {
            MYASSERT( dbproCache[pmv2DNoElems+i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DBkgPrbs+i]))[dbpronr2], "Inconsistency.");
            MYASSERT( dbposCache[pmv2DNoElems+i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DTrgFrqs+i]))[db1pos2], "Inconsistency.");
        }
        atomicAdd(&d_testvar,1UL);
    }
    atomicAdd(&d_testvar,1UL);
    if( querposoffset < 1)
    if((size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset) <= d_testvar )
        printf(" Counter %lu col2= %u qcnt= %u(%lu) dbcnt= %u(%lu)   ss= %lu\n", 
               d_testvar, col2, nqyposs, querposoffset, ndb1poss + ndbCposs, bdb1posoffset + bdbCposoffset,
               (size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset)
        );
#endif

    //perform coalescent write of scores
    MYASSERT(score1>0.0f,"Invalid score1.");
    score1 = __logf(score1) - CONSTINITSHIFT;
    row = row * (ndb1poss + ndbCposs + dbxpad);
    outscores[row + col] = score1;
    //outmodscores[row + col] = score1;

    if( col2 < (ndb1poss + ndbCposs)) {
        MYASSERT(score2>0.0f,"Invalid score2.");
        score2 = __logf(score2) - CONSTINITSHIFT;
        outscores[row + col2] = score2;
        //outmodscores[row + col2] = score2;
    }
}


// #define TEST_CalcSMInitSMEMUnroll2
// -------------------------------------------------------------------------
// CalcSMInitSMEMUnroll2: device code for calculating initial score 
// matrix using shared memory and twofold unrolling;
// each block processes two blocks actually; 
// db profile data remains the same for these two spatially prallel blocks
// 
__global__ void CalcSMInitSMEMUnroll2(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* /*outmodscores*/ )
{
    __shared__ FPTYPE qrproCache[TIMES2(pmv2DNoElems)],//cache for query background probabilities
            dbproCache[pmv2DNoElems][SMINIT_2DCACHE_DIM];//cache for background probabilities of db profiles
    __shared__ FPTYPE 
            qrposCache[TIMES2(pmv2DNoElems)][SMINIT_2DCACHE_DIM],//cache for a tile of query positions
            dbposCache[pmv2DNoElems][SMINIT_2DCACHE_DIM];//cache for a tile of db profile positions
    //
    uint blockbeg_y = blockIdx.y * blockDim.y * 2;
    uint row = blockbeg_y + threadIdx.y;
    uint row2 = row + blockDim.y;
    uint col = blockIdx.x * blockDim.x + threadIdx.x;
    //size_t querypos = row + querposoffset;
    //size_t querypos2 = row2 + querposoffset;
    uint qpronr, qpronr2, dbpronr;
    //physical indices:
    size_t db1pos;// = col < ndb1poss ? col + bdb1posoffset: (col - ndb1poss + bdbCposoffset)/*jump to section ndbCposs*/;
    uint dbfldsndx;// = col < ndb1poss ? pmv2DTotFlds: TIMES2(pmv2DTotFlds);//db index for dc_pm2dvfields_
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    //
    //for comments, see CalcSMInitSMEM
    if( blockbeg_y + threadIdx.x < nqyposs ) 
    {
        if( threadIdx.y < pmv2DNoElems ) {
            qrposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            //read query background probabilities
            //NOTE: valid when processing one query at a time
            if( threadIdx.x < 1 ) {
                qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                //read only one element and once per block line (thread.y)
                qrproCache[threadIdx.y] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+threadIdx.y]))[qpronr];
            }
        }
    }
    if( blockbeg_y + threadIdx.x + blockDim.y < nqyposs ) 
    {
        if( threadIdx.y < pmv2DNoElems ) {
            qrposCache[pmv2DNoElems+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x + blockDim.y 
                ];
            //read query background probabilities
            //NOTE: valid when processing one query at a time
            if( threadIdx.x < 1 ) {
                qpronr2 = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+blockDim.y+querposoffset];
                //read only one element and once per block line (thread.y)
                qrproCache[pmv2DNoElems+threadIdx.y] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+threadIdx.y]))[qpronr2];
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+threadIdx.y]))[db1pos];
            dbproCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+threadIdx.y]))[dbpronr];
        }
    }
    //the warp reads data written by other warps, hence sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    float score1 = 0.0f, score2 = 0.0f;
    //
    //
    //#pragma unroll
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        //exchanged access for qrposCache
        score1 += __fdividef( 
            ( qrposCache[i][threadIdx.y] * dbposCache[i][threadIdx.x]), 
            ( 0.5f * (qrproCache[i] + dbproCache[i][threadIdx.x]))
        );
        if( row2 < nqyposs ) {
            score2 += __fdividef( 
                ( qrposCache[pmv2DNoElems+i][threadIdx.y] * dbposCache[i][threadIdx.x]), 
                ( 0.5f * (qrproCache[pmv2DNoElems+i] + dbproCache[i][threadIdx.x]))
            );
        }
    }
#ifdef TEST_CalcSMInitSMEMUnroll2
    size_t querypos = row + querposoffset;
    size_t querypos2 = row2 + querposoffset;
    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    qpronr2 = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos2];
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        MYASSERT( qrproCache[i] == ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr], "Inconsistency.");
        MYASSERT( qrposCache[i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos], "Inconsistency.");
        MYASSERT( dbproCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr], "Inconsistency.");
        MYASSERT( dbposCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+i]))[db1pos], "Inconsistency.");
    }
    if( row2 < nqyposs ) {
        for( uint i = 0; i < pmv2DNoElems; i++ ) {
            MYASSERT( qrproCache[pmv2DNoElems+i] == ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr2], "Inconsistency.");
            MYASSERT( qrposCache[pmv2DNoElems+i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos2], "Inconsistency.");
        }
        atomicAdd(&d_testvar,1UL);
    }
    atomicAdd(&d_testvar,1UL);
    if( querposoffset < 1)
    if((size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset) <= d_testvar )
        printf(" Counter %lu row2= %u qcnt= %u(%lu) dbcnt= %u(%lu)   ss= %lu\n", 
               d_testvar, row2, nqyposs, querposoffset, ndb1poss + ndbCposs, bdb1posoffset + bdbCposoffset,
               (size_t)(nqyposs + querposoffset) * (size_t)(ndb1poss + ndbCposs + bdb1posoffset + bdbCposoffset)
        );
#endif

    //perform coalescent write of scores
    MYASSERT(score1>0.0f,"Invalid score1.");
    score1 = __logf(score1) - CONSTINITSHIFT;
    row = row * (ndb1poss + ndbCposs);
    outscores[row + col] = score1;
    //outmodscores[row + col] = score1;

    if( row2 < nqyposs ) {
        MYASSERT(score2>0.0f,"Invalid score2.");
        score2 = __logf(score2) - CONSTINITSHIFT;
        row2 = row2 * (ndb1poss + ndbCposs);
        outscores[row2 + col] = score2;
        //outmodscores[row2 + col] = score2;
    }
}


// #define TEST_CalcSMInitSMEM
// -------------------------------------------------------------------------
// CalcSMInitSMEM: device code for computing initial score 
// matrix using shared memory
// 
__global__ void CalcSMInitSMEM(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* /*outmodscores*/ )
{
    __shared__ FPTYPE qrproCache[pmv2DNoElems],//cache for query background probabilities
            dbproCache[pmv2DNoElems][SMINIT_2DCACHE_DIM];//cache for background probabilities of db profiles
    __shared__ FPTYPE 
            qrposCache[pmv2DNoElems][SMINIT_2DCACHE_DIM],//cache for a tile of query positions
            dbposCache[pmv2DNoElems][SMINIT_2DCACHE_DIM];//cache for a tile of db profile positions
    //
    uint blockbeg_y = blockIdx.y * blockDim.y;
    uint row = blockbeg_y + threadIdx.y;
    uint col = blockIdx.x * blockDim.x + threadIdx.x;
    //size_t querypos = row + querposoffset;
    uint qpronr, dbpronr;
    //physical indices:
    size_t db1pos;// = col < ndb1poss ? col + bdb1posoffset: (col - ndb1poss + bdbCposoffset)/*jump to section ndbCposs*/;
    uint dbfldsndx;// = col < ndb1poss ? pmv2DTotFlds: TIMES2(pmv2DTotFlds);//db index for dc_pm2dvfields_
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    //
    //if( row < nqyposs ) {
        //size_t querypos = row + querposoffset;
        //qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
        //strided and therefore hugely inefficient access for
        // __shared__ FPTYPE qrposCache[SMINIT_2DCACHE_DIM][pmv2DNoElems];
        // (access type qrposCache[threadIdx.y][i]):
        //if( threadIdx.x < pmv2DNoElems )
        //    qrposCache[threadIdx.y][threadIdx.x] = 
        //        ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+threadIdx.x]))[querypos];
    //}
    if( blockbeg_y + threadIdx.x < nqyposs ) {
        //on the contrary, this implementation guarantees coalescent memory access
        //and shared memory write with no bank conflicts
        //NOTE: block x and y dimensions should be equal;
        //NOTE also: calculating matrix requires accessing query and db profile 
        // positions multiple times, and this reduces effective bandwidth; 
        // # transactions depends on the cache line (32 or 128B)
        if( threadIdx.y < pmv2DNoElems ) {
            qrposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+threadIdx.y]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            //read query background probabilities;
            //NOTE: valid when processing one query at a time
            if( threadIdx.x < 1 ) {
                qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                //read only one element and once per block line (thread.y)
                qrproCache[threadIdx.y] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+threadIdx.y]))[qpronr];
            }
        }
    }
    if( col < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+threadIdx.y]))[db1pos];
            //read db profile background probabilities; 
            //usually, warps will access the same profile vector of probabilities, but it is
            //possible that these vectors will represent different profiles;
            //when accessing the same vector, L2 cache will be exploited and the performance 
            //penalty will be likely diminished much;
            //note that this is not large loss in bandwidth, as there are no strided accesses
            dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
            dbproCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+threadIdx.y]))[dbpronr];
        }
    }
    //the warp reads data written by other warps, hence sync
    __syncthreads();
    //
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    //
    float score = 0.0f;
    //
    //pragma unroll won't work because of dependencies in the warp
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        //exchanged access for qrposCache
        score += __fdividef( 
            ( qrposCache[i][threadIdx.y] * dbposCache[i][threadIdx.x]), 
            ( 0.5f * (qrproCache[i] + dbproCache[i][threadIdx.x]))
        );
    }
#ifdef TEST_CalcSMInitSMEM
    size_t querypos = row + querposoffset;
    qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        MYASSERT( qrproCache[i] == ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr], "Inconsistency.");
        MYASSERT( qrposCache[i][threadIdx.y] == ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos], "Inconsistency.");
        MYASSERT( dbproCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr], "Inconsistency.");
        MYASSERT( dbposCache[i][threadIdx.x] == ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+i]))[db1pos], "Inconsistency.");
    }
#endif

//     if( blockIdx.x == 4 && blockIdx.y == 5 ) {
//         printf(" block (%u,%u) thread (%u,%u) lane %u   row %u col %u   db1pos %lu(%lu)\n\n",
//             blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, lane_id(), row, col, db1pos, bdb1posoffset );
//     }
//     if( row == nqyposs-1 && col == ndb1poss-1 ) {
//         size_t querypos = row + querposoffset;
//         printf("CalcSMInitSMEM: block (%u,%u) thread (%u,%u)  "
//             "nqyposs %u ndb1poss %u ndbCposs %u  "
//             "row %u col %u  querypos %lu(%lu) db1pos %lu(%lu)\n\n", 
//             blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, 
//             nqyposs, ndb1poss, ndbCposs,  row, col,  querypos, querposoffset, db1pos, bdb1posoffset );
//         qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
//         dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr]);
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", qrproCache[i] );
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", ((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos]);
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", qrposCache[i][threadIdx.y] );
//         printf("\n\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr]);
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", dbproCache[i][threadIdx.x] );
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+i]))[db1pos]);
//         printf("\n");
//         for( uint i = 0; i < pmv2DNoElems; i++ )
//             printf( " %f", dbposCache[i][threadIdx.x] );
//         printf("\n\n");
//     }

    //perform coalescent write of scores
    MYASSERT(score>0.0f,"Invalid score.");
    score = __logf(score) - CONSTINITSHIFT;
    row = row * (ndb1poss + ndbCposs);
    outscores[row + col] = score;
    //outmodscores[row + col] = score;
}


// -------------------------------------------------------------------------
// CalcSMInit: naive implementation
//
__global__ void CalcSMInit(
    uint nqyposs, uint ndb1poss, uint ndbCposs, 
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores,
    CUBSM_TYPE* /*outmodscores*/ )
{
    uint ndx = blockIdx.x * blockDim.x + threadIdx.x;
    uint row = ndx / (ndb1poss + ndbCposs);
    uint col = ndx % (ndb1poss + ndbCposs);
    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        return;
    size_t querypos = row + querposoffset;
    //physical indices:
    size_t db1pos;// = col < ndb1poss ? col + bdb1posoffset: (col - ndb1poss + bdbCposoffset)/*jump to section ndbCposs*/;
    uint dbfldsndx;// = col < ndb1poss ? pmv2DTotFlds: TIMES2(pmv2DTotFlds);//db index for dc_pm2dvfields_
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    //
    uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    uint dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    //
    float score = 0.0f;
    //
    #pragma unroll
    for( uint i = 0; i < pmv2DNoElems; i++ ) {
        score +=
        (((FPTYPE*)(dc_pm2dvfields_[pmv2DTrgFrqs+i]))[querypos]) * 
        (((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+i]))[db1pos]) /
        ( 0.5f *
            (
                ( ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr] ) +
                ( ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr] )
//                 __ldg( &((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+i]))[qpronr] ) +
//                 __ldg( &((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+i]))[dbpronr] )
            )
        );
    }

    //perform coalescent write of scores
    MYASSERT(score>0.0f,"Invalid score.");
    score = __logf(score) - CONSTINITSHIFT;
    row = row * (ndb1poss + ndbCposs);
    outscores[row + col] = score;
    //outmodscores[row + col] = score;
    //
}


// =========================================================================
// -------------------------------------------------------------------------
// KernelStreamedIterativeTest: simple test kernel
//
__global__ void KernelStreamedIterativeTest(
    uint serqposnr,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores )
{
    uint blockbeg_y = blockIdx.y * blockDim.y;
    uint row = blockbeg_y + threadIdx.y;
    const uint col = blockIdx.x * blockDim.x * 2 + threadIdx.x;//logical column
    const uint col2 = col + blockDim.x;//logical column
    //physical indices:
//     size_t db1pos;
//     size_t db1pos2;
//     uint dbfldsndx;
//     uint dbfldsndx2;
//     if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
//                             dbfldsndx = pmv2DTotFlds;
//     } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
//                             dbfldsndx = TIMES2(pmv2DTotFlds);
//     }
//     if( col2 < ndb1poss ) { db1pos2 = col2 + bdb1posoffset;
//                             dbfldsndx2 = pmv2DTotFlds;
//     } else {                db1pos2 = col2 - ndb1poss + bdbCposoffset;//jump to section ndbCposs
//                             dbfldsndx2 = TIMES2(pmv2DTotFlds);
//     }
    //
    float score1, score2;
    score1 = serqposnr * row + nqyposs + col + querposoffset; 
    score2 = serqposnr * row + nqyposs + col2 + querposoffset;
    //
    for(int i=0;i<32;i++) {
        score1 *= serqposnr;
        score2 *= serqposnr;
    }
    //perform coalescent write of scores
    row = row * (ndb1poss + ndbCposs + dbxpad);
    outscores[row + col] = score1;
    if( col2 < (ndb1poss + ndbCposs)) {
        outscores[row + col2] = score2;
    }
}

// -------------------------------------------------------------------------
// KernelStreamedIterativeTest: simple test kernel
//
__global__ void KernelParentStreamedIterativeTest(
    uint serqposnr,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    size_t querposoffset, size_t bdb1posoffset, size_t bdbCposoffset,
    CUBSM_TYPE* outscores )
{
    dim3 nthrds(SMINIT_2DCACHE_DIM,SMINIT_2DCACHE_DIM,1);
    dim3 nblcks = dim3(100/*4900*/,1,1);
    for( uint j = 0; j < 469; j++ ) {
        KernelStreamedIterativeTest<<<nblcks,nthrds>>>( 
            serqposnr,
            (uint)nqyposs, (uint)ndb1poss, (uint)ndbCposs, (uint)dbxpad,
            querposoffset, bdb1posoffset, bdbCposoffset,
            outscores );
        MYCUDACHECKLAST;
    }
}
