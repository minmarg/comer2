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
#include "CuBatchScoreMatrix_hdp1s.cuh"
#include "CuBatchScoreMatrix_Init_hdp1s.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// hdp1sTexo, texture object encapsulating serialized HDP scores
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

// -------------------------------------------------------------------------
// CalcSM_Init_HDP1S_SMEMUnroll2x: device code for calculating initial and 
// HDP scores using shared memory and twofold unrolling along the x axis;
// each block processes two blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses more shared memory but allows for greater occupancy and 
// reduced load bank conflicts!
// NOTE: For HDP scores, it also uses texture memory (for big data), which 
// due to random acces provides greater efficiency in comparison to other 
// types of memory;
// NOTE: output pointers should be aligned!
// NOTE: results are not written to outmodscores;
// 
__global__ void CalcSM_Init_HDP1S_SMEMUnroll2x(
    cudaTextureObject_t hdp1sTexo,
    SerializedScoresAttr h1attr,
    float hdp1swgt,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    const CUBSM_TYPE CONSTINITSHIFT,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ /*outmodscores*/ )
{
    __shared__ FPTYPE qrproCache[pmv2DNoElems],//cache for query background probabilities
            dbproCache[TIMES2(pmv2DNoElems)][SMINIT_2DCACHE_DIM];//cache for background probabilities of db profiles
    __shared__ FPTYPE 
            qrposCache[pmv2DNoElems][SMINIT_2DCACHE_DIM],//cache for a tile of query positions
            dbposCache[TIMES2(pmv2DNoElems)][SMINIT_2DCACHE_DIM];//cache for a tile of db profile positions
    //
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
    uint pronr;
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
                pronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
                //read only one element and once per block line (thread.y)
                qrproCache[threadIdx.y] = 
                    ((FPTYPE*)(dc_pm2dvfields_[pps2DBkgPrbs+threadIdx.y]))[pronr];
                if( threadIdx.y < 1 ) {
                    //read query ENO
                    //read only one element per block (blockDim.y x blockDim.x)
                    qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[pronr];
                }
            }
            if( threadIdx.y < 1 ) {
                qrhdp1prbCache[threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[pmv2DHDP1prb]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
                qrhdp1ndxCache[threadIdx.x] = 
                ((INTYPE*)(dc_pm2dvfields_[pmv2DHDP1ind]))[
                    blockbeg_y + querposoffset + threadIdx.x
                ];
            }
        }
    }
    //
    if( col < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DTrgFrqs+threadIdx.y]))[db1pos];
            pronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
            dbproCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DBkgPrbs+threadIdx.y]))[pronr];
            if( threadIdx.y < 1 ) {
                dbhdp1prbCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1prb/*+threadIdx.y*/]))[db1pos];
                dbhdp1ndxCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DHDP1ind/*+threadIdx.y*/]))[db1pos];
                dbenoCache[0/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx+pps2DENO/*+threadIdx.y*/]))[pronr];
            }
        }
    }
    if( col2 < (ndb1poss + ndbCposs)) {
        if( threadIdx.y < pmv2DNoElems ) {
            dbposCache[pmv2DNoElems+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DTrgFrqs+threadIdx.y]))[db1pos2];
            pronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
            dbproCache[pmv2DNoElems+threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DBkgPrbs+threadIdx.y]))[pronr];
            if( threadIdx.y < 1 ) {
                dbhdp1prbCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1prb/*+threadIdx.y*/]))[db1pos2];
                dbhdp1ndxCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DHDP1ind/*+threadIdx.y*/]))[db1pos2];
                dbenoCache[1/*+threadIdx.y*/][threadIdx.x] = 
                    ((FPTYPE*)(dc_pm2dvfields_[dbfldsndx2+pps2DENO/*+threadIdx.y*/]))[pronr];
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
    pronr = col2 < (ndb1poss + ndbCposs);
    //
    float score1, score2;
    float s1, s2;

// score1=0.0f;score2=0.0f;
    CalcInitScoreSMEMUnroll2x( 
        qrproCache,
        dbproCache,
        qrposCache,
        dbposCache,
        CONSTINITSHIFT,
        pronr/*col2 < (ndb1poss + ndbCposs)*/,
        score1, score2
    );

    CalcHDP1ScoreSMEMUnroll2x( 
        h1attr,
        hdp1swgt,
        hdp1sTexo,
        qrenoCache,
        dbenoCache,
        qrhdp1prbCache,
        dbhdp1prbCache,
        qrhdp1ndxCache,
        dbhdp1ndxCache,
        pronr/*col2 < (ndb1poss + ndbCposs)*/,
        s1, s2
    );

    /*if( s1 ) */score1 += s1;

    row = row * (ndb1poss + ndbCposs + dbxpad);

    //perform coalescent write of scores
    outscores[row + col] = score1;
    if( pronr/*col2 < (ndb1poss + ndbCposs)*/) {
        /*if( s2 ) */score2 += s2;
        outscores[row + col2] = score2;
    }
}
