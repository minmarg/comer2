/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_init_h__
#define __CuBatchDP_init_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"

//inner dimension shift for shared memory arrays
#define CUDPIDS 1

//device functions for executing dynamic programming

__global__ void ExecDPUnroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores, 
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer 
);


// =========================================================================
// -------------------------------------------------------------------------
__device__ __forceinline__
CUBDP_TYPE DPLocGetLog0()
{
    return CUBDP_Q(-32768);
}

// -------------------------------------------------------------------------
// GetNoBlockDiagonalsProfile: get the number of block diagonals given the 
// lengths of query and db profile (see also CuBatchDP.cu)
// 
__device__ __forceinline__
int GetNoBlockDiagonalsProfile(
    int nqyposs,
    int dbpro1len,
    int DIM_X,
    int DIM_D )
{
    //return ((dbpro1len + nqyposs - 1) + DIM_X-1)/DIM_X + (nqyposs - 1)/DIM_D;
    //REVISION: due to the positioning of the first block, the first 1-position diagonal of the 
    // first diagonal block is out of bounds: remove -1
    //TODO: make this calculation unique over the whole code
    return ((dbpro1len + nqyposs) + DIM_X-1)/DIM_X + (nqyposs - 1)/DIM_D;
}

// -------------------------------------------------------------------------
// GetMaqxNoIterations: get the maximum number of iterations for the 
// block to perform 
// 
__device__ __forceinline__
int GetMaqxNoIterations(
    int x, int y, int quelen, int dbprolen, int blklen = CUDP_2DCACHE_DIM_X )
{
    x = dbprolen - x;
    if( x <= 0 )
        return 0;
    y -= quelen-1;//NOTE
    if( 0 < y )
        x -= y;
    if( x <= 0 )
        return 0;
    if( blklen <= x )
        return blklen;
    return x;
}

// -------------------------------------------------------------------------
// DPLocinitTrnProbs: initialize transition probabilities
//
__device__ __forceinline__
void DPLocInitTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_D])
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         smtrpCache[i][threadIdx.x] = CUBDP_Q(-32768);
//     }
    smtrpCache[dptrMMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMIc][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMDp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIIc][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDMp][threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDDp][threadIdx.x] = DPLocGetLog0();
}

__device__ __forceinline__
void DPLocInitTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_DpX],
    int ndxoff )
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         smtrpCache[i][ndxoff+threadIdx.x] = CUBDP_Q(-32768);
//     }
    smtrpCache[dptrMMp][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMIc][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrMDp][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIMp][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrIIc][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDMp][ndxoff+threadIdx.x] = DPLocGetLog0();
    smtrpCache[dptrDDp][ndxoff+threadIdx.x] = DPLocGetLog0();
}

// -------------------------------------------------------------------------
// DPLocGetLogTrnProb: read transition probability value and take its 
// logarithm
__device__ __forceinline__
FPTYPE DPLocGetLogTrnProb( int fldsndx, int trn, int pos, FPTYPE fct = 1.0f )
{
    return ((FPTYPE*)(dc_pm2dvfields_[fldsndx+ptr2DTrnPrbs+trn]))[pos];
//     FPTYPE val = ((FPTYPE*)(dc_pm2dvfields_[fldsndx+ptr2DTrnPrbs+trn]))[pos];
//     if( val )
//         return __logf(val) * fct;
//     else
//         return DPLocGetLog0();
}

// -------------------------------------------------------------------------
// DPLocCacheTrnProbs: cache transition probabilities to SMEM at position 
// pos
//
__device__ __forceinline__
void DPLocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_D],
    int fldsndx,
    int pos )
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         if( i == dptrMIc || i == dptrIIc )
//             smtrpCache[i][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i+1, pos );
//         else
//             smtrpCache[i][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i, pos );
//     }
    smtrpCache[dptrMMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMMp, pos );
    smtrpCache[dptrMIc][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMIc, pos+1 );
    smtrpCache[dptrMDp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrMDp, pos );
    smtrpCache[dptrIMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrIMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrIIc][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrIIc, pos+1, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDMp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrDMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDDp][threadIdx.x] = DPLocGetLogTrnProb( fldsndx, dptrDDp, pos, CUDP_TRNPRB_LOGFCT );
}

// DPLocCacheTrnProbs: overloaded for caching CUDP_2DCACHE_DIMx2 positions
//
__device__ __forceinline__
void DPLocCacheTrnProbs( 
    FPTYPE (* __restrict__ smtrpCache)[CUDP_2DCACHE_DIM_DpX],
    int ndxoff,
    int fldsndx,
    int pos )
{
//     #pragma unroll
//     for( int i = 0; i < nTDPUsedProTrn; i++ ) {
//         if( i == dptrMIc || i == dptrIIc )
//             smtrpCache[i][ndxoff+threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i+1, pos );
//         else
//             smtrpCache[i][ndxoff+threadIdx.x] = DPLocGetLogTrnProb( fldsndx, i, pos );
//     }
    int ndx = ndxoff+threadIdx.x;
    smtrpCache[dptrMMp][ndx] = DPLocGetLogTrnProb( fldsndx, dptrMMp, pos );
    smtrpCache[dptrMIc][ndx] = DPLocGetLogTrnProb( fldsndx, dptrMIc, pos+1 );
    smtrpCache[dptrMDp][ndx] = DPLocGetLogTrnProb( fldsndx, dptrMDp, pos );
    smtrpCache[dptrIMp][ndx] = DPLocGetLogTrnProb( fldsndx, dptrIMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrIIc][ndx] = DPLocGetLogTrnProb( fldsndx, dptrIIc, pos+1, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDMp][ndx] = DPLocGetLogTrnProb( fldsndx, dptrDMp, pos, CUDP_TRNPRB_LOGFCT );
    smtrpCache[dptrDDp][ndx] = DPLocGetLogTrnProb( fldsndx, dptrDDp, pos, CUDP_TRNPRB_LOGFCT );
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// DPLocInitCache: initialize cache of a DP buffer
// shft, shift of inner-most index for writing values
//
template<unsigned int _2DCACHE_DIM_D>
__device__ __forceinline__
void DPLocInitCache( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_D+1],
    int shft = 0,
    CUBDP_TYPE initval = CUBDP_Q(0),
    CUBDP_TYPE initmmval = CUBDP_Q(0))
{
    smdpbufCache[dpdsssStateMM][threadIdx.x+shft] = initmmval;
    #pragma unroll
    for( int i = dpdsssStateMM+1; i < nTDPDiagScoreSubsections; i++ )
        smdpbufCache[i][threadIdx.x+shft] = initval;
}

template<unsigned int _2DCACHE_DIM_X>
__device__ __forceinline__
void DPLocInitCache( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_X],
    CUBDP_TYPE initval = CUBDP_Q(0),
    CUBDP_TYPE initmmval = CUBDP_Q(0))
{
    smdpbufCache[dpdsssStateMM][threadIdx.x] = initmmval;
    #pragma unroll
    for( int i = dpdsssStateMM+1; i < nTDPDiagScoreSubsections; i++ )
        smdpbufCache[i][threadIdx.x] = initval;
}

// -------------------------------------------------------------------------
// DPLocCacheBuffer: cache one of the DP buffers;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to read data from;
// x, x position at which the buffer has to be read;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index for writing values
//
template<unsigned int _2DCACHE_DIM_D>
__device__ __forceinline__
void DPLocCacheBuffer( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_D+1],
    const CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        smdpbufCache[i][threadIdx.x+shft] = gmemtmpbuffer[y + x];
    }
}

template<unsigned int _2DCACHE_DIM_X>
__device__ __forceinline__
void DPLocCacheBuffer( 
    CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_X],
    const CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        smdpbufCache[i][threadIdx.x] = gmemtmpbuffer[y + x];
    }
}

// -------------------------------------------------------------------------
// DPLocWriteBuffer: write one of the cahced DP buffers back to GMEM;
// smdpbufCache, SMEM cache;
// gmemtmpbuffer, address of the buffer to write data to;
// x, x position at which the data has to be written;
// y, initial y position in the buffer;
// stride, stride to refer to the same x position in the buffer;
// shft, shift of inner-most index of the cache
//
template<unsigned int _2DCACHE_DIM_D>
__device__ __forceinline__
void DPLocWriteBuffer( 
    const CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_D+1],
    CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride,
    int shft = 0 )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        gmemtmpbuffer[y + x] = smdpbufCache[i][threadIdx.x+shft];
    }
}

template<unsigned int _2DCACHE_DIM_X>
__device__ __forceinline__
void DPLocWriteBuffer( 
    const CUBDP_TYPE (* __restrict__ smdpbufCache)[_2DCACHE_DIM_X],
    CUBDP_TYPE* __restrict__ gmemtmpbuffer,
    int x,
    int y,
    int stride )
{
    #pragma unroll
    for( int i = 0; i < nTDPDiagScoreSubsections; i++, y += stride ) {
        gmemtmpbuffer[y + x] = smdpbufCache[i][threadIdx.x];
    }
}

#endif//__CuBatchDP_init_h__
