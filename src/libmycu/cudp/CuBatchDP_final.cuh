/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_final_h__
#define __CuBatchDP_final_h__

#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"

//make an unsorted list of profile attributes
#define UNSRT_PASSED_PROFILES

//attributes of profiles passed through the DP phase
enum TPassedPros {
    dpppNPassedPros,//number of profiles passed to a further stage of processing
    dpppNPosits,//total number of positions of passed profiles
    dpppMaxProLen,//maximum profile length over passed profiles
    nTPassedPros
};

//device variable: attributes of passed profiles
extern __device__ unsigned int d_gDPPassedPros[nTPassedPros];

//device functions for finalizing dynamic programming

__global__ void FinalizeDP(
    CUBSM_TYPE scorethld,
    uint ndb1pros, uint ndbCpros,
    uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores,
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata
);

// =========================================================================
// -------------------------------------------------------------------------
// dpfinmaxndx: find maximum value between two values and save its index
template<typename T>
__device__ __forceinline__ 
void dpfinmaxndx(
    T* __restrict__ maxscCache,
    int* __restrict__ maxscndxCache,
    int ndx1,
    int ndx2 )
{
    //to avoid race conditions
    T valatndx2 = maxscCache[ndx2];
    __syncwarp();
    if( maxscCache[ndx1] < valatndx2 ) { 
        maxscCache[ndx1] = valatndx2;
        maxscndxCache[ndx1] = maxscndxCache[ndx2];
    }
//     else {
//         maxscndxCache[ndx1] = ndx1;
//     }
}

template<typename T>
__device__ __forceinline__ 
void dpfinmaxndxinit(
    T* __restrict__ maxscCache,
    int* __restrict__ maxscndxCache,
    int ndx1,
    int ndx2 )
{
    T valatndx2 = maxscCache[ndx2];
    __syncwarp();
    //NOTE: not using sync as long as ndx2-ndx1>=warp_size
    if( maxscCache[ndx1] < valatndx2 ) { 
        maxscCache[ndx1] = valatndx2;
        maxscndxCache[ndx1] = ndx2;
    } else {
        maxscndxCache[ndx1] = ndx1;
    }
}

// -------------------------------------------------------------------------

#endif//__CuBatchDP_final_h__
