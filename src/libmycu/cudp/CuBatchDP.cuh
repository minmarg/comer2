/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_h__
#define __CuBatchDP_h__

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/msg.h"
#include "libpro/srcpro/MOptions.h"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchDP
// Batch computation of dynamic programming for multiple profile-profile 
// pairs using the CUDA architecture
//
class CuBatchDP
{
public:
    CuBatchDP();
    virtual ~CuBatchDP();

    void PerformCompleteDynProgDevice(
        cudaStream_t streamproc,
        CUBSM_TYPE scorethld,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t querprosOmtd,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t dbpro1len,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        CUBSM_TYPE* scores,
        CUBDP_TYPE* tmpdpdiagbuffers, 
        CUBDP_TYPE* tmpdpbotbuffer,
        unsigned int* maxcoordsbuf,
        char* btckdata,
        //output values for host:
        unsigned int* globvarsbuf
    );

    void PerformDynProgDevice(
        cudaStream_t streamproc,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t querprosOmtd,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t dbpro1len,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        CUBSM_TYPE* scores,
        CUBDP_TYPE* tmpdpdiagbuffers,
        CUBDP_TYPE* tmpdpbotbuffer,
        unsigned int* maxcoordsbuf,
        char* btckdata
    );

    void FinalizeDynProgDevice(
        cudaStream_t streamproc,
        CUBSM_TYPE scorethld,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t dbpro1len,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        CUBSM_TYPE* scores,
        CUBDP_TYPE* tmpdpdiagbuffers,
        CUBDP_TYPE* tmpdpbotbuffer,
        unsigned int* maxcoordsbuf,
        char* btckdata,
        unsigned int* globvarsbuf
    );

    void FinalizeALNDynProgDevice(
        cudaStream_t streamproc,
        float logevthld,
        size_t ndb1pros,
        size_t querprosOmtd,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        CUBSM_TYPE* scores,//[in]
        CUBDP_TYPE* tmpdpdiagbuffers,//[in]
        float* dp2alndatbuffers,//[in/out]
        unsigned int* attrpassed,
        size_t dbxpad2,
        size_t cbdbalnlen2,
        unsigned int* maxcoordsbuf,//[in]
        char* btckdata,//[in]
        char* outalns//[out]
    );

    static unsigned int GetMaxRegBlockDiagonalElems( 
        size_t maxdiagelems, size_t querylen , size_t blockwidth );

    static unsigned int GetMaxBlockDiagonalElems( 
        size_t dblen, size_t querylen, size_t blockwidth, size_t blocklen );

    static unsigned int GetMaxBlockDiagonalElemsBlock2( 
        size_t dblen, size_t querylen, size_t blockwidth, size_t blocklen );

protected:
//     explicit CuBatchDP();

private:
};

// -------------------------------------------------------------------------
// INLINES ...
//
// GetMaxRegBlockDiagonalElems: get maximum number of elements of a regular 
// block diagonal;
// NOTE: we consider oblique block diagonals as opposed to rectangular block 
// diagonals; neighbour blocks in a diagonal share an edge (top or bottom);
// It equals the minimum of [l/w] and [(N-g)/w]+1 (>=[N/w]), where N is 
// maxdiagelems, w is blockwidth, g=1 if l>=w and g=w-l+1 otherwise; l, 
// query length; [], ceil rounding;
// this follows from the correspondence between a block row and a oblique 
// block diagonal;
// maxdiagelems, maximum number of elements in a regular diagonal 
// (of length 1);
// querylen, query length;
// blockwidth, oblique block's width (e.g., 32);
//
inline
unsigned int CuBatchDP::GetMaxRegBlockDiagonalElems( 
    size_t maxdiagelems, size_t querylen, size_t blockwidth )
{
    const size_t blkdiagsalongquery = (querylen + blockwidth - 1) / blockwidth;
    size_t g = 1;
    if( querylen < blockwidth)
        g = blockwidth - querylen + 1;
    size_t maxregblkdiagelems = 1;
    if( g < maxdiagelems )
        maxregblkdiagelems = ((maxdiagelems-g) + blockwidth - 1) / blockwidth + 1;
    return (unsigned int)SLC_MIN( blkdiagsalongquery, maxregblkdiagelems );
}

// -------------------------------------------------------------------------
// GetMaxBlockDiagonalElems: get maximum number of elements (blocks) a block 
// diagonal can contain, where blocks in a diagonal share a ppoint 
// (top-right or bottom-left);
// NOTE: consider oblique block diagonals as opposed to rectangular block 
// diagonals;
// It equals the minimum of [l/w] and [(x+w)/(b+w)], where x is dblen, 
// l is the query length, w is blockwidth (the length of an oblique edge),
// b is the block height; [], ceil rounding;
// dblen, number of db profile positions in processing (x coord.);
// querylen, query length (y coord.);
// blockwidth, oblique block's width (e.g., 32);
// blocklen, block's length (e.g., 16, 32);
inline
unsigned int CuBatchDP::GetMaxBlockDiagonalElems( 
    size_t dblen, size_t querylen, size_t blockwidth, size_t blocklen )
{
    const size_t blocksalongquery = (querylen + blockwidth - 1) / blockwidth;
    const size_t bpw = blockwidth + blocklen;
    size_t maxblkdiagelems = ((dblen + bpw) + bpw - 1) / bpw;
    return (unsigned int)SLC_MIN( blocksalongquery, maxblkdiagelems );
}

// -------------------------------------------------------------------------
// GetMaxBlockDiagonalElemsBlock2: get maximum number of RECTANGULAR 
// blocks in one diagonal given the score matrix dimensions
// dblen, number of db profile positions in processing (x coord.);
// querylen, query length (y coord.);
// blockwidth, rectangular block's width, or height (e.g., 32);
// blocklen, block's length (e.g., 16, 32);
inline
unsigned int CuBatchDP::GetMaxBlockDiagonalElemsBlock2( 
        size_t dblen, size_t querylen, size_t blockwidth, size_t blocklen )
{
    const size_t blocksalongquery = (querylen + blockwidth - 1) / blockwidth;
    const size_t blocksalongdbpro = (dblen + blocklen - 1) / blocklen;
    return (unsigned int)SLC_MIN( blocksalongquery, blocksalongdbpro );
}

#endif//__CuBatchDP_h__
