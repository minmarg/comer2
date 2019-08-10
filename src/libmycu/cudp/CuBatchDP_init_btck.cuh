/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchDP_init_btck_h__
#define __CuBatchDP_init_btck_h__

#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/btckcoords.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "CuBatchDP_com.h"

//device functions for executing dynamic programming

__global__ void ExecDP_Btck_Unroll32x(
    uint blkdiagnum,
    uint lastydiagnum,
    uint ndb1pros,
    uint querprosOmtd, uint ndb1prosOmtd, uint ndbCprosOmtd,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ scores,
    CUBDP_TYPE* __restrict__ tmpdpdiagbuffers,
    CUBDP_TYPE* __restrict__ tmpdpbotbuffer,
    uint* __restrict__ maxscoordsbuf,
    char* __restrict__ btckdata
);

// =========================================================================
// -------------------------------------------------------------------------
// dpmaxandcoords: update maximum value along with the coordinates of this 
// value
template<typename T>
__device__ __forceinline__ 
void dpmaxandcoords( T& a, T b, uint& xy, uint x, uint y )
{ 
    if( a < b ) { 
        a = b;
        xy = CombineCoords(x,y);
    }
}

// =========================================================================
// -------------------------------------------------------------------------
// GetDbProfileCorners: get left and right position of a right and left 
// triangular corner respectively in the DP matrix between the query and a 
// db profile;
// nqyposs, query length;
// dbpro1len, db profile length;
// y, query position under processing;
// left, left position of the db profile for given y;
// right, right position of the db profile for given y;
//
__device__ __forceinline__ 
void GetDbProfileCorners( 
    int nqyposs,
    int dbpro1len,
    int y,
    int* left,
    int* right )
{
    int maxlen = (int)rintf(CUDP_CORNERLENPRC * (float)myhdmax(nqyposs, dbpro1len));
    if( CUDP_MAXCORNERLEN < maxlen )
        maxlen = CUDP_MAXCORNERLEN;

    //set boundaries of DP matrix triangles
    *left = 0; *right = dbpro1len-1;
    if( nqyposs - 1 < y + maxlen )
        *left += maxlen + y - nqyposs + 1;
    if( y < maxlen )
        *right -= maxlen - y;
}

// -------------------------------------------------------------------------
// GetDbProfileCorners: get left and right position of a right and left 
// triangular corner respectively in the DP matrix between the query and a 
// db profile;
// nqyposs, query length;
// dbpro1len, db profile length;
// y1, upper query position under processing;
// y2, lower query position under processing;
// left, left position of the db profile for given y;
// right, right position of the db profile for given y;
//
__device__ __forceinline__ 
void GetDbProfileCorners( 
    int nqyposs,
    int dbpro1len,
    int y1,
    int y2,
    int* left,
    int* right )
{
    int maxlen = (int)rintf(CUDP_CORNERLENPRC * (float)myhdmax(nqyposs, dbpro1len));
    if( CUDP_MAXCORNERLEN < maxlen )
        maxlen = CUDP_MAXCORNERLEN;
    //set boundaries of DP matrix triangles
    *left = 0; *right = dbpro1len-1;
    if( nqyposs - 1 < y1 + maxlen )
        *left += maxlen + y1 - nqyposs + 1;
    if( y2 < maxlen )
        *right -= maxlen - y2;
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// GetCornerMaxlen: get the maximum length of the distant corners to be 
// ignored;
// nqyposs, query length;
// dbpro1len, db profile length;
//
__device__ __forceinline__ 
int GetCornerMaxlen( 
    int nqyposs,
    int dbpro1len )
{
    int edgelen = (int)rintf(CUDP_CORNERLENPRC * (float)myhdmax(nqyposs, dbpro1len));
    if( CUDP_MAXCORNERLEN < edgelen )
        edgelen = CUDP_MAXCORNERLEN;
    return edgelen;
}

// -------------------------------------------------------------------------
// CellYXinCorners: check whether the DP cell is in the triangular 
// area of the distant corners;
// nqyposs, query length;
// dbpro1len, db profile length;
// edgelen, length of the calculated triangular corner edge;
// y, query position in the DP matrix;
// x, db profile position in the DP matrix;
//
__device__ __forceinline__ 
int CellYXinCorners( 
    int nqyposs,
    int dbpro1len,
    int edgelen,
    int y,
    int x )
{
    return 
        (nqyposs <= y + edgelen && x <= edgelen + y - nqyposs) ||
        (dbpro1len - edgelen + y <= x && y < edgelen);
}

// -------------------------------------------------------------------------
// CellYXinValidArea: check whether the DP cell is in valid area of the 
// DP matrix;
// nqyposs, query length;
// dbpro1len, db profile length;
// edgelen, length of the calculated triangular corner edge;
// y, query position in the DP matrix;
// x, db profile position in the DP matrix;
//
__device__ __forceinline__ 
int CellYXinValidArea( 
    int nqyposs,
    int dbpro1len,
    int edgelen,
    int y,
    int x )
{
    return 
        !CellYXinCorners(nqyposs, dbpro1len, edgelen, y, x) &&
        0 <= y && y < nqyposs && 0 <= x && x < dbpro1len;
}

#endif//__CuBatchDP_init_btck_h__
