/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __btckcoords_h__
#define __btckcoords_h__

// -------------------------------------------------------------------------
// CombineCoords: combine (x,y) coordinates (profile-profile positions) into 
// one integer value;
// NOTE: the arguments x and y are supposed to contain 16-bit values!
__host__ __device__ __forceinline__
uint CombineCoords(uint x, uint y)
{
    return (x << 16) | y;//(y & 0xffff);
}
// GetCoordX and GetCoordY extract x and y coordinates from their 
// combination
__host__ __device__ __forceinline__
uint GetCoordX(uint xy)
{
    return (xy >> 16) & 0xffff;
}
__host__ __device__ __forceinline__
uint GetCoordY(uint xy)
{
    return xy & 0xffff;
}

#endif//__btckcoords_h__
