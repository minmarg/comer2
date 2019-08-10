/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __warpscan_h__
#define __warpscan_h__

// mywarpprefixsum: exclusive prefix sum perfomed in a warp
template <typename T> 
__device__ __forceinline__ 
T mywarpexcprefixsum( T reg )
{
    T tmp;
    int lid = threadIdx.x & 0x1f;//lane id
    tmp = __shfl_up_sync(0xffffffff, reg, 1); if (1 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 2); if (2 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 4); if (4 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 8); if (8 <= lid) reg += tmp;
    tmp = __shfl_up_sync(0xffffffff, reg, 16); if (16 <= lid) reg += tmp;
    return reg;
}

// mywarprevexcprefixsum: exclusive prefix sum accumulated in the reversed order
template <typename T> 
__device__ __forceinline__ 
T mywarprevexcprefixsum( T reg )
{
    T tmp;
    int lid = threadIdx.x & 0x1f;//lane id
    tmp = __shfl_down_sync(0xffffffff, reg, 1); if (lid <= 1) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 2); if (lid <= 2) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 4); if (lid <= 4) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 8); if (lid <= 8) reg += tmp;
    tmp = __shfl_down_sync(0xffffffff, reg, 16); if (lid <= 16) reg += tmp;
    return reg;
}

#endif//__warpscan_h__
