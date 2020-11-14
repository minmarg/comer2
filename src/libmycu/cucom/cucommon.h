/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __cucommon_h__
#define __cucommon_h__

#include "liblib/mybase.h"

#include <assert.h>

#include "liblib/msg.h"
#include "myassert.h"
// #include "warpscan.cuh"

// CUDA L2 Cache line size
#define CUL2CLINESIZE 32

// evaluate err anyway; if it does not evaluate to cudaSuccess, reset the 
// error by calling cudaGetLastError() and return false:
#define MYCUDACHECKPASS(err) ((err), cudaGetLastError() == cudaSuccess)
#define MYCUDACHECKLAST mycudaCheck(cudaGetLastError(),__FILE__,__LINE__,__func__)
#define MYCUDACHECK(err) mycudaCheck(err,__FILE__,__LINE__,__func__)
inline __host__ __device__ void mycudaCheck( 
    cudaError_t err,
    const char* file, unsigned int line, const char* func )
{
    if( err != cudaSuccess ) {
        const char* sub = DIRSEPSTR "src" DIRSEPSTR;
        const char* p = file;
        int i = 0;
        for(; p && *p; p += (i? i: 1)) {
            for(i = 0; p[i] && sub[i] && p[i] == sub[i]; i++);
            if(!sub[i] || !p[i])
                break;
        }
#if defined(__CUDA_ARCH__)
        printf( "CUDA error %s \"%s\"%s"
                "    (File: %s; Line: %u; Function: %s)%s",
                cudaGetErrorName(err), cudaGetErrorString(err), NL, 
                sub[i]? file: p+1, line, func, NL );
        //__syncthreads();
        asm("trap;");
#else
        fprintf( stderr, "CUDA error %s \"%s\"%s"
                "    (File: %s; Line: %u; Function: %s)%s",
                cudaGetErrorName(err), cudaGetErrorString(err), NL, 
                sub[i]? file: p+1, line, func, NL );
        throw myruntime_error(
            mystring("CUDA error ")+cudaGetErrorName(err)+" \""+cudaGetErrorString(err)+"\"", 
                sub[i]? file: p+1, line, func );
#endif
    }
}

//simple swap template function
template <typename T> __host__ __device__ __forceinline__ 
void myhdswap( T& a, T& b ){ T c(a); a=b; b=c; }

template <typename T> __host__ __device__ __forceinline__ 
T myhdmax( T a, T b ){ return a<b?b:a; }

template <typename T> __host__ __device__ __forceinline__ 
T myhdmin( T a, T b ){ return a<b?a:b; }

template <typename T, typename T2> __host__ __device__ __forceinline__ 
void myhdmaxassgn( T& a, T b, T2& c, T2 d ){ if(a<b){a=b; c=d;} }

// template<typename T>
// __device__ __forceinline__ T ldg(const T* ptr) {
// #if __CUDA_ARCH__ >= 350
//     return __ldg(ptr);
// #else
//     return *ptr;
// #endif
// }

//__forceinline__ __device__ unsigned int lane_id()
//{
//    unsigned int ret; 
//    asm volatile ("mov.u32 %0, %laneid;" : "=r"(ret));
//    return ret;
//}

#endif//__cucommon_h__
