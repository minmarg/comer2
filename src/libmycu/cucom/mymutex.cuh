/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mymutex_h__
#define __mymutex_h__

__device__ __forceinline__ 
void LOCK(volatile int* mutex)
{
    while(atomicCAS((int*)mutex, 0, 1) != 0);
}

__device__ __forceinline__ 
void UNLOCK(volatile int* mutex)
{
    *mutex = 0;
    __threadfence();
}

// lock class;
// NOTE: uses an additional register
class UNIQUE_LOCK
{
public:
    __device__ __forceinline__ 
    UNIQUE_LOCK(volatile int* mutex): mutex_(mutex) {lock();}

    __device__ __forceinline__ 
    ~UNIQUE_LOCK() {unlock();}

    __device__ __forceinline__ void lock() {LOCK(mutex_);}
    __device__ __forceinline__ void unlock() {UNLOCK(mutex_);}
private:
    volatile int* mutex_;
};

#endif//__mymutex_h__
