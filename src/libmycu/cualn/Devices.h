/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __Devices_h__
#define __Devices_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include "liblib/BinarySearchStructure.h"

#define MAXNDEVS 32

class Devices;

extern Devices DEVPROPs;

// -------------------------------------------------------------------------
//
struct DeviceProperties {
    enum {
        DEVMINMEMORYRESERVE = 200 * ONEM
    };
    DeviceProperties():
        devid_(-1),
        ccmajor_(0),ccminor_(0),
        totmem_(0),reqmem_(0),
        textureAlignment_(0),maxTexture1DLinear_(0),
        deviceOverlap_(0),asyncEngineCount_(0),
        computeMode_(0)
    {};
    DeviceProperties( const DeviceProperties& dprop ):
        devid_(dprop.devid_),
        ccmajor_(dprop.ccmajor_),
        ccminor_(dprop.ccminor_),
        totmem_(dprop.totmem_),
        reqmem_(dprop.reqmem_),
        textureAlignment_(dprop.textureAlignment_),
        maxTexture1DLinear_(dprop.maxTexture1DLinear_),
        deviceOverlap_(dprop.deviceOverlap_),
        asyncEngineCount_(dprop.asyncEngineCount_),
        computeMode_(dprop.computeMode_)
    {};
    int devid_;//device id
    int ccmajor_;//major compute capability
    int ccminor_;//minor compute capability
    size_t totmem_;//total global memory
    size_t reqmem_;//requested amount of global memory
    size_t textureAlignment_;//alignment size for textures
    int maxTexture1DLinear_;//1D texture size
    int deviceOverlap_;//concurrent memory transfer and kernel execution
    int asyncEngineCount_;//number of asynchronous memory transfer engines
    int computeMode_;//enum:
    //cudaComputeModeDefault = 0
    //cudaComputeModeExclusive = 1
    //cudaComputeModeProhibited = 2
    //cudaComputeModeExclusiveProcess = 3
};

// _________________________________________________________________________
// Class Devices
//
// Implementation of queries for device properties
//
class Devices
{
public:
    Devices( int maxdvs = MAXNDEVS ): 
        maxdevs_(maxdvs),
        maxmem_(-1L),
        devs_(NULL)
    {};
    ~Devices() {
        DestroyDevs();
    };

    void SetMaxMemoryAmount( ssize_t mem_in_mb ) { maxmem_ = mem_in_mb * ONEM; };

    void PrintDevices( FILE* );
    void PrettyPrintUsedDevices();

    void RegisterDevices();

    int GetNDevices() const;
    const DeviceProperties* GetDevicePropertiesAt(int n) const;

protected:
    void ReadDevices();
    void SortDevices();
    bool RegisterDeviceProperties( int devid, ssize_t maxmem, bool checkduplicates );
    void PruneRegisteredDevices( int ndevs );
    void DestroyDevs();
    void NewDevs();

private:
    int maxdevs_;//maximum number of devices
    ssize_t maxmem_;//maximum memory amount for all devices
    SimpleVector* devs_;//list of devices
};

////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
//
inline
void Devices::DestroyDevs()
{
    MYMSG( "Devices::DestroyDevs", 6 );
    if( devs_ == NULL )
        return;
    for( int i = 0; i < (int)devs_->GetSize(); i++ ) {
        if( devs_->GetValueAt(i))
            delete (DeviceProperties*)devs_->GetValueAt(i);
    }
    delete devs_;
    devs_ = NULL;
}

// -------------------------------------------------------------------------
//
inline
void Devices::NewDevs()
{
    MYMSG( "Devices::NewDevs", 6 );
    DestroyDevs();
    devs_ = new SimpleVector(MAXNDEVS);
    if( devs_ == NULL )
        throw MYRUNTIME_ERROR("Devices::NewDevs: Not enough memory.");
}

// -------------------------------------------------------------------------
//
inline
int Devices::GetNDevices() const
{
    MYMSG( "Devices::GetNDevices", 7 );
    if( devs_ == NULL )
        return 0;
    return (int)devs_->GetSize();
}

// -------------------------------------------------------------------------
//
inline
const DeviceProperties* Devices::GetDevicePropertiesAt(int n) const
{
    MYMSG( "Devices::GetDevicePropertiesAt", 6 );
    if( devs_ == NULL || (int)devs_->GetSize() <= n )
        return NULL;
    return (const DeviceProperties*)devs_->GetValueAt(n);
}

#endif//__Devices_h__
