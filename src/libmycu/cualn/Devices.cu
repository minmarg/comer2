/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime_api.h>

#include "extsp/psl.h"
#include "liblib/mysort.h"
#include "libpro/srcpro/CLOptions.h"
#include "libmycu/cucom/cucommon.h"

#include "Devices.h"

Devices DEVPROPs;

// _________________________________________________________________________
// Class DeviceProperties
//
void Devices::PrintDevices( FILE* fp )
{
    MYMSG( "Devices::PrintDevices", 6 );

    int deviceCount = 0;
    MYCUDACHECK( cudaGetDeviceCount( &deviceCount ));

    if( deviceCount < 1 ) {
        fprintf( fp, "%sThere are no available device(s) that support CUDA.%s", NL,NL);
        return;
    }

    fprintf( fp, "%sDetected %d CUDA-capable device(s):%s%s", NL, deviceCount, NL,NL);

    for( int dev = 0; dev < deviceCount; dev++ ) {
        MYCUDACHECK( cudaSetDevice( dev ));
        cudaDeviceProp deviceProp;
        MYCUDACHECK( cudaGetDeviceProperties( &deviceProp, dev ));
        fprintf( fp, "Device id %d: \"%s\"%s%s", 
            dev, deviceProp.name,
            (deviceProp.computeMode == cudaComputeModeProhibited)? " [Prohibited!]": "",
            NL);
    }

    fprintf( fp, "%s",NL);
}

// -------------------------------------------------------------------------
// RegisterDevices: register devices scheduled for use
//
void Devices::RegisterDevices()
{
    MYMSG( "Devices::RegisterDevices", 6 );

    SetMaxMemoryAmount( CLOptions::GetDEV_MEM());

    mystring devN = CLOptions::GetDEV_N();
    size_t pos;
    char* p = NULL;
    long int did, ndevs = -1;

    //parse command line option here to get the list of devices
    if( !devN.empty() && devN[0] == ',') {
        NewDevs();
        devN = devN.substr(1);
        //device ids are listed
        for( pos = devN.find(','); ; pos = devN.find(',')) {
            mystring strval = devN.substr( 0, pos );
            errno = 0;
            did = strtol( strval.c_str(), &p, 10 );
            if( errno || *p || strval.empty())
                throw MYRUNTIME_ERROR("Invalid ID specified by option dev-N.");
            if( did < 0 || 1000000 < did )
                throw MYRUNTIME_ERROR("Invalid ID specified by command-line option dev-N.");
            RegisterDeviceProperties( did, maxmem_, true/*checkduplicates*/);
            if( pos == mystring::npos )
                break;
            devN = devN.substr(pos+1);
        }
    }
    else {
        //certain number of devices is given
        errno = 0;
        ndevs = strtol( devN.c_str(), &p, 10 );
        if( errno || *p )
            throw MYRUNTIME_ERROR("Invalid number of devices given: "
            "Begin option dev-N with ',' if a list of IDs is intended.");
        if( ndevs <= 0 || 1000000 < ndevs )
            throw MYRUNTIME_ERROR("Invalid number of devices given by option dev-N.");
        ReadDevices();
        SortDevices();
    }

    PruneRegisteredDevices( ndevs );

    PrettyPrintUsedDevices();

    //NOTE:reset errno as it may have been previously set by CUDA calls!
    errno = 0;
}

// -------------------------------------------------------------------------
// ReadDevices: read all CUDA-capable devices available on the system
//
void Devices::ReadDevices()
{
    MYMSG( "Devices::ReadDevices", 6 );

    int deviceCount = 0;
    MYCUDACHECK( cudaGetDeviceCount( &deviceCount ));

    if( deviceCount < 1 )
        return;

    NewDevs();

    for( int dev = 0; dev < deviceCount; dev++ ) {
        RegisterDeviceProperties( dev, maxmem_, false/*checkduplicates*/);
    }
}

// -------------------------------------------------------------------------
// SortDevices: sort CUDA-capable devices that have been saved in the list
//

//DevCmp: comparison function for sorting devices
static inline
int DevCmp( const void* vector, size_t n1, size_t n2 )
{
    const SimpleVector* vec = (const SimpleVector*)vector;
    if( vec == NULL )
        return 0;
    const DeviceProperties* dev1 = (const DeviceProperties*)vec->GetValueAt(n1);
    const DeviceProperties* dev2 = (const DeviceProperties*)vec->GetValueAt(n2);
    if( dev1 == NULL || dev2 == NULL )
        return 0;
    //[n2]-[n1], to sort in descending order
    return (int)
        (dev2->ccmajor_ == dev1->ccmajor_)
        ?   ((dev2->ccminor_ == dev1->ccminor_)
             ?  (int)((dev2->totmem_>>20) - (dev1->totmem_>>20))
             :  (dev2->ccminor_ - dev1->ccminor_))
        :   (dev2->ccmajor_ - dev1->ccmajor_)
    ;
}

//DevSwp: swap function for sorting devices
static inline
int DevSwp( void* vector, size_t n1, size_t n2 )
{
    SimpleVector* vec = (SimpleVector*)vector;
    if( vec == NULL )
        return 0;
    const DeviceProperties* dev1 = (const DeviceProperties*)vec->GetValueAt(n1);
    const DeviceProperties* dev2 = (const DeviceProperties*)vec->GetValueAt(n2);
    if( dev1 == NULL || dev2 == NULL )
        return 0;
    vec->SetValueAt(n1, dev2);
    vec->SetValueAt(n2, dev1);
    return 0;
}

void Devices::SortDevices()
{
    MYMSG( "Devices::SortDevices", 6 );
    if( devs_ == NULL )
        return;
    HeapSort( devs_, devs_->GetSize(), DevCmp, DevSwp );
}

// -------------------------------------------------------------------------
// RegisterDeviceProperties: register the properties of device with id devid
//

//DevIdCmptor: comparison of registered device ids
static inline
int DevIdCmptor( const void* key1, const void* key2 )
{
    const DeviceProperties* dev1 = (const DeviceProperties*)key1;
    const DeviceProperties* dev2 = (const DeviceProperties*)key2;
    if( dev1 == NULL || dev2 == NULL )
        return -1;
    return dev1->devid_ - dev2->devid_;
}

bool Devices::RegisterDeviceProperties( int devid, ssize_t maxmem, bool checkduplicates )
{
    MYMSG( "Devices::RegisterDeviceProperties", 6 );

#ifdef __DEBUG__
    if( devs_ == NULL )
        throw MYRUNTIME_ERROR("Devices::GetDeviceProperties: Null devices.");
#endif

    if( checkduplicates ) {
        //if checking for duplicate ids is on, 
        // ignore devices that have been already registered (same id)
        DeviceProperties dp;
        dp.devid_ = devid;
        if( devs_->SimpleFind( &dp, DevIdCmptor ))
            return false;
    }

    cudaDeviceProp deviceProp;
    MYCUDACHECK( cudaSetDevice( devid ));
    MYCUDACHECK( cudaGetDeviceProperties( &deviceProp, devid ));
    //
    if( deviceProp.computeMode == cudaComputeModeProhibited )
        return false;
    //
    DeviceProperties* dprop = new DeviceProperties();
    if( dprop == NULL )
        throw MYRUNTIME_ERROR("Devices::GetDeviceProperties: Not enough memory.");
    dprop->devid_ = devid;
    dprop->ccmajor_ = deviceProp.major;
    dprop->ccminor_ = deviceProp.minor;
    dprop->totmem_ = deviceProp.totalGlobalMem;
    dprop->reqmem_ = dprop->totmem_;
    if( 0 < maxmem && maxmem < (ssize_t)dprop->totmem_ )
        dprop->reqmem_ = maxmem;
    if( dprop->reqmem_ + DeviceProperties::DEVMINMEMORYRESERVE > dprop->totmem_ )
        dprop->reqmem_ = (dprop->totmem_ <= DeviceProperties::DEVMINMEMORYRESERVE)?
            dprop->totmem_: dprop->totmem_ - DeviceProperties::DEVMINMEMORYRESERVE;
    dprop->textureAlignment_ = deviceProp.textureAlignment;
    dprop->maxTexture1DLinear_ = deviceProp.maxTexture1DLinear;
    dprop->deviceOverlap_ = deviceProp.deviceOverlap;
    dprop->asyncEngineCount_ = deviceProp.asyncEngineCount;
    dprop->computeMode_ = deviceProp.computeMode;
    //
    MYCUDACHECK( cudaDeviceReset());
    //
    devs_->Push( dprop );
    return true;
}

// -------------------------------------------------------------------------
// PruneRegisteredDevices: prune registered devices so that their number 
// does not exceed the allowed number
//
void Devices::PruneRegisteredDevices( int ndevs )
{
    MYMSG( "Devices::PruneRegisteredDevices", 6 );

    if( devs_ == NULL )
        return;

    if( ndevs < 0 )
        ndevs = maxdevs_;
    
    for( int i = ndevs; i < (int)devs_->GetSize(); ) {
        const DeviceProperties* dprop = (const DeviceProperties*)devs_->RemoveLast();
        if( dprop )
            delete dprop;
    }
}

// -------------------------------------------------------------------------
// PrintSavedDevicesTest: print devices that are currently in the list
//
void Devices::PrettyPrintUsedDevices()
{
    if( devs_ == NULL )
        return;
    MYMSGBEGl(1)
        char msgbuf[KBYTE];
        if( devs_->GetSize())
            MYMSG("Devices scheduled for use:",1);
        for( int i = 0; i < (int)devs_->GetSize(); i++ ) {
            const DeviceProperties* dev = (const DeviceProperties*)devs_->GetValueAt(i);
            if( dev == NULL )
                continue;
            sprintf( msgbuf, 
                "Device id %d: Compute_capability %d.%d Tot_Mem %zu(MB) Used_Mem %zu(MB) "
                " TextAln %zu Text1D %d(MB) Ovlp %d Engs %d Mode %d", 
                dev->devid_, dev->ccmajor_, dev->ccminor_, dev->totmem_>>20, dev->reqmem_>>20,
                dev->textureAlignment_, dev->maxTexture1DLinear_>>20, dev->deviceOverlap_,
                dev->asyncEngineCount_, dev->computeMode_ );
            MYMSG(msgbuf,1);
        }
    MYMSGENDl
}
