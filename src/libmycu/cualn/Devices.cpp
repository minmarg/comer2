/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime_api.h>

#include "liblib/mybase.h"
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
// SetUsedDevices: set devices scheduled for use
//
void Devices::SetUsedDevices()
{
    MYMSG( "Devices::SetUsedDevices", 6 );

    SetMaxMemoryAmount( CLOptions::GetDEV_MEM());

    mystring devN = CLOptions::GetDEV_N();
    size_t pos;
    char* p;
    int val;
    long int did, ndevs;

    //parse command line option here to get the list of devices
    if((pos = devN.find(',')) == mystring::npos ) {
        //certain number is given
        ndevs = strtol( devN.c_str(), &p, 10 );
        if( errno || *p )
            throw MYRUNTIME_ERROR("Invalid argument of option dev-N.");
        if( ndevs < 0 || 1000000 < ndevs )
            throw MYRUNTIME_ERROR("Invalid command-line option dev-N.");
        ReadDevices();
        SortDevices();
    }
    else {
        //device ids are listed
        for( ; pos != mystring::npos; pos = devN.find(',')) {
            mystring strval = devN.substr( 0, pos );
            did = strtol( strval.c_str(), &p, 10 );
            if( errno || *p )
                throw MYRUNTIME_ERROR("Invalid id specified by option dev-N.");
            if( did < 0 || 1000000 < did )
                throw MYRUNTIME_ERROR("Invalid id specified by command-line option dev-N.");
            devN = devN.substr(pos+1);
            MYCUDACHECK( cudaSetDevice( did ));
            GetDeviceProperties( did, maxmem_ );
        }
    }


    PrettyPrintUsedDevices();
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
        MYCUDACHECK( cudaSetDevice( dev ));
        GetDeviceProperties( dev, maxmem_ );
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
             ?  ((dev2->totmem_>>20) - (dev1->totmem_>>20))
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
// GetDeviceProperties: get the properties of device with id devid
//
void Devices::GetDeviceProperties( int devid, ssize_t maxmem )
{
    MYMSG( "Devices::GetDeviceProperties", 6 );

#ifdef __DEBUG__
    if( devs_ == NULL )
        throw MYRUNTIME_ERROR("Devices::GetDeviceProperties: Null devices.");
#endif

    cudaDeviceProp deviceProp;
    MYCUDACHECK( cudaGetDeviceProperties( &deviceProp, devid ));
    //
    if( deviceProp.computeMode == cudaComputeModeProhibited )
        return;
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
    dprop->textureAlignment_ = deviceProp.textureAlignment;
    dprop->maxTexture1DLinear_ = deviceProp.maxTexture1DLinear;
    dprop->deviceOverlap_ = deviceProp.deviceOverlap;
    dprop->asyncEngineCount_ = deviceProp.asyncEngineCount;
    dprop->computeMode_ = deviceProp.computeMode;
    //
    devs_->Push( dprop );
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
        MYMSG("Devices scheduled for use:",1);
        for( int i = 0; i < (int)devs_->GetSize(); i++ ) {
            const DeviceProperties* dev = (const DeviceProperties*)devs_->GetValueAt(i);
            if( dev == NULL )
                continue;
            sprintf( msgbuf, 
                "Device id %d: Compute_capability %d.%d Tot_Mem %lu(MB) Used_Mem %ld(MB) "
                " TextAln %lu Text1D %d(MB) Ovlp %d Engs %d Mode %d", 
                dev->devid_, dev->ccmajor_, dev->ccminor_, dev->totmem_>>20, dev->reqmem_>>20,
                dev->textureAlignment_, dev->maxTexture1DLinear_>>20, dev->deviceOverlap_,
                dev->asyncEngineCount_, dev->computeMode_ );
            MYMSG(msgbuf,1);
        }
    MYMSGENDl
}
