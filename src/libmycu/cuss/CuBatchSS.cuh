/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchSS_h__
#define __CuBatchSS_h__

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/MOptions.h"
#include "libmycu/cusco/CuBatchScoreMatrix_com.h"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "CuBatchSS_com.h"

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchSS
// Batch computation of statistical significance for multiple 
// profile-profile pairs using the CUDA architecture
//
class CuBatchSS
{
public:
    CuBatchSS();
    virtual ~CuBatchSS();

    void CalculateAlnStatisticsDevice(
        cudaStream_t streamproc,
        bool Xuninf,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        //
        int ssemodel,
        float qyeno,
        float searchspace,
        float reflambda, float refK,
        float expgappedlambda, float expgappedK,
        //
        CUBSM_TYPE* scores,//in:substitution scores
        CUBDP_TYPE* tmpdpbuffer,//in:temporary information of passed profiles
        float* tmpss2datbuffers,//[out:profile-specific]
        float* dp2alndatbuffers,//[out:profile-specific]
        //summary of profiles passed to phase 2
        unsigned int* attrpassed
    );

    void CalculateStatisticsDeviceHelper(
        cudaStream_t streamproc,
        bool Xuninf,
        size_t ndb1pros,
        size_t ndbCpros,
        size_t ndb1prosOmtd,
        size_t ndbCprosOmtd,
        size_t nqyposs,
        size_t ndb1poss,
        size_t ndbCposs,
        size_t dbxpad,
        size_t querposoffset,
        size_t bdb1posoffset,
        size_t bdbCposoffset,
        //
        int ssemodel,
        float qyeno,
        float searchspace,
        float reflambda, float refK,
        float expgappedlambda, float expgappedK,
        //
        CUBSM_TYPE* scores,//in:substitution scores
        CUBDP_TYPE* tmpdpbuffer,//in:temporary information of passed profiles
        float* tmpss2datbuffers,//[out:profile-specific]
        float* dp2alndatbuffers,//[out:profile-specific]
        //summary of profiles passed to phase 2
        unsigned int* attrpassed
    );

protected:
//     explicit CuBatchSS();

private:
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__CuBatchSS_h__
