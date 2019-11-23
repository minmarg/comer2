/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchMAPDP_h__
#define __CuBatchMAPDP_h__

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
#include "CuBatchMAPDP_com.h"

////////////////////////////////////////////////////////////////////////////
// CLASS CuBatchMAPDP
// Batch computation of MAP dynamic programming for multiple profile-profile 
// pairs using the CUDA architecture
//
class CuBatchMAPDP
{
public:
    CuBatchMAPDP();
    virtual ~CuBatchMAPDP();

    void PerformCompleteMAPDynProgDevice(
        cudaStream_t streamproc,
        float logevthld,
        bool modscoresinuse,
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
        //
        CUBSM_TYPE* scores,
        CUBSM_TYPE* modscores,
        float* tmpdpdiagbuffers, 
        float* tmpdpbotbuffer,
//float* tmpblockprobscales, 
//float* tmpprobscales,
        unsigned int* maxcoordsbuf,
        char* btckdata,
        float* tmpfwdprobs,
//float* tmpbwdprobs,
        float* tmpss2datbuffers,//[tmp] for total probability
        float* dp2alndatbuffers,//[in]
        char* outalns,//[out]
        //summary of profiles passed to phase 2
        unsigned int* attrpassed,
        size_t dbxpad2,
        size_t cbdbalnlen2
    );

    void CalculateFwdBwdProbsMAPDynProgDevice(
        cudaStream_t streamproc,
        float logevthld,
        bool modscoresinuse,
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
        CUBSM_TYPE* modscores,
        float* tmpdpdiagbuffers,
        float* tmpdpbotbuffer,
//float* tmpblockprobscales,
//float* tmpprobscales,
        float* tmpfwdprobs,
//float* tmpbwdprobs,
        float* tmpss2datbuffers,//[tmp]
        float* dp2alndatbuffers,//[in]
        unsigned int* attrpassed,
        size_t dbxpad2,
        unsigned int* maxcoordsbuf,
        char* btckdata
    );

    void FinalizeMAPDynProgDevice(
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
        float* tmpdpdiagbuffers,//[in]
        float* dp2alndatbuffers,//[in/out]
        unsigned int* attrpassed,
        size_t dbxpad2,
        size_t cbdbalnlen2,
        unsigned int* maxcoordsbuf,//[in]
        char* btckdata,//[in]
        char* outalns//[out]
    );

protected:
//     explicit CuBatchMAPDP();

private:
};

// -------------------------------------------------------------------------
// INLINES ...
//

#endif//__CuBatchMAPDP_h__
