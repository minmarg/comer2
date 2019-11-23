/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PM2DVectorFields_h__
#define __PM2DVectorFields_h__

#include "liblib/alpha.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/TRANSPROBS.h"

// for the case when FPTYPE==INTYPE
// #define FPTYPEisINTYPE
// #define LNTYPEisINTYPE
// #define CHTYPEisINTYPE

//NOTE: for alignment and proper data type conversion in a device by CUDA
// make all types of the same size

#define FPTYPE      float
#define SZFPTYPE    (sizeof(FPTYPE))

#define LNTYPE      unsigned int
#define SZLNTYPE    (sizeof(LNTYPE))

#define INTYPE      int //short
#define SZINTYPE    (sizeof(INTYPE))

#ifdef CHTYPEisINTYPE
#   define CHTYPE   int
#else
#   define CHTYPE   char
#endif
#define SZCHTYPE    (sizeof(CHTYPE))

//constants
enum TPSVectorConsts { 
        pmv2DNoElems = NUMAA,//effective number of frequency values
        pmv2DNoCVEls = NUMAA-1,//number of context vector entries
        pmv2DNoSSSps = SS_NSTATES,//number of SS state probabilities
        ptr2DNoElems = P_NSTATES-2//number of transition states; two transitions are omitted [see CuBatchDP_com.h]
};

//NOTE: data alignment is NOT required as a separate buffer is 
// allocated for each field

enum TPM2DVectorFields {
//fields and their indices of the 2D vector, representing general data for a 
//given profile model:
        pps2DBkgPrbs = 0,//profile background amino acid probabilities; pmv2DNoElems fields in total
        pps2DENO = pps2DBkgPrbs + pmv2DNoElems,//effective number of observations
        pps2DLen = pps2DENO + 1,//profile length; [SHORT INT]
        pps2DDist = pps2DLen + 1,//the total number of positions up to this profile; [UINT]

//fields and their indices of 2D vector representing transitions of 
//one profile model position;
//note: each profile has a vector of beginning transitions
//         ptr2DTrnPrbsExp = pps2DDist + 1,//log transition probabilities exponentiated; ptr2DNoElems fields
//         ptr2DTrnPrbs = ptr2DTrnPrbsExp + ptr2DNoElems,//log transition probabilities; ptr2DNoElems fields in total
        ptr2DTrnPrbs = pps2DDist + 1,//log transition probabilities; ptr2DNoElems fields in total

//{{Fields and their indices of 2D vector representing one position of 
//profile model for scoring;
//NOTE: profile (and CV and SS vector) address of type unsigned short is safe enough, as the 
// max number of profiles processed in parallel is unlikely currently (on both CPU and GPU)
        pmv2DTrgFrqs = ptr2DTrnPrbs + ptr2DNoElems,//target frequencies; pmv2DNoElems in total
        pmv2DCVentrs = pmv2DTrgFrqs + pmv2DNoElems,//context vector; pmv2DNoCVEls entries; [*OPTIONAL*]
        pmv2DCVprior = pmv2DCVentrs + pmv2DNoCVEls,//the prior probability of context vector; [*OPTIONAL*]
        pmv2DCVnorm2 = pmv2DCVprior + 1,//the squared norm of context vector; [*OPTIONAL*]
        pmv2DSSsprbs = pmv2DCVnorm2 + 1,//SS state probabilities; pmv2DNoSSSps entries; [*OPTIONAL*]
        pmv2DHDP1prb = pmv2DSSsprbs + pmv2DNoSSSps,//the largest probability for pos vector to belong to a cluster
        pmv2DHDP1ind = pmv2DHDP1prb + 1,//the index of the most probable cluster; [SHORT INT]
        pmv2DAddrPro = pmv2DHDP1ind + 1,//the address of profile in the vector of profile-specific data; [SHORT INT]
        pmv2DAddrCV = pmv2DAddrPro + 1,//the address of CV in vectors of position-specific profile model; [INT]
        pmv2DAddrSS = pmv2DAddrCV + 1,//the address of SS data in vectors of position-specific profile model; [INT]
        pmv2Daa = pmv2DAddrSS + 1,//amino acid; [CHAR]
        pmv2DSSstate = pmv2Daa + 1,//SS state indicator; [CHAR]; [*OPTIONAL*]
//}}
        pmv2DTotFlds//total number of fields
};

#endif//__PM2DVectorFields_h__
