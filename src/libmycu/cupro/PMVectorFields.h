/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMVectorFields_h__
#define __PMVectorFields_h__

#include "liblib/alpha.h"
#include "liblib/mybase.h"
#include "libpro/srcpro/TRANSPROBS.h"

// for the case when FPTYPE==INTYPE
// #define FPTYPEisINTYPE

#define FPTYPE      float
#define SZFPTYPE    (sizeof(FPTYPE))

#define INTYPE      short
#define SZINTYPE    (sizeof(INTYPE))

#define CHTYPE      char
#define SZCHTYPE    (sizeof(CHTYPE))

//fields and their addresses of vector representing general data for a 
//given profile model
enum TPSVector { 
        ppsNoElems = NUMAA,//effective number of frequency values
        //
        //SORT the fields by the size of their elements and ALIGN them
        //(sorted fields help avoid padding in between to align them)
        ppsBkgPrbs = 0,//profile background amino acid probabilities
        ppsENO = ppsBkgPrbs + ppsNoElems*SZFPTYPE,//effective number of observations
        ppsLen = ppsENO + SZFPTYPE,//profile length
        //
        ppsTotal = ppsLen + SZINTYPE,//total size of all fields
        ppsEnd = ALIGN_UP( ppsTotal, SZFPTYPE )//end of the data vector icludes padding wrt SZFPTYPE
};

//{{Fields and their addresses of vector representing one position of 
//profile model for scoring.
//Codes describing information present in the vector:
enum TPosVecInfCodes { 
        pmvicFull = 0,//00b, full information
        pmvicNoSS,//01b, SS data is absent
        pmvicNoCV,//10b, CV data is absent
        pmvicNoSSCV,//11b, SS and CV data is absent
        pmvicN
};
enum TPosVectorPreamb { 
        pmvNoElems = NUMAA,//effective number of frequency values
        pmvNoCVEls = NUMAA-1,//number of context vector entries
        pmvNoSSSps = SS_NSTATES,//number of SS state probabilities
};
//addresses of fields: the size of each field is represented by the added term in the next field;
//SORT the fields by the size of their elements and ALIGN them
//(sorted fields help avoid padding in between to align them)
//NOTE: profile address of type unsigned short is safe enough, as the max number of 
// profiles processed in parallel is unlikely currently (on both CPU and GPU)
enum TPosVectorFull {
        pmvfullInfCode = 0,//indicator of information present in the vector; type char
        pmvfullTrgFrqs = ALIGN_UP( pmvfullInfCode+SZCHTYPE, SZFPTYPE ),//target frequencies
        pmvfullCVentrs = pmvfullTrgFrqs + pmvNoElems*SZFPTYPE,//context vector
        pmvfullCVprior = pmvfullCVentrs + pmvNoCVEls*SZFPTYPE,//the prior probability of context vector
        pmvfullCVnorm2 = pmvfullCVprior + SZFPTYPE,//the squared norm of context vector
        pmvfullSSsprbs = pmvfullCVnorm2 + SZFPTYPE,//SS state probabilities
        pmvfullHDP1prb = pmvfullSSsprbs + pmvNoSSSps*SZFPTYPE,//the largest probability for pos vector to belong to a cluster
        pmvfullHDP1ind = pmvfullHDP1prb + SZFPTYPE,//the index of the most probable cluster
        pmvfullAddrPro = pmvfullHDP1ind + SZINTYPE,//the address of profile in the vector of profile-specific data
        pmvfullSSstate = pmvfullAddrPro + SZINTYPE,//SS state indicator
        //
        pmvfullTotal = pmvfullSSstate + SZCHTYPE,//total size of all fields
        //end of the data vector icludes padding: assume SZFPTYPE to be of the greatest size
        pmvfullEnd = ALIGN_UP( pmvfullTotal, SZFPTYPE )
};
enum TPosVectorNoSS {
        pmvnossInfCode = 0,//indicator of information present in the vector; type char
        pmvnossTrgFrqs = ALIGN_UP( pmvnossInfCode+SZCHTYPE, SZFPTYPE ),//target frequencies
        pmvnossCVentrs = pmvnossTrgFrqs + pmvNoElems*SZFPTYPE,//context vector
        pmvnossCVprior = pmvnossCVentrs + pmvNoCVEls*SZFPTYPE,//the prior probability of context vector
        pmvnossCVnorm2 = pmvnossCVprior + SZFPTYPE,//the squared norm of context vector
        pmvnossSSsprbs = 0,//no SS data
        pmvnossHDP1prb = pmvnossCVnorm2 + SZFPTYPE,//the largest probability for pos vector to belong to a cluster
        pmvnossHDP1ind = pmvnossHDP1prb + SZFPTYPE,//the index of the most probable cluster
        pmvnossAddrPro = pmvnossHDP1ind + SZINTYPE,//the address of profile in the vector of profile-specific data
        pmvnossSSstate = 0,//no SS data
        //
        pmvnossTotal = pmvnossAddrPro + SZINTYPE,//total size of all fields
        pmvnossEnd = ALIGN_UP( pmvnossTotal, SZFPTYPE )//end of the data vector icludes padding
};
enum TPosVectorNoCV {
        pmvnocvInfCode = 0,//indicator of information present in the vector; type char
        pmvnocvTrgFrqs = ALIGN_UP( pmvnocvInfCode+SZCHTYPE, SZFPTYPE ),//target frequencies
        pmvnocvCVentrs = 0,//no context vector data
        pmvnocvCVprior = 0,//no CV data
        pmvnocvCVnorm2 = 0,//no CV data
        pmvnocvSSsprbs = pmvnocvTrgFrqs + pmvNoElems*SZFPTYPE,//SS state probabilities
        pmvnocvHDP1prb = pmvnocvSSsprbs + pmvNoSSSps*SZFPTYPE,//the largest probability for pos vector to belong to a cluster
        pmvnocvHDP1ind = pmvnocvHDP1prb + SZFPTYPE,//the index of the most probable cluster
        pmvnocvAddrPro = pmvnocvHDP1ind + SZINTYPE,//the address of profile in the vector of profile-specific data
        pmvnocvSSstate = pmvnocvAddrPro + SZINTYPE,//SS state indicator
        //
        pmvnocvTotal = pmvnocvSSstate + SZCHTYPE,//total size of all fields
        pmvnocvEnd = ALIGN_UP( pmvnocvTotal, SZFPTYPE )//end of the data vector icludes padding
};
enum TPosVectorNoSSCV {
        pmvnosscvInfCode = 0,//indicator of information present in the vector; type char
        pmvnosscvTrgFrqs = ALIGN_UP( pmvnosscvInfCode+SZCHTYPE, SZFPTYPE ),//target frequencies
        pmvnosscvCVentrs = 0,//no context vector data
        pmvnosscvCVprior = 0,//no CV data
        pmvnosscvCVnorm2 = 0,//no CV data
        pmvnosscvSSsprbs = 0,//no SS data
        pmvnosscvHDP1prb = pmvnosscvTrgFrqs + pmvNoElems*SZFPTYPE,//the largest probability for pos vector to belong to a cluster
        pmvnosscvHDP1ind = pmvnosscvHDP1prb + SZFPTYPE,//the index of the most probable cluster
        pmvnosscvAddrPro = pmvnosscvHDP1ind + SZINTYPE,//the address of profile in the vector of profile-specific data
        pmvnosscvSSstate = 0,//no SS data
        //
        pmvnosscvTotal = pmvnosscvAddrPro + SZINTYPE,//total size of all fields
        pmvnosscvEnd = ALIGN_UP( pmvnosscvTotal, SZFPTYPE )//end of the data vector icludes padding
};
//}}

//fields and their addresses of vector representing transitions of 
//one profile model position;
//note: each profile has a vector of beginning transitions
enum TPosTrVector { 
        ptrNoElems = P_NSTATES,//number of transition states
        //
        ptrTrnPrbs = 0,//transition probabilities
        //
        ptrTotal = ptrTrnPrbs + ptrNoElems*SZFPTYPE,//total size of fields
        ptrEnd = ALIGN_UP( ptrTotal, SZFPTYPE )////end of the data vector icludes padding wrt SZFPTYPE
};

#endif//__PMVectorFields_h__
