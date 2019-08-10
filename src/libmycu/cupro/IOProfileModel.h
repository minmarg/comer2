/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __IOProfileModel_h__
#define __IOProfileModel_h__

#include "liblib/mybase.h"
#include "liblib/myfiler.h"

extern "C" {
extern const char* patstrDATBINVER;
extern const char* patstrDATVER;
extern const char* patstrDESC;
extern const char* patstrFILE;
extern const char* patstrCMD;
extern const char* patstrLEN;
extern const char* patstrNOS;
extern const char* patstrEFFNOS;
extern const char* patstrSCALE;
extern const char* patstrNULL;
extern const char* patstrPOST;
extern const char* patstrCT;
extern const char* patstrCV;
extern const char* patstrSS;
extern const char* patstrSPCOMP;
extern const char* patstrSPREFR;
extern const char* patstrENTROPY;
extern const char* patstrEXPECTED;
extern const char* patstrEXPNN;
extern const char* patstrEND;
extern const size_t lenstrEND;
extern const int   INTSCALE;
}

// -------------------------------------------------------------------------

void BinaryReadProfileHeader(TCharStream* fp, 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata);

void BinaryReadProfileData(TCharStream* fp, 
    unsigned int distn,
    int profnr,
    int prolen,
    int /*scale*/, 
    char** pmorgaddr, 
    char** pmdata );

// -------------------------------------------------------------------------

template<typename T>
void TextReadProfileHeader(T* fp, 
                    mystring* desc, mystring* file,
                    int* prolen, int* scale, char** pmdata);
template<typename T>
void TextReadProfileData(T* fp,
                    unsigned int distn,
                    int profnr, int prolen, int scale,
                    char** pmorgaddr, char** pmdata );

void TextWriteProfileHD(FILE* fp, 
                    const mystring& desc, const mystring& file, 
                    int scale, 
                    char** pmorgaddr, 
                    char** pmdata );

template<typename T>
void TextReadProfileHeaderBufferless(T* fp, 
                    mystring* desc, mystring* file,
                    int* prolen, int* scale, char** pmdata);
template<typename T>
void TextReadProfileDataBufferless(T* fp,
                    unsigned int distn,
                    int profnr, int prolen, int scale,
                    char** pmorgaddr, char** pmdata );

// -------------------------------------------------------------------------

#endif//__IOProfileModel_h__
