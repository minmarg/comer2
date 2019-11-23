/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "liblib/mybase.h"
// #include "liblib/myfiler.h"
#include "liblib/ConfigFile.h"
#include "MOptions.h"

using namespace MOptions;

const char* MOptions::SECOPTIONS = "OPTIONS";

// global file-private filename
static const char* gsMOptions_filename;

// private function definitions
static inline const char* GetFilename() { return gsMOptions_filename; }

// private functional interface
static void PrivRead();

// =========================================================================
// MACROs
//
#define MOPTPRIV_int_KEY( KEY1 ) \
    int         value; \
    const char* section = SECOPTIONS; \
    ConfigFile  config( section, GetFilename()); \
    if( !config.GetInt( KEY1, &value )) \
        return;//not found, return

#define MOPTPRIV_float_KEY( KEY1 ) \
    float       value; \
    const char* section = SECOPTIONS; \
    ConfigFile  config( section, GetFilename()); \
    if( !config.GetFloat( KEY1, &value )) \
        return;//not found, return

#define MOPTPRIV_double_KEY( KEY1 ) \
    double      value; \
    const char* section = SECOPTIONS; \
    ConfigFile  config( section, GetFilename()); \
    if( !config.GetDouble( KEY1, &value )) \
        return;//not found, return

#define MOPTPRIV_mystring_KEY( KEY1 ) \
    char        cvalue[BUF_MAX]; \
    mystring    value; \
    const char* section = SECOPTIONS; \
    ConfigFile  config( section, GetFilename()); \
    if( !config.GetString( KEY1, cvalue, BUF_MAX )) \
        return;/*not found, return*/ \
    value = cvalue;

// -------------------------------------------------------------------------

#define DEFINEOPTION( NAME, TYPE, DEFVALUE, COND ) \
  TYPE MOptions::val##NAME##_ = DEFVALUE; \
  static void Read##NAME() { \
    MOPTPRIV_##TYPE##_KEY( TOSTR(NAME) ); \
    if(!( COND )) \
        throw MYRUNTIME_ERROR(mystring(GetFilename())+": Invalid option " TOSTR(NAME)); \
    val##NAME##_ = value; \
  }

// -------------------------------------------------------------------------
// Output options
DEFINEOPTION( EVAL, float, 10.0f, value>=0.0f )
DEFINEOPTION( NOHITS, int, 700, value>=1 )
DEFINEOPTION( NOALNS, int, 700, value>=1 )
DEFINEOPTION( SHOW, int, 1, value==0 || value==1 )
DEFINEOPTION( DSCLEN, int, 1024, value>=0 );
DEFINEOPTION( DSCWIDTH, int, 60, value>=40 && value<=100000);
DEFINEOPTION( ALNWIDTH, int, 60, value>=20 && value<=100000);

// Profile construction options
DEFINEOPTION( IDENTITY, int, 90, value>=1 && value<=100 )
DEFINEOPTION( SHOWCMD, int, 0, value==0 || value==1 );

DEFINEOPTION( DELSTATE, int, 1, 1 )

// ADVANCED OPTIONS:
// Profile construction options
DEFINEOPTION( PCFWEIGHT, float, 1.5f, value>=0.0f )
DEFINEOPTION( MINALNFRN, int, 5, value>=1 && value<=100 )
DEFINEOPTION( MINALNPOS, int, 10, value>=1 )

DEFINEOPTION( SUBMAT, mystring, "Gonnet", 1 )

DEFINEOPTION( TFRMIX, mystring, "no", value=="no" || value=="hdpctx" )
DEFINEOPTION( MIXWGT, float, 0.1f, value>0.0f && value<=1.0f )
DEFINEOPTION( SCOADJ, mystring, "hdpsco", value=="no" || value=="hdpctx" || value=="hdpsco" )
DEFINEOPTION( SUPCLT, int, 5, value>=1 || value==-1 )
DEFINEOPTION( ADJWGT, float, 0.33f, value>0.0f && value<1.0f )
DEFINEOPTION( cADJWGT, float, 0.33f, value>0.0f && value<1.0f )

DEFINEOPTION( CVSWGT, float, 0.15f, value>=0.0f && value<1.0f )

DEFINEOPTION( SSSWGT, float, 0.12f, value>=0.0f && value<1.0f )

// model for the estimation of statistical significance
DEFINEOPTION( SSEMODEL, int, 1, value>=0 && value<=2 )

// SEG options; HC Filtering
DEFINEOPTION( HCFILTER, int, 0, 1 )
DEFINEOPTION( HCWINDOW, int, 12, value>1 )
DEFINEOPTION( HCLOWENT, float, 3.3f, value>=0.0f )
DEFINEOPTION( HCHIGHENT, float, 3.4f, value>=0.0f )

// LC Filtering
DEFINEOPTION( INVLCFILTER, int, 0, 1 )
DEFINEOPTION( LCFILTEREACH, int, 1, 1 )
DEFINEOPTION( LCWINDOW, int, 8, value>1 )
DEFINEOPTION( LCLOWENT, float, 1.6f, value>=0.0f )
DEFINEOPTION( LCHIGHENT, float, 1.6f, value>=0.0f )
DEFINEOPTION( DISTANCE, float, 12.96f, value>=0.0f )

DEFINEOPTION( SCHEME, mystring, "psLSO", 1 )
DEFINEOPTION( MAPALN, int, 1, value==0 || value==1 );
DEFINEOPTION( MINPP, float, 0.3f, value>=0.0f && value<1.0f )

// High-scoring segment pairs
DEFINEOPTION( HSPLEN, int, 3, value>2 && value<=50 )
DEFINEOPTION( HSPMINSCORE, int, 0, value>=0 )
DEFINEOPTION( HSPMAXDIST, int, 60, value>0 )
DEFINEOPTION( NOHSPS, int, 3, value>0 && value<=10 )

// =========================================================================
// Read: read options from file
//
void MOptions::Read( const char* filename )
{
    gsMOptions_filename = filename;
//     if( !filename || !file_exists( filename ))
//         throw MYRUNTIME_ERROR( mystring("MOptions::Read: Options file not found: ") + filename );
    PrivRead();
}

// -------------------------------------------------------------------------
// Read: read options from file
//
static void PrivRead()
{
    ReadEVAL();
    ReadNOHITS();
    ReadNOALNS();
    ReadSHOW();
    ReadDSCLEN();
    ReadDSCWIDTH();
    ReadALNWIDTH();

    ReadIDENTITY();
    ReadSHOWCMD();

    ReadDELSTATE();

    ReadPCFWEIGHT();
    ReadMINALNFRN();
    ReadMINALNPOS();
    ReadSUBMAT();

    ReadTFRMIX();
    ReadMIXWGT();
    ReadSCOADJ();
    ReadSUPCLT();
    ReadADJWGT();
    ReadcADJWGT();

    ReadCVSWGT();

    ReadSSSWGT();

    ReadSSEMODEL();

    ReadHCFILTER();
    ReadHCWINDOW();
    ReadHCLOWENT();
    ReadHCHIGHENT();

    ReadINVLCFILTER();
    ReadLCFILTEREACH();
    ReadLCWINDOW();
    ReadLCLOWENT();
    ReadLCHIGHENT();
    ReadDISTANCE();

    ReadSCHEME();
    ReadMAPALN();
    ReadMINPP();

    ReadHSPLEN();
    ReadHSPMINSCORE();
    ReadHSPMAXDIST();
    ReadNOHSPS();
}
