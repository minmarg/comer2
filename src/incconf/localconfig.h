/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __localconfig_h__
#define __localconfig_h__

#include "liblib/preproc.h"
#include "liblib/platform.h"

#define CONCATSTRSEP( arg1 )        CONCATSTRA( arg1, DIRSEPSTR )

// #if !defined( COMER_INSTALL_DIR )
// #   error   "Unable to compile: Undefined installation directory (COMER_INSTALL_DIR)."
// #endif

#if !defined( COMERSTATEDIR )
#   error   "Unable to compile: Undefined configuration directory (COMERSTATEDIR)."
#endif

#define VARDIR var
#define PARAM_DIRNAME  COMERSTATEDIR

#define PARAM_FILENAME comer.par
#define PARAM_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( PARAM_FILENAME )
//
#define HDP1_FILENAME hdp1.par
#define HDP1_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( HDP1_FILENAME )
#define HDPsco_FILENAME hdps
#define HDPsco_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( HDPsco_FILENAME )
//
// #define HDPCtx_FILENAME hdpctx.par
// #define HDPCtx_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( HDPCtx_FILENAME )
// #define HDPcts_FILENAME hdpcts
// #define HDPcts_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( HDPcts_FILENAME )
//
#define CVS2S_FILENAME cvs2s
#define CVS2S_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( CVS2S_FILENAME )
//
#define SSSS_FILENAME ssss
#define SSSS_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( SSSS_FILENAME )
//
// #define HDPSSSS_FILENAME hdpssss
// #define HDPSSSS_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( HDPSSSS_FILENAME )
// #define iHDPSSSS_FILENAME ihdpssss
// #define iHDPSSSS_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( iHDPSSSS_FILENAME )
//
#define OPTIONS_FILENAME options.txt
#define OPTIONS_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( OPTIONS_FILENAME )


static const char*  var_param_DIR = TOSTR( VARDIR );
static const char*  var_param_DIRNAME = TOSTR( PARAM_DIRNAME );

static const char*  var_param_FILENAME = TOSTR( PARAM_FILENAME );
static const char*  var_param_FULLPATHNAME = PARAM_FULLNAME;
//
static const char*  var_hdp1_FILENAME = TOSTR( HDP1_FILENAME );
static const char*  var_hdp1_FULLPATHNAME = HDP1_FULLNAME;
static const char*  var_hdpsco_FILENAME = TOSTR( HDPsco_FILENAME );
static const char*  var_hdpsco_FULLPATHNAME = HDPsco_FULLNAME;
//
// static const char*  var_hdpctx_FILENAME = TOSTR( HDPCtx_FILENAME );
// static const char*  var_hdpctx_FULLPATHNAME = HDPCtx_FULLNAME;
// static const char*  var_hdpcts_FILENAME = TOSTR( HDPcts_FILENAME );
// static const char*  var_hdpcts_FULLPATHNAME = HDPcts_FULLNAME;
//
static const char*  var_cvs2s_FILENAME = TOSTR( CVS2S_FILENAME );
static const char*  var_cvs2s_FULLPATHNAME = CVS2S_FULLNAME;
//
static const char*  var_ssss_FILENAME = TOSTR( SSSS_FILENAME );
static const char*  var_ssss_FULLPATHNAME = SSSS_FULLNAME;
//
// static const char*  var_hdpssss_FILENAME = TOSTR( HDPSSSS_FILENAME );
// static const char*  var_hdpssss_FULLPATHNAME = HDPSSSS_FULLNAME;
// static const char*  var_ihdpssss_FILENAME = TOSTR( iHDPSSSS_FILENAME );
// static const char*  var_ihdpssss_FULLPATHNAME = iHDPSSSS_FULLNAME;
//
static const char*  var_options_FILENAME = TOSTR( OPTIONS_FILENAME );
static const char*  var_options_FULLPATHNAME = OPTIONS_FULLNAME;


inline const char* GetParamDirectory()      {   return var_param_DIR;  }
inline const char* GetFullParamDirname()    {   return var_param_DIRNAME;  }

inline const char* GetParamFilename()       {   return var_param_FILENAME;  }
inline const char* GetFullParamFilename()   {   return var_param_FULLPATHNAME;  }
//
inline const char* GetHDP1Filename()        {   return var_hdp1_FILENAME;  }
inline const char* GetFullHDP1Filename()    {   return var_hdp1_FULLPATHNAME;  }
inline const char* GetHDPscoFilename()      {   return var_hdpsco_FILENAME;  }
inline const char* GetFullHDPscoFilename()  {   return var_hdpsco_FULLPATHNAME;  }
//
// inline const bool  GetHDPctsUsed()          {   return false; }
// inline const char* GetHDPCtxFilename()      {   return var_hdpctx_FILENAME;  }
// inline const char* GetFullHDPCtxFilename()  {   return var_hdpctx_FULLPATHNAME;  }
// inline const char* GetHDPctsFilename()      {   return var_hdpcts_FILENAME;  }
// inline const char* GetFullHDPctsFilename()  {   return var_hdpcts_FULLPATHNAME;  }
//
inline const char* GetCVS2SFilename()       {   return var_cvs2s_FILENAME;  }
inline const char* GetFullCVS2SFilename()   {   return var_cvs2s_FULLPATHNAME;  }
//
inline const char* GetSSSSFilename()        {   return var_ssss_FILENAME;  }
inline const char* GetFullSSSSFilename()    {   return var_ssss_FULLPATHNAME;  }
//
// inline const char* GetHDPSSSSFilename()     {   return var_hdpssss_FILENAME;  }
// inline const char* GetFullHDPSSSSFilename() {   return var_hdpssss_FULLPATHNAME;  }
// inline const char* GetiHDPSSSSFilename()    {   return var_ihdpssss_FILENAME;  }
// inline const char* GetFulliHDPSSSSFilename(){   return var_ihdpssss_FULLPATHNAME;  }
//
inline const char* GetOptionsFilename()     {   return var_options_FILENAME;  }
inline const char* GetFullOptionsFilename() {   return var_options_FULLPATHNAME;  }

#endif//__localconfig_h__
