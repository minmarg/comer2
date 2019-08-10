/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __logitnormal_h__
#define __logitnormal_h__

// interface
//
void LogitNormalErrCorr( double* f, size_t fsz, double tol );
void LogitNormal2Normal( double* f, size_t fsz, double tol, bool correct = true );
void Normal2LogitNormal( const double* nv, size_t nvsz, double* lnv, size_t lnvsz );

// equivalents for data of the float type
void LogitNormalErrCorr_f( float* f, size_t fsz, float tol );
void LogitNormal2Normal_f( float* f, size_t fsz, float tol, bool correct = true );
void Normal2LogitNormal_f( const float* nv, size_t nvsz, float* lnv, size_t lnvsz );

#endif//__logitnormal_h__
