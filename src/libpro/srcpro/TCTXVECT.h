/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __TCTXVECT_h__
#define __TCTXVECT_h__

#include "extsp/pslvector.h"
#include "liblib/CtxtCoefficients.h"

// class CtxtCoefficients;
// class Pslvector;

//parameters of context vector distributed as mv. t-distribution
struct _TCTXVECT
{
    _TCTXVECT();
    //
    const int   AVGLEN, hAVGLEN;//length for averaging
    const float AVGCWGT;//central weight in averaging
    const CtxtCoefficients AVGCOEFFS;//coefficients in averaging
    //
    const int   CTXLEN;//context length
    const float CWGT;//central weight of context coefficients
    const bool  MIX;//mix normal vectors within the context
    const bool  AUG;//augment vector 3-fold: into left, center, right partitions
    const bool  MEAN;//mean vector in use
    const int   DIM;//dimensions of context vector
    const float KAPPA0;//strictly positive
    const float NU0;//must be greater than dim-1
    const extspsl::Pslvector MU0;//mean vector
    const float loKAPPA0;
    const float PowerNU0;
    const float CTERM;//constant term in calculating log probability of observations
};
extern _TCTXVECT CVS;

#endif//__TCTXVECT_h__
