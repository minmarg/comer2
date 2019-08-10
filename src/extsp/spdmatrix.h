/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_spdmatrix__
#define __extspsl_spdmatrix__

#include "nsmatrix.h"

#ifndef SPDmatrixTESTPRINT
// #define SPDmatrixTESTPRINT
#endif

namespace extspsl {

// -------------------------------------------------------------------------
// Symmetric, Positive Definite square matrix
//
class SPDmatrix: public NSmatrix
{
public:
    SPDmatrix( int nrc );
    SPDmatrix( const SPDmatrix& );
    virtual ~SPDmatrix();

    virtual SPDmatrix&  operator=( const SPDmatrix& );
    virtual SPDmatrix&  operator=( const NSmatrix& );
    virtual SPDmatrix&  operator=( const Pslmatrix& );

    int     CholeskyDecompose();
    int     CDedSolve( Pslvector& xb ) const;
    int     CDedInvert();
    int     CDedDet( float* ) const;
    int     CDedLogDet( float* ) const;

protected:
    explicit    SPDmatrix();
};


// -------------------------------------------------------------------------

}//namespace extspsl

#endif//__extspsl_spdmatrix__
