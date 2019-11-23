/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <ctype.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/alpha.h"
#include "liblib/logitnormal.h"
#include "liblib/CtxtCoefficients.h"
#include "SUBSTABLE.h"
#include "TCTXVECT.h"


_TCTXVECT CVS;

// =========================================================================
// _TCTXVECT: set constants required for processing context vectors
//
_TCTXVECT::_TCTXVECT()
:
    AVGLEN( 1/*7,5*/),//length for averaging
    hAVGLEN( AVGLEN >> 1 ),
    AVGCWGT( 1.f/*0.20,0.40*/),//central weight in averaging
    AVGCOEFFS( AVGLEN, AVGCWGT ),//coefficients in averaging
    //
    CTXLEN( 15/*9*/),//context length
    CWGT( 0.40f/*0.30*/),//central weight of context coefficients
    MIX( true ),//mix normal vectors within the context
    AUG( false ),//augment vector 3-fold: into left, center, right partitions
    MEAN( true ),//mean vector in use
    //
    DIM( MIX? (NUMAA-1)*(AUG? 3: 1): CTXLEN*(NUMAA-1)),//dimensions of context vector
    KAPPA0( 0.1f ),//strictly positive
    NU0( (float)DIM ),//must be greater than dim-1
    MU0( NUMAA-1 ),
    loKAPPA0( 1.f/ KAPPA0 ),
    PowerNU0( -0.5f*(NU0+2.f)),
    CTERM(//constant term in calculating log probability of observations
      (0.5f*DIM+PowerNU0)*(logf(KAPPA0)-logf(KAPPA0+2.f)) )
{
    int ii;
    const int noeffress = NUMAA;
    float probs[noeffress];
    float sum = 0.0;
    for( ii = 1; ii <= DIM; ii++ )
        sum += logf(NU0+1.0f-ii) - SLC_LN2;
    const_cast<float&>(CTERM) += sum;
    //
    if( MEAN ) {
        for( ii = 0; ii < noeffress; ii++ )
            probs[ii] = Robinson_PROBS[ii];
        //transform to normal
        ::LogitNormal2Normal_f( probs, noeffress, 1.e-1f, false );
        for( ii = 0; ii < noeffress-1; ii++ )
            const_cast<extspsl::Pslvector&>(MU0).SetValueAt( ii, probs[ii]);
    }
    //
    const_cast<CtxtCoefficients&>(AVGCOEFFS).FindCoefficients();
}
