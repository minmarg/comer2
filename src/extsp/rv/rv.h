/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_rv__
#define __extspsl_rv__

#include "extsp/rng.h"

namespace extspsl {

// -------------------------------------------------------------------------
// Interface for random variable
//
class RVar
{
public:
    RVar( Rng& rng ): rng_( rng ) {};
    virtual ~RVar() {};

    virtual int     Gen( float* ) = 0;
    Rng&            GetRng() { return rng_; }

protected:
    Rng&    rng_;
};

}//namespace extspsl

#endif//__extspsl_rv__
