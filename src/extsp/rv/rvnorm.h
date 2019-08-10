/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_rvnorm__
#define __extspsl_rvnorm__

#include <stddef.h>
#include "pslerror.h"
#include "rng.h"
#include "rv.h"

namespace extspsl {

// -------------------------------------------------------------------------
// Normal random variable
//
class RVNorm: public RVar
{
public:
    enum TRVN {
        TRVN_Polar,
        TRVN_Ratio,
        TRVN_Ziggurat_M00
    };
    RVNorm( Rng& rng, TRVN met = TRVN_Ziggurat_M00 );
    virtual ~RVNorm();

    virtual int     Gen( float* res ) { return Gen( res, (float*)0 ); }
    int             Gen( float*, float* );
    int             GenPolar( float*, float* = (float*)0 );
    int             GenRatio( float* );
    int             GenZigMsg00( float* );

    TRVN    GetMet() const { return met_; }
    void    SetMet( TRVN value ) { met_ = value; }

    float   GetStd() const { return std_; }
    void    SetStd( float value ) { std_ = value; }

    float   GetMean() const { return mean_; }
    void    SetMean( float value ) { mean_ = value; }

private:
    TRVN    met_;
    float   std_;
    float   mean_;
};

// -------------------------------------------------------------------------
//
inline
int RVNorm::Gen( float* rv1, float* rv2 )
{
    if( rv1 == NULL )
        return PSL_ERR_ADDRESS;
    int err;
    float   std = GetStd();
    float   mean = GetMean();
    bool    brv2 = GetMet() == TRVN_Polar && rv2 != NULL;

    switch( GetMet()) {
      case TRVN_Polar: err = GenPolar( rv1, rv2 ); break;
      case TRVN_Ratio: err = GenRatio( rv1 ); break;
      case TRVN_Ziggurat_M00: err = GenZigMsg00( rv1 ); break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    if( std != 1.0f ) {
        *rv1 *= std;
        if( brv2 )
            *rv2 *= std;
    }
    if( mean ) {
        *rv1 += mean;
        if( brv2 )
            *rv2 += mean;
    }
    return err;
}

}//namespace extspsl

#endif//__extspsl_rvnorm__
