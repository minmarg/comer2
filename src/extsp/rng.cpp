/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pslerror.h"
#include "rng.h"

namespace extspsl {

// -------------------------------------------------------------------------
//
Rng::Rng()
{
}

Rng::~Rng()
{
}

// Get: get random number in interval from min to min+interval(<max)
//
unsigned int Rng::GetI( unsigned int interval )
{
    unsigned int    offset = GetMin();
    unsigned int    range = GetMax() - offset;
    unsigned int    scale;
    unsigned int    rn = interval;
    const int maxn = 10;
    int n;

    if( interval == 0 )
        return 0;
    if( range < interval )
        return Get();

    scale = range / interval;

    for( n = 0; interval <= rn && n < maxn; n++ ) {
        rn = ( Get() - offset )/ scale;
    }
    return rn;
}

// GetSingle0: get random number in interval ]0,1]
//
float Rng::GetSingle0()
{
    float u = 0.0;
    const int maxn = 10;
    int n;
    for( n = 0; u == 0.0f && n < maxn; n++ )
        u = GetSingle();
    if( u == 0.0f )
        PRIVERROR("Rng::GetSingle0: Failed to generate RN.");
    return u;
}

// GetSingle1: get random number in interval [0,1[
//
float Rng::GetSingle1()
{
    float u = 1.0f;
    const int maxn = 10;
    int n;
    for( n = 0; u == 1.0f && n < maxn; n++ )
        u = GetSingle();
    if( u == 1.0f )
        PRIVERROR("Rng::GetSingle1: Failed to generate RN.");
    return u;
}

// GetSingle01: get random number in interval ]0,1[
//
float Rng::GetSingle01()
{
    float u = 0.0f;
    const int maxn = 10;
    int n;
    for( n = 0; ( u == 0.0f || u == 1.0f ) && n < maxn; n++ )
        u = GetSingle();
    if( u == 0.0f || u == 1.0f )
        PRIVERROR("Rng::GetSingle01: Failed to generate RN.");
    return u;
}

// =========================================================================
// =========================================================================
// MTRng
const unsigned int MTRng::sM_UPPER = 0x80000000U;
const unsigned int MTRng::sM_LOWER = 0x7fffffffU;
const unsigned int MTRng::sD_seed = 4357;
// =========================================================================
// MTRng
//
MTRng::MTRng()
{
    SetMin( 0 );
    SetMax( 0xffffffffU );
    SetSize( sizeof(mti_) + sizeof(*mt_)*P1 );
    Set();
}

MTRng::~MTRng()
{
}

// -------------------------------------------------------------------------
// Get: generate random number
//
unsigned int MTRng::Get()
{
    unsigned int k, y;
    int kk;

    if( P1 <= mti_ ) 
    {   //generate P1 words at one time
        for( kk = 0; kk < P1 - P2; kk++ ) {
            y = ( mt_[kk] & sM_UPPER ) | ( mt_[kk+1] & sM_LOWER );
            mt_[kk] = mt_[kk+P2] ^ ( y >> 1 ) ^ MAGIC( y );
        }
        for (; kk < P1 - 1; kk++ ) {
            y = ( mt_[kk] & sM_UPPER ) | ( mt_[kk+1] & sM_LOWER );
            mt_[kk] = mt_[kk+(P2-P1)] ^ ( y >> 1 ) ^ MAGIC( y );
        }
        y = ( mt_[P1-1] & sM_UPPER ) | ( mt_[0] & sM_LOWER );
        mt_[P1-1] = mt_[P2-1] ^ ( y >> 1 ) ^ MAGIC( y );

        mti_ = 0;
    }
    //tempering
    k = mt_[mti_];
    k ^= ( k >> 11 );
    k ^= ( k <<  7 ) & 0x9d2c5680U;
    k ^= ( k << 15 ) & 0xefc60000U;
    k ^= ( k >> 18 );

    mti_++;
    return k;
}

// -------------------------------------------------------------------------
// GetSingle: get single-precision value in interval [0;1)
//
float MTRng::GetSingle()
{
    if( !GetMax())
        return 0.0f;
    return Get() / (float)GetMax();
}

// -------------------------------------------------------------------------
// Set: initialize rng with seed
//
void MTRng::Set( unsigned int seed )
{
    int n;

    if( !seed )
        seed = sD_seed;//default seed

    mt_[0] = seed & 0xffffffffU;

    for( n = 1; n < P1; n++ ) {
        //GSL: cf. Knuth "Art of Computer Programming" Vol. 2, 3rd ed. p.106 for multiplier
        mt_[n]  = 1812433253U * ( mt_[n-1] ^ ( mt_[n-1] >> 30 )) + n;
        mt_[n] &= 0xffffffffU;
    }
    mti_ = n;
}

// =========================================================================
// =========================================================================
// GFSR4Rng
const unsigned int GFSR4Rng::sD_seed = 4357;
// =========================================================================
// GFSR4Rng
//
GFSR4Rng::GFSR4Rng()
{
    SetMin( 0 );
    SetMax( 0xffffffffU );
    SetSize( sizeof(ri_) + sizeof(*rs_)*P1 );
    Set();
}

GFSR4Rng::~GFSR4Rng()
{
}

// -------------------------------------------------------------------------
// Get: generate random number
//
unsigned int GFSR4Rng::Get()
{
    ri_ = ( ri_ + 1 ) & P;
    return rs_[ri_] =
        rs_[( ri_ + ( P1 - R1 )) & P] ^
        rs_[( ri_ + ( P1 - R2 )) & P] ^
        rs_[( ri_ + ( P1 - R3 )) & P] ^
        rs_[( ri_ + ( P1 - R4 )) & P];
}

// -------------------------------------------------------------------------
// GetSingle: get single-precision value in interval [0;1)
//
float GFSR4Rng::GetSingle()
{
    if( !GetMax())
        return 0.0f;
    return Get() / (float)GetMax();
}

// -------------------------------------------------------------------------
// Set: initialize rng with seed
//
void GFSR4Rng::Set( unsigned int seed )
{
    int n, m, k;
    unsigned int msku = 0x80000000U;
    unsigned int mask = 0xffffffffU;
    unsigned int sp, b;

    if( !seed )
        seed = sD_seed;

    //to avoid low-order bit correlations
    for( n = 0; n < P1; n++ ) {
        sp = 0 ;
        b = msku;
        for( m = 0; m < 32; m++ ) {
            seed = LCG( seed );
            if( seed & msku ) 
                sp |= b;
            b >>= 1;
        }
        rs_[n] = sp;
    }

    //matrix orthogonalisation
    for( n = 0; n < 32; n++ ) {
        k = n * 3 + 7;
        rs_[k] &= mask;
        rs_[k] |= msku;
        mask >>= 1;
        msku >>= 1;
    }
    ri_ = n;
}

// =========================================================================
// =========================================================================
// TausRng
const unsigned int TausRng::sD_seed = 1;
// =========================================================================
// TausRng
//
TausRng::TausRng()
{
    SetMin( 0 );
    SetMax( 0xffffffffU );
    SetSize( sizeof(z1_) + sizeof(z2_) + sizeof(z3_) + sizeof(z4_));
    Set();
}

TausRng::~TausRng()
{
}

// -------------------------------------------------------------------------
// Get: generate random number
//
unsigned int TausRng::Get()
{
    unsigned int b1, b2, b3, b4;
    unsigned int mask = 0xffffffffU;

    b1 = ((( z1_ << 6 ) & mask ) ^ z1_ ) >> 13;
    z1_ = ((( z1_ & 4294967294U ) << 18 ) & mask ) ^ b1;

    b2 = ((( z2_ << 2 ) & mask ) ^ z2_ ) >> 27;
    z2_ = ((( z2_ & 4294967288U ) << 2 ) & mask ) ^ b2;

    b3 = ((( z3_ << 13 ) & mask ) ^ z3_ ) >> 21;
    z3_ = ((( z3_ & 4294967280U ) << 7 ) & mask ) ^ b3;

    b4 = ((( z4_ << 3 ) & mask ) ^ z4_ ) >> 12;
    z4_ = ((( z4_ & 4294967168U ) << 13 ) & mask ) ^ b4;

    return z1_ ^ z2_ ^ z3_ ^ z4_;
}

// -------------------------------------------------------------------------
// GetSingle: get single-precision value in interval [0;1)
//
float TausRng::GetSingle()
{
    if( !GetMax())
        return 0.0f;
    return Get() / (float)GetMax();
}

// -------------------------------------------------------------------------
// Set: initialize rng with seed
//
void TausRng::Set( unsigned int seed )
{
    int n;
    if( !seed )
        seed = sD_seed;

    z1_ = LCG( seed );
    if( z1_ < 2 )
        z1_ += 2;
    z2_ = LCG( z1_ );
    if( z2_ < 8 )
        z2_ += 8;
    z3_ = LCG( z2_ );
    if( z3_ < 16 )
        z3_ += 16;
    z4_ = LCG( z3_ );
    if( z4_ < 128 )
        z4_ += 128;

    //to satify recurrence condition
    for( n = 0; n < 10; n++ ) Get();
}

// =========================================================================
// =========================================================================
// RANLUXdRng
const float         RANLUXdRng::s_mprec = 1.0f / 281474976710656.0f;//1/2^48
const unsigned int  RANLUXdRng::sD_seed = 1;
// =========================================================================
// RANLUXdRng
//
RANLUXdRng::RANLUXdRng()
{
    int n;
    SetMin( 0 );
    SetMax( 0xffffffffU );
    SetSize( sizeof(*reg_)*P1 + sizeof( carrier_ ) +
      sizeof(iti_) + sizeof(itj_) + sizeof(itiprev_) + sizeof(itmax_) + 
      sizeof(*nindxs_)*P1 );
    for( n = 0; n < P1 - 1; n++ )
        nindxs_[n] = n + 1;
    nindxs_[n] = 0;
    Set();
}

RANLUXdRng::~RANLUXdRng()
{
}

// -------------------------------------------------------------------------
// Increment: progress state
//
void RANLUXdRng::Increment()
{
    int     n, nmax;
    float   r1, r2, r3;

    for( n = 0; 0 < iti_; n++ ) {
        r1 = reg_[itj_] - reg_[iti_];
        r2 = r1 - carrier_;
        if( r2 < 0 ) {
            carrier_ = s_mprec;
            r2 += 1;
        }
        else
            carrier_ = 0;
        reg_[iti_] = r2;
        iti_ = nindxs_[iti_];
        itj_ = nindxs_[itj_];
    }

    nmax = itmax_ - 12;

    for( ; n <= nmax; n += 12 ) {
        r1 = reg_[7] - reg_[0];
        r1 -= carrier_;

        RANLUXstep( r2, r1, 8, 1, 0 );
        RANLUXstep( r3, r2, 9, 2, 1 );
        RANLUXstep( r1, r3, 10, 3, 2 );
        RANLUXstep( r2, r1, 11, 4, 3 );
        RANLUXstep( r3, r2, 0, 5, 4 );
        RANLUXstep( r1, r3, 1, 6, 5 );
        RANLUXstep( r2, r1, 2, 7, 6 );
        RANLUXstep( r3, r2, 3, 8, 7 );
        RANLUXstep( r1, r3, 4, 9, 8 );
        RANLUXstep( r2, r1, 5, 10, 9 );
        RANLUXstep( r3, r2, 6, 11, 10 );

        if( r3 < 0 ) {
            carrier_ = s_mprec;
            r3 += 1;
        }
        else
            carrier_ = 0;
        reg_[11] = r3;
    }

    nmax = itmax_;

    for( ; n < nmax; n++ ) {
        r1 = reg_[itj_] - reg_[iti_];
        r2 = r1 - carrier_;
        if( r2 < 0 ) {
            carrier_ = s_mprec;
            r2 += 1;
        }
        else
            carrier_ = 0;
        reg_[iti_] = r2;
        iti_ = nindxs_[iti_];
        itj_ = nindxs_[itj_];
    }
    itiprev_ = iti_;
}

// -------------------------------------------------------------------------
// Get: generate random number
//
unsigned int RANLUXdRng::Get()
{
    return (unsigned int)( GetSingle() * (float)GetMax());
}

// -------------------------------------------------------------------------
// GetSingle: get single-precision value in interval [0;1)
//
float RANLUXdRng::GetSingle()
{
    iti_ = nindxs_[iti_];
    if( iti_ == itiprev_ )
        Increment();
    return reg_[iti_];
}

// -------------------------------------------------------------------------
// SetLux: initialize rng with the given luxury level and seed
//
void RANLUXdRng::SetLux( unsigned int seed, unsigned int luxury )
{
    int ibit, jbit, s, n, m, xbit[31];
    float r, rr;

    if( !seed )
        seed = sD_seed;
    s = seed & 0xffffffffU;

    for( n = 0; n < 31; n++ ) {
        xbit[n] = s % 2;
        s /= 2;
    }
    ibit = 0;
    jbit = 18;

    for( n = 0; n < 12; n++ ) {
        r = 0;
        for( m = 1; m <= 48; m++ ) {
            rr = (float)(( xbit[ibit] + 1 ) % 2 );
            r += r + rr;
            xbit[ibit] = ( xbit[ibit] + xbit[jbit]) % 2;
            ibit = ( ibit + 1 ) % 31;
            jbit = ( jbit + 1 ) % 31;
        }
        reg_[n] = s_mprec * r;
    }

    carrier_ = 0.0f;
    iti_ = 11;
    itj_ = 7;
    itiprev_ = 0;
    itmax_ = luxury;
}

}//namespace extspsl
