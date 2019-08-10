/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_rng__
#define __extspsl_rng__

namespace extspsl {

// -------------------------------------------------------------------------
// GNU SL implementations of the Random number generators
//
// -------------------------------------------------------------------------
// Interface for Random number generator
//
class Rng
{
public:
    Rng();
    virtual ~Rng();

    virtual unsigned int    Get() = 0;
    virtual float           GetSingle() = 0;
    virtual void            Set( unsigned int seed = 0 ) = 0;

    unsigned int            GetI( unsigned int interval );
    float                   GetSingle0();
    float                   GetSingle1();
    float                   GetSingle01();

    unsigned int            GetMin() const { return min_; }
    unsigned int            GetMax() const { return max_; }
    unsigned int            GetSize() const { return size_; }

protected:
    void    SetMin( unsigned int value ) { min_ = value; }
    void    SetMax( unsigned int value ) { max_ = value; }
    void    SetSize( unsigned int value ) { size_ = value; }

protected:
    unsigned int    min_;
    unsigned int    max_;
    unsigned int    size_;
};

// -------------------------------------------------------------------------
// Mersenne Twister RNG with seeding of 2002 release; period ~2^{19937} (GSL)
//  M.Matsumoto, T.Nishimura. (1998) Mersenne Twister: A 623-
//  dimensionally equidistributerd uniform pseudorandom number generator. 
//  ACM Transactions on Modeling and Computer Simulation, 8(1), 3-30.
//
class MTRng: public Rng
{
    enum {
        P1 = 624,
        P2 = 397
    };
public:
    MTRng();
    virtual ~MTRng();

    virtual unsigned int    Get();
    virtual float           GetSingle();
    virtual void            Set( unsigned int seed = 0 );

    static unsigned int     GetDefSeed() { return sD_seed; }

protected:
    static unsigned int     MAGIC( unsigned int y ) { return ( y & 1 )? 0x9908b0dfU: 0; }
    static unsigned int     GetMaskUPPER() { return sM_UPPER; }
    static unsigned int     GetMaskLOWER() { return sM_LOWER; }

private:
    unsigned int    mt_[P1];
    int             mti_;
    static const unsigned int   sM_UPPER;
    static const unsigned int   sM_LOWER;
    static const unsigned int   sD_seed;
};

// -------------------------------------------------------------------------
// Generalized Feedback Shift-Register (GFSR) based on xor-sum of particular 
//  past lagged values; period ~2^{9689} (GSL)
//  R.M.Ziff. (1998) Four-tap shift-register-sequence random-number
//  generators. Computers in Physics, 12(4), 385-392. 
//
class GFSR4Rng: public Rng
{
    enum {
        R1 = 471,
        R2 = 1586,
        R3 = 6988,
        R4 = 9689,
        P  = 16383, //2^14-1
        P1 = 16384
    };
public:
    GFSR4Rng();
    virtual ~GFSR4Rng();

    virtual unsigned int    Get();
    virtual float           GetSingle();
    virtual void            Set( unsigned int seed = 0 );

    static unsigned int     GetDefSeed() { return sD_seed; }

protected:
    static unsigned int     LCG( unsigned int s ) { return ( 69069 * s )& 0xffffffffU; }

private:
    unsigned int    rs_[P1];//register sequence
    int             ri_;
    static const unsigned int sD_seed;
};

// -------------------------------------------------------------------------
// Maximally equidistributed combined, collision free Tausworthe generator;
//  period ~2^{113} (GSL)
//  P.L'Ecuyer. (1999) Tables of Maximally-Equidistributed Combined LFSR 
//  Generators. Mathematics of Computation, 68(225), 261-269.
//  P.L'Ecuyer. (1996) Maximally Equidistributed Combined Tausworthe 
//  Generators. Mathematics of Computation, 65(213), 203-213.
//
class TausRng: public Rng
{
public:
    TausRng();
    virtual ~TausRng();

    virtual unsigned int    Get();
    virtual float           GetSingle();
    virtual void            Set( unsigned int seed = 0 );

    static unsigned int     GetDefSeed() { return sD_seed; }

protected:
    static unsigned int     LCG( unsigned int s ) { return ( 69069 * s )& 0xffffffffU; }

private:
    unsigned int z1_, z2_, z3_, z4_;
    static const unsigned int sD_seed;
};

// -------------------------------------------------------------------------
// Martin Luescher's second generation double-precision (48-bit) 
//  version of the RANLUX (lagged fibonacci) generator (GSL)
//  M.Luescher. (1994) A portable high-quality random number generator
//  for lattice field theory calculations. Computer Physics Communications, 
//  79, 100-110.
//
class RANLUXdRng: public Rng
{
    enum {
        P1 = 12
    };
public:
    RANLUXdRng();
    virtual ~RANLUXdRng();

    virtual unsigned int    Get();
    virtual float           GetSingle();
    virtual void            Set( unsigned int seed = 0 ) { Set1( seed ); }
    void    Set1( unsigned int seed = 0 ) { SetLux( seed, 202 ); }
    void    Set2( unsigned int seed = 0 ) { SetLux( seed, 397 ); }

    static unsigned int     GetDefSeed() { return sD_seed; }

protected:
    void    SetLux( unsigned int seed, unsigned int luxury );
    void    Increment();
    void    RANLUXstep( float&, float&, unsigned int, unsigned int, unsigned int );

private:
    float           reg_[P1];
    float           carrier_;
    unsigned int    iti_;
    unsigned int    itj_;
    unsigned int    itiprev_;
    unsigned int    itmax_;
    //
    int             nindxs_[P1];
    static const float          s_mprec;
    static const unsigned int   sD_seed;
};

// -------------------------------------------------------------------------
//
inline
void RANLUXdRng::RANLUXstep( float& r1, float& r2, 
        unsigned int it1, unsigned int it2, unsigned int it3 )
{
    r1 = reg_[it1] - reg_[it2];
    if( r2 < 0 ) {
        r1 -= s_mprec;
        r2 += 1;
    }
    reg_[it3] = r2;
}

}//namespace extspsl

#endif//__extspsl_rng__
