/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SUBSTABLE_h__
#define __SUBSTABLE_h__

#include <stddef.h>
#include "liblib/alpha.h"

class mystring;

// AB Robinson and LR Robinson background probabilities
// (PNAS USA 88, 1991, 8880-4)
extern const float Robinson_PROBS[NUMAA];

// Background probabilities by pscores
extern const float Pscores_PROBS[NUMAA];

// _________________________________________________________________________
// CLASS SUBSTABLE (substitution table)
//
class SUBSTABLE
{
public:
    enum TScoreClass {
        Class80,
        Class62,
        Class45,
        ClassGn,
        ClassPS
    };
    SUBSTABLE();

    TScoreClass GetClass() const         { return class_; }
    void SetClass( TScoreClass value )   { class_ = value; }

    float   PROBABility( int a );
    float   LogPROBABility( int a );

    float   PROBABILITY_1( int a );
    float   LogPROBABILITY_1( int a );

    float   PROBABILITY_2( int a );
    float   LogPROBABILITY_2( int a );

    float   FreqRatio( int a, int b );
    float   PrecomputedEntry( int a, int b );
    int     Entry( int a, int b );
    float   StatisParam( int scheme, int field );

    void    StoreProbabilities_1( const float* probs );
    void    RestoreProbabilities_1();

    void    StoreProbabilities_2( const float* probs );
    void    RestoreProbabilities_2();

    void    StoreProbabilities( const float* probs );
    void    RestoreProbabilities();

protected:
    const float* GetProbabilities_1() { return probs_1_; }
    const float* GetProbabilities_2() { return probs_2_; }
    const float* GetProbabilities() { return auxprobabs_; }

private:
    void    ComputeLogProbabilities();
    void    ComputeLogProbabilities( const float* probs, float* logprobs );
    void    InitLogProbabs( float* logprobs );

private:
    TScoreClass     class_;
    float           logprobs_[NUMAA];
    const float*    probs_1_;
    float           logprobs_1_[NUMAA];
    const float*    probs_2_;
    float           logprobs_2_[NUMAA];
    const float*    auxprobabs_;
    float           auxlogps_[NUMAA];
};
extern SUBSTABLE    STABLE;
void SetSTABLE( const mystring&, const char* filename = NULL );
// //


// Frequency ratios for BLOSUM80
extern const float BLOSUM80_FREQRATIOS[NUMAA][NUMAA];

// Frequency ratios for BLOSUM62 as determined by Stephen Altschul;
//  Stephen and Jorja Henikoff used different number for B, Z, X.
//  Each entry in the table equals to substitution frequency qij
//  devided by the product of background probabilities pi * pj
// extern const double BLOSUM62_FREQRATIOS[ NUMALPH ][ NUMALPH ];

// Frequency ratios for BLOSUM62 computed purely from BLOSUM62 scores;
//  Each entry in the table equals to target probability qij
//  devided by the product of marginal probabilities pi * pj
//  (computed originaly)
extern const float BLOSUM62_FREQRATIOS[NUMAA][NUMAA];

// Frequency ratios for BLOSUM45
extern const float BLOSUM45_FREQRATIOS[NUMAA][NUMAA];


// Frequency ratios obtained by program pscores from this (or previous version of) software package
extern const float PSCORES_FREQRATIOS[NUMAA][NUMAA];


// -------------------------------------------------------------------------
// (Gonnet marginal probabilities solved from O p = 1, where 
//  O is odds matrix, and p is probability vector to solve for)
extern const float Gonnet_PROBS[NUMAA];

// Gonnet frequency ratios;
extern const float GONNET_FREQRATIOS[NUMAA][NUMAA];


// -------------------------------------------------------------------------
// Computed BLOSUM80 table given frequency ratios
struct _TCOMPUTED_BLOSUM80
{
    _TCOMPUTED_BLOSUM80();
    float operator()( int a, int b ) { return data[a][b]; }
public:
    float data[NUMAA][NUMAA];
};
extern _TCOMPUTED_BLOSUM80  COMPUTED_BLOSUM80;
//
// Computed BLOSUM62 table given frequency ratios
struct _TCOMPUTED_BLOSUM62
{
    _TCOMPUTED_BLOSUM62();
    float operator()( int a, int b ) { return data[a][b]; }
public:
    float data[NUMAA][NUMAA];
};
extern _TCOMPUTED_BLOSUM62  COMPUTED_BLOSUM62;
//
// Computed BLOSUM45 table given frequency ratios
struct _TCOMPUTED_BLOSUM45
{
    _TCOMPUTED_BLOSUM45();
    float operator()( int a, int b ) { return data[a][b]; }
public:
    float data[NUMAA][NUMAA];
};
extern _TCOMPUTED_BLOSUM45  COMPUTED_BLOSUM45;
//
// Computed PSCORES_ table given frequency ratios
struct _TCOMPUTED_PSCORES_
{
    _TCOMPUTED_PSCORES_();
    float operator()( int a, int b ) { return data[a][b]; }
    void Compute();
public:
    float data[NUMAA][NUMAA];
};
extern _TCOMPUTED_PSCORES_  COMPUTED_PSCORES_;

//
// Computed GONNET table given frequency ratios
struct _TCOMPUTED_GONNET_
{
    _TCOMPUTED_GONNET_();
    float operator()( int a, int b ) { return data[a][b]; }
    void Compute();
private:
public:
    float data[NUMAA][NUMAA];
};
extern _TCOMPUTED_GONNET_  COMPUTED_GONNET;


//Substitution table: 'BLOSUM80'
extern const int BLOSUM80[NUMAA][NUMAA];

//Substitution table: 'BLOSUM62'
extern const int BLOSUM62[NUMAA][NUMAA];

//Substitution table: 'BLOSUM45'
extern const int BLOSUM45[NUMAA][NUMAA];

//Substitution table produced by pscores
extern const int PSCORES_[NUMAA][NUMAA];



extern const float Blosum80ScalingConstant;// = 2.0;
extern const float Blosum62ScalingConstant;// = 2.0;
extern const float Blosum45ScalingConstant;// = 3.0;
extern const float PScores_ScalingConstant;// = 3.0;
// extern const float Gonnet_ScalingConstant;// = 4.0;
extern const float Gonnet_ScalingConstant;// = 3.0;



enum TBlosumFields {
    Open,       //gap existence penalty
    Extend,     //gap extension penalty
    Lambda,     //statistical parameter Lambda
    K,          //statistical parameter K
    H,          //statistical parameter H
    alpha,      //statistical parameter alpha (slope of edge-effect correction with linear regression)
    beta,       //statistical parameter beta (intercept of edge-effect correction with linear regression)
    NumFields
};

enum TBlosum80Entries {
    Ungapped80,
    g25_2_80,
    g13_2_80,
    g9_2_80,
    g8_2_80,
    g7_2_80,
    g6_2_80,
    g11_1_80,
    g10_1_80,
    g9_1_80,
    Num80Entries
};

// (Data from the NCBI BLAST toolkit)
extern const float BLOSUM80_VALUES[Num80Entries][NumFields];

enum TBlosumEntries {
    Ungapped,
    g11_2,
    g10_2,
    g9_2,
    g8_2,
    g7_2,
    g6_2,
    g13_1,
    g12_1,
    g11_1,
    g10_1,
    g9_1,
    NumEntries
};

// (Data from the NCBI BLAST toolkit)
extern const float BLOSUM62_VALUES[NumEntries][NumFields];

enum TBlosum45Entries {
    Ungapped45,
    g13_3_45,
    g12_3_45,
    g11_3_45,
    g10_3_45,
    g16_2_45,
    g15_2_45,
    g14_2_45,
    g13_2_45,
    g12_2_45,
    g19_1_45,
    g18_1_45,
    g17_1_45,
    g16_1_45,
    Num45Entries
};

// (Data from the NCBI BLAST toolkit)
extern const float BLOSUM45_VALUES[Num45Entries][NumFields];

enum TPscores_Entries {
    UngappedPS,
    NumPSEntries
};

extern const float PSCORES_VALUES[NumPSEntries][NumFields];

enum TGonnet_Entries {
    UngappedGn,
    NumGnEntries
};

extern const float GONNET_VALUES[NumGnEntries][NumFields];

#endif//__SUBSTABLE_h__
