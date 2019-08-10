/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __HDPbase__
#define __HDPbase__

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>

#include "extsp/psl.h"
#include "extsp/ivector.h"
#include "hdplib.h"
#include "basin.h"
#include "dish.h"
#include "menu.h"
#include "rntchain.h"
#include "HDPscores.h"
// #include "HDP_SSSScores.h"
// #include "iHDP_SSSScores.h"

// -------------------------------------------------------------------------
// Choices for the application of the HDP mixture model --------------------
// mixing of profile target frequencies
enum TTFrMix {
    tfrmixNo,       //no mixing
    tfrmixHDPCtx    //mix target frequencies under the HDP framework
};
//score adjustment
enum TScoAdj {
    scoadjNo,       //no adjustment
    scoadjHDPCtx,   //score adjustment under the HDP framework
    scoadjHDPsco    //apply scores derived from HDP cluster membership 
                    //  distribution along structure alignments
};
// End of choices ----------------------------------------------------------


//global object for application only
class HDPbase;
extern HDPbase  HDPBASE;
// extern HDPbase  HDPctBASE;
void SetHDPBASE( const char* filename );
// void SetHDPctBASE( const char* filename );

//define deg.of freedom nu of student's t distribution
// #define NU_t_o2 ( nu_t_o2 - do2 ) //original deg. of freedom
#define NU_t_o2 ( nu_t_o2 + do2 * GetDegFAdjustment()) //adjusted deg.of freedom

// -------------------------------------------------------------------------
// class HDPbase: base class for a Hierarchical Dirichlet Process 
// mixture model
//
class HDPbase
{
public:
    enum TResType {
        TRT_MHUpdate,//results obtained from MH update
        TRT_GibbsUpdate,//results obtained from Gibbs sampling update
        TRT_novalues
    };

public:
    HDPbase();
    virtual ~HDPbase();

    bool        GetUninfPrior() const { return uninfprior_; }
    void        SetUninfPrior( bool value ) { uninfprior_ = value; }

    bool        GetCluster4Each() const { return clst4each; }
    void        SetCluster4Each( bool value ) { clst4each = value; }

    float       GetS0ScaleFac() const { return S0scalefac_; }
    void        SetS0ScaleFac( float value ) { S0scalefac_ = value; }

    float       GetDegFAdjustment() const { return adjdegf_; }
    void        SetDegFAdjustment( float value ) { adjdegf_ = value; }


    float       GetMixWeight() const { return mixweight_; }
    void        SetMixWeight( float value ) { mixweight_ = value; }

    const HDPscores* GetScores() const {return scores_; }
    void        SetScores( const HDPscores* value );


    void        SetPriorProbs();

    int         GetReadNoGroups() const { return rdnorsts_; }
    void        SetReadNoGroups( int value ) { rdnorsts_ = value; }

    int         GetNoSupClusters() const { return nosupclsts_; }
    void        SetNoSupClusters( int value ) { nosupclsts_ = value + 1; }

    float       GetAdjWeight() const { return adjweight_; }
    void        SetAdjWeight( float value ) { adjweight_ = value; }


    void        ReserveMenu( int );
    void        ReserveBasin( int );
    int         AddToBasin( extspsl::Pslvector* nv );

    void        CalcPPProbs( extspsl::Pslvector& lnvar, float* bppr, 
                             extspsl::Pslvector* ppprobs, extspsl::Ivector* cindcs, 
                             float lpfact = 0.02f, bool usepriors = true, bool tonormal = true ) const;
    void        MixCLNVar( extspsl::Pslvector& lnvar, extspsl::Pslvector* lnmixed ) const;

    bool        CalcPriorParams( bool addnoise = false, float std = 1.0f );
    int         CalcPriorParams( const extspsl::Pslvector& stds );
    int         SetUninfPriorParamsS0I( int dim, int ctx, 
                                        const extspsl::Pslvector& mvec, float = -1.0f, float = -1.0f );
    void        CalcPriorProbFact();
    void        RecalcMenuParams();
    void        RecalcDishParams( int k );
    void        CalcProbOfDish( int k, float* lprob );
    void        CalcProbOfData( float* lprob );

    void            ReadGroups( const char*, extspsl::Ivector* dids );
    void            PrintGroups( const char*, float* = NULL, int* = NULL );

    void            ReadParameters( const char*, const extspsl::Ivector* dids = NULL );
    virtual void    PrintParameters( const char* );

    void            PrintDishes( const char*, float* = NULL, int* = NULL );


    const Basin*    GetBasin() const { return basin_; }
    Basin*          GetBasin() { return basin_; }

    const Menu*     GetMenu() const { return menu_; }
    Menu*           GetMenu() { return menu_; }
    void            SetMenu( Menu* value ) { DestroyMenu(); menu_ = value; }
    void            ReplaceMenu( Menu* value ) { menu_ = value; }

    const RntChain* GetChain() const { return chain_; }
    RntChain*       GetChain() { return chain_; }

    int         GetCtxtSize() const { return ctxtsize_; }
    float       GetDPMTau() const { return tau_; }
    float       GetDPMGamma() const { return gamma_; }

public:
    ///{{M-VARIATE PROBABILITIES
    //prior
    void        PriorProbVec( const extspsl::Pslvector*, float* lprob ) const;
    void        PriorProbVecHlp( const extspsl::Pslvector*, float* lprob ) const;
    void        PriorProbVecHlpObs( const extspsl::Pslvector*, float* lprob ) const;
    //posterior predictive
    void        ProbVecOfDish( const extspsl::Pslvector*, int k, float* lprob ) const;
    void        ProbVecOfDishHlp( const extspsl::Pslvector*, int k, float* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbVecOfDishHlpObs( const extspsl::Pslvector*, int k, float* lprob ) const;
#endif
    //exclusive posterior predictive
    void        ProbVecOfDishExc( const extspsl::Pslvector*, int k, float* lprob ) const;
    void        ProbVecOfDishExcHlp( const extspsl::Pslvector*, int k, float* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbVecOfDishExcHlpObs( const extspsl::Pslvector*, int k, float* lprob ) const;
#endif
    ///}}

    ///{{MATRIX PROBABILITIES
    //prior matrix probabilities
    void        PriorProbMtx( const Table*, float* lprob ) const;
    void        PriorProbMtx( const extspsl::Pslmatrix*, float* lprob ) const;
    void        PriorProbMtxHlp( const Table*, float* lprob ) const;
    void        PriorProbMtxHlp( const extspsl::Pslmatrix*, float* lprob ) const;
    void        PriorProbMtxHlpObs( const extspsl::Pslmatrix*, float* lprob ) const;
    //posterior predictive matrix probabilities
    void        ProbMtxOfDish( const Table*, int k, float* lprob ) const;
    void        ProbMtxOfDish( const extspsl::Pslmatrix*, int k, float* lprob ) const;
    void        ProbMtxOfDishHlp( const Table*, int k, float* lprob ) const;
    void        ProbMtxOfDishHlp( const extspsl::Pslmatrix*, int k, float* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbMtxOfDishHlpObs( const extspsl::Pslmatrix*, int k, float* lprob ) const;
#endif
    //exclusive posterior predictive matrix probabilities (not including given data)
    void        ProbMtxOfDishExc( const Table*, int k, float* lprob ) const;
    void        ProbMtxOfDishExc( const extspsl::Pslmatrix*, int k, float* lprob ) const;
    void        ProbMtxOfDishExcHlp( const Table*, int k, float* lprob ) const;
    void        ProbMtxOfDishExcHlp( const extspsl::Pslmatrix*, int k, float* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbMtxOfDishExcHlpObs( const extspsl::Pslmatrix*, int k, float* lprob ) const;
#endif
    ///}}

protected:
    int         GetIterationRead() const { return iterread_; }
    void        SetIterationRead( int value ) { iterread_ = value; }

    float       GetMaxLProbData() const { return mlpd_; }
    void        SetMaxLProbData( float value ) { mlpd_ = value; }

    float       GetLastLProbData() const { return lastlpd_; }
    void        SetLastLProbData( float value ) { lastlpd_ = value; }

    TResType    GetResType() const { return restype_; }
    void        SetResType( TResType value ) { restype_ = value; }

    int         GetMHUpdateNum() const { return mhupdate_; }
    void        SetMHUpdateNum( int value ) { mhupdate_ = value; }

    int         GetGibbsIt() const { return gibbsit_; }
    void        SetGibbsIt( int value ) { gibbsit_ = value; }

    bool        GetRestarted() const { return restarted_; }
    void        SetRestarted( bool value ) { restarted_ = value; }

    void        SetDPMTau( float value ) { tau_ = value; }
    void        SetDPMGamma( float value ) { gamma_ = value; }

    void            ReadParameters( FILE*, const extspsl::Ivector* );
    bool            ReadSample( FILE*, extspsl::Pslvector**, int dim = 0 );
    virtual void    PrintParameters( FILE* );

    void            ReadPriorParams( FILE* );
    void            PrintPriorParams( FILE* );

    void            ReadDishParams( FILE*, const extspsl::Ivector* );
    void            PrintDishParams( FILE* );

    void            PrintGroups( FILE*, float* = NULL, int* = NULL );
    void            PrintDishes( FILE*, float* = NULL, int* = NULL );
    void            PrintSummary( const char*, float* = NULL, int* = NULL );
    void            PrintSummary( FILE*, float* = NULL, int* = NULL );
    void            PrintBasin( FILE* );//TEST

    void        Table2Mtx( const Table*, extspsl::Pslmatrix* ) const;

    void        DestroyBasin();
    void        DestroyMenu();
    void        DestroyChain();

    void        InitBasin( int size );
    void        InitMenu( int size );
    void        InitChain( int size );

    bool        TestMenu() const;

    void        ResetTotalNoSamples()       { totalnos_ = 0; }
    void        IncTotalNoSamples()         { totalnos_++; }
    int         GetTotalNoSamples() const   { return totalnos_; }

    void        SetCtxtSize( int value ) { ctxtsize_ = value; }

    static float GetDefDPMTau() { return s_defdpmtau_; }
    static float GetDefDPMGamma() { return s_defdpmgamma_; }

    static int  GetDefSampleDim() { return s_defsampledim_; }
    static int  GetDefDishSize() { return s_defdishsize_; }
    static int  GetDefTableSize() { return s_deftablesize_; }
    static int  GetDefBasinSize() { return s_defbasinsize_; }

    static float GetDefKappa_pp_a() {return s_defkappa_pp_a_; }
    static float GetDefKappa_pp_b() {return s_defkappa_pp_b_; }
    static float GetDefKappaParam() {return s_defkappa_; }
    static float GetDefNuParam() { return s_defnu_; }

private:
    bool        uninfprior_; //uninformative prior
    bool        clst4each;  //assign each sample to new cluster
    Basin*      basin_; //basin of entire set of vectors
    Menu*       menu_;  //menu of dishes (global clusters)
    RntChain*   chain_; //chain of restaurants (groups)
    int         ctxtsize_;  //size of context
    int         rdnorsts_; //proposal (read) number of restaurants (groups)
    int         totalnos_;  //total number of samples
    TResType    restype_;  //type of results
    int         mhupdate_;  //MH update number
    int         gibbsit_;  //Gibbs iteration
    bool        restarted_;  //just restarted

    //{{application parameters
    float       mixweight_; //PME weight of mixing
    int         nosupclsts_;//number of support clusters (dishes)
    float       adjweight_; //weight of score adjustment under HDP framework
    const HDPscores* scores_;//dish id-id scores
    //}}

    float       tau_;       //DPM concentration parameter
    float       gamma_;     //DPM concentration parameter

    static float s_defdpmtau_;//default value of the DPM concentration parameter tau_
    static float s_defdpmgamma_;//default value of the DPM concentration parameter gamma_

    static int s_defsampledim_;//default dimensionality of samples
    static int s_defdishsize_;//default dish size
    static int s_deftablesize_;//default table size
    static int s_defbasinsize_;//default basin size

    static float s_defkappa_pp_a_;//default value of the prior param. a for kappa
    static float s_defkappa_pp_b_;//default value of the prior param. b for kappa
    static float s_defkappa_;//default value of the NIW kappa parameter (scale parameter)
    static float s_defnu_;//default value of the NIW nu parameter (deg. of freedom)

    int         iterread_;//iteration number read
    float       mlpd_;//maximum value of log probability of data
    float       lastlpd_;//last obtained value of log probability of data

    float       S0scalefac_;//S0 scale factor
    float       adjdegf_;//adjustment to deg. of freedom in terms of dim. over 2
};


// -------------------------------------------------------------------------
// DestroyBasin: Destroy basin of vectors
//
inline
void HDPbase::DestroyBasin()
{
    if( basin_ )
        delete basin_;
    basin_ = NULL;
}

// -------------------------------------------------------------------------
// DestroyMenu: Destroy menu of dishes
//
inline
void HDPbase::DestroyMenu()
{
    if( menu_ )
        delete menu_;
    menu_ = NULL;
}

// -------------------------------------------------------------------------
// DestroyChain: Destroy chain of restaurants
//
inline
void HDPbase::DestroyChain()
{
    if( chain_ )
        delete chain_;
    chain_ = NULL;
}

// -------------------------------------------------------------------------
// InitBasin: Initialize basin of vectors
//
inline
void HDPbase::InitBasin( int size )
{
    DestroyBasin();

    basin_ = new Basin( size );
    if( basin_ == NULL )
        throw MYRUNTIME_ERROR("HDPbase::InitBasin: Not enough memory.");
    basin_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitMenu: Initialize menu of dishes
//
inline
void HDPbase::InitMenu( int size )
{
    DestroyMenu();

    menu_ = new Menu( size );
    if( menu_ == NULL )
        throw MYRUNTIME_ERROR("HDPbase::InitMenu: Not enough memory.");
    menu_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitChain: Initialize chain of restaurants
//
inline
void HDPbase::InitChain( int size )
{
    DestroyChain();

    chain_ = new RntChain( size );
    if( chain_ == NULL )
        throw MYRUNTIME_ERROR("HDPbase::InitChain: Not enough memory.");
    chain_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// ReserveBasin: reserve basin size
//
inline
void HDPbase::ReserveMenu( int size )
{
    InitMenu( size );
}

// -------------------------------------------------------------------------
// ReserveBasin: reserve basin size
//
inline
void HDPbase::ReserveBasin( int size )
{
    InitBasin( size );
}

// =========================================================================
// AddToBasin: add normal multivariate variable to basin
//
inline
int HDPbase::AddToBasin( extspsl::Pslvector* nv )
{
    if( GetBasin() == NULL )
        InitBasin( GetDefBasinSize());
    return GetBasin()->NewValue( nv );
}

// =========================================================================
// SetScores: set HDP scores
//
inline
void HDPbase::SetScores( const HDPscores* value )
{
    scores_ = NULL;
    if( value == NULL )
        return;
    if( !GetMenu() || GetMenu()->GetSize() != value->GetCardinality())
        throw MYRUNTIME_ERROR("HDPbase::SetScores: Inconsistent score table.");
    scores_ = value;
}

#endif//__HDPbase__
