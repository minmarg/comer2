/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __menu__
#define __menu__

#include <stdio.h>
#include <stdlib.h>
#include "hdplib.h"
#include "liblib/mybase.h"
#include "extsp/ivector.h"
#include "extsp/spdmatrix.h"
#include "dish.h"

// -------------------------------------------------------------------------
// class Menu: reusable vector of dishes (global clusters)
//
class Menu
{
public:
    Menu( int size );
    Menu( const Menu& );
    explicit Menu();
    ~Menu();

    Menu&  operator=( const Menu& );

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    Dish*       GetDishAt( int n ) const;
    void        SetDishAt( int n, Dish* value );
    int         NewDish( Dish* value );//push dish
    void        RemDishAt( int loc, Dish* value );//remove dish

    int         NewTmpDish( Dish* value );//push temp. dish
    void        RemTmpDishAt( int loc, Dish* value );//remove temp. dish

    float       GetProbAt( int n ) const;
    void        SetProbAt( int n, float value );

    float       GetPriorProbAt( int n ) const;
    void        SetPriorProbAt( int n, float value );

    float       GetPriorProbNewDish() const { return npprob_; }
    void        SetPriorProbNewDish( float value ) { npprob_ = value; }

    extspsl::Pslvector* GetMuVectorAt( int n ) const;
    void        SetMuVectorAt( int n, extspsl::Pslvector* value );

    extspsl::SPDmatrix* GetSMatrixAt( int n ) const;
    void        SetSMatrixAt( int n, extspsl::SPDmatrix* value );

#ifdef PROBMTXINVSM
    extspsl::SPDmatrix* GetInvSMatrixAt( int n ) const;
    void        SetInvSMatrixAt( int n, extspsl::SPDmatrix* value );
#endif
    float       GetLDetSMAt( int n ) const;
    void        SetLDetSMAt( int n, float value );

    //{{Prior variables
    float       GetLogPriorProbFact() const { return lmvrprobfact_; }
    void        SetLogPriorProbFact( float value ) { lmvrprobfact_ = value; }

    float       GetLogMtxPriorProbFact() const { return lmtxprobfact_; }
    void        SetLogMtxPriorProbFact( float value ) { lmtxprobfact_ = value; }

    float       GetS0ScaleFac() const { return S0scalefac_; }
    void        SetS0ScaleFac( float value ) { S0scalefac_ = value; }

    extspsl::SPDmatrix* GetS0() const { return S0_; }
    void        SetS0( extspsl::SPDmatrix* value ) { if( S0_ && GetDestroyPriors()) delete S0_; S0_ = value; }

    extspsl::SPDmatrix* GetInvS0() const { return S0inv_; }
    void        SetInvS0( extspsl::SPDmatrix* value ) { if( S0inv_ && GetDestroyPriors()) delete S0inv_; S0inv_ = value; }

    float       GetLDetS0() const { return ldetS0_; }
    void        SetLDetS0( float value ) { ldetS0_ = value; }

    extspsl::Pslvector* GetMu0() const { return m0_; }
    void        SetMu0( extspsl::Pslvector* value ) { if( m0_ && GetDestroyPriors()) delete m0_; m0_ = value; }

    float       GetKappa0_pp_a() const { return a_k0_; }
    void        SetKappa0_pp_a( float value ) { a_k0_ = value; }

    float       GetKappa0_pp_b() const { return b_k0_; }
    void        SetKappa0_pp_b( float value ) { b_k0_ = value; }

    float       GetKappa0() const { return k0_; }
    void        SetKappa0( float value ) { k0_ = value; }

    float       GetNu0() const { return v0_; }
    void        SetNu0( float value ) { v0_ = value; }

    int         GetDim() const { return dim_; }
    void        SetDim( int value ) { dim_ = value; }

    int         GetCtx() const { return ctx_; }
    void        SetCtx( int value ) { ctx_ = value; }
    //}}

    int         GetNoTables() const { return notbls_; }
    void        SetNoTables( int value ) { notbls_ = value; }
    void        CalcNoTables();


    void        SoftCopy( const Menu& vector ); //copy addresses
    void        HardCopy( const Menu& vector ); //copy elements
    void        SoftCopyPriors( const Menu& vector );//copy addresses
    void        HardCopyPriors( const Menu& vector );//copy elements
    void        Clear();                        //clear all elements

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    bool        GetDestroyPriors() const { return destroypriors_; }
    void        SetDestroyPriors( bool value ) { destroypriors_ = value; }

    void        Reserve( int size );
    void        ReserveDishes( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyDishes();
    int         GetCapacity() const { return capacity_; }

    Dish**      GetDishes() const { return values_; }

    void        ReallocVacans( int newcap );
    void        DestroyVacans();
    int         GetCapVacans() const { return capvac_; }

    const int*  GetVacancies() const { return vacancies_; }
    int         GetVacantAt( int n ) const;
    void        SetVacantAt( int n, int index );
    void        PushVacant( int index );

    int         GetNoVacans() const { return novacs_; }
    void        SetNoVacans( int value ) { novacs_ = value; }

    void        SetSize( int value ) { length_ = value; }
    void        SetActualSize( int value ) { actlen_ = value; }

private:
    Dish**      values_;        //dishes
    float*      probs_;         //probabilities of dishes
    float*      pprobs_;        //prior probabilities of dishes
    float       npprob_;        //prior probability of new dish
    extspsl::Pslvector** mu_;   //mean vectors of NIW dist.
    extspsl::SPDmatrix** lambda_;//scale matrices of NIW dist.
#ifdef PROBMTXINVSM
    extspsl::SPDmatrix** lambdainv_//inverse scale matrices of NIW dist.
#endif
    float*      ldetlambda_;    //log determinants of scale matrices

    float       lmvrprobfact_;  //log of constant factor of probability of new vector
    float       lmtxprobfact_;  //log of constant factor of probability of new matrix (set of vecs)
    float       S0scalefac_;    //S0 scale factor
    extspsl::SPDmatrix* S0_;    //prior scale matrix of NIW distribution
    extspsl::SPDmatrix* S0inv_; //inverse of prior scale matrix of NIW distribution
    float       ldetS0_;        //log determinant of prior scale matrix
    extspsl::Pslvector* m0_;    //prior mean of NIW distribution
    float       a_k0_, b_k0_;   //prior parameters to kappa0
    float       k0_, v0_;       //no. measurements (kappa) and deg. of freedom (nu) of prior dist.
    int         dim_;           //number of dimensions of data
    int         ctx_;           //context length
    int         notbls_;        //number of tables

    int         actlen_;        //actual number of values
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the vector
    bool        destroy_;       //flag of destroying object 
    bool        destroypriors_; //flag of destroying prior parameters
    int*        vacancies_;     //vector of indices of vacancies
    int         novacs_;        //number of vacant elements
    int         capvac_;        //capcity of vacancies_
};


// -------------------------------------------------------------------------
// Reserve: Reserve space for vectors
//
inline
void Menu::Reserve( int size )
{
    ReserveDishes( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveDishes: Reserve space for vector of dishes and accompanying data
//
inline
void Menu::ReserveDishes( int size )
{
    if( 0 < size ) {
        Realloc( size );
//         SetSize( size );
    }
}

// -------------------------------------------------------------------------
// ReserveVacans: Reserve space for vector of vacancies
//
inline
void Menu::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void Menu::Destroy()
{
    DestroyDishes();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyDishes: Destroy vector of dishes and related data
//
inline
void Menu::DestroyDishes()
{
    int n;
    if( GetDestroy()) {
        if( values_ ) {
            for( n = 0; n < length_; n++ )
                if( values_[n] )
                    delete values_[n];
            free( values_ );
        }
        if( mu_ ) {
            for( n = 0; n < length_; n++ )
                if( mu_[n] )
                    delete mu_[n];
            free( mu_ );
        }
        if( lambda_ ) {
            for( n = 0; n < length_; n++ )
                if( lambda_[n] )
                    delete lambda_[n];
            free( lambda_ );
        }
#ifdef PROBMTXINVSM
        if( lambdainv_ ) {
            for( n = 0; n < length_; n++ )
                if( lambdainv_[n] )
                    delete lambdainv_[n];
            free( lambdainv_ );
        }
#endif
        if( probs_ )
            free( probs_ );
        if( pprobs_ )
            free( pprobs_ );
        if( ldetlambda_ )
            free( ldetlambda_ );
    }
    values_ = NULL;
    probs_ = NULL;
    pprobs_ = NULL;
    mu_ = NULL;
    lambda_ = NULL;
#ifdef PROBMTXINVSM
    lambdainv_ = NULL;
#endif
    ldetlambda_ = NULL;
    actlen_ = 0;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// DestroyVacans: Destroy vector of vacancies
//
inline
void Menu::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}

// -------------------------------------------------------------------------
// CalcNoTables: calculate number of tables over all dishes
//
inline
void Menu::CalcNoTables()
{
    Dish*   dish;
    int     ntbls = 0;
    int     n;
    for( n = 0; n < GetSize(); n++ ) {
        dish = GetDishAt( n );
        if( dish == NULL )
            continue;
        ntbls += dish->GetNoTables();
    }
    SetNoTables( ntbls );
}

// -------------------------------------------------------------------------
// GetDishAt: get dish at position loc
//
inline
Dish* Menu::GetDishAt( int loc ) const
{
// #ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetDishAt: Memory access error.");
// #endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetDishAt: set dish at position loc
//
inline
void Menu::SetDishAt( int loc, Dish* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetDishAt: Memory access error.");
#endif
    values_[loc] = value;
}

// -------------------------------------------------------------------------
// NewTmpDish: push temporary dish to the end of queue
//
inline
int Menu::NewTmpDish( Dish* value )
{
    int loc;
    if( capacity_ <= ( loc = length_ ))
        Realloc( TIMES2( capacity_ + 1 ));
    values_[loc] = value;
    length_++;
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemTmpDishAt: remove temporary dish at position loc
//
inline
void Menu::RemTmpDishAt( int loc, Dish* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::RemTmpDishAt: Memory access error.");
#endif
    if( loc+1 != length_ )
        throw MYRUNTIME_ERROR("Menu::RemTmpDishAt: Invalid index.");
    if( values_[loc] != value )
        throw MYRUNTIME_ERROR("Menu::RemTmpDishAt: Memory access error.");
    if( GetDestroy()) {
        if( values_[loc]) delete values_[loc];
        if( mu_[loc]) delete mu_[loc];
        if( lambda_[loc]) delete lambda_[loc];
#ifdef PROBMTXINVSM
        if( lambdainv_[loc]) delete lambdainv_[loc];
#endif
    }
    values_[loc] = NULL;
    probs_[loc] = 0.0f;
    pprobs_[loc] = 0.0f;
    mu_[loc] = NULL;
    lambda_[loc] = NULL;
#ifdef PROBMTXINVSM
    lambdainv_[loc] = NULL;
#endif
    ldetlambda_[loc] = 0.0f;
    length_--;
    actlen_--;
}

// -------------------------------------------------------------------------
// NewDish: push dish
//
inline
int Menu::NewDish( Dish* value )
{
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw MYRUNTIME_ERROR("Menu::NewDish: Memory access error.");
        values_[loc] = value;
        novacs_--;
        actlen_++;
        return loc;
    }
    if( capacity_ <= ( loc = length_ ))
        Realloc( TIMES2( capacity_ + 1 ));
    values_[loc] = value;
    length_++;
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemDishAt: remove dish at position loc
//
inline
void Menu::RemDishAt( int loc, Dish* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::RemDishAt: Memory access error.");
#endif
    if( values_[loc] != value )
        throw MYRUNTIME_ERROR("Menu::RemDishAt: Memory access error.");
    if( GetDestroy()) {
        if( values_[loc]) delete values_[loc];
        if( mu_[loc]) delete mu_[loc];
        if( lambda_[loc]) delete lambda_[loc];
#ifdef PROBMTXINVSM
        if( lambdainv_[loc]) delete lambdainv_[loc];
#endif
    }
    values_[loc] = NULL;
    probs_[loc] = 0.0f;
    pprobs_[loc] = 0.0f;
    mu_[loc] = NULL;
    lambda_[loc] = NULL;
#ifdef PROBMTXINVSM
    lambdainv_[loc] = NULL;
#endif
    ldetlambda_[loc] = 0.0f;
    PushVacant( loc );
    actlen_--;
}



// -------------------------------------------------------------------------
// GetProbAt: get probability of dish
//
inline
float Menu::GetProbAt( int loc ) const
{
#ifdef __DEBUG__
    if( !probs_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetProbAt: Memory access error.");
#endif
    return probs_[loc];
}

// -------------------------------------------------------------------------
// SetProbAt: set probability of dish
//
inline
void Menu::SetProbAt( int loc, float value )
{
#ifdef __DEBUG__
    if( !probs_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetProbAt: Memory access error.");
#endif
    probs_[loc] = value;
}

// -------------------------------------------------------------------------
// GetPriorProbAt: get prior probability of dish
//
inline
float Menu::GetPriorProbAt( int loc ) const
{
#ifdef __DEBUG__
    if( !pprobs_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetPriorProbAt: Memory access error.");
#endif
    return pprobs_[loc];
}

// -------------------------------------------------------------------------
// SetPriorProbAt: set prior probability of dish
//
inline
void Menu::SetPriorProbAt( int loc, float value )
{
#ifdef __DEBUG__
    if( !pprobs_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetPriorProbAt: Memory access error.");
#endif
    pprobs_[loc] = value;
}

// -------------------------------------------------------------------------
// GetMuVectorAt: get mean vector associated with the dish at loc
//
inline
extspsl::Pslvector* Menu::GetMuVectorAt( int loc ) const
{
#ifdef __DEBUG__
    if( !mu_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetMuVectorAt: Memory access error.");
#endif
    return mu_[loc];
}

// -------------------------------------------------------------------------
// SetMuVectorAt: set mean vector for the dish at loc
//
inline
void Menu::SetMuVectorAt( int loc, extspsl::Pslvector* value )
{
#ifdef __DEBUG__
    if( !mu_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetMuVectorAt: Memory access error.");
#endif
    if( mu_[loc])
        delete mu_[loc];
    mu_[loc] = value;
}

// -------------------------------------------------------------------------
// GetSMatrixAt: get scale matrix of the dish at loc
//
inline
extspsl::SPDmatrix* Menu::GetSMatrixAt( int loc ) const
{
#ifdef __DEBUG__
    if( !lambda_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetSMatrixAt: Memory access error.");
#endif
    return lambda_[loc];
}

// -------------------------------------------------------------------------
// SetSMatrixAt: set scale matrix for the dish at loc
//
inline
void Menu::SetSMatrixAt( int loc, extspsl::SPDmatrix* value )
{
#ifdef __DEBUG__
    if( !lambda_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetSMatrixAt: Memory access error.");
#endif
    if( lambda_[loc])
        delete lambda_[loc];
    lambda_[loc] = value;
}

#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// GetInvSMatrixAt: get inverse scale matrix of the dish at loc
//
inline
extspsl::SPDmatrix* Menu::GetInvSMatrixAt( int loc ) const
{
#ifdef __DEBUG__
    if( !lambdainv_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetInvSMatrixAt: Memory access error.");
#endif
    return lambdainv_[loc];
}

// -------------------------------------------------------------------------
// SetInvSMatrixAt: set inverse scale matrix for the dish at loc
//
inline
void Menu::SetInvSMatrixAt( int loc, extspsl::SPDmatrix* value )
{
#ifdef __DEBUG__
    if( !lambdainv_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetInvSMatrixAt: Memory access error.");
#endif
    if( lambdainv_[loc])
        delete lambdainv_[loc];
    lambdainv_[loc] = value;
}

#endif//PROBMTXINVSM

// -------------------------------------------------------------------------
// GetLDetSMAt: get log determinant of scale matrix of the dish at loc
//
inline
float Menu::GetLDetSMAt( int loc ) const
{
#ifdef __DEBUG__
    if( !ldetlambda_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetLDetSMAt: Memory access error.");
#endif
    return ldetlambda_[loc];
}

// -------------------------------------------------------------------------
// SetLDetSMAt: set log determinant of scale matrix for the dish at loc
//
inline
void Menu::SetLDetSMAt( int loc, float value )
{
#ifdef __DEBUG__
    if( !ldetlambda_ || loc < 0 || length_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetLDetSMAt: Memory access error.");
#endif
    ldetlambda_[loc] = value;
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int Menu::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Menu::GetVacantAt: Memory access error.");
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void Menu::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw MYRUNTIME_ERROR("Menu::SetVacantAt: Memory access error.");
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index in the vector of vacancies
//
inline
void Menu::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__menu__
