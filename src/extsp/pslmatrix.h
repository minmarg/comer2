/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_pslmatrix__
#define __extspsl_pslmatrix__

#include <stdio.h>
#include "pslerror.h"
#include "pslvector.h"

namespace extspsl {

class Pslmatrix
{
public:
    Pslmatrix( int nr, int nc );
    Pslmatrix( const Pslmatrix& );
    explicit Pslmatrix();
    virtual ~Pslmatrix();

    virtual Pslmatrix&  operator=( const Pslmatrix& );

    int         GetNoRows() const { return norows_; }
    int         GetNoCols() const { return nocols_; }

    void        SetIdentity();
    int         SetCentering( float nrm );

    //{{MATRIX OPERATORS
    int         Transpose();
    int         Transpose( Pslmatrix& tr ) const;
    int         Add( const Pslmatrix& mt );//addition
    int         Add( float );
    int         AddToDiag( float value );
    int         AddIdentity() { return AddToDiag( 1.0f ); }
    int         Sub( const Pslmatrix& mt );//subtraction
    int         SubFrom( const Pslmatrix& mt );//subtraction: mt-this
    int         SumColumns( Pslvector& vv ) const;//sum of all columns
    int         Scale( float );
    int         Mul( const Pslmatrix& m1, const Pslmatrix& m2 );//multiplication
    int         Mul( const Pslvector& v1, const Pslvector& v2 );//multiplication
    int         KronProduct( const Pslmatrix& m1, const Pslmatrix& m2 );//Kronecker product
    //}}

    //{{SOLVE methods
    int         SolveLower( Pslvector& xb, bool unitdiag = false ) const;
    int         SolveUpper( Pslvector& xb, bool unitdiag = false ) const;

    int         MultiplyLower( Pslvector& bx ) const;
    int         MultiplyUpper( Pslvector& bx ) const;
    //}}

    //{{SWAP operations
    int         SwapRows( int i1, int i2 );
    int         SwapCols( int j1, int j2 );
    //}}

    //{{SUB- vector, matrix operations
    const Pslmatrix SubMatrix( int i, int j, int nr, int nc ) const;

    const Pslvector Stack() const;
    const Pslvector RowVector( int i ) const;
    const Pslvector SubRowVector( int i, int offset, int n ) const;

    const Pslvector ColVector( int j ) const;
    const Pslvector SubColVector( int j, int offset, int n ) const;

    const Pslvector DiagVector() const;
    const Pslvector SubDiagVector( int k ) const;
    const Pslvector SuperDiagVector( int k ) const;
    //}}


    float       GetValueAt( int n, int m ) const;   //get value at the given position
    void        SetValueAt( int n, int m, float value );//set value at the given position
    void        AddValueAt( int n, int m, float value );//add value at the given position
    void        MulValueAt( int n, int m, float value );//multiply by value at the given position
    void        ChangeSignAt( int n, int m );//change of sign

    void        Copy( const Pslmatrix& mtx );       //copy elements
    void        Zero();                             //assign all elements to zero
    void        Clear();                            //clear all elements
    void        Print( FILE*, const char* format = NULL ) const;//print matrix

    void        Reserve( int nr, int nc );

protected:
    void        Realloc( int newcap );
    int         GetCapacity() const { return capacity_; }

    void        DestroyValues();

    const float* GetValues() const { return values_; }

    bool        GetMaster() const { return master_; }
    void        SetMaster() { master_ = true; }
    void        SetNotMaster() { master_ = false; }

    int         GetPhysRowDim() const { return rowdim_; }
    void        SetPhysRowDim( int value ) { rowdim_ = value; }

    void        SetDims( int nr, int nc ) { if( !GetMaster()) return; norows_ = nr; nocols_ = nc; }

protected:
    float*      values_;        //single-precision values of matrix
    int         norows_;        //number of rows
    int         nocols_;        //number of columns
    int         capacity_;      //current capacity of data
    int         rowdim_;        //physical row dimension (occupied in memory)
    bool        master_;        //flag of mastering an object
};

// -------------------------------------------------------------------------
// DestroyValues: destroy matrix
//
inline
void Pslmatrix::DestroyValues()
{ 
    if( master_ ) 
        if( values_ ) 
            free( values_ ); 
    values_ = NULL;
    norows_ = 0;
    nocols_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// Reserve: Reserve space for matrix
//
inline
void Pslmatrix::Reserve( int nr, int nc )
{
    if( !GetMaster())
        PRIVERROR("Pslmatrix::Reserve: Only master is allowed to allocate memory.");
    if( 0 < nr && 0 < nc  ) {
        SetPhysRowDim( nc );
        Realloc( nr * nc );
        SetDims( nr, nc );
    }
}

// -------------------------------------------------------------------------
// GetValueAt: get value at the position
//
inline
float Pslmatrix::GetValueAt( int n, int m ) const
{
#ifdef __DEBUG__
    if( !values_ || norows_ <= n || nocols_ <= m || rowdim_ < 1 )
        PRIVERROR("Pslmatrix::GetValueAt: Memory access error.");
#endif
    return values_[ n*rowdim_ + m ];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at the position
//
inline
void Pslmatrix::SetValueAt( int n, int m, float value )
{
#ifdef __DEBUG__
    if( !values_ || norows_ <= n || nocols_ <= m || rowdim_ < 1 )
        PRIVERROR("Pslmatrix::SetValueAt: Memory access error.");
#endif
    values_[ n*rowdim_ + m ] = value;
}

// -------------------------------------------------------------------------
// AddValueAt: add value at the position
//
inline
void Pslmatrix::AddValueAt( int n, int m, float value )
{
#ifdef __DEBUG__
    if( !values_ || norows_ <= n || nocols_ <= m || rowdim_ < 1 )
        PRIVERROR("Pslmatrix::AddValueAt: Memory access error.");
#endif
    values_[ n*rowdim_ + m ] += value;
}

// -------------------------------------------------------------------------
// MulValueAt: multiply by value at the position
//
inline
void Pslmatrix::MulValueAt( int n, int m, float value )
{
#ifdef __DEBUG__
    if( !values_ || norows_ <= n || nocols_ <= m || rowdim_ < 1 )
        PRIVERROR("Pslmatrix::MulValueAt: Memory access error.");
#endif
    values_[ n*rowdim_ + m ] *= value;
}

// -------------------------------------------------------------------------
// ChangeSignAt: change sign of the value at the position
//
inline
void Pslmatrix::ChangeSignAt( int n, int m )
{
#ifdef __DEBUG__
    if( !values_ || norows_ <= n || nocols_ <= m || rowdim_ < 1 )
        PRIVERROR("Pslmatrix::ChangeSignAt: Memory access error.");
#endif
    values_[ n*rowdim_ + m ] = -values_[ n*rowdim_ + m ];
}

}//namespace extspsl

#endif//__extspsl_pslmatrix__
