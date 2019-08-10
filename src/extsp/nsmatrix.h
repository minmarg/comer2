/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __extspsl_nsmatrix__
#define __extspsl_nsmatrix__

#include <stdlib.h>
#include <string.h>
#include "pslerror.h"
#include "pslmatrix.h"

#ifndef NSmatrixTESTPRINT
#define NSmatrixTESTPRINT
#endif

namespace extspsl {

// -------------------------------------------------------------------------
// Invertibale square matrix
//
class NSmatrix: public Pslmatrix
{
public:
    NSmatrix( int nrc );
    NSmatrix( const NSmatrix& );
    virtual ~NSmatrix();

    virtual NSmatrix&   operator=( const NSmatrix& );
    virtual NSmatrix&   operator=( const Pslmatrix& );

    int     LUDecompose();
    int     LUDedSolve( Pslvector& xb ) const;
    int     LUDedInvert( Pslmatrix& inv ) const;
    int     LUDedDet( float* ) const;
    int     LUDedLogDet( float* ) const;

protected:
    explicit NSmatrix();

    void    InitPermuts( int size );
    void    CopyPermuts( const NSmatrix& );
    void    DestroyPermuts();
    int     GetPermutsAt( int pos ) const;
    int     SwapPermPoss( int i1, int i2 );
    int     PermuteVector( Pslvector& v ) const;
    int     GetSizeOfRowPerms() const {return szoperm_; }

private:
    bool    LUDedSingular() const;

private:
    int*    rowperm_;//row permutation vector containing indices of permuted positions
    int     szoperm_;//length of rowperm_
    int     moneton_;//(-1)^n; positive if number of interchanges is even; negative orherwise
};


// -------------------------------------------------------------------------
// InitPermuts: initialize vector of row permutations
//
inline
void NSmatrix::InitPermuts( int size )
{
    if( !GetMaster()) {
        PRIVERROR("NSmatrix::InitPermuts: Slave cannot initialize vector of permutations.");
        return;
    }
    if( size < 1 )
        return;
    if( size != szoperm_ ) {
        DestroyPermuts();
        rowperm_ = ( int* )malloc( size * sizeof( int ));
        if( rowperm_ == NULL ) {
            PRIVERROR("NSmatrix::InitPermuts: Not enough memory.");
            return;
        }
        szoperm_ = size;
    }
    for( int n = 0; n < szoperm_; n++ )
        rowperm_[n] = n;
}

// -------------------------------------------------------------------------
// CopyPermuts: copy vector of row permutations
//
inline
void NSmatrix::CopyPermuts( const NSmatrix& right )
{
    DestroyPermuts();
    if( GetMaster()) {
        rowperm_ = NULL;
        szoperm_ = right.szoperm_;
        if( szoperm_ < 1 || right.rowperm_ == NULL )
            return;
        rowperm_ = ( int* )malloc( szoperm_ * sizeof( int ));
        if( rowperm_ == NULL ) {
            PRIVERROR("NSmatrix::CopyPermuts: Not enough memory.");
            return;
        }
        memcpy( rowperm_, right.rowperm_, szoperm_ * sizeof( int ));
        return;
    }
    else {
        rowperm_ = right.rowperm_;
        szoperm_ = right.szoperm_;
    }
}

// -------------------------------------------------------------------------
// GetPermutsAt: get value of row permutations
//
inline
int NSmatrix::GetPermutsAt( int pos ) const
{
    if( rowperm_ == NULL || GetSizeOfRowPerms() <= pos ) {
        PRIVERROR("NSmatrix::GetPermutsAt: Memory access error.");
        return pos;
    }
    return rowperm_[pos];
}

// -------------------------------------------------------------------------
// SwapPermPoss: swap position values in the permutation vector
//
inline
int NSmatrix::SwapPermPoss( int i1, int i2 )
{
    if( !GetMaster()) {
        PRIVERROR("NSmatrix::SwapPermPoss: Slave cannot swap permutation vector values.");
        return PSL_ERR_ILLEGAL;
    }
    int size = GetSizeOfRowPerms();
    int tmp;
    if( size <= i1 || size <= i2 || rowperm_ == NULL ) {
        PRIVERROR("NSmatrix::SwapPermPoss: Memory access error.");
        return PSL_ERR_ADDRESS;
    }
    if( i1 == i2 )
        return PSL_SUCCESS;
    tmp = rowperm_[i1];
    rowperm_[i1] = rowperm_[i2];
    rowperm_[i2] = tmp;
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// PermuteVector: permute vector using information of row permutations
//
inline
int NSmatrix::PermuteVector( Pslvector& v ) const
{
    int size = GetSizeOfRowPerms();
    int n, k, pk;
    float tmp;

    if( rowperm_ == NULL ) {
        PRIVERROR("NSmatrix::PermuteVector: Memory access error.");
        return PSL_ERR_ADDRESS;
    }
    if( v.GetSize() != GetSizeOfRowPerms()) {
        PRIVERROR("NSmatrix::PermuteVector: Inconsistent dimensions.");
        return PSL_ERR_DIM;
    }
    for( n =  0; n < size; n++ ) {
        for( k = GetPermutsAt( n ); n < k; k = GetPermutsAt( k ));
        if( k < n )
            continue;
        //k == n
        pk = GetPermutsAt( k );
        if( pk == n )
            continue;

        tmp = v.GetValueAt( n );
        while( pk != n ) {
            v.SetValueAt( k, v.GetValueAt( pk ));
            k = pk;
            pk = GetPermutsAt( k );
        }
        v.SetValueAt( k, tmp );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// DestroyPermuts: destroy vector of row permutations
//
inline
void NSmatrix::DestroyPermuts()
{
    if( !GetMaster())
        return;
    if( rowperm_ ) {
        free( rowperm_ );
        rowperm_ = NULL;
    }
    szoperm_ = 0;
}

}//namespace extspsl

#endif//__extspsl_nsmatrix__
