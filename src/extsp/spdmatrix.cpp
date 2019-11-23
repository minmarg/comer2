/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psl.h"
#include "pslerror.h"
#include "spdmatrix.h"

using namespace extspsl;

// -------------------------------------------------------------------------
// constructor: initialization
//
SPDmatrix::SPDmatrix( int nrc )
:   NSmatrix( nrc )
{
}

// -------------------------------------------------------------------------
// constructor: unmastered copy
//
SPDmatrix::SPDmatrix( const SPDmatrix& right )
:   NSmatrix()
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
SPDmatrix::SPDmatrix()
:   NSmatrix()
{
}

// -------------------------------------------------------------------------
// destructor:
//
SPDmatrix::~SPDmatrix()
{
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
SPDmatrix& SPDmatrix::operator=( const SPDmatrix& right )
{
    NSmatrix::operator=( right );
    return *this;
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
SPDmatrix& SPDmatrix::operator=( const NSmatrix& right )
{
    NSmatrix::operator=( right );
    return *this;
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
SPDmatrix& SPDmatrix::operator=( const Pslmatrix& right )
{
    NSmatrix::operator=( right );
    return *this;
}

// =========================================================================
// CholeskyDecompose: factorize the matrix into Cholesky Decomposition:
//  A = LL^T. Only diagonal and lower triangular elements are used. On output
//  diagonal and lower triangular part elements are changed with L, upper
//  triangular part is overwritten with L^T. If the matrix is not positive 
//  definite, an error is returned.
//
int SPDmatrix::CholeskyDecompose()
{
    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();

#ifdef SPDmatrixTESTPRINT
    const char* errdom = "Matrix must be positive definite";
#endif

    if( norows != nocols ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CholeskyDecompose: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }

    float sum;
    float A_00, L_00;
    float A_10, L_10;
    float A_11, L_11;
    float A_ii, A_jj, A_ij;
    float L_ii;//diagonal element
    int i, j;

    // do the first 2 rows explicitly; it's faster
    A_00 = GetValueAt( 0, 0 );

    if( A_00 <= 0.0f ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CholeskyDecompose: PSL_ERR_DOMAIN: %s.\n", errdom );
#endif
        return PSL_ERR_DOMAIN;
    }
    L_00 = sqrtf( A_00 );
    SetValueAt( 0, 0, L_00 );

    if( 1 < norows ) {
        A_10 = GetValueAt( 1, 0 );
        A_11 = GetValueAt( 1, 1 );

        L_10 = A_10 / L_00;
        L_ii = A_11 - L_10 * L_10;

        if( L_ii <= 0.0f ) {
#ifdef SPDmatrixTESTPRINT
            fprintf( stderr, "SPDmatrix::CholeskyDecompose: PSL_ERR_DOMAIN: %s.\n", errdom );
#endif
            return PSL_ERR_DOMAIN;
        }
        L_11 = sqrtf( L_ii );
        SetValueAt( 1, 0, L_10 );        
        SetValueAt( 1, 1, L_11 );
    }

    for( i = 2; i < norows; i++ ) {
        A_ii = GetValueAt( i, i );

        for( j = 0; j < i; j++ ) {
            sum = 0.0f;

            A_ij = GetValueAt( i, j );
            A_jj = GetValueAt( j, j );

            const Pslvector ci = RowVector( i );
            const Pslvector cj = RowVector( j );

            if( 0 < j ) {
                const Pslvector di = ci.SubVector( 0, j );
                const Pslvector dj = cj.SubVector( 0, j );
                if( dj.DotProduct( di, &sum) != PSL_SUCCESS ) {
#ifdef SPDmatrixTESTPRINT
                    fprintf( stderr, "SPDmatrix::CholeskyDecompose: PSL_ERR_DIM: "
                                     "Dot product failed.\n" );
#endif
                    return PSL_ERR_DIM;
                }
            }
            //A_jj always positive: A_jj=L_jj for j<i
            A_ij = ( A_ij - sum ) / A_jj;
            SetValueAt( i, j, A_ij );
        }

        const Pslvector cci = RowVector( i );
        const Pslvector ddi = cci.SubVector( 0, i );

        sum = ddi.Norm2();
        L_ii = A_ii - sum * sum;

        if( L_ii <= 0.0f ) {
#ifdef SPDmatrixTESTPRINT
            fprintf( stderr, "SPDmatrix::CholeskyDecompose: PSL_ERR_DOMAIN: %s.\n", errdom );
#endif
            return PSL_ERR_DOMAIN;
        }
        L_ii = sqrtf( L_ii );
        SetValueAt( i, i, L_ii );
    }

    //copy the transposed lower triangle to the upper triangle;
    //diagonal is shared
    for( i = 1; i < norows; i++ )
        for( j = 0; j < i; j++ )
            SetValueAt( j, i, GetValueAt( i, j ));

    return PSL_SUCCESS;
}

// =========================================================================
// CDedSolve: solve system of equations: Ax = b, where A = LL^T 
//  decomposed into two matrices L and L^T written in the lower and upper
//  triangular parts of A, respectively. On input xb is b vector, on output
//  it holds the solution.
//  First, Ly = b is solved for y; second, L^T x = y is solved for x.
//
int SPDmatrix::CDedSolve( Pslvector& xb ) const
{
    int status;
    if( GetNoRows() != GetNoCols()) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedSolve: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != xb.GetSize()) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedSolve: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }

    if(( status = SolveLower( xb )) != PSL_SUCCESS ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedSolve: %s\n", TranslatePSLError( status ));
#endif
        return status;
    }
    if(( status = SolveUpper( xb )) != PSL_SUCCESS ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedSolve: %s\n", TranslatePSLError( status ));
#endif
        return status;
    }
    return PSL_SUCCESS;
}

// =========================================================================
// CDedInvert: invert the matrix factorized by Cholesky decomposition.
//  On input this matrix is assumed to have already been decomposed: A=LL^T.
//  On output this matrix is overwritten with inverted one which is
//  A^{-1} = L^{-T} L^{-1}
//
int SPDmatrix::CDedInvert()
{
    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();

    if( norows != nocols ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedInvert: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }

    float A_jj, sum;
    int i, j, n;

    //invert the lower triangle (which is L) by solving
    //X L = I, where X is L^{-1}
    for( i = 0; i < norows; i++ ) {
        j = norows - i - 1;

        A_jj = GetValueAt( j, j );
        if( A_jj == 0.0f ) {
#ifdef SPDmatrixTESTPRINT
            fprintf( stderr, "SPDmatrix::CDedInvert: PSL_ERR_ILLEGAL.\n");
#endif
            return PSL_ERR_ILLEGAL;
        }
        A_jj = 1.0f / A_jj;
        SetValueAt( j, j, A_jj );

        if( 0 < i ) {//j < norows - 1
            n = i;//norows - j - 1;
            const Pslmatrix subm = SubMatrix( j+1, j+1, n, n );//right-most triangular sub matrix
            Pslvector   cl = SubColVector( j, j+1, n );//lower column vector

            subm.MultiplyLower( cl );
            cl.MultiplyBy( -A_jj );//divide by -A_jj
        }
    }

    //the lower triangle now contains L^{-1};
    //compute A^{-1} = L^{-T} L^{-1}
    //A^{-1}_{ij} is dot product of column i of L^{-1} and column j of L^{-1}
    for( i = 0; i < norows; i++ ) {
        for( j = i; j < nocols; j++ ) {
            const Pslvector v1 = SubColVector( i, j, norows - j );
            const Pslvector v2 = SubColVector( j, j, norows - j );

            v1.DotProduct( v2, &sum );
            //store in the upper triangle
            SetValueAt( i, j, sum );
        }
    }

    //copy the transposed upper triangle to the lower triangle
    for( j = 1; j < nocols; j++ )
        for( i = 0; i < j; i++ )
            SetValueAt( j, i, GetValueAt( i, j ));

    return PSL_SUCCESS;
}

// =========================================================================
// CDedDet: calculate determinant of the matrix factorized by 
//  Cholesky decomposition. Since the matrix is: A = L L^T, determinant
//  is equal to the squared product of the diagonal elements.
//
int SPDmatrix::CDedDet( float* det ) const
{
    if( !det )
        return 1;

    const int norows = GetNoRows();
    const int nocols = GetNoCols();
    float res = 1.0f;
    int i;

    if( norows != nocols ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedDet: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }
    for( i = 0; i < norows; i++ )
        res *= GetValueAt( i, i );

    *det = res * res;
    return PSL_SUCCESS;
}
// -------------------------------------------------------------------------
// CDedLogDet: calculate log determinant of the matrix factorized by 
//  Cholesky decomposition. Since the matrix is: A = L L^T, determinant
//  is equal to the squared product of the diagonal elements.
//
int SPDmatrix::CDedLogDet( float* ldet ) const
{
    if( !ldet )
        return 1;

    const int norows = GetNoRows();
    const int nocols = GetNoCols();
    float val, lres = 0.0f;
    int i;

    if( norows != nocols ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "SPDmatrix::CDedLogDet: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }
    for( i = 0; i < norows; i++ ) {
        val = fabsf( GetValueAt( i, i ));
        if( val == 0.0f ) {
#ifdef SPDmatrixTESTPRINT
            fprintf( stderr, "SPDmatrix::CDedLogDet: PSL_ERR_ILLEGAL.\n");
#endif
            return PSL_ERR_ILLEGAL;
        }
        val = logf( val );
        lres += val + val;
    }
    *ldet = lres;
    return PSL_SUCCESS;
}
