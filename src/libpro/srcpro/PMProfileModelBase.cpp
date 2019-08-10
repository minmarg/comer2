/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "liblib/alpha.h"
// #include "PMTransModel.h"
#include "PMProfileModelBase.h"

// // #include "logitnormal.h"
// // #include "CtxtCoefficients.h"

// // #include "Serializer.h"

namespace pmodel {

// -------------------------------------------------------------------------
// default constructor: 
//
PMProfileModelBase::PMProfileModelBase()
:   values_( NULL ),
    residues_( NULL ),
    length_( 0 ),
    allocated_( 0 )
{
    init();
}

// -------------------------------------------------------------------------
// destructor: deallocate memory
//
PMProfileModelBase::~PMProfileModelBase()
{
    destroy();
}

// -------------------------------------------------------------------------
// init: initialize members
//
void PMProfileModelBase::init()
{
    values_ = NULL;
    residues_ = NULL;
    length_ = 0;
    allocated_ = 0;
}

// -------------------------------------------------------------------------
// IsConsistentWith: verify the consistency of two objects of this class
//
bool PMProfileModelBase::IsConsistentWith( const PMProfileModelBase& one ) const
{
    if( GetSize() != one.GetSize())
        return false;

    for( int n = 0; n < GetSize(); n++ )
        if( residues_[n] != one.residues_[n] )
            return false;

    return true;
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset members
//
void PMProfileModelBase::destroy()
{
    if( values_ ) { free( values_ ); values_ = NULL; }
    if( residues_ ) { free( residues_ ); residues_ = NULL; }
    length_ = 0;
    allocated_ = 0;
}

// -------------------------------------------------------------------------
// Clear: erase information but leave the space allocated
//
void PMProfileModelBase::Clear()
{
    if( allocated_ ) {
        memset( values_, 0, sizeof(float) * PVDIM * allocated_ );
        memset( residues_, 0, sizeof(char) * allocated_ );
    }

    length_ = 0;
}

// -------------------------------------------------------------------------
// reallocate: (re)allocate memory
//
void PMProfileModelBase::reallocate( int size )
{
    float (*tmp_values)[PVDIM] = NULL;
    char* tmp_ress = NULL;

    if( size <= allocated_ )
        return;

    if( allocated_ <= 0 ) {
        tmp_values = ( float(*)[PVDIM] )malloc( sizeof(float) * PVDIM * size );
        tmp_ress = ( char* )malloc( sizeof(char) * size );

    } else {
        tmp_values = ( float(*)[PVDIM] )realloc( values_, sizeof(float) * PVDIM * size );
        tmp_ress = ( char* )realloc( residues_, sizeof(char) * size );
    }

    if( !tmp_values || !tmp_ress ) {
        if( tmp_values ) { free( tmp_values ); tmp_values = NULL; }
        if( tmp_ress ) { free( tmp_ress ); tmp_ress = NULL; }
        throw MYRUNTIME_ERROR( "PMProfileModelBase::reallocate: Not enough memory." );
    }

    values_ = tmp_values;
    residues_ = tmp_ress;

    // fill uninitialized memory with zeros
    memset( values_ + allocated_, 0, sizeof(float) * PVDIM * ( size - allocated_ ));
    memset( residues_ + allocated_, 0, sizeof(char) * ( size - allocated_ ));

    allocated_ = size;
}

// -------------------------------------------------------------------------
// Push: push a vector of values
//
void PMProfileModelBase::Push( const float values[PVDIM], char rs )
{
    if( allocated_ <= length_ ) {
        int newcap = TIMES2( allocated_ );
        if( newcap <= length_ )
            newcap = length_ + 1;
        reallocate( newcap );
    }
    PushAt( values, rs, GetSize());
}

// -------------------------------------------------------------------------
// PushAt: assign a vector of values to the given position
//
void PMProfileModelBase::PushAt( const float values[PVDIM], char rs, int pos )
{
    if( allocated_ <= pos )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::PushAt: Memory access error." );

    for( int i = 0; i < PVDIM; i++ )
        values_[pos][i] = values[i];

    residues_[pos] = rs;

    if( length_ <= pos )
        length_ = pos + 1;
}

// -------------------------------------------------------------------------
// Print: print information
//
void PMProfileModelBase::Print( const char* filename ) const
{
    int nores = PVDIM;//effective number of residues
    //
    myruntime_error mre;
    FILE* fp = stdout;
    char r;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::Print: Failed to open file for writing." );

    try {
        fprintf( fp, "%43c Target values%s", 32, NL );
        fprintf( fp, "%9c", 32 );

        for( r = 0; r < nores; r++ )
            fprintf( fp, "%4c", DehashCode(r));

        for( int p = 0; p < GetSize(); p++ ) {
            fprintf( fp, "%s%5d %c   ", NL, p+1, DehashCode( GetResidueAt(p)));

            for( r = 0; r < nores; r++ )
                fprintf( fp, "%4d", (int)rintf( 100.0f* GetValueAt(p,r)) );

        }
        fprintf( fp, "%s", NL );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( fp != stdout )
        fclose( fp );

    if( mre.isset())
        throw mre;
}

// // // -------------------------------------------------------------------------
// // // CheckIntegrity: checks integrity of the structure; teh amino acid symbols
// // //     must be valid
// // // -------------------------------------------------------------------------
// // 
// // void PMProfileModelBase::CheckIntegrity() const
// // {
// //     for( int n = 0; n < GetColumns(); n++ )
// //         if( NUMALPH <= ( *this )[n] )
// //             throw myruntime_error(
// //                 mystring( "PMProfileModelBase: Data corruption." ));
// // }

// // // -------------------------------------------------------------------------
// // // CheckForAllZeros: verifies wether there extists positions with all
// // //     values of zero; if so, the appropriate value is set to 1, indicating
// // //     100 percent amiino acid frequency in that position.
// // //     This is necessary for scoring matrix made of two profiles to take its
// // //     effect.
// // //     The method is applicable for frequency matrix only
// // // -------------------------------------------------------------------------
// // 
// // void PMProfileModelBase::CheckForAllZeros()
// // {
// //     const double    zval = 0.001;
// //     int             r;
// // 
// //     const int       efective_number = NUMAA;    // effective number of residues
// //     const double    eqpart = rint(( double )FREQUENCY_SUM / efective_number ) / FREQUENCY_SUM;
// // 
// //     static int      symB = HashAlphSymbol('B');
// //     static int      symZ = HashAlphSymbol('Z');
// //     static int      resN = HashAlphSymbol('N');
// //     static int      resD = HashAlphSymbol('D');
// //     static int      resQ = HashAlphSymbol('Q');
// //     static int      resE = HashAlphSymbol('E');
// // 
// //     static double   Bprob = ( LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD ));
// //     static double   Zprob = ( LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE ));
// // 
// //     for( int n = 0; n < GetColumns(); n++ ) {
// //         for( r = 0; r < NUMALPH; r++ )
// //             if( zval < ( *this )( n, r ))
// //                 break;
// //         if( r == NUMALPH ) {
// //             // If all values are zero
// // //             for( r = 0; r < NUMALPH; r++ )
// //                 //set weighted frequencies to the background frequencies
// // //                 ( *this )( n, r ) = LOSCORES.PROBABility( r );
// //             r = GetResidueAt( n );
// //             if( r == X ) {
// //                 //X is at this position; set all amino acids equally probable
// //                 for( int e = 0; e < efective_number; e++ )
// //                     ( *this )( n, e ) = LOSCORES.PROBABility( e );//eqpart;
// //             } else
// //             if( r == symB ) {
// //                 //B is at the position; make appropriate amino acids available and
// //                 //round a floating point number to precision of 2
// //                 ( *this )( n, resN ) = rint( LOSCORES.PROBABility( resN ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
// //                 ( *this )( n, resD ) = rint( LOSCORES.PROBABility( resD ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
// //             } else
// //             if( r == symZ ) {
// //                 //Z is at the position; make appropriate amino acids available and
// //                 //round a floating point number to precision of 2
// //                 ( *this )( n, resQ ) = rint( LOSCORES.PROBABility( resQ ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
// //                 ( *this )( n, resE ) = rint( LOSCORES.PROBABility( resE ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
// //             } else
// //                 //set corresponding amino acid fully conserved else
// //                 ( *this )( n, r ) = 1.0; //set 1 so that scores for one profile repeat the LOSCORES scores
// //         }
// //     }
// // }

// // // -------------------------------------------------------------------------
// // // Serialize: write the class data to file for reading them later
// // // -------------------------------------------------------------------------
// // 
// // void PMProfileModelBase::Serialize( Serializer& serializer ) const
// // {
// //     serializer.Write(( char* )&columns, sizeof( columns ), 1 );
// // 
// //     for( int n = 0; n < columns; n++ )
// //         serializer.Write(( char* )values[n], sizeof( double ), NUMALPH );
// // 
// //     if( columns > 0 )
// //         serializer.Write( aacids, sizeof( char ), columns );
// // }

// // // -------------------------------------------------------------------------
// // // Deserialize: read data into the class members
// // // -------------------------------------------------------------------------
// // 
// // void PMProfileModelBase::Deserialize( Serializer& serializer )
// // {
// //     //there's no need for destroying variables since memory allocated previously is reused
// // //     destroy();
// // 
// //     serializer.Read(( char* )&columns, sizeof( columns ), 1 );
// // 
// //     if( columns > MAXCOLUMNS )
// //         throw myruntime_error(
// //             mystring( "PMProfileModelBase: Number of positions read from file is larger than the maximum allowed." ));
// //         
// //     if( columns <= 0 )
// //         throw myruntime_error(
// //             mystring( "PMProfileModelBase: Invalid number of positions read from file." ));
// // 
// //     Reserve( columns ); // memory allocation
// // 
// //     for( int n = 0; n < columns; n++ )
// //         serializer.Read(( char* )values[n], sizeof( double ), NUMALPH );
// // 
// //     serializer.Read( aacids, sizeof( char ), columns );
// //     //
// //     CheckIntegrity();
// // }

}//namespace pmodel
