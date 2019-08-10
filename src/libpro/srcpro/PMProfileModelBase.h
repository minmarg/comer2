/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __PMProfileModelBase_h__
#define __PMProfileModelBase_h__

#include "liblib/mybase.h"
#include "liblib/alpha.h"
#include "PMTransModel.h"

// // class Serializer;
// // class GapScheme;

// profile model namespace
namespace pmodel {

// _________________________________________________________________________
// Class PMProfileModelBase
// base class for profile model
//
class PMProfileModelBase
{
public:
    enum {
        //dimensions
        PVDIM = NUMALPH
    };
public:
    PMProfileModelBase();
    virtual ~PMProfileModelBase();

    int GetSize/*GetColumns*/() const { return length_; }

// //     double          operator()( int m, int a ) const    { return GetValueAt( m, a ); }
// //     char            operator[]( int m ) const           { return GetResidueAt( m );  }

// //     double&         operator()( int m, int a );         //modification of values
// //     char&           operator[]( int m );                //modification of aacids

    char GetResidueAt( int m ) const;
    float GetValueAt( int m, int a ) const;//get value at specified position for specified residue type
    const float ( *GetVectorAt( int m ) const )[PVDIM];
    const float ( *GetVector() const )[PVDIM] { return values_; }

    void            Push( const float values[PVDIM], char );//push a vector of values
    virtual void    PushAt( const float values[PVDIM], char, int pos );//push vector at a position

// //     virtual void    Serialize( Serializer& ) const;
// //     virtual void    Deserialize( Serializer& );

    // check consistency
    bool        IsConsistentWith/*IsCompatible*/( const PMProfileModelBase& pm ) const;
    bool        IsConsistentWith/*IsCompatible*/( const PMTransModel& tm ) const;

    virtual void    Print/*OutputMatrix*/( const char* = NULL ) const;

    void        Reserve( int size ) { reallocate( size ); }

// //     void            CheckForAllZeros();//verify wether extists positions with all values of zero

    virtual void    Clear();//clear all information

    const char*     GetResidues() const { return residues_; }

protected:
    float           (*GetVector())[PVDIM] { return values_; }

    virtual void    destroy();//memory deallocation and reset of values
    virtual void    reallocate( int size );//memory allocation
    virtual void    init();//initialization of the members

// //     void            SetColumns( int col )   { columns = col; }
// //     void            CheckIntegrity() const;                 //check integrity of the structure

protected:
    float       (*values_)[PVDIM];//target values: mixed frequencies, target probability, log-odds, etc.
    char*       residues_/*aacids*/;//profile's representative residue sequence
    int         length_/*columns*/; //profile length
    int         allocated_;         //positions allocated

};

// /////////////////////////////////////////////////////////////////////////
// INLINES -----------------------------------------------------------------
//
// GetValue: get value for the specified residue at the given position 
//
inline
float PMProfileModelBase::GetValueAt( int m, int a ) const
{
#ifdef __DEBUG__
    if( length_ <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::GetValueAt: Memory access error." );

    if( PVDIM <= a || a < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::GetValueAt: Memory access error." );
#endif
    return values_[m][a];
}

// -------------------------------------------------------------------------
// GetVectorAt: get the vector of values at the given position
//
inline
const float ( *PMProfileModelBase::GetVectorAt( int m ) const )[PVDIM]
{
#ifdef __DEBUG__
    if( !values_ || length_ <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::GetVectorAt: Memory access error." );
#endif
    return values_ + m;
}

// // // -------------------------------------------------------------------------
// // // operator(): used to modify score value at the specified position and for
// // //     the specified amino acid
// // // -------------------------------------------------------------------------
// // 
// // inline
// // double& PMProfileModelBase::operator()( int m, int a )
// // {
// // #ifdef __DEBUG__
// //     if( columns <= m || m < 0 )
// //         throw myruntime_error(
// //             mystring( "PMProfileModelBase: Memory access error." ));
// // 
// //     if( NUMALPH <= a || a < 0 )
// //         throw myruntime_error(
// //             mystring( "PMProfileModelBase: Memory access error." ));
// // #endif
// //     return values[m][a];
// // }

// -------------------------------------------------------------------------
// GetResidueAt: get residue at the specified position
//
inline
char PMProfileModelBase::GetResidueAt( int m ) const
{
#ifdef __DEBUG__
    if( length_ <= m || m < 0 )
        throw MYRUNTIME_ERROR( "PMProfileModelBase::GetResidueAt: Memory access error." );
#endif
    return residues_[m];
}

// // // -------------------------------------------------------------------------
// // // operator[]: returns amino acid to be modified
// // // -------------------------------------------------------------------------
// // 
// // inline
// // char& PMProfileModelBase::operator[]( int m )
// // {
// // #ifdef __DEBUG__
// //     if( columns <= m || m < 0 )
// //         throw myruntime_error(
// //             mystring( "PMProfileModelBase: Memory access error." ));
// // #endif
// //     return aacids[m];
// // }

// -------------------------------------------------------------------------
// IsConsistentWith: verify the consistency of this and PMTransModel class 
// object
//
inline
bool PMProfileModelBase::IsConsistentWith( const PMTransModel& goc ) const
{
    if( GetSize() != goc.GetSize())
        return false;

    return true;
}

}//namespace pmodel

#endif//__PMProfileModelBase_h__
