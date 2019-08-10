/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __BMSequence_h__
#define __BMSequence_h__

#include "liblib/mybase.h"

// #define UNUSED  0
// #define USED    1
// #define USED_IN_EXTENT  3

// _________________________________________________________________________
// Class BMSequence
// characterization of biomolecular sequence in the MSA
//
class BMSequence {
    enum TBMSUsage {
        BMS_UNUSED = 0,
        BMS_USED = 1,
        BMS_USED_IN_EXTENT = 3
    };
public:
    BMSequence( size_t reservation );
    virtual ~BMSequence();

//     unsigned char&      operator[]( int n ) { return GetResidueAt(n); };
    unsigned char   operator[]( int n ) const { return GetResidueAt(n); };

    size_t          GetSize() const { return length_; }
    size_t          GetCapacity() const { return capacity_; }

    unsigned int    GetCluster() const { return cluster_; }
    void            SetCluster( unsigned int value ) { cluster_ = value; }
    void            IncCluster() { cluster_++; }
    void            ResetCluster() { cluster_ = (unsigned int)-1; }

    size_t          GetEffectiveSize() const { return efflength_; }
    void            SetEffectiveSize( size_t value ) { efflength_ = value; }

    const char*     GetDescription() const { return description_.c_str(); }
    const unsigned char* GetResidues() const { return residues_; }

    unsigned char   GetResidueAt( size_t n ) const;
    void            SetResidueAt( unsigned char r, size_t );

    bool            IsUsedAt ( size_t n ) const;
    void            SetUnusedAt ( size_t );
    void            SetUsedAt   ( size_t );
    bool            IsUsedInExtentAt ( size_t n ) const;
    void            SetUsedInExtentAt( size_t );
    void            UnsetUsedInExtentAt( size_t );

    bool            IsUsedAndInfAt( size_t p ) const;

    float           GetWeightAt( size_t n ) const;
    void            SetWeightAt ( float w, size_t );

    float           GetGlbWeight() const { return glbwght_; }
    void            SetGlbWeight( float value ) { glbwght_ = value; }

    void AppendDescription( const char* desc, size_t pos, size_t len ) { description_.append( desc, pos, len ); }

    bool            GetUsed() const { return used_; }
    void            SetUsed( bool u ) { used_ = u; }

    size_t          GetFirstUsed() const { return firstused_; }
    void            SetFirstUsed( size_t value ) { firstused_ = value; }
    size_t          GetLastUsed() const { return lastused_; }
    void            SetLastUsed( size_t value ) { lastused_ = value; }

    virtual void    push( unsigned char r, unsigned char f = BMS_USED, float w = 0.0 );
    virtual void    clear();//clear sequence

protected:
    explicit        BMSequence( const BMSequence& );
    explicit        BMSequence();

    virtual BMSequence& operator=( const BMSequence& one );

    virtual void    Realloc( size_t newcap );
    virtual void    Init();

protected:
    mystring            description_;   //description of sequence
    unsigned char*      residues_;      //residues (sequence)
    unsigned char*      flags_;         //flags of residue usage for each position
    float*              weights_;       //residue weights for each positions
    float               glbwght_;       //global weight
    size_t              length_;        //sequence length
    size_t              efflength_;     //effective length of sequence (gaps excluded)
    size_t              capacity_;      //current capacity of sequence
    bool                used_;          //whether or not the sequence is used
    size_t              firstused_;     //index of the first used position
    size_t              lastused_;      //index of the last used position
    //
    unsigned int        cluster_;       //cluster number if clustering is in effect
};

// =========================================================================
// INLINES
//
// get residue at the given position
//
inline
unsigned char BMSequence::GetResidueAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !residues_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::GetResidueAt: Memory access error.");
#endif
    return residues_[n];
}

// set residue at position n
//
inline
void BMSequence::SetResidueAt( unsigned char r, size_t n )
{
#ifdef __DEBUG__
    if( residues_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::SetResidueAt: Memory access error.");
#endif
    residues_[n] = r;
}

// -------------------------------------------------------------------------
// get `used' flag at the given position
//
inline
bool BMSequence::IsUsedAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !flags_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::IsUsedAt: Memory access error.");
#endif
    return GetUsed() && flags_[n] != BMS_UNUSED;
}

// Set the flag `unused' at position n
//
inline
void BMSequence::SetUnusedAt( size_t n )
{
#ifdef __DEBUG__
    if( flags_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::SetUnusedAt: Memory access error.");
#endif
    flags_[n] = BMS_UNUSED;
}

// set the flag `used' at position n
//
inline
void BMSequence::SetUsedAt( size_t n )
{
#ifdef __DEBUG__
    if( flags_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::SetUsedAt: Memory access error.");
#endif
    flags_[n] = BMS_USED;
}

// check whether the sequence is used and informative at position n
//
inline
bool BMSequence::IsUsedAndInfAt( size_t n ) const
{
    if( !IsUsedAt( n ))
        return false;
    return GetFirstUsed() <= n && n <= GetLastUsed();
}

// -------------------------------------------------------------------------
// IsUsedInExtentAt: get a flag of whether the sequence is used in the 
// extent computed for the given position
//
inline
bool BMSequence::IsUsedInExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !flags_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::IsUsedInExtentAt: Memory access error.");
#endif
    return flags_[n] == BMS_USED_IN_EXTENT;
}

// SetUsedInExtentAt: set the flag of being used in extent for the given 
// position 
//
inline
void BMSequence::SetUsedInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( flags_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMSequence::SetUsedInExtentAt: Memory access error." );
#endif
    flags_[n] = BMS_USED_IN_EXTENT;
}

// UnsetUsedInExtentAt: unset the flag of being used in extent for the 
// given position
//
inline
void BMSequence::UnsetUsedInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( flags_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR( "BMSequence::UnsetUsedInExtentAt: Memory access error." );
#endif
    flags_[n] &= BMS_USED;
}

// -------------------------------------------------------------------------
// GetWeightAt: get weight at position n
//
inline
float BMSequence::GetWeightAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !weights_ || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::GetWeightAt: Memory access error.");
#endif
    return weights_[n];
}

// SetWeightAt: set weight at position n
//
inline
void BMSequence::SetWeightAt( float w, size_t n )
{
#ifdef __DEBUG__
    if( weights_ == NULL || length_ <= n )
        throw MYRUNTIME_ERROR("BMSequence::SetWeightAt: Memory access error.");
#endif
    weights_[n] = w;
}

#endif//__BMSequence_h__
