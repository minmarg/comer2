/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdlib.h>

#include "liblib/alpha.h"
#include "BMSequence.h"

// -------------------------------------------------------------------------
// Default constructor
//
BMSequence::BMSequence()
{
    Init();
}

// Initialization constructor
//
BMSequence::BMSequence( size_t reservation )
{
    Init();
    Realloc( reservation );
}

// Copy constructor
//
BMSequence::BMSequence( const BMSequence& one )
{
    Init();
    if( !one.GetCapacity())
        throw MYRUNTIME_ERROR("BMSequence::BMSequence: Uninialized argument.");

    Realloc( one.GetCapacity());
    *this = one;
}

// -------------------------------------------------------------------------
// Destructor
//
BMSequence::~BMSequence()
{
    if( residues_ ) free( residues_ );
    if( flags_ ) free( flags_ );
    if( weights_ ) free( weights_ );
}

// -------------------------------------------------------------------------
// Init: initialization
//
void BMSequence::Init()
{
    residues_ = NULL;
    flags_ = NULL;
    weights_ = NULL;
    glbwght_ = 0.0f;
    length_ = 0;
    efflength_ = 0;
    capacity_ = 0;
    used_ = true;
    firstused_ = MYSIZE_MAX;
    lastused_ = MYSIZE_MAX;
    ResetCluster();
}

// -------------------------------------------------------------------------
// Assignment
//
BMSequence& BMSequence::operator=( const BMSequence& one )
{
    if( capacity_ < one.capacity_ )
        Realloc( one.capacity_ );

#ifdef __DEBUG__
    if( !residues_ || !flags_ || !weights_ )
        throw MYRUNTIME_ERROR("BMSequence::operator=: Memory access error.");
#endif

    memcpy( residues_,  one.residues_,  sizeof( unsigned char ) * one.length_ );
    memcpy( flags_,     one.flags_,     sizeof( unsigned char ) * one.length_ );
    memcpy( weights_,   one.weights_,   sizeof( float ) * one.length_ );

    glbwght_ = one.glbwght_;
    used_ = one.used_;
    length_ = one.length_;
    efflength_ = one.efflength_;
    firstused_ = one.firstused_;
    lastused_ = one.lastused_;
    cluster_ = one.cluster_;

    memset( residues_ + length_, 0, sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( flags_ + length_,    0, sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( weights_ + length_,  0, sizeof( float ) * ( capacity_ - length_ ));

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: (re)allocate memory
//
void BMSequence::Realloc( size_t newcap )
{
    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        residues_ = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        flags_ = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        weights_ = ( float* )malloc( sizeof( float ) * newcap );
    } else {
        residues_ = ( unsigned char* )realloc( residues_, sizeof( unsigned char ) * newcap );
        flags_ = ( unsigned char* )realloc( flags_, sizeof( unsigned char ) * newcap );
        weights_ = ( float* )realloc( weights_, sizeof( float ) * newcap );
    }

    if( !residues_ || !flags_ || !weights_ )
        throw MYRUNTIME_ERROR("BMSequence::Realloc: Not enough memory.");

    unsigned char*  tress = residues_;
    unsigned char*  tflgs = flags_;
    float*          tweis = weights_;

    if( capacity_ != 0 ) {
        tress = residues_ + capacity_;
        tflgs = flags_ + capacity_;
        tweis = weights_ + capacity_;
    }

    memset( tress, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tflgs, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tweis, 0, sizeof( float ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: push information describing one position of the sequence
//
void BMSequence::push( unsigned char r, unsigned char f, float w )
{
    if( capacity_ <= length_ ) {
        size_t newcap = TIMES2( capacity_ );
        if( newcap <= length_ )
            newcap = length_ + 1;
        Realloc( newcap );
    }

    residues_[length_] = r;
    flags_[length_] = f;
    weights_[length_] = w;

    length_++;

    if( r != GAP )
        efflength_++;
}

// -------------------------------------------------------------------------
// clear: clear all buffers allocated for the sequence
//
void BMSequence::clear()
{
    memset( residues_, 0, sizeof( unsigned char ) * capacity_ );
    memset( flags_, 0, sizeof( unsigned char ) * capacity_ );
    memset( weights_, 0, sizeof( float ) * capacity_ );

    length_ = 0;
    efflength_ = 0;
//
    // TODO: modify appropriately MSA sources, where they use the method
    description_.erase();
    glbwght_ = 0.0f;
    used_ = true;
    firstused_ = MYSIZE_MAX;
    lastused_ = MYSIZE_MAX;
    ResetCluster();
}
