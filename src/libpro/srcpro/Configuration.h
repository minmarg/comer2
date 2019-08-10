/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __Configuration_h__
#define __Configuration_h__

#include <stdlib.h>
#include "liblib/mybase.h"

class Configuration;

// // Version history:
// //
// static const char* configversion = "0.1";

void SetUngappedParams( Configuration& config );

// Known gap penalty schemes
enum TConfigType {
    CTUngapped,
    CTGapped,
    NoCTypes
};

////////////////////////////////////////////////////////////////////////////
// Class Configuration
//
// Interface for accessing reference statistical parameters
//
class Configuration
{
public:
    enum { NOGAPVAL = 0 };

    static const char*  s_section_ungapped_;
    static const char*  s_section_gapped_;
    static const char*  s_space_;

    static const char*  s_gap_open_cost_key_;
    static const char*  s_gap_extn_cost_key_;
    static const char*  s_lambda_key_;
    static const char*  s_K_key_;
    static const char*  s_H_key_;
    static const char*  s_alpha_key_;
    static const char*  s_beta_key_;
    static const char*  s_scalef_key_;

public:
    Configuration();
    Configuration( const char* fullname );
    Configuration( const char* fullname, int, int );
    Configuration( const Configuration& );
    ~Configuration();

    const Configuration& operator=( const Configuration& );

    void            ReadUngapped();//read parameters for ungapped configuration 
    void            WriteUngapped();//write parameters for ungapped configuration
    void            Read();//read configuration from file
    void            Write() const;//write configuration to file

    bool            GetAutoGapOpenCost() const          { return auto_gap_cost_;    }
    int             GetGapOpenCost() const              { return gap_open_cost_;    }
    int             GetGapExtendCost() const            { return gap_extend_cost_;  }
    float           GetLambda() const                   { return lambda_;           }
    float           GetK() const                        { return K_;                }
    float           GetH() const                        { return H_;                }
    float           GetAlpha() const                    { return alpha_;            }
    float           GetBeta() const                     { return beta_;             }
    float           GetScaleFactor() const              { return scale_;            }

    void            SetAutoGapOpenCost( bool value )    { auto_gap_cost_    = value; }
    void            SetGapOpenCost( int value )         { gap_open_cost_    = abs(value); }
    void            SetGapExtendCost( int value )       { gap_extend_cost_  = abs(value); }
    void            SetLambda( float value )            { lambda_           = value; }
    void            SetK( float value )                 { K_                = value; }
    void            SetH( float value )                 { H_                = value; }
    void            SetAlpha( float value )             { alpha_            = value; }
    void            SetBeta( float value )              { beta_             = value; }
    void            SetScaleFactor( float value )       { scale_            = value; }

    const char*     GetFilename() const                 { return filename_; }
    void            SetFilename( const char* name )     { filename_ = name; }

    bool            IsUngapped() const { return GetGapOpenCost() == NOGAPVAL && GetGapExtendCost() == NOGAPVAL; }

protected:
    void            ReadSection( const char* section, bool read_costs, bool read_scale );
    void            WriteSection( const char* section, bool write_costs, bool write_scale ) const;

    static void     CheckCost( int value );
    static void     CheckLambda( float value );
    static void     CheckK( float value );
    static void     CheckH( float value );
    static void     CheckAlpha( float value );
    static void     CheckBeta( float value );
    static void     CheckScaleFactor( float value );

    static int      GetDefault_gap_open_cost()          { return -1; }
    static int      GetDefault_gap_extend_cost()        { return -1; }
    static float    GetDefault_lambda()                 { return -1.0f; }
    static float    GetDefault_K()                      { return -1.0f; }
    static float    GetDefault_H()                      { return -1.0f; }
    static float    GetDefault_alpha()                  { return 0.0f; }
    static float    GetDefault_beta()                   { return 0.0f; }
    static float    GetDefault_scale()                  { return 1.0f; }

private:
//     enum {
//             BEGTEXT,
//             CONFVER,
//             BINSIGN,
//             END
//     };

    const char* filename_;

    bool    auto_gap_cost_;//calculated gap costs
    int     gap_open_cost_;//gap open penalty
    int     gap_extend_cost_;//gap extension penalty
    float   lambda_;//statistical parameter Lambda
    float   K_;//statistical parameter K
    float   H_;//statistical parameter H (entropy)
    float   alpha_;//statistical parameter alpha (slope of edge-effect correction by linear regression)
    float   beta_;//statistical parameter beta (intercept of edge-effect correction by linear regression)
    float   scale_;//scale factor to scale score matrix

    static const char*  s_autocostsym_;//symbol indicating the calculation of gap costs
//     static const char*  signature[];//signature to put/get it in/from the configuration file
};

// -------------------------------------------------------------------------
// INLINES
//
// CheckCost: check gap cost value
//
inline
void Configuration::CheckCost( int value )
{
    if( value < 0 )
        throw MYRUNTIME_ERROR( "Configuration::CheckCost: Invalid gap cost." );
}

// CheckLambda: check lambda value
//
inline
void Configuration::CheckLambda( float value )
{
    if( value <= 0.0f /*|| 1.0f <= value */)
        throw MYRUNTIME_ERROR( "Configuration::CheckLambda: Invalid Lambda." );
}

// CheckK: check K value
//
inline
void Configuration::CheckK( float value )
{
    if( value <= 0.0f || 1.0f <= value )
        throw MYRUNTIME_ERROR( "Configuration::CheckK: Invalid parameter K." );
}

// CheckH: check H value
//
inline
void Configuration::CheckH( float value )
{
    if( value <= 0.0f || 10.0f <= value )
        throw MYRUNTIME_ERROR( "Configuration::CheckH: Invalid parameter H." );
}

// CheckAlpha: check alpha value
//
inline
void Configuration::CheckAlpha( float value )
{
    if( value <= 0.0f || 10.0f <= value )
        throw MYRUNTIME_ERROR( "Configuration::CheckAlpha: Invalid parameter Alpha." );
}

// CheckBeta: check beta value
//
inline
void Configuration::CheckBeta( float value )
{
    if( 0.0f < value && value < 0.0f )//the second included intentionally
        throw MYRUNTIME_ERROR( "Configuration::CheckBeta: Invalid parameter Beta." );
}

// CheckScaleFactor: check scale factor value
//
inline
void Configuration::CheckScaleFactor( float value )
{
    if( value <= 0.0f || 10.0f <= value )
        throw MYRUNTIME_ERROR( "Configuration::CheckScaleFactor: Invalid scale factor." );
}

#endif//__Configuration_h__
