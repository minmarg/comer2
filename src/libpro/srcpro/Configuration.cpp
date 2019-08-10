/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdlib.h>

#include "liblib/mybase.h"
#include "liblib/ConfigFile.h"
#include "libpro/srcpro/SUBSTABLE.h"
#include "Configuration.h"

//designates calculated gap costs
const char* Configuration::s_autocostsym_ = "A";

// //configformat1.0| -0x2a
// const char* Configuration::signature[] = {
//     "\nCONFIGURATION VERSION ",
//     configversion,
//     "\x39\x45\x44\x3c\x3f\x3d\x3c\x45\x48\x43\x37\x4a\x07\x04\x06\x00",
//     NULL
// };

//key names:
const char* Configuration::s_section_ungapped_  = "UNGAPPED";
const char* Configuration::s_section_gapped_    = "GAPPED";
const char* Configuration::s_space_             = "_";
//
const char* Configuration::s_gap_open_cost_key_ = "Gap-open-cost";
const char* Configuration::s_gap_extn_cost_key_ = "Gap-extension-cost";
const char* Configuration::s_lambda_key_        = "Lambda";
const char* Configuration::s_K_key_             = "K";
const char* Configuration::s_H_key_             = "H";
const char* Configuration::s_alpha_key_         = "Alpha";
const char* Configuration::s_beta_key_          = "Beta";
const char* Configuration::s_scalef_key_        = "Scale-factor";

// -------------------------------------------------------------------------
// SetUngappedParams: modifies Configuration to contain ungapped
// parameters
//
void SetUngappedParams( Configuration& config )
{
    config.SetGapOpenCost( Configuration::NOGAPVAL );
    config.SetGapExtendCost( Configuration::NOGAPVAL );
    config.SetLambda(   STABLE.StatisParam( Ungapped, Lambda ));
    config.SetK(        STABLE.StatisParam( Ungapped, K ));
    config.SetH(        STABLE.StatisParam( Ungapped, H ));
    config.SetAlpha(    STABLE.StatisParam( Ungapped, alpha ));
    config.SetBeta(     STABLE.StatisParam( Ungapped, beta ));
    config.SetScaleFactor( 1.0f );
}

// =========================================================================
// Constructor
//
Configuration::Configuration( const char* fullname )
:   filename_( fullname )
{
    auto_gap_cost_  = false;
    gap_open_cost_  = NOGAPVAL;
    gap_extend_cost_= NOGAPVAL;
    lambda_         = GetDefault_lambda();
    K_              = GetDefault_K();
    H_              = GetDefault_H();
    alpha_          = GetDefault_alpha();
    beta_           = GetDefault_beta();
    scale_          = GetDefault_scale();
}

// Constructor
//
Configuration::Configuration( const char* fullname, int cost_open, int cost_extend )
:   filename_( fullname )
{
    auto_gap_cost_  = false;
    gap_open_cost_  = cost_open;
    gap_extend_cost_= cost_extend;
    lambda_         = GetDefault_lambda();
    K_              = GetDefault_K();
    H_              = GetDefault_H();
    alpha_          = GetDefault_alpha();
    beta_           = GetDefault_beta();
    scale_          = GetDefault_scale();
}

// -------------------------------------------------------------------------
// Default construction
//
Configuration::Configuration()
:   filename_( NULL )
{
    auto_gap_cost_  = false;
    gap_open_cost_  = GetDefault_gap_open_cost();
    gap_extend_cost_= GetDefault_gap_extend_cost();
    lambda_         = GetDefault_lambda();
    K_              = GetDefault_K();
    H_              = GetDefault_H();
    alpha_          = GetDefault_alpha();
    beta_           = GetDefault_beta();
    scale_          = GetDefault_scale();
}

// Copy constructor
//
Configuration::Configuration( const Configuration& one )
{
    operator=( one );
}

// -------------------------------------------------------------------------
// Destructor
//
Configuration::~Configuration()
{
}

// -------------------------------------------------------------------------
// Assignment operator
//
const Configuration& Configuration::operator=( const Configuration& one )
{
    SetFilename( one.GetFilename());

    SetAutoGapOpenCost( one.GetAutoGapOpenCost());
    SetGapOpenCost( one.GetGapOpenCost());
    SetGapExtendCost( one.GetGapExtendCost());
    SetLambda( one.GetLambda());
    SetK( one.GetK());
    SetH( one.GetH());
    SetAlpha( one.GetAlpha());
    SetBeta( one.GetBeta());
    SetScaleFactor( one.GetScaleFactor());

    return *this;
}

// -------------------------------------------------------------------------
// ReadUngapped: read parameters for ungapped configuration from file
//
void Configuration::ReadUngapped()
{
    SetGapOpenCost( NOGAPVAL );
    SetGapExtendCost( NOGAPVAL );
    Read();
}

// -------------------------------------------------------------------------
// WriteUngapped: write parameters of ungapped configuration to file
//
void Configuration::WriteUngapped()
{
    SetGapOpenCost( NOGAPVAL );
    SetGapExtendCost( NOGAPVAL );
    Write();
}

// -------------------------------------------------------------------------
// Read: reads configuration from file;
// NOTE: gap open cost and gap extension cost must be specified before
// calling this method, because section name is determined by them
//
void Configuration::Read()
{
    char    section[BUF_MAX];
    bool    auto_open = GetAutoGapOpenCost();
    int     cost_open = GetGapOpenCost();
    int     cost_extend = GetGapExtendCost();
    bool    read_costs = false;//true;
    bool    read_scale = false;

    if(( !auto_open &&
        cost_open == GetDefault_gap_open_cost()) ||
        cost_extend == GetDefault_gap_extend_cost())
        throw MYRUNTIME_ERROR( "Configuration::Read: Gap scheme must be specified before reading parameters." );

    if( cost_open == NOGAPVAL && cost_extend == NOGAPVAL ) {
        strcpy( section, s_section_ungapped_ );
        read_costs = false;//for ungapped configuration, do not read gap cost information
        read_scale = true;//for ungapped configuration, scale factor must be read
    } else {
        if( cost_open == NOGAPVAL || cost_extend == NOGAPVAL )
            throw MYRUNTIME_ERROR( "Configuration::Read: Invalid gap costs specified." );
        else {
            if( auto_open )
                sprintf( section, "%s%s%s%s%d", s_section_gapped_, s_space_, s_autocostsym_, s_space_, cost_extend );
            else
                sprintf( section, "%s%s%d%s%d", s_section_gapped_, s_space_, cost_open, s_space_, cost_extend );
        }
    }

    ReadSection( section, read_costs, read_scale );

    if(( !auto_open &&
        cost_open != GetGapOpenCost()) ||
        cost_extend != GetGapExtendCost())
        throw MYRUNTIME_ERROR( "Configuration::Read: Section for given gap scheme not found." );
}

// -------------------------------------------------------------------------
// Write: write configuration to file; at the same time, check whether the 
// costs specified are valid 
//
void Configuration::Write() const
{
    char    section[BUF_MAX];
    bool    auto_open = GetAutoGapOpenCost();
    int     cost_open = GetGapOpenCost();
    int     cost_extend = GetGapExtendCost();
    bool    write_costs = true;
    bool    write_scale = false;

    if(( !auto_open &&
        cost_open == GetDefault_gap_open_cost()) ||
        cost_extend == GetDefault_gap_extend_cost())
        throw MYRUNTIME_ERROR( "Configuration::Write: Unspecified gap scheme." );

    if( cost_open == NOGAPVAL && cost_extend == NOGAPVAL ) {
        strcpy( section, s_section_ungapped_ );
        write_costs = false;//for ungapped configuration, do not write gap cost information
        write_scale = true;//for ungapped configuration, write scale factor
    } else {
        if( cost_open == NOGAPVAL || cost_extend == NOGAPVAL )
            throw MYRUNTIME_ERROR( "Configuration::Write: Invalid gap costs specified." );
        else {
            if( auto_open )
                sprintf( section, "%s%s%s%s%d", s_section_gapped_, s_space_, s_autocostsym_, s_space_, cost_extend );
            else
                sprintf( section, "%s%s%d%s%d", s_section_gapped_, s_space_, cost_open, s_space_, cost_extend );
        }
    }

    WriteSection( section, write_costs, write_scale );
}

// -------------------------------------------------------------------------
// ReadSection: read data from section
//
void Configuration::ReadSection( const char* section, bool read_costs, bool read_scale )
{
    const mystring preamb = "Configuration::ReadSection: ";
    const int tmpsz = 10;
    char    tmpstrval[tmpsz];
    int     tmpintval;
    float   tmpfpval;
    bool    auto_open = GetAutoGapOpenCost();
    bool    found = false;

    ConfigFile config( section, GetFilename());

    if( read_costs ) {
        // Gap open cost
        if( auto_open )
            found = config.GetString( s_gap_open_cost_key_, tmpstrval, tmpsz );
        else
            found = config.GetInt( s_gap_open_cost_key_, &tmpintval );

        if( !found )
            throw MYRUNTIME_ERROR( preamb + "Key " + s_gap_open_cost_key_ + " not found." );

        if( auto_open ) {
            if( strlen( tmpstrval ) != strlen( s_autocostsym_ ) &&
                strcmp( tmpstrval, s_autocostsym_ ))
                throw MYRUNTIME_ERROR( preamb + "Inconcsistent value of " + s_gap_open_cost_key_ );
        } else
            SetGapOpenCost(( tmpintval < 0 )? -tmpintval: tmpintval );

        // Gap extension cost
        found = config.GetInt( s_gap_extn_cost_key_, &tmpintval );

        if( !found )
            throw MYRUNTIME_ERROR( preamb + "Key " + s_gap_extn_cost_key_ + " not found." );

        SetGapExtendCost(( tmpintval < 0 )? -tmpintval: tmpintval );
    }

    // parameter Lambda
    found = config.GetFloat( s_lambda_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_lambda_key_ + " not found." );

    CheckLambda( tmpfpval );
    SetLambda( tmpfpval );

    // parameter K
    found = config.GetFloat( s_K_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_K_key_ + " not found." );

    CheckK( tmpfpval );
    SetK( tmpfpval );

    // parameter H 
    found = config.GetFloat( s_H_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_H_key_ + " not found." );

    CheckH( tmpfpval );
    SetH( tmpfpval );

    // parameter alpha 
    found = config.GetFloat( s_alpha_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_alpha_key_ + " not found." );

    CheckAlpha( tmpfpval );
    SetAlpha( tmpfpval );

    // parameter beta
    found = config.GetFloat( s_beta_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_beta_key_ + " not found." );

    CheckBeta( tmpfpval );
    SetBeta( tmpfpval );

    if( ! read_scale )
        return;

    // scale factor
    found = config.GetFloat( s_scalef_key_, &tmpfpval );

    if( !found )
        throw MYRUNTIME_ERROR( preamb + "Key " + s_scalef_key_ + " not found." );

    CheckScaleFactor( tmpfpval );
    SetScaleFactor( tmpfpval );
}

// -------------------------------------------------------------------------
// WriteSection: write section data to file
//
void Configuration::WriteSection( const char* section, bool write_costs, bool write_scale ) const
{
    ConfigFile  config( section, GetFilename());
    bool        auto_open = GetAutoGapOpenCost();

    if( write_costs ) {
        // Gap open cost
        if( auto_open )
            config.WriteString( s_gap_open_cost_key_, s_autocostsym_ );
        else {
            CheckCost( GetGapOpenCost());
            config.WriteInt( s_gap_open_cost_key_, GetGapOpenCost());
        }

        // Gap extension cost
        CheckCost( GetGapExtendCost());
        config.WriteInt( s_gap_extn_cost_key_, GetGapExtendCost());
    }

    // parameter Lambda
    CheckLambda( GetLambda());
    config.WriteFloat( s_lambda_key_, GetLambda());

    // parameter K
    CheckK( GetK());
    config.WriteFloat( s_K_key_, GetK());

    // parameter H 
    CheckH( GetH());
    config.WriteFloat( s_H_key_, GetH());

    // parameter alpha 
    CheckAlpha( GetAlpha());
    config.WriteFloat( s_alpha_key_, GetAlpha());

    // parameter beta
    CheckBeta( GetBeta());
    config.WriteFloat( s_beta_key_, GetBeta());

    if( ! write_scale )
        return;

    // scale factor
    CheckScaleFactor( GetScaleFactor());
    config.WriteFloat( s_scalef_key_, GetScaleFactor());
}
