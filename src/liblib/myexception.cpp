/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "platform.h"
#include "mylimits.h"
#include "myexception.h"
#include <stdio.h>

// -------------------------------------------------------------------------
// CLASS myruntime_error
//
// pretty_format: place exception information on mystring object
//
mystring myruntime_error::pretty_format( mystring preamb ) const throw()
{
    char buf[BUF_MAX];
    bool emp = _errdsc.empty();
    if( emp )
        preamb += myexception::what();
    else
        preamb += _errdsc;
    int bfile = 0, bline = 0, bfunc = 0;
    if( file())
        bfile = 1;
    if( line())
        bline = 2;
    if( function())
        bfunc = 4;
    if( bfile | bline | bfunc ) {
        preamb += NL;
        preamb += "    (";
        if( file()) {
            preamb += "File: ";
            preamb += file();
            if( bline | bfunc ) {
                preamb += "; ";
            }
        }
        if( line()) {
            preamb += "Line: ";
            sprintf( buf, "%u", line());
            preamb += buf;
            if( bfunc ) {
                preamb += "; ";
            }
        }
        if( function()) {
            preamb += "Function: ";
            preamb += function();
        }
        preamb += ")";
    }
    return preamb;
}

