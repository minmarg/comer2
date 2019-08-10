/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <ctype.h>
#include "msg.h"
#include "mybase.h"
#include "alpha.h"

const char* gAAcids = "ARNDCQEGHILKMFPSTWYVBZJX";
const char* gAlphabet = "ARNDCQEGHILKMFPSTWYVBZJX-";
const char* gSSAlphabet = "CEH";
const char* gSSlcAlphabet = " eh";

// -------------------------------------------------------------------------
// alphabet size
//
int NumAlphabet()
{
    return (int)strlen(gAlphabet);
}

// -------------------------------------------------------------------------
// DehashCode: return ascii amino acid symbol given its code
//
char DehashCode( unsigned char in )
{
    if( in < NUMALPH )
        return gAlphabet[in];
    throw MYRUNTIME_ERROR("Unrecognized amino acid code.");
}

// -------------------------------------------------------------------------
// HashResidue: get code for amino acid symbol
//
int HashAlphSymbol( char in )
{
    char upper = toupper( in );

    for( int n = 0; gAlphabet[n]; n++ )
        if( upper == gAlphabet[n] )
            return n;

    char instr[] = { in, 0 };
    if( upper == 'O' || upper == 'U' ) {
        warning(( mystring( "Amino acid '" ) + instr + mystring( "' replaced by X" )).c_str(), false );
        return X;
    }
    throw MYRUNTIME_ERROR( mystring("Unrecognized amino acid '") + instr + mystring("'."));
}

// -------------------------------------------------------------------------
// IsValidRes: check if `code' represents valid amino acid
//
bool IsValidResSym( unsigned char code )
{
    if( code <= X )
        return true;
    return false;
}

// =========================================================================
// DehashSSCode: returns ascii SS state symbol given its code
//
char DehashSSCode( unsigned char in )
{
    if( in < SS_NSTATES )
        return gSSAlphabet[in];
    throw MYRUNTIME_ERROR("Unrecognized SS state.");
}

// =========================================================================
// DehashSSCode: returns ascii SS state symbol for output given its code
//
char DehashSSCodeLc( unsigned char in )
{
    if( in < SS_NSTATES )
        return gSSlcAlphabet[in];
    throw MYRUNTIME_ERROR("Unrecognized SS state.");
}

// -------------------------------------------------------------------------
// HashSSState: get code for SS state symbol
//
int HashSSState( char in )
{
    char upper = toupper( in );

    for( int n = 0; gSSAlphabet[n]; n++ )
        if( upper == gSSAlphabet[n])
            return n;

    char instr[] = { in, 0 };
    throw MYRUNTIME_ERROR( mystring("Unrecognized SS state '") + instr + mystring("'."));
}

