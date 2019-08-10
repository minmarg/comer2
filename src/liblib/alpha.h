/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __alpha_h__
#define __alpha_h__

// effective number of amino acids
#define NUMAA 20
#define NUMALPH 25

#define X           23
// #define ASTERISK    23
#define GAP         24

// interface of functions
//
char DehashCode( unsigned char );
int HashAlphSymbol( char );
int NumAlphabet();
bool IsValidResSym( unsigned char );
//
char DehashSSCode( unsigned char );
char DehashSSCodeLc( unsigned char );
int HashSSState( char );
//

//                             //    012345678901234567890123
extern const char*     gAAcids;// = "ARNDCQEGHILKMFPSTWYVBZJX";
// static const char*     gAAcids = "ARNDCQEGHILKMFPSTWYVBZX*";

//                               //    0123456789012345678901234
extern const char*     gAlphabet;// = "ARNDCQEGHILKMFPSTWYVBZJX-";
// static const char*     gAlphabet = "ARNDCQEGHILKMFPSTWYVBZX*-";

//{{SS states
enum SSSTATES {
    SS_C,
    SS_E,
    SS_H,
    SS_NSTATES
};
extern const char*  gSSAlphabet;// = "CEH";
//}}

#endif//__alpha_h__
