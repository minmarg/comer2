/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include <thread>
#include <chrono>

#include "platform.h"
#include "mytype.h"
#include "mystring.h"
#include "myexception.h"
#include "mylimits.h"
#include "msg.h"

std::chrono::high_resolution_clock::time_point gtSTART = 
std::chrono::high_resolution_clock::now();

const char* PROGDIR = NULL;
const char* PROGNAME = NULL;
const char* PROGVERSION = NULL;
const char* PROGREFERENCES[] = {
    "Margelevicius, M. (2020). "
    "COMER2: GPU-accelerated sensitive and specific homology searches. "
    "Bioinformatics 36, 3570-3572.",

    "Margelevicius, M. (2016). "
    "Bayesian nonparametrics in protein remote homology search. "
    "Bioinformatics 32, 2744-2752.",

    NULL
};

int*        __PARGC__ = NULL;
char***     __PARGV__ = NULL;

int         VERBOSE = 0;
bool        WARNINGSRECORDED = false;
bool        ERRORSRECORDED = false;


// -------------------------------------------------------------------------
// set global variables
//
void SetProgramName( const char* name, const char* version )
{
    PROGDIR = my_dirname( name );
    PROGNAME = my_basename( name );
    if( version )
        PROGVERSION = version;
}

void SetArguments( int* pargc, char*** pargv )
{
    __PARGC__ = pargc;
    __PARGV__ = pargv;
}

void SetVerboseMode( int value )
{
    VERBOSE = value;
}

// -------------------------------------------------------------------------
// messaging routines of the toolset
//
void error( const char* str, bool putnl )
{
    ERRORSRECORDED = true;
    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] ERROR: %s%s%s", PROGNAME, str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] ERROR: %s%s", PROGNAME, str, NL);
    else if(str && putnl)
        fprintf( stderr, "ERROR: %s%s%s", str, NL, NL);
    else if(str)
        fprintf( stderr, "ERROR: %s%s", str, NL);
}

void warning( const char* str, bool putnl, int/* minlevel*/)
{
//     if( VERBOSE < minlevel )
//         return;
    WARNINGSRECORDED = true;
    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] WARNING: %s%s%s", PROGNAME, str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] WARNING: %s%s", PROGNAME, str, NL);
    else if(str && putnl)
        fprintf( stderr, "WARNING: %s%s%s", str, NL, NL);
    else if(str)
        fprintf( stderr, "WARNING: %s%s", str, NL);
}

void checkforwarnings()
{
    if( !WARNINGSRECORDED && !ERRORSRECORDED)
        return;
    if(PROGNAME)
        fprintf( stderr, "%s[%s] There are %s.%s", NL, PROGNAME,
                ERRORSRECORDED? "ERRORS": "WARNINGS", NL);
    else 
        fprintf( stderr, "%sThere are %s.%s", NL,
                ERRORSRECORDED? "ERRORS": "WARNINGS", NL);
}

void message( const char* str, bool putnl, int minlevel )
{
    if( VERBOSE < minlevel )
        return;

    std::chrono::high_resolution_clock::time_point tnow = 
    std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elpsd = tnow - gtSTART;

    if(str && PROGNAME && putnl)
        fprintf( stderr, "[%s] {%.6fs} (%zu) %s%s%s", PROGNAME, elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL, NL);
    else if(str && PROGNAME)
        fprintf( stderr, "[%s] {%.6fs} (%zu) %s%s", PROGNAME, elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL);
    else if(str && putnl)
        fprintf( stderr, "{%.6fs} (%zu) %s%s%s", elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL, NL);
    else if(str)
        fprintf( stderr, "{%.6fs} (%zu) %s%s", elpsd.count(),
                std::hash<std::thread::id>{}(std::this_thread::get_id()), 
                str, NL);
}

void messagec( char sym, bool indent, int minlevel )
{
    if( VERBOSE < minlevel )
        return;
    if( indent )
        fprintf( stderr, "%20c", sym );
    else
        fprintf( stderr, "%c", sym );
}

void progname_and_version( FILE* fp )
{
    if( PROGNAME )
        fprintf( fp, "%s", PROGNAME );
    if( PROGVERSION )
        fprintf( fp, " %s", PROGVERSION );
    if( PROGNAME || PROGVERSION )
        fprintf( fp, "%s%s", NL, NL );
}

// -------------------------------------------------------------------------
// print_cmdline: print command line
//
void print_cmdline( TPrintFunction pfunc, void* vpn )
{
    int n;
    if( !vpn )
        return;
    if( __PARGV__ && __PARGC__  && *__PARGV__ ) {
        for( n = 0; n < *__PARGC__; n++ )
            pfunc( vpn, "%s ", (*__PARGV__)[n] );
    }
    pfunc( vpn, "%s", NL );
}

// -------------------------------------------------------------------------
// print_dtime: print date and local time
//
void print_dtime(int minlevel)
{
    if( VERBOSE < minlevel )
        return;
    char tmstr[BUF_MAX];
    time_t t = time(NULL);
    struct tm* ctm = localtime(&t);
    if( ctm ) {
        if( strftime(tmstr, sizeof(tmstr), "%c", ctm) != 0)
            fprintf(stderr,"[%s] %s%s%s", PROGNAME, tmstr,NL,NL);
    }
}

// -------------------------------------------------------------------------
// path-processing functions
//
const char* my_basename( const char* name )
{
    const char* bp = NULL;
    if( name )
        for(bp  = name + strlen( name );
            bp != name && bp[-1] != DIRSEP;
            bp-- );
    return  bp;
}

// my_dirname: return directory name given pathname; if directory is found
//     in the given path, return the current directory (`.')
//
const char* my_dirname( const char* name )
{
    static const size_t size = KBYTE;
    static char dir[size];
    const char* bp = name;

    dir[0] = '.'; dir[1] = 0;

    if( !bp )
        return dir;

    for( bp += strlen( name ); bp != name && bp[-1] != DIRSEP; bp-- );

    if( bp != name && *bp == 0 ) {
        //the last character is dir separator
        for( ; bp != name && bp[-1] == DIRSEP; bp-- );

        if( bp != name )
            //find separator before the last one
            for( ; bp != name && bp[-1] != DIRSEP; bp-- );
    }

    if( bp == name )
        //no directory separator has been found
        return dir;

    for( ; bp != name && bp[-1] == DIRSEP; bp-- );

    if( bp == name ) {
        bp++;

        if( *bp == DIRSEP && bp[1] != DIRSEP )
            //there are exactly two separators at the beginning; save them
            bp++;
    }

    if( bp - name < (ssize_t)size ) {
        //if enough memory allocated
        memcpy( dir, name, bp - name );
        dir[ bp - name ] = 0;
    }

    return dir;
}

// -------------------------------------------------------------------------
// usage: returns <instructions> translated by appropriately inserting
//     program name, version, and date information
//
mystring usage( const char* path, const char* instructions, const char* version, const char* date )
{
    size_t      pos;
    mystring    instr = instructions;
    mystring    prog = path;
    mystring    full;
    bool        first = true;

    if( !path || !instructions )
        return instr;

    if(( pos = prog.rfind( DIRSEP )) != mystring::npos ){
        prog = prog.substr( pos+1 );
    }

    full = prog + mystring(" ");

    if( version && strlen( version ))
        full += version;

    if( date && strlen( date ))
        full += mystring( " (" ) + date + mystring( ")" );

    while(( pos = instr.find( "<>" )) != mystring::npos ){
        instr.erase( pos, 2 );
        instr.insert( pos, first? full: prog );
        if( first ) 
            first = false;
    }
    while(( pos = instr.find("-*")) != mystring::npos ){
        instr.erase( pos, 2 );
        for( size_t len = full.length(); len; len-- )
            instr.insert( pos++, "-");
    }
    return instr;
}

// -------------------------------------------------------------------------
// file_print: redirect formatted output string to file; file
//     pointer (vpn) must be preintialized and must be valid
//
int file_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    FILE*   fp = ( FILE* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vfprintf( fp, format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error("file_print: Formatted print to file failed.", __EXCPOINT__ );

    return ret;
}

// -------------------------------------------------------------------------
// string_print: same as file_print except that redirection is to chr. 
//     string, which is assumed to have enough space to contain information;
//     string pointer (vpn) must be preallocated and must be valid
//
int string_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    char*   sp = ( char* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vsprintf( sp + strlen( sp ), format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error("string_print: Formatted print to string failed.", __EXCPOINT__ );

    return ret;
}

