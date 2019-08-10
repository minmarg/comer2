/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

// #include <ctype.h>
// #include <stdlib.h>
// #include <string.h>
#include "BinarySearchStructure.h"
#include "mylimits.h"
#include "myexception.h"
#include "mystring.h"
#include "mygetopt.h"

// file globals:

const size_t ggonMaxopts = BUF_MAX;
const char* ggosTer = "--";

// -------------------------------------------------------------------------
// OptionNameCmp: compare two option names
//
int OptionNameCmp( const void* key1, const void* key2, void* /*pars*/ )
{
//     const MyGetopt* pthis = (const MyGetopt*)pars;
//     if( pthis == NULL )
//         return 0;
    const myoption* opts1 = (const myoption*)key1;
    const myoption* opts2 = (const myoption*)key2;
    if( opts1 == NULL || opts2 == NULL )
        return 0;
    if( opts1->sLongname_ == NULL || opts2->sLongname_ == NULL )
        return 0;
    return strcmp( opts1->sLongname_, opts2->sLongname_ );
}

// -------------------------------------------------------------------------
// CLASS MyGetopt
//
// constructor
//
MyGetopt::MyGetopt( myoption* opts, const char *argv[], int argc )
:   ncnt_( 0 ),
    stopped_( false ),
    argv_( NULL ),
    argc_( 0 ),
    options_( NULL ),
    srtopts_( NULL )
{
    Init( opts, argv, argc );
}

// destructor
//
MyGetopt::~MyGetopt()
{
    argv_ = NULL;
    options_ = NULL;
    DestroySrtOpts();
}

// -------------------------------------------------------------------------
// DestroySrtOpts: desroy member srtopts_
//
void MyGetopt::DestroySrtOpts()
{
    if( srtopts_ ) {
        delete srtopts_;
        srtopts_ = (BinarySearchStructure*)(void*)0;
    }
}

// -------------------------------------------------------------------------
// Init: initialize object
//
void MyGetopt::Init( myoption* opts, const char *argv[], int argc )
{
    size_t  n;
    if( opts == NULL || argv == NULL )
        throw myruntime_error("MyGetopt::Init: Memory access error.", __EXCPOINT__ );
    options_ = opts;
    argv_ = argv;
    argc_ = argc;
    DestroySrtOpts();
    srtopts_ = new BinarySearchStructure( 
        OptionNameCmp, 10/*initial size*/, false/*duplicates*/, this/*priv. params*/);
    if( srtopts_ == NULL )
        throw myruntime_error("MyGetopt::Init: Not enough memory.", __EXCPOINT__ );
    for( n = 0; ; n++) {
        if( ggonMaxopts <= n ) {
            DestroySrtOpts();
            throw myruntime_error("MyGetopt::Init: Too long list of options.", __EXCPOINT__ );
        }
        if( options_[n].sLongname_ == NULL )
            break;
        srtopts_->Push((const myoption*)( options_ + n ));
    }
}

// -------------------------------------------------------------------------
// GetNextOption: get the next option from processing the command line;
//      return value:
//          option's value if the option has been found;
//          '?' if option has not been found;
//          -1 if the processing of the command line has to be 
//      stopped; the processing is stopped on the end of the command line or 
//      options terminator `--'.
//
int MyGetopt::GetNextOption( mystring* argument )
{
    if( srtopts_ == NULL )
        throw myruntime_error("MyGetopt::GetNextOption: Memory access error.", __EXCPOINT__ );

    if( stopped_ )
        return my_gorv_term;
    if( ncnt_ < 0 || ncnt_ > argc_ )
        return my_gorv_term;
    const char* carg = argv_[ncnt_];
    ncnt_++;
    if( carg == NULL || *carg == 0 )
        return my_gorv_term;
    if( strcmp( carg, ggosTer ) == 0 ) {
        stopped_ = true;
        return my_gorv_term;
    }
    if( *carg != '-' ) {
        if( argument )
            *argument = carg;
        return my_gorv_value;
    }
    carg++;
    if( *carg == 0 )
        return my_gorv_illformed;
    if( *carg == '-' ) {
        ++carg;
        if( *carg == 0 || *carg == '-' )
            return my_gorv_illformed;
    }

    int loc = -1;
    const myoption* fndopt = NULL;//found option
    myoption privopt = {NULL,my_n_targflags,0};//option read from the command line
    mystring arg = carg;
    mystring key = arg;
    size_t p = arg.find ('=');

    if( p != mystring::npos ) {
        //command line option is provided with '='
        key = arg.substr( 0, p );
        if( argument )
            *argument = arg.substr( p+1 );
    }

    privopt.sLongname_ = key.c_str();

    if( srtopts_->Find((const void*)&privopt, &loc ) == false || 
        loc < 0 ) {
        return my_gorv_notfound;
    }

    fndopt = ( const myoption* )srtopts_->GetValueAt(( size_t )loc );
    if( fndopt == NULL ) {
        return my_gorv_notfound;
    }

    while(( fndopt->eFlag_ == my_required_argument || 
            fndopt->eFlag_ == my_optional_argument ) && p == mystring::npos ) {
        //command line option found; next, read argument (value) if required
        if( ncnt_ > argc_ ) {
            if( fndopt->eFlag_ == my_optional_argument ) {
                if( argument )
                    argument->erase();
                break;
            }
            return my_gorv_noarg;
        }
        carg = argv_[ncnt_];
        ncnt_++;
        if( carg == NULL || *carg == 0 || *carg == '-') {
            if( fndopt->eFlag_ == my_optional_argument ) {
                if( argument )
                    argument->erase();
                ncnt_--;
                break;
            }
            return my_gorv_noarg;
        }
        if( argument )
            *argument = carg;
        break;
    }

    return fndopt->nRet_;
}
