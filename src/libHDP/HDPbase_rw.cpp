/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

// #include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "HDPbase.h"

//default menu size
static const int gcsz_defmenusz = 100;

using namespace extspsl;

// =========================================================================
// Read: read set of vectors from file; vectors are supposed to be divided 
//  into groups;
//  dids, dish indices
//
void HDPbase::ReadGroups( const char* filename, Ivector* dids )
{
    mystring        preamb = "HDPbase::ReadGroups: ";
    FILE*           fp = NULL;
    size_t          length, rbts, c;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    int             emsg;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
//     const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";

    const size_t    lenstrendofsmp = strlen( patstrendofsmp );
    const size_t    lenstrendoftbl = strlen( patstrendoftbl );
    const size_t    lenstrendofgrp = strlen( patstrendofgrp );
    float           lprob;
    int             iter;
    int             ctxtsize = 0, dim = 0;
    int             totsamples;
    int             t, notables;
    int             n, nosamples;
    int             nogroups;
    int             nodishes, nocdshs = 0;
    int             dishid;
    int             npos, kpos = -1/*to make compiler happy*/, dpos;
    myruntime_error mre;
    Pslvector*      frq = NULL;
    Dish*           dish = NULL;
    Table*          tble = NULL;
    Restaurant*     rest = NULL;

    if( dids == NULL )
        throw MYRUNTIME_ERROR( preamb + "Memory access error (dish indices).");

    if( !filename ||
      ( fp = fopen( filename, "r" )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Failed to open input file.");

    try{
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        //read iteration number
        if(( p = strstr( locbuffer, patstrgss )) != NULL ) {
            p += strlen( patstrgss );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &iter, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( iter < 0 )
                throw MYRUNTIME_ERROR( preamb + "Iteration number is invalid." );

            SetIterationRead( iter );

            //read next line
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );
        }

        //read log probability
        if(( p = strstr( locbuffer, patstrlpd )) != NULL ) {
            p += strlen( patstrlpd );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &lprob, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            SetMaxLProbData( lprob );

            //read next line
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );
        }

        //read number of groups
        if(( p = strstr( locbuffer, patstrnog )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        p += strlen( patstrnog );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &nogroups, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( nogroups < 1 )
            throw MYRUNTIME_ERROR( preamb + "Number of groups is invalid." );

        //read number of dishes
        if(( p = strstr( locbuffer, patstrnds )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        p += strlen( patstrnds );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &nodishes, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( nodishes < 1 )
            throw MYRUNTIME_ERROR( preamb + "Number of dishes is invalid." );

        //read total number of samples
        if(( p = strstr( locbuffer, patstrtns )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        p += strlen( patstrtns );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &totsamples, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( totsamples < nogroups || totsamples < nodishes )
            throw MYRUNTIME_ERROR( preamb + "Total number of samples is invalid." );

        //allocate space for dish ids
        dids->Reserve( TIMES2( totsamples ));
        dids->SetAllToValue( -1 );

        //read size of contexts
        if(( p = strstr( locbuffer, patstrssz )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        p += strlen( patstrssz );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &ctxtsize, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( ctxtsize < 1 )
            throw MYRUNTIME_ERROR( preamb + "Sample size is invalid." );
        SetCtxtSize( ctxtsize );
//
        InitBasin( totsamples );//initialize the main storage of samples
        InitMenu( totsamples );//initially, allocate a separate dish for each of the samples
        InitChain( nogroups );//allocate chain of restaurants
        ResetTotalNoSamples();
//
        while( !feof( fp )) {
            //read next group
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ))
                break;

            if( !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            if(( p = strstr( locbuffer, patstrgrp )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            //read number of tables
            if(( p = strstr( locbuffer, patstrnot )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            p += strlen( patstrnot );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

            if(( emsg = 
                read_integer( p, length - (size_t)( p - locbuffer ), &notables, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( notables < 1 )
                throw MYRUNTIME_ERROR( preamb + "Number of tables is invalid." );
//
            rest = new Restaurant( nosamples );
            if( rest == NULL )
                throw MYRUNTIME_ERROR( preamb + "Read: Not enough memory." );
            rest->SetDestroy( true );
            GetChain()->NewRestaurant( rest );
//
            for( t = 0; t < notables && !feof( fp ); t++ )
            {
                //read next table
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( !length )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                if(( p = strstr( locbuffer, patstrtbl )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                //read dish id
                if(( p = strstr( locbuffer, patstrdid )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                p += strlen( patstrdid );

                if( length <= (size_t)( p - locbuffer ))
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                if(( emsg = 
                    read_integer( p, length - (size_t)( p - locbuffer ), &dishid, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( dishid < 0 )
                    throw MYRUNTIME_ERROR( preamb + "Dish id is invalid." );
                if( TIMES2( totsamples ) <= dishid )
                    throw MYRUNTIME_ERROR( preamb + "Too large value of Dish id." );

                dish = NULL;
                if( !GetCluster4Each()) {
                    //create new dish or use existing
                    if( dids->GetValueAt( dishid ) < 0 ) {
                        dish = new Dish( GetDefDishSize());
                        if( dish == NULL )
                            throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
                        dish->SetBasin( GetBasin());
                        kpos = GetMenu()->NewDish( dish );
                        dids->SetValueAt( dishid, kpos );
                        nocdshs++;
                    }
                    else {
                        dish = GetMenu()->GetDishAt( kpos = dids->GetValueAt( dishid ));
                    }
                    if( dish == NULL )
                        throw MYRUNTIME_ERROR( preamb + "Null dish." );
                }
                //read number of samples
                if(( p = strstr( locbuffer, patstrnos )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                p += strlen( patstrnos );

                if( length <= (size_t)( p - locbuffer ))
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format." );

                if(( emsg = 
                    read_integer( p, length - (size_t)( p - locbuffer ), &nosamples, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( nosamples < 1 )
                    throw MYRUNTIME_ERROR( preamb + "Number of samples is invalid." );

                tble = NULL;
                if( !GetCluster4Each()) {
                    //create new table
                    tble = new Table( GetDefTableSize());
                    if( tble == NULL )
                        throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
                    tble->SetBasin( GetBasin());
                    tble->SetMenu( GetMenu());
                    tble->SetDishIndex( kpos );//kpos is from dish assignment above
                    rest->NewTable( tble );
                }
//
                for( n = 0; n < nosamples && !feof( fp ); n++ )
                {
                    if( !ReadSample( fp, &frq, dim ))
                        throw MYRUNTIME_ERROR( preamb + "Failed to read data." );

                    if( frq == NULL )
                        throw MYRUNTIME_ERROR( preamb + "Failed to read data." );

                    if( !dim ) {
                        dim = frq->GetSize();
                        GetMenu()->SetDim( dim );
                    }
                    else if( dim != frq->GetSize())
                        throw MYRUNTIME_ERROR( preamb + "Invalid dimensionality of samples." );
//
                    npos = GetBasin()->NewValue( frq );
//
                    if( GetCluster4Each() || dish == NULL ) {
                        //dish for each sample
                        dish = new Dish( GetDefDishSize());
                        if( dish == NULL )
                            throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
                        dish->SetBasin( GetBasin());
                        kpos = GetMenu()->NewDish( dish );
                        nocdshs++;
                    }
                    dpos = dish->NewVectorNInd( npos );
//
                    if( GetCluster4Each() || tble == NULL ) {
                        //table for each sample
                        tble = new Table( GetDefTableSize());
                        if( tble == NULL )
                            throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
                        tble->SetBasin( GetBasin());
                        tble->SetMenu( GetMenu());
                        tble->SetDishIndex( kpos );
                        rest->NewTable( tble );
                    }
                    tble->NewVectorNInd( npos, dpos );
//
                    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                    if( length < lenstrendofsmp )
                        throw MYRUNTIME_ERROR( preamb + "No end of sample." );

                    for( c = 0; c < lenstrendofsmp; c++ )
                        if( locbuffer[c] != patstrendofsmp[c] )
                            throw MYRUNTIME_ERROR( "Invalid end of sample." );

                    IncTotalNoSamples();
                }//end of one table
//
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( length < lenstrendoftbl )
                    throw MYRUNTIME_ERROR( preamb + "No end of table." );

                for( c = 0; c < lenstrendoftbl; c++ )
                    if( locbuffer[c] != patstrendoftbl[c] )
                        throw MYRUNTIME_ERROR( preamb + "Invalid end of table." );
            }//end of one group

            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( length < lenstrendofgrp )
                throw MYRUNTIME_ERROR( preamb + "No end of group." );

            for( c = 0; c < lenstrendofgrp; c++ )
                if( locbuffer[c] != patstrendofgrp[c] )
                    throw MYRUNTIME_ERROR( preamb + "Invalid end of group." );
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;

    if(( GetCluster4Each() && nocdshs != totsamples )||
      ( !GetCluster4Each() && nocdshs != nodishes ))
        throw MYRUNTIME_ERROR( preamb + "Inconsistent number of dishes." );
}

// -------------------------------------------------------------------------
// ReadSample: read observed frequencies from file
//
bool HDPbase::ReadSample( FILE* fp, Pslvector** pfv, int dim )
{
    if( fp == NULL || pfv == NULL )
        return false;

    mystring        preamb = "HDPbase::ReadSample: ";
    myruntime_error mre;
    int             emsg, n;
    size_t          len, rbts, read = 0;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize] = {0};
    int             ctxtlen = GetCtxtSize();
    float           val;
    Pslvector*      f;

    *pfv = NULL;

    if( dim <= 0 )
        dim = GetDefSampleDim() * ctxtlen;
    f = new Pslvector();
    if( f == NULL )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory." );
    f->Allocate( dim );

    try {
        for( n = 0; n < ctxtlen; read = 0, n++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format: End of file.");

            while( emsg == 0 ) {
                //read frequency
                emsg = read_float( locbuffer + read, len - read, &val, &rbts );
                if( emsg == ERR_RD_NOVL )
                    break;
                if( emsg != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                read += rbts;

                f->Push( val );
            }
        }
        if( f->GetSize() <= 0 )
            throw MYRUNTIME_ERROR( preamb + "No sample values read." );
        *pfv = f;

    } catch( myexception const& ex ) {
        mre = ex;
        if( f ) {
            delete f;
            f = NULL;
        }
    }

    if( mre.isset())
        throw mre;

    return true;
}





// =========================================================================
// ReadParameters: read HDP parameters
//
void HDPbase::ReadParameters( const char* filename, const Ivector* dids )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL )
        throw MYRUNTIME_ERROR( mystring("HDPbase::ReadParameters: Failed to open file ") + 
                  filename );

    try {
        ReadParameters( fp, dids );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// ReadParameters: read HDP parameters given file descriptor
//
void HDPbase::ReadParameters( FILE* fp, const Ivector* dids )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        InitMenu( gcsz_defmenusz );

    ReadPriorParams( fp );
    ReadDishParams( fp, dids );
    if( GetCtxtSize() && GetCtxtSize() != GetMenu()->GetCtx())
        throw MYRUNTIME_ERROR("HDPbase::ReadParameters: Inconsistent context length.");
    SetCtxtSize( GetMenu()->GetCtx());
    //probability factors
    CalcPriorProbFact();
}





// =========================================================================
// ReadPriorParams: read prior parameters from file
//
void HDPbase::ReadPriorParams( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    const mystring  preamb = "HDPbase::ReadPriorParams: ";
    myruntime_error mre;
    char*           p;
    int             emsg;

    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrmu0 = "mu0 =";
    const char*     patstrs0m = "S0 =";
    const char*     patstrs0i = "S0^{-1} =";
    const char*     patstrds0 = "ln|S0| =";
    const char*     patstrendofdsh = "////";

    Pslvector*      f = NULL;
    SPDmatrix*      s = NULL;
    float   fpval;
    int     dim, ctx;
    int     i, j;

    try {
        if( GetMenu() == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null menu.");

        //read tau
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrtau )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrtau );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid `tau' value.");

        SetDPMTau( fpval );

        //read gamma
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrgam )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrgam );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid `gamma' value.");

        SetDPMGamma( fpval );

        //read dim
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrdim )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrdim );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &dim, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( dim < 1 || 10000 < dim )
            throw MYRUNTIME_ERROR( preamb + "Invalid `dim' value.");

        if( GetMenu()->GetDim() && dim != GetMenu()->GetDim())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent dimensions.");
        GetMenu()->SetDim( dim );

        //read ctx
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrctx )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrctx );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &ctx, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( ctx < 1 || 100 < ctx )
            throw MYRUNTIME_ERROR( preamb + "Invalid `ctx' value.");

        GetMenu()->SetCtx( ctx );

        //read prior parameter a for kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0_a )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrkp0_a );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid prior parameter `a' of kappa0.");

        GetMenu()->SetKappa0_pp_a( fpval );

        //read prior parameter b for kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0_b )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrkp0_b );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid prior parameter `b' of kappa0.");

        GetMenu()->SetKappa0_pp_b( fpval );

        //read kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0 )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrkp0 );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid `kappa0' value.");

        GetMenu()->SetKappa0( fpval );

        //read nu_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrnu0 )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrnu0 );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( fpval <= 0.0f )
            throw MYRUNTIME_ERROR( preamb + "Invalid `nu0' value.");

        GetMenu()->SetNu0( fpval );

        //read mu_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrmu0 )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        f = new Pslvector();
        if( f == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
        f->Allocate( dim );

        for( i = 0, read = 0; i < dim; i++ ) {
            if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
            read += rbts;
            f->Push( fpval );
        }

        GetMenu()->SetMu0( f );
        f = NULL;

        //read S_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrs0m )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        s = new SPDmatrix( dim );
        if( s == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        for( i = 0; i < dim; i++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            for( j = 0, read = 0; j < dim; j++ ) {
                if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                read += rbts;
                s->SetValueAt( i, j, fpval );
            }
        }

        GetMenu()->SetS0( s );
        s = NULL;

        //read S_0^{-1}
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrs0i )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        s = new SPDmatrix( dim );
        if( s == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

        for( i = 0; i < dim; i++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            for( j = 0, read = 0; j < dim; j++ ) {
                if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                read += rbts;
                s->SetValueAt( i, j, fpval );
            }
        }

        GetMenu()->SetInvS0( s );
        s = NULL;

        //read ln|S_0|
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrds0 )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrds0 );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        GetMenu()->SetLDetS0( fpval );

        //read footer
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrendofdsh )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( f ) { delete f; f = NULL; }
    if( s ) { delete s; s = NULL; }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// PrintPriorParams: print prior parameters to file
//
void HDPbase::PrintPriorParams( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR("HDPbase::PrintPriorParams: Null Menu.");

    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrmu0 = "mu0 =";
    const char*     patstrs0m = "S0 =";
    const char*     patstrs0i = "S0^{-1} =";
    const char*     patstrds0 = "ln|S0| =";
    const char*     patstrendofdsh = "////";

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
#ifdef OS_MS_WINDOWS
	localtime_s(&t_ttm, &t_tt);
	asctime_s(t_ttma, KBYTE-1, &t_ttm);
#else
    localtime_r( &t_tt, &t_ttm );
    asctime_r( &t_ttm, t_ttma );
#endif

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## -- Concentration parameters --%s##%s", NL, NL );
    fprintf( fp, "%s %.21g%s", patstrtau, GetDPMTau(), NL );
    fprintf( fp, "%s %.21g%s%s", patstrgam, GetDPMGamma(), NL, NL );

    fprintf( fp, "## -- Prior parameters --%s##%s", NL, NL );
    fprintf( fp, "%s %d%s%s %d%s%s %.5g%s%s %.5g%s%s %.21g%s%s %.21g%s", 
          patstrdim, GetMenu()->GetDim(), NL, patstrctx, GetMenu()->GetCtx(), NL, 
          patstrkp0_a, GetMenu()->GetKappa0_pp_a(), NL, patstrkp0_b, GetMenu()->GetKappa0_pp_b(), NL, 
          patstrkp0, GetMenu()->GetKappa0(), NL, patstrnu0, GetMenu()->GetNu0(), NL );
    fprintf( fp, "%s%s", patstrmu0, NL );
    if( GetMenu()->GetMu0())
        GetMenu()->GetMu0()->Print( fp, " %.21g");
    fprintf( fp, "%s%s", patstrs0m, NL );
    if( GetMenu()->GetS0())
        GetMenu()->GetS0()->Print( fp, " %.21g");
    fprintf( fp, "%s%s", patstrs0i, NL );
    if( GetMenu()->GetInvS0())
        GetMenu()->GetInvS0()->Print( fp, " %.21g");
    fprintf( fp, "%s %.21g%s", patstrds0, GetMenu()->GetLDetS0(), NL );
    fprintf( fp, "%s%s%s", patstrendofdsh, NL, NL );
}





// =========================================================================
// ReadDishParams: read dish parameters from file
//
void HDPbase::ReadDishParams( FILE* fp, const Ivector* dids )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    const mystring  preamb = "HDPbase::ReadDishParams: ";
    myruntime_error mre;
    char*           p;
    int             emsg;

    const char*     patstrgps = "Number of groups =";
    const char*     patstrnod = "Number of dishes =";
    const char*     patstrdsh = "Dish";
    const char*     patstrdid = "Id=";
    const char*     patstrtm = "time=";
    const char*     patstrnk = "nk =";
    const char*     patstrmk = "mk =";
    const char*     patstrlp = "ln p =";
    const char*     patstrmu = "mu =";
    const char*     patstrsm = "S =";
#ifdef PROBMTXINVSM
    const char*     patstrsi = "S^{-1} =";
#endif
    const char*     patstrds = "ln|S| =";
    const char*     patstrendofdsh = "////";
    Dish*           dish;
    Pslvector*      f = NULL;
    SPDmatrix*      s = NULL;
    long long int   time;
    float   fpval;
    int     intval;
    int     nodshs, nk, mk, dim, ctx;
    int     i, j, d, dishid;

    try {
        if( GetMenu() == NULL )
            throw MYRUNTIME_ERROR( preamb + "Null Menu.");

        if( dids == NULL && ( GetMenu()->GetActualSize() || GetMenu()->GetSize()))
            throw MYRUNTIME_ERROR( preamb + "Menu configuration already exists.");

        dim = GetMenu()->GetDim();
        ctx = GetMenu()->GetCtx();

        if( dim < 1 || ctx < 1 )
            throw MYRUNTIME_ERROR( preamb + "Invalid dimensions and/or context length.");

        //read Number of groups
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrgps )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrgps );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( intval < 0 )
            throw MYRUNTIME_ERROR( preamb + "Invalid number of groups.");

        SetReadNoGroups( intval );

        //read Number of dishes
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( p = strstr( locbuffer, patstrnod )) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        p += strlen( patstrnod );

        if( length <= (size_t)( p - locbuffer ))
            throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

        if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        if( intval < 0 )
            throw MYRUNTIME_ERROR( preamb + "Invalid number of dishes.");

        if( dids && GetMenu()->GetActualSize() != intval )
            throw MYRUNTIME_ERROR( preamb + "Inconsistent number of dishes.");

        nodshs = intval;

        for( d = 0; d < nodshs; d++ )
        {
            //read Dish index
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrdsh )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrdsh );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &intval, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( intval != d )
                throw MYRUNTIME_ERROR( preamb + "Invalid dish index.");

            //read dish id
            if(( p = strstr( locbuffer, patstrdid )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrdid );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = 
                read_integer( p, length - (size_t)( p - locbuffer ), &dishid, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( dishid < 0 )
                throw MYRUNTIME_ERROR( preamb + "Invalid Dish id.");

            if( dids ) {
                if( dids->GetSize() <= dishid )
                    throw MYRUNTIME_ERROR( preamb + "Wrong dish id.");
                dishid = dids->GetValueAt( dishid );
                if(( dish = GetMenu()->GetDishAt( dishid )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Wrong dish id.");
            }
            else {
                dish = new Dish( HDPbase::GetDefDishSize());
                if( dish == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
                dish->SetBasin( NULL );//no basin
                dishid = GetMenu()->NewDish( dish );
                if( dishid != d )
                    throw MYRUNTIME_ERROR( preamb + "Invalid index of new dish obtained.");
            }

            //read dish time
            if(( p = strstr( locbuffer, patstrtm )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrtm );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = 
                read_llinteger( p, length - (size_t)( p - locbuffer ), &time, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( time < 0 )
                throw MYRUNTIME_ERROR( preamb + "Invalid Dish time.");

            dish->SetTime(( time_t )time );

            //read nk
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrnk )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrnk );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &nk, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( nk < 1 )
                throw MYRUNTIME_ERROR( preamb + "Invalid `nk' value.");

            if( dids && nk != dish->GetActualSize())
                throw MYRUNTIME_ERROR( preamb + "Wrong `nk' value.");

            dish->SetReadSize( nk );

            //read mk
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrmk )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrmk );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = read_integer( p, length - (size_t)( p - locbuffer ), &mk, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( mk < 1 || nk < mk )
                throw MYRUNTIME_ERROR( preamb + "Invalid `mk' value.");

            dish->SetReadNoTables( mk );

            //read log probability
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrlp )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrlp );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            GetMenu()->SetProbAt( dishid, fpval );

            //read mu
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrmu )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            f = new Pslvector();
            if( f == NULL )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
            f->Allocate( dim );

            for( i = 0, read = 0; i < dim; i++ ) {
                if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                read += rbts;
                f->Push( fpval );
            }

            GetMenu()->SetMuVectorAt( dishid, f );//set mu
            f = NULL;

            //read S_0
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrsm )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            s = new SPDmatrix( dim );
            if( s == NULL )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

            for( i = 0; i < dim; i++ ) {
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( feof( fp ) || !length )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

                for( j = 0, read = 0; j < dim; j++ ) {
                    if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                    read += rbts;
                    s->SetValueAt( i, j, fpval );
                }
            }

            GetMenu()->SetSMatrixAt( dishid, s );
            s = NULL;

#ifdef PROBMTXINVSM
            //read S^{-1}
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrsi )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            s = new SPDmatrix( dim );
            if( s == NULL )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

            for( i = 0; i < dim; i++ ) {
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

                if( feof( fp ) || !length )
                    throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

                for( j = 0, read = 0; j < dim; j++ ) {
                    if(( emsg = read_float( locbuffer + read, length - read, &fpval, &rbts )) != 0 )
                        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                    read += rbts;
                    s->SetValueAt( i, j, fpval );
                }
            }

            GetMenu()->SetInvSMatrixAt( dishid, s );
            s = NULL;
#endif

            //read ln|S|
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrds )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            p += strlen( patstrds );

            if( length <= (size_t)( p - locbuffer ))
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( emsg = read_float( p, length - (size_t)( p - locbuffer ), &fpval, &rbts )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            GetMenu()->SetLDetSMAt( dishid, fpval );

            //read footer
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

            if(( p = strstr( locbuffer, patstrendofdsh )) == NULL )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format.");
        }

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( f ) { delete f; f = NULL; }
    if( s ) { delete s; s = NULL; }

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// PrintDishParams: print dish parameters to file
//
void HDPbase::PrintDishParams( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        throw MYRUNTIME_ERROR("HDPbase::PrintDishParams: Null Menu.");
    if( GetChain() == NULL )
        throw MYRUNTIME_ERROR("HDPbase::PrintDishParams: Null Chain.");

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgps = "Number of groups =";
    const char*     patstrnod = "Number of dishes =";
    const char*     patstrdsh = "Dish";
    const char*     patstrdid = "Id=";
    const char*     patstrtm = "time=";
    const char*     patstrnk = "nk =";
    const char*     patstrmk = "mk =";
    const char*     patstrlp = "ln p =";
    const char*     patstrmu = "mu =";
    const char*     patstrsm = "S =";
#ifdef PROBMTXINVSM
    const char*     patstrsi = "S^{-1} =";
#endif
    const char*     patstrds = "ln|S| =";
    const char*     patstrendofdsh = "////";
    const Dish*     dish;
    int c, d;

    fprintf( fp, "## -- Dish parameters --%s##%s", NL, NL );
    fprintf( fp, "## %s %g%s", patstrlpd, GetLastLProbData(), NL );//TODO:get last probability
    fprintf( fp, "%s %d%s", patstrgps, GetChain()->GetSize(), NL );
    fprintf( fp, "%s %d%s", patstrnod, GetMenu()->GetActualSize(), NL );
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s%d %s%lu)%s", patstrdsh, c++, patstrdid, d, 
                 patstrtm, ( unsigned long int )dish->GetTime(), NL );
        fprintf( fp, "%s %d%s", patstrnk, dish->GetActualSize(), NL );
        fprintf( fp, "%s %d%s", patstrmk, dish->GetNoTables(), NL );
        fprintf( fp, "%s %g%s", patstrlp, GetMenu()->GetProbAt(d), NL );
        fprintf( fp, "%s%s", patstrmu, NL );
        if( GetMenu()->GetMuVectorAt( d ))
            GetMenu()->GetMuVectorAt( d )->Print( fp, " %.21g");
        fprintf( fp, "%s%s", patstrsm, NL );
        if( GetMenu()->GetSMatrixAt( d ))
            GetMenu()->GetSMatrixAt( d )->Print( fp, " %.21g");
#ifdef PROBMTXINVSM
        fprintf( fp, "%s%s", patstrsi, NL );
        if( GetMenu()->GetInvSMatrixAt( d ))
            GetMenu()->GetInvSMatrixAt( d )->Print( fp, " %.21g");
#endif
        fprintf( fp, "%s %.21g%s", patstrds, GetMenu()->GetLDetSMAt( d ), NL );
        fprintf( fp, "%s%s%s", patstrendofdsh, NL, NL );
    }
}

// =========================================================================




// =========================================================================
// PrintGroups: print HDP structure
//
void HDPbase::PrintGroups( const char* filename, float* lprob, int* iter )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw MYRUNTIME_ERROR( "HDPbase::PrintGroups: Failed to open output file." );

    try {
        PrintGroups( fp, lprob, iter );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// Print: print HDP structure to file
//
void HDPbase::PrintGroups( FILE* fp, float* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetChain() == NULL || GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
    const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";
    const int       dim = GetMenu()->GetDim();
    const int       ctx = GetCtxtSize();
    const int       nvals = dim / ctx;
    const Restaurant*   rest;
    const Table*        tble;
    const Pslvector*    frq;
    int r, tt, t, n, v;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
#ifdef OS_MS_WINDOWS
	localtime_s(&t_ttm, &t_tt);
	asctime_s(t_ttma, KBYTE - 1, &t_ttm);
#else
	localtime_r(&t_tt, &t_ttm);
	asctime_r(&t_ttm, t_ttma);
#endif

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## " ); 
    print_cmdline( &file_print, fp );
    fprintf( fp, "##%s", NL );

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "## MH update %d-%d%s", GetGibbsIt(), GetMHUpdateNum(), NL );
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "## Gibbs sampling update %d%s", GetGibbsIt(), NL );
    else
        fprintf( fp, "## N/A update%s", NL );

    if( iter )
        fprintf( fp, "%s %d%s", patstrgss, *iter, NL );
    if( lprob )
        fprintf( fp, "%s %g%s", patstrlpd, *lprob, NL );
    fprintf( fp, "##%s", NL ); 
    fprintf( fp, "%s %d; %s %d; %s %d; %s %d%s", 
              patstrnog, GetChain()->GetSize(), 
              patstrnds, GetMenu()->GetActualSize(),
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize(), NL );
    for( r = 0; r < GetChain()->GetSize(); r++ ) {
        rest = GetChain()->GetRestaurantAt( r );
        if( rest == NULL )
            continue;
        fprintf( fp, "%s %d; %s %d%s", patstrgrp, r, patstrnot, rest->GetActualSize(), NL );
        for( tt = t = 0; t < rest->GetSize(); t++ )
        {
            tble = rest->GetTableAt( t );
            if( tble == NULL )
                continue;
            fprintf( fp, "%s %d (%s %d; %s %d); %s %d%s", 
                     patstrtbl, tt++, patstrtid, t, patstrdid, tble->GetDishIndex(),
                     patstrnos, tble->GetActualSize(), NL );
            for( n = 0; n < tble->GetSize(); n++ )
            {
                if( tble->GetVectorNIndAt( n ) < 0 )
                    continue;
                frq = tble->GetVectorNAt( n );
                if( frq == NULL )
                    throw MYRUNTIME_ERROR("HDPbase::PrintGroups: Null vector.");
                for( v = 0; v < frq->GetSize(); v++ ) {
                    fprintf( fp, " %12.6g", frq->GetValueAt( v ));
                    if(( v+1 ) % nvals == 0 && ( v+1 ) < frq->GetSize())
                        fprintf( fp, "%s", NL );
                }
                fprintf( fp, "%s%s%s", NL, patstrendofsmp, NL );
            }
            fprintf( fp, "%s%s", patstrendoftbl, NL );//end of table
        }
        fprintf( fp, "%s%s%s", patstrendofgrp, NL, NL );
    }
}

// =========================================================================
// PrintDishes: print HDP dish structure
//
void HDPbase::PrintDishes( const char* filename, float* lprob, int* iter )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw MYRUNTIME_ERROR( mystring("HDPbase::PrintDishes: Failed to open file ") +
                filename );
    try {
        PrintDishes( fp, lprob, iter );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// Print: print HDP structure to file
//
void HDPbase::PrintDishes( FILE* fp, float* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of clusters =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Cluster";
    const char*     patstrcid = "Id=";
    const char*     patstrtm  = "time=";
    const char*     patstrnot = "Number of samples =";
    const char*     patstrlp  = "Log Prob =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendofgrp = "////";
    const int       dim = GetMenu()->GetDim();
    const int       ctx = GetCtxtSize();
    const int       nvals = dim / ctx;
    const Dish*         dish;
    const Pslvector*    frq;
    int c, d, v, a;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
#ifdef OS_MS_WINDOWS
	localtime_s(&t_ttm, &t_tt);
	asctime_s(t_ttma, KBYTE - 1, &t_ttm);
#else
	localtime_r(&t_tt, &t_ttm);
	asctime_r(&t_ttm, t_ttma);
#endif

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## " ); 
    print_cmdline( &file_print, fp );
    fprintf( fp, "##%s", NL ); 

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "## MH update %d-%d%s", GetGibbsIt(), GetMHUpdateNum(), NL );
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "## Gibbs sampling update %d%s", GetGibbsIt(), NL );
    else
        fprintf( fp, "## N/A update%s", NL );

    if( iter )
        fprintf( fp, "%s %d%s", patstrgss, *iter, NL );
    if( lprob )
        fprintf( fp, "%s %g%s", patstrlpd, *lprob, NL );
    fprintf( fp, "##%s", NL ); 
    fprintf( fp, "%s %d; %s %d; %s %d%s", 
              patstrnog, GetMenu()->GetActualSize(), 
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize(), NL );
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s%d %s%lu); %s %d; %s %g%s", 
                 patstrgrp, c++, patstrcid, d, 
                 patstrtm, ( unsigned long int )dish->GetTime(), 
                 patstrnot, dish->GetActualSize(), 
                 patstrlp, GetMenu()->GetProbAt(d), NL );
        for( v = 0; v < dish->GetSize(); v++ )
        {
            if( dish->GetVectorNIndAt( v ) < 0 )
                continue;
            frq = dish->GetVectorNAt( v );
            if( frq == NULL )
                throw MYRUNTIME_ERROR("HDPbase::PrintDishes: Null vector.");
            for( a = 0; a < frq->GetSize(); a++ ) {
                fprintf( fp, " %12.6g", frq->GetValueAt( a ));
                if(( a+1 ) % nvals == 0 && ( a+1 ) < frq->GetSize())
                    fprintf( fp, "%s", NL );
            }
            fprintf( fp, "%s%s%s", NL, patstrendofsmp, NL );
        }
        fprintf( fp, "%s%s%s", patstrendofgrp, NL, NL );
    }
}

// -------------------------------------------------------------------------
// PrintParameters: print parameters of the HDP structure
//
void HDPbase::PrintParameters( const char* filename )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw MYRUNTIME_ERROR( mystring("HDPbase::PrintParameters: Failed to open file ") +
                filename );
    try {
        PrintParameters( fp );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// PrintParameters: print parameters of the HDP structure
//
void HDPbase::PrintParameters( FILE* fp )
{
    if( fp == NULL )
        return;
    PrintPriorParams( fp );
    PrintDishParams( fp );
}





// =========================================================================
// PrintSummary: print summary information at the end of file
//
void HDPbase::PrintSummary( const char* filename, float* lprob, int* iter )
{
    FILE* fp = NULL;
    myruntime_error mre;

    if( !filename )
        return;

    if(( fp = fopen( filename, "a" )) == NULL )
        throw MYRUNTIME_ERROR("HDPbase::PrintSummary: Failed to open output report file.");

    try {
        PrintSummary( fp, lprob, iter );
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
//
void HDPbase::PrintSummary( FILE* fp, float* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrnog = "Number of clusters =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrgrp = "Cluster";
    const char*     patstrcid = "Id =";
    const char*     patstrnot = "nk =";
    const char*     patstrmk  = "mk =";
    const char*     patstrlp  = "ln p =";
    const char*     patstrtm  = "time =";
    const char*     patstrendofgrp = "////";
    const Dish*     dish;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];
    int c, d;

    time( &t_tt );
#ifdef OS_MS_WINDOWS
	localtime_s(&t_ttm, &t_tt);
	asctime_s(t_ttma, KBYTE - 1, &t_ttm);
#else
	localtime_r(&t_tt, &t_ttm);
	asctime_r(&t_ttm, t_ttma);
#endif

    fprintf( fp, "## %s", t_ttma );

    if( iter )
        fprintf( fp, "%s %d %s%s", patstrgss, *iter, GetRestarted()? "restart": "", NL );

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "MH update %d-%d%s", GetGibbsIt(), GetMHUpdateNum(), NL );
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "Gibbs sampling update %d%s", GetGibbsIt(), NL );
    else
        fprintf( fp, "N/A update%s", NL );

    if( lprob )
        fprintf( fp, "%s %g%s", patstrlpd, *lprob, NL );
    fprintf( fp, "%s %.5g%s", patstrtau, GetDPMTau(), NL );
    fprintf( fp, "%s %.5g%s", patstrgam, GetDPMGamma(), NL );
    fprintf( fp, "%s %.5g%s%s %.5g%s%s %.5g%s%s %.5g%s%s %d%s%s %d%s", 
          patstrkp0_a, GetMenu()->GetKappa0_pp_a(), NL, patstrkp0_b, GetMenu()->GetKappa0_pp_b(), NL, 
          patstrkp0, GetMenu()->GetKappa0(), NL, patstrnu0, GetMenu()->GetNu0(), NL, 
          patstrdim, GetMenu()->GetDim(), NL, patstrctx, GetMenu()->GetCtx(), NL );
    fprintf( fp, "%s %d; %s %d%s", 
              patstrnog, GetMenu()->GetActualSize(), 
              patstrtns, GetTotalNoSamples(), NL );
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s %d); %s %d; %s %d; %s %g; %s %lu%s", 
              patstrgrp, c++, patstrcid, d, 
              patstrnot, dish->GetActualSize(), patstrmk, dish->GetNoTables(), 
              patstrlp, GetMenu()->GetProbAt(d), patstrtm, ( unsigned long int )dish->GetTime(), NL );
    }
    fprintf( fp, "%s%s%s", patstrendofgrp, NL, NL );
}





// =========================================================================
// PrintBasin: print basin of samples
//
void HDPbase::PrintBasin( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetBasin() == NULL )
        return;

    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Sample size =";
    const char*     patstrendofsmp = "*";
    const Pslvector* frq;
    int n, v;

    fprintf( fp, "%s %d; %s %d%s",
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize(), NL );
    for( n = 0; n < GetBasin()->GetSize(); n++ ) {
        frq = GetBasin()->GetValueAt( n );
        if( frq == NULL )
            continue;
        for( v = 0; v < frq->GetSize(); v++ ) {
            fprintf( fp, " %12.6g", frq->GetValueAt( v ));
        }
        fprintf( fp, "%s%s%s", NL, patstrendofsmp, NL );
    }
}
