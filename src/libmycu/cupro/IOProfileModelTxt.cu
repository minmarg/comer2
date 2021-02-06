/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "liblib/myfiler.h"

#include <stdio.h>
#include <string.h>

#include "extsp/psl.h"
#include "extsp/pslvector.h"
#include "liblib/msg.h"
#include "liblib/alpha.h"
#include "libpro/srcpro/TRANSPROBS.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PM2DVector.cuh"
#include "libmycu/cudp/CuBatchDP_com.h"
#include "IOProfileModel.h"

// //note: constant arrays in host are ok; 
// //note: however, device array constants are located in global memory; better solutions follow;
// const int iopm_infcode[pmvicN] = {pmvfullInfCode,pmvnossInfCode,pmvnocvInfCode,pmvnosscvInfCode};
// const int iopm_trgfrqs[pmvicN] = {pmvfullTrgFrqs,pmvnossTrgFrqs,pmvnocvTrgFrqs,pmvnosscvTrgFrqs};
// const int iopm_cventrs[pmvicN] = {pmvfullCVentrs,pmvnossCVentrs,pmvnocvCVentrs,pmvnosscvCVentrs};
// const int iopm_cvprior[pmvicN] = {pmvfullCVprior,pmvnossCVprior,pmvnocvCVprior,pmvnosscvCVprior};
// const int iopm_cvnorm2[pmvicN] = {pmvfullCVnorm2,pmvnossCVnorm2,pmvnocvCVnorm2,pmvnosscvCVnorm2};
// const int iopm_sssprbs[pmvicN] = {pmvfullSSsprbs,pmvnossSSsprbs,pmvnocvSSsprbs,pmvnosscvSSsprbs};
// const int iopm_hdp1prb[pmvicN] = {pmvfullHDP1prb,pmvnossHDP1prb,pmvnocvHDP1prb,pmvnosscvHDP1prb};
// const int iopm_hdp1ind[pmvicN] = {pmvfullHDP1ind,pmvnossHDP1ind,pmvnocvHDP1ind,pmvnosscvHDP1ind};
// const int iopm_addrpro[pmvicN] = {pmvfullAddrPro,pmvnossAddrPro,pmvnocvAddrPro,pmvnosscvAddrPro};
// const int iopm_ssstate[pmvicN] = {pmvfullSSstate,pmvnossSSstate,pmvnocvSSstate,pmvnosscvSSstate};
// // const int iopm_total[pmvicN] = {pmvfullTotal,pmvnossTotal,pmvnocvTotal,pmvnosscvTotal};
// const int iopm_end[pmvicN] = {pmvfullEnd,pmvnossEnd,pmvnocvEnd,pmvnosscvEnd};

// -------------------------------------------------------------------------
// TextReadProfileHeader: read profile header data from file;
// desc, profile description; required for host;
// file, profile filename; required for host;
// prolen, profile length;
// scale, scale factor read from file;
// pmdata, profile-specific data required for both host and device;
// NOTE: on return, desc, file, and pmdata do not point to the end of data;
// 
template<typename T>
bool TextReadProfileHeader(T* fp, 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata)
{
    MYMSG( "TextReadProfileHeader", 5 );
    const mystring preamb = "TextReadProfileHeader: ";

    if( !fp )
        throw MYRUNTIME_ERROR( preamb + "Null file descriptor.");
    if( !prolen || !scale || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");

    const float accuracy = 1.0e-3f;
    const int noress = NUMAA;
    mystring buffer;
    const char* p;
    size_t rbts;
    int emsg;
    const char* descp, *filep;
    float effnos = 0.0f;
    int intval;
    float value, consv;
    int r, n = 0;

    PM2DVector pm2dvec( pmdata );

    //read version number
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if(feof(fp) && buffer.empty())
        return false;

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( buffer.c_str(), patstrDATVER )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrDATVER );

    if( buffer.length() <= (size_t)(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( p, pmodel::PMProfileModel::dataversion )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Inappropriate profile version number." );


    //read description line
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( buffer.c_str(), patstrDESC )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    descp = p + strlen( patstrDESC );
    for( ; *descp == ' ' || *descp == '\t'; descp++ );
    for( n = (int)buffer.length() - 1; 0 <= n && ( buffer[n]=='\n' || buffer[n]=='\r' ); n-- )
        buffer[n] = 0;
    *desc++ = descp;


    //read filename
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( buffer.c_str(), patstrFILE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    filep = p + strlen( patstrFILE );
    for( ; *filep == ' ' || *filep == '\t' ; filep++ );
    for( n = (int)buffer.length() - 1; 0 <= n && ( buffer[n]=='\n' || buffer[n]=='\r' ); n-- )
        buffer[n] = 0;
    *file++ = filep;


    //read command line without analysis
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    //read profile length
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( buffer.c_str(), patstrLEN )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrLEN );

    if( buffer.length() <= (size_t)(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), prolen, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( *prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length." );

    if( MAXCOLUMNS < *prolen )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Too large profile length." );


    //read number of sequences w/o analysis
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    //read effective number of sequences
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( buffer.c_str(), patstrEFFNOS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrEFFNOS );

    if( buffer.length() <= (size_t)(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_float( p, buffer.length() - (size_t)(p-buffer.c_str()), &effnos, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( effnos <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid effective number of observations." );


    //read scale factor
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( buffer.c_str(), patstrSCALE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrSCALE );

    if( buffer.length() <= (size_t)(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), scale, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( *scale < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale factor." );


    //read residue letters; omit checking the order
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));


    //{{read background probabilities
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( buffer.c_str(), patstrNULL )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No background probabilities." );

    p += strlen( patstrNULL );

    consv = 0.0f;
    for( r = 0; r < noress; r++ ) {
        if( buffer.length() <= (size_t)(p-buffer.c_str()))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete background probability data." );

        if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || *scale < intval )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Invalid background probability value." );

        value = (float)intval / (float)*scale;
        consv += value;

        //NOTE: SET profile-specific data
        pm2dvec.SetField( pps2DBkgPrbs + r, 0, (FPTYPE)value );
    }

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid profile data: Invalid background probabilities." );

    //NOTE: SET profile-specific data
    pm2dvec.SetField( pps2DENO, 0, (FPTYPE)effnos );
    pm2dvec.SetField( pps2DLen, 0, (INTYPE)*prolen );
    //}}


    //read posterior (generalized target) probabilities, but do not analyze them
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    return true;
}

// -------------------------------------------------------------------------
// TextReadProfileData: read profile data from file;
// distn, distance in positions to profile to write to buffer;
// profnr, profile serial number;
// prolen, profile length;
// scale, scale factor for profile data;
// pmorgaddr, 2D vector of original addresses for a batch of profile models,
//          used to calculate a shift for CV and SS data from the beginning;
// pmdata is a 2D vector representing a batch of profile models, each of 
// which includes: 
// seqn, amino acid sequence the profile has been coonstructed for; 
//       required for host;
// prodata, position-specific profile model data required for device and 
//          host;
// ptrdata, position-specific profile model transition data required for 
//          device and host;
// NOTE: on return, the individual vectors of pmdata point to the end of 
// data;
//
template<typename T>
void TextReadProfileData(T* fp, 
    unsigned int distn,
    int profnr,
    int prolen,
    int scale, 
    char** pmorgaddr, 
    char** pmdata )
{
    MYMSG( "TextReadProfileData", 5 );
    const mystring preamb = "TextReadProfileData: ";

    if( !fp )
        throw MYRUNTIME_ERROR( preamb + "Null file descriptor.");
    if( !pmorgaddr || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");

    //const float accuracy = 1.0e-3f;
    const int maxnocls = 1000;
    const int noress = NUMAA;
    char msgbuf[BUF_MAX];
    mystring buffer;
    const char* p;
    size_t rbts;
    bool lineread;
    int emsg;
    int intval;
    char res, sss;
    float value/*, consv*/;
    size_t noppps;
    float HDP1prb;
    int HDP1ind;
    bool HDPctset = false;//HDP ctx data's present
    bool CtxVecset = false;//CV's present
    extspsl::Pslvector ctxvec;//context vector
    float lpprb = 0.0f, norm2 = 0.0f;
    int vsize = 0;//cv size
    bool SSSset = false;//SS information's present
    bool SSSP3 = false;//SS information format
    float sssprob[SS_NSTATES];
    int m, r, tt, t, t0 = 0, tn = 0, n = 0;

    PM2DVector pm2dvec( pmdata );

    //{{beginning transition probabilities
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p = buffer.c_str();

    for( t = tt = 0; t < P_NSTATES; t++ )
    {
        if( buffer.length() <= (size_t)(p-buffer.c_str()))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete transition probability data." );

        if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        value = (float)intval / (float)scale;

        //NOTE: SET the beginning transitions of profile model and 
        // MOVE the pointer to the next position
//         tt = t;
//         if( P_DI < t )
//             tt -= 2;
//         else if( P_ID < t && t < P_DI )
//             tt -= 1;
//         else if( P_ID <= t )
//             continue;
        if( P_ID == t || t == P_DI )
            continue;
        value = value? (P_MD<t? logf(value)*MAIN_TRNPRB_LOGFCT: logf(value)): -32768.0f;
        pm2dvec.SetFieldNext( ptr2DTrnPrbs + tt, 0, (FPTYPE)value );
//         pm2dvec.SetFieldNext( ptr2DTrnPrbsExp + tt, 0, (FPTYPE)((value<SLC_LOG_SP_MIN)? 0.0f: expf(value)));
        tt++;
    }

    //NOTE: SET profile-specific data
    pm2dvec.SetFieldNext( pps2DDist, 0, (LNTYPE)distn );
    //}}


    //beginning MID state expected observations w/o analysis
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    lineread = false;

    for( m = 0; m < prolen; m++ )
    {
        //target probabilities
        if( !lineread ) {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }
        }

        lineread = false;

        if(( emsg = read_integer( p = buffer.c_str(), buffer.length(), &n, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( n != m + 1 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: Wrong numbering.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p += rbts;

        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
            sprintf( msgbuf, "Wrong profile format at pos %d: No residue.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_symbol( p, buffer.length() - (size_t)(p-buffer.c_str()), &res, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        p += rbts;

        //NOTE: SET amino acid and MOVE to the next position
        //NOTE:symbol is not hashed!
        //res = HashAlphSymbol( res );
        pm2dvec.SetFieldNext( pmv2Daa, 0, (CHTYPE)res );

        for( r = 0; r < noress; r++ )
        {
            if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete target probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            value = (float)intval / (float)scale;

            //NOTE: SET profile model position-specific data and 
            // MOVE the pointer to the next position; 
            //target frequencies are always present and occupy the same address in the vector
            pm2dvec.SetFieldNext( pmv2DTrgFrqs + r, 0, (FPTYPE)value );
        }


        //observed frequencies w/o analysis
        //NOTE: condensed format does not include this section
//         if(( emsg = skip_comments( fp, buffer )) != 0 ) {
//             sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
//         }
// 
//         if( feof( fp ) || buffer.empty()) {
//             sprintf( msgbuf, "Wrong profile format at pos %d.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }


        //transition probabilities
        if(( emsg = skip_comments( fp, buffer )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || buffer.empty()) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p = buffer.c_str();

        for( t = tt = 0; t < P_NSTATES; t++ )
        {
            if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete transition probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            value = (float)intval / (float)scale;

            //NOTE: SET the transitions of profile model at a given position and
            // MOVE to the next
//             tt = t;
//             if( P_DI < t )
//                 tt -= 2;
//             else if( P_ID < t && t < P_DI )
//                 tt -= 1;
//             else if( P_ID <= t )
//                 continue;
            if( P_ID == t || t == P_DI )
                continue;
            value = value? (P_MD<t? logf(value)*MAIN_TRNPRB_LOGFCT: logf(value)): -32768.0f;
            pm2dvec.SetFieldNext( ptr2DTrnPrbs + tt, 0, (FPTYPE)value );
//             pm2dvec.SetFieldNext( ptr2DTrnPrbsExp + tt, 0, (FPTYPE)((value<SLC_LOG_SP_MIN)? 0.0f: expf(value)));
            tt++;
        }


        //MID state observations w/o analysis
        if(( emsg = skip_comments( fp, buffer )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || buffer.empty()) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }


        //{{CLUSTERS: HDP1
        //bck posterior probability, no. posterior probability values
        if(( emsg = skip_comments( fp, buffer )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || buffer.empty()) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p = buffer.c_str();

        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        //bppprob = (float)intval / (float)scale;

        p += rbts;

        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }
        if( maxnocls < intval ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Too large number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        noppps = (size_t)intval;

        if( noppps < 1 ) {
            HDP1prb = 0.0f;
            HDP1ind = -1;
        }
        else {
            //posterior predictive probabilities
            if(( emsg = skip_comments( fp, buffer )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p = buffer.c_str();

            //read only the first probability value
            if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Incomplete HDP1 posterior probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( intval < 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Negative HDP1 posterior probability value.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p += rbts;

            HDP1prb = (float)intval / (float)scale;

            //indices of posteriors
            if(( emsg = skip_comments( fp, buffer )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                  "No HDP1 cluster indices.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p = buffer.c_str();

            //read only the index of the first cluster
            if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Incomplete HDP1 cluster index data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( intval != -1 && intval < 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Negative HDP1 cluster index.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p += rbts;

            HDP1ind = intval;
        }//if( noppps )
        //}}HDP1


        //{{HDP ctx
        if( m < 1 )
            HDPctset = true;//probe for HDP ctx information
        if( HDPctset ) {
            //bck posterior probability, no. posterior probability values
            if(( emsg = skip_comments( fp, buffer )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( buffer.c_str(), patstrCT )) == NULL ) {
                if( m < 1 )
                    HDPctset = false;//no HDP ctx information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No HDP ctx probabilities.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCT );
                //do not read and check data futher in this line dedicated to HDP ctx information
            }//else
        }//if( HDPctset )
        //}}HDP ctx


        //{{Context vector
        if( m < 1 )
            CtxVecset = true;//probe for context vector data
        //Initialize even in the case of no use of CV scoring;
        //init each entry to ensure CVS2S score==0 when CVS information is not to be in use
        ctxvec.Reserve(pmv2DNoCVEls);
        ctxvec.AssignTo(0.17302947f);
        lpprb = 0.0f;
        norm2 = 0.0f;
        if( CtxVecset ) {
            if( !lineread )
                if(( emsg = skip_comments( fp, buffer )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( buffer.c_str(), patstrCV )) == NULL ) {
                if( m < 1 )
                    CtxVecset = false;//no context vector data
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCV );

                if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No prior for context vector.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                lpprb = (float)intval / (float)scale;

                p += rbts;

                if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                norm2 = (float)intval / (float)scale;

                p += rbts;

                if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( 1000 < intval ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Too large context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( !vsize && intval != pmv2DNoCVEls ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Only context vectors of size %d are supported.", m, pmv2DNoCVEls );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( vsize && vsize != intval && vsize != pmv2DNoCVEls ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Inconsistent context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                vsize = intval;

                p += rbts;

                if( vsize ) {
                    ctxvec.Allocate( vsize );
                    ctxvec.Clear();

                    //vector elements
                    for( t = 0; t < vsize; t++ )
                    {
                        if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "Incomplete context vector data.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        if(( emsg = read_integer( p, 
                                    buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                        }

                        p += rbts;

                        ctxvec.Push( (float)intval / (float)scale );
                    }
                }//if( vsize )
            }//else
        }//if( CtxVecset )
        //}}Context vector


        //{{SS state
        if( m < 1 ) {
            SSSset = true;//probe SS information existance
            SSSP3 = true;
            t0 = 0; tn = SS_NSTATES-1;
        }
        //Initialize even in the case of no use of SS information
        memset( sssprob, 0, SS_NSTATES*sizeof(float));
        sss = 0;
        if( SSSset ) {
            if( !lineread )
                if(( emsg = skip_comments( fp, buffer )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || buffer.empty()) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( buffer.c_str(), patstrSS )) == NULL ) {
                if( m < 1 )
                    SSSset = false;//no SS information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: No SS state.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrSS );

                if(( emsg = read_symbol( p, buffer.length() - (size_t)(p-buffer.c_str()), &sss, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                p += rbts;

                sss = HashSSState( sss );

                if( !SSSP3 ) {
                    t0 = sss; tn = sss;
                }

                for( t = t0; t <= tn; t++ )
                {
                    if( buffer.length() <= (size_t)(p-buffer.c_str())) {
                        if( m < 1 && t == t0+1 ) {
                            SSSP3 = false;
                            sssprob[(int)sss] = (float)intval / (float)scale;
                            break;
                        }
                        else {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "No SS state probability.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }
                    }

                    if(( emsg = read_integer( p, 
                                buffer.length() - (size_t)(p-buffer.c_str()), &intval, &rbts )) != 0 ) {
                        sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                        throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                    }

                    p += rbts;
                    for( ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n'; p++ );

                    sssprob[t] = (float)intval / (float)scale;
                }
            }
        }//if( SSSset )
        //}}SS state


        //{{NOTE: SET profile model position-specific data and 
        // MOVE all pointers to the next position; 
        //target frequencies have been set already
        if( CtxVecset ) {
            //first, set the address
            pm2dvec.SetFieldNext( pmv2DAddrCV, 0, 
                (LNTYPE)(pm2dvec.GetFieldVector(pmv2DCVentrs)-pmorgaddr[pmv2DCVentrs]) / (LNTYPE)SZFPTYPE
            );
        }
        else
            pm2dvec.SetFieldNext( pmv2DAddrCV, 0, (LNTYPE)-1 );
        //set CV data
        for( r = 0; r < pmv2DNoCVEls; r++ )
            pm2dvec.SetFieldNext( pmv2DCVentrs + r, 0, (FPTYPE)ctxvec.GetValueAt(r));
        pm2dvec.SetFieldNext( pmv2DCVprior, 0, (FPTYPE)lpprb );
        pm2dvec.SetFieldNext( pmv2DCVnorm2, 0, (FPTYPE)norm2 );
        //
        if( SSSset ) {
            //first, set the address
            pm2dvec.SetFieldNext( pmv2DAddrSS, 0, 
                (LNTYPE)(pm2dvec.GetFieldVector(pmv2DSSsprbs)-pmorgaddr[pmv2DSSsprbs]) / (LNTYPE)SZFPTYPE
            );
        }
        else
            pm2dvec.SetFieldNext( pmv2DAddrSS, 0, (LNTYPE)-1 );
        //set SS data
        for( r = 0; r < pmv2DNoSSSps; r++ )
            pm2dvec.SetFieldNext( pmv2DSSsprbs + r, 0, (FPTYPE)sssprob[r]);
        pm2dvec.SetFieldNext( pmv2DSSstate, 0, (CHTYPE)DehashSSCodeLc(sss));
        //
        pm2dvec.SetFieldNext( pmv2DHDP1prb, 0, (FPTYPE)HDP1prb );
        pm2dvec.SetFieldNext( pmv2DHDP1ind, 0, (INTYPE)HDP1ind );
        pm2dvec.SetFieldNext( pmv2DAddrPro, 0, (INTYPE)profnr );
        //}}SET

    }//for(;m<prolen;)



    //NOTE: Initialize() for transitions



    //{{statistical parameters
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: No statistical parameters." );

    if(( p = strstr( buffer.c_str(), patstrEXPNN )) != NULL ) {
        p += strlen( patstrEXPNN );
        //do not do further analysis
    }
    else {
        //computed
        if(( emsg = skip_comments( fp, buffer )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || buffer.empty())
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameters." );

        //omit further analysis


        //reference
        if(( emsg = skip_comments( fp, buffer )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || buffer.empty())
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameters." );

        //omit further analysis


        //entropy, expected
        if(( emsg = skip_comments( fp, buffer )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || buffer.empty())
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No entropy." );

        //omit further analysis
    }//else (strstr( locbuffer, patstrEXPNN )) == NULL)
    //}}statistical parameters


    //footer
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( buffer.length() < lenstrEND )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No ending." );

    for( n = 0; n < (int)lenstrEND; n++ )
        if( buffer[n] != patstrEND[n] )
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Invalid ending." );
}

// -------------------------------------------------------------------------
// TextWriteProfileHD: write full profile data to file
//
void TextWriteProfileHD(FILE* fp, 
                    const mystring& desc, const mystring& file, 
                    int scale, 
                    char** pmorgaddr, 
                    char** pmdata )
{
    MYMSG("TextWriteProfileHD",3);
    const mystring preamb = "TextWriteProfileHD: ";

    if( !fp )
        throw MYRUNTIME_ERROR( preamb + "Null file descriptor.");
    if( !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null arguments.");

    if( scale < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale factor." );

    const int noress = NUMAA;
    const int vsize = pmv2DNoCVEls;//size of context vector
//     char msgbuf[BUF_MAX];
    char res;
    int m, z, r, t, tt, n = 0;

    PM2DVector pm2dvec( pmdata );


    int prolen = (int)*(INTYPE*)pm2dvec.GetField( pps2DLen, 0 );
    if( prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length." );


    fprintf( fp, "%s%s%s", patstrDATVER, pmodel::PMProfileModel::dataversion, NL );
    fprintf( fp, "%-9s%s%s", patstrDESC, desc.empty()? "": desc.c_str(), NL );
    fprintf( fp, "%-9s%s%s", patstrFILE, file.empty()? "": file.c_str(), NL );

    fprintf( fp, "%-9s", patstrCMD );   print_cmdline( &file_print, fp );
    fprintf( fp, "%-9s%d%s", patstrLEN, prolen, NL );
    fprintf( fp, "%-9s%d%s", patstrNOS, 0/*nodata*/, NL );
    fprintf( fp, "%-9s%.1f%s", patstrEFFNOS, (float)*(FPTYPE*)pm2dvec.GetField(pps2DENO,0), NL );
    fprintf( fp, "%-9s%d%s", patstrSCALE, scale, NL );

    fprintf( fp, "# Target / Observed / ");
    for( t = 0; t < P_NSTATES; t++ ) fprintf( fp, "%s ", gTPTRANS_NAMES[t]);
    fprintf( fp, "/ NexpM NexpI NexpD; Weight; Information /%s", NL );
    fprintf( fp, "# BPProb NoPosts / PProbs / PPNdxs / CT:NoPosts BPProb PProbs PPNdxs ");
    fprintf( fp, "/ CV: Prior Norm2 N Vector / SS:Pred Prob (C,E,H)%s", NL );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < noress; r++ )
        fprintf( fp, " %7c", DehashCode( r ) );


    fprintf( fp, "%s%7s   ", NL, patstrNULL );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", 
            (int)rintf( scale * (float)*(FPTYPE*)pm2dvec.GetField( pps2DBkgPrbs + r, 0 )));

    fprintf( fp, "%s%7s   ", NL, patstrPOST );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", 0/*nodata*/ );

    fprintf( fp, "%s%10c", NL, 32 );
    for( t = 0; t < P_NSTATES; t++ ) {
        tt = t;
        if( P_DI < t )
            tt -= 2;
        else if( P_ID < t && t < P_DI )
            tt -= 1;
        else if( P_ID <= t )
            continue;
        fprintf( fp, "%7d ", 
            (int)rintf( scale * (float)*(FPTYPE*)pm2dvec.GetField( ptr2DTrnPrbs + tt, 0 )));
    }

    fprintf( fp, "%s%10c", NL, 32 );
    for( t = 0; t < PS_NSTATES; t++ )
        fprintf( fp, "%7d ", 0/*nodata*/);


    for( m = 0; m < prolen; m++ ) {
        res = *pm2dvec.GetField( pmv2Daa, m );//get and...
        // omit unused positions and gaps in query
        if( res == '-' )
            continue;

        fprintf( fp, "%s%5d %c   ", NL, ++n, res );

//         if( pmvicN <= z || z < 0 ) {
//             MYMSGBEG
//                 sprintf( msgbuf, "%s: pos.: %d provec:", preamb.c_str(), m );
//                 for( int i = 0; i < 10; i++ )
//                     sprintf( msgbuf+strlen(msgbuf), " %x", provec[i]);
//                 sprintf( msgbuf+strlen(msgbuf), "..." );
//                 MYMSG( msgbuf,4 );
//             MYMSGEND
//             throw MYRUNTIME_ERROR( preamb + "Invalid the information indicator field.");
//         }

        for( r = 0; r < noress; r++ )
            fprintf( fp, "%7d ", 
                (int)rintf( scale * (float)*(FPTYPE*)pm2dvec.GetField( pmv2DTrgFrqs + r, m*SZFPTYPE )));

        fprintf( fp, "%s%10c", NL, 32 );

        for( r = 0; r < noress; r++ )//NOTE:!
            fprintf( fp, "%7d ", 0/*nodata*/);

        fprintf( fp, "%s%10c", NL, 32 );

        for( t = 0; t < P_NSTATES; t++ ) {
            tt = t;
            if( P_DI < t )
                tt -= 2;
            else if( P_ID < t && t < P_DI )
                tt -= 1;
            else if( P_ID <= t )
                continue;
            fprintf( fp, "%7d ", 
                (int)rintf( scale * (float)*(FPTYPE*)pm2dvec.GetField( ptr2DTrnPrbs + tt, (m+1)*SZFPTYPE )));
        }

        fprintf( fp, "%s%10c", NL, 32 );

        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, "%7d ", 0/*nodata*/);

        fprintf( fp, "%7d %7d ", 0, 0/*nodata*/);

        //{{HDP1
        fprintf( fp, "%s%10c", NL, 32 );
        fprintf( fp, "%7d %7d ", 0/*nodata*/, 1 );

        fprintf( fp, "%s%10c", NL, 32 );

        fprintf( fp, "%7d ", 
            (int)rintf( scale * (float)*(FPTYPE*)pm2dvec.GetField( pmv2DHDP1prb, m*SZFPTYPE )));

        fprintf( fp, "%s%10c", NL, 32 );

        fprintf( fp, "%7d ", (int)*(INTYPE*)pm2dvec.GetField( pmv2DHDP1ind, m*SZINTYPE ));
        //}}

        //{{HDP ctx
        //}}

        //{{context vector
        z = (int)*(LNTYPE*)pm2dvec.GetField( pmv2DAddrCV, m*SZLNTYPE );
        if( 0 <= z ) {
            fprintf( fp, "%s%14c%s %7d %7d %7d ", NL, 32, patstrCV,
                (int)rintf( scale * (float)*(FPTYPE*)(pmorgaddr[pmv2DCVprior] + z*SZFPTYPE)),
                (int)rintf( scale * (float)*(FPTYPE*)(pmorgaddr[pmv2DCVnorm2] + z*SZFPTYPE)),
                vsize );

            for( t = 0; t < vsize; t++ )
                fprintf( fp, "%7d ", 
                    (int)rintf( scale * (float)*(FPTYPE*)(pmorgaddr[pmv2DCVentrs+t] + z*SZFPTYPE)));
        }
        //}}

        //{{SS data
        z = (int)*(LNTYPE*)pm2dvec.GetField( pmv2DAddrSS, m*SZLNTYPE );
        if( 0 <= z ) {
            //use floor in rounding SS state probability
            fprintf( fp, "%s%13c%s%c", NL, 32, 
                patstrSS,
                /*DehashSSCode( */(char)*(CHTYPE*)(pmorgaddr[pmv2DSSstate] + z*SZCHTYPE ))/*)*/;
            for( t = 0; t < SS_NSTATES; t++ )
                fprintf( fp, " %7d", 
                    (int)( scale * (float)*(FPTYPE*)(pmorgaddr[pmv2DSSsprbs+t] + z*SZFPTYPE)));
        }
        //}}

    }//for(;m<prolen;)

    //{{NOTE:THIS BLOCK TO BE REMOVED
    fprintf( fp, "%s%10c%7d %7d", NL, 32,
                ( int )rintf( scale * (-7.f)),
                ( int )rintf( scale * (-1.f)));
    //}}
    fprintf( fp, "%s", NL );

    //{{statistical parameters
    fprintf( fp, "%-25s  %-6s   %-6s%s", " ", "K", "Lambda", NL );
    fprintf( fp, "%-25s  %6.4f   %6.4f%s", patstrSPCOMP, 0.0f, 0.0f/*nodata*/, NL );
    fprintf( fp, "%-25s  %6.4f   %6.4f%s", patstrSPREFR, 0.0f, 0.0f/*nodata*/, NL );
    fprintf( fp, "%s %6.4f; %s %6.4f%s", patstrENTROPY, 0.0f, 
             patstrEXPECTED, 0.0f/*nodata*/, NL );
    //}}

    fprintf( fp, "%s%s", patstrEND, NL );
}

// =========================================================================
// =========================================================================
// Instantiations

template
bool TextReadProfileHeader<FILE>(FILE* fp, 
                    mystring* desc, mystring* file,
                    int* prolen, int* scale, char** pmdata);
template
bool TextReadProfileHeader<TCharStream>(TCharStream* fp, 
                    mystring* desc, mystring* file,
                    int* prolen, int* scale, char** pmdata);

template
void TextReadProfileData<FILE>(FILE* fp,
                    unsigned int distn,
                    int profnr, int prolen, int scale,
                    char** pmorgaddr, char** pmdata );
template
void TextReadProfileData<TCharStream>(TCharStream* fp,
                    unsigned int distn,
                    int profnr, int prolen, int scale,
                    char** pmorgaddr, char** pmdata );

// =========================================================================
// =========================================================================



// -------------------------------------------------------------------------
// TextReadProfileHeaderBufferless: read profile header data from file 
// without using an additional buffer;
// desc, profile description; required for host;
// file, profile filename; required for host;
// prolen, profile length;
// scale, scale factor read from file;
// pmdata, profile-specific data required for both host and device;
// NOTE: on return, desc, file, and pmdata do not point to the end of data;
//
template<typename T>
bool TextReadProfileHeaderBufferless(T* fp, 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata)
{
    MYMSG( "TextReadProfileHeaderBufferless", 5 );
    const mystring preamb = "TextReadProfileHeaderBufferless: ";

    if( !fp )
        throw MYRUNTIME_ERROR( preamb + "Null file descriptor.");
    if( !prolen || !scale || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");

    const float accuracy = 1.0e-3f;
    const int noress = NUMAA;
    char* pdata = NULL;
    const char* p;
    size_t rbts, datlen = 0;
    int emsg;
    const char* descp, *filep;
    float effnos = 0.0f;
    int intval;
    float value, consv;
    int r;

    PM2DVector pm2dvec( pmdata );

    //read version number
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if(feof(fp) && datlen < 1)
        return false;

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( pdata, patstrDATVER )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrDATVER );

    if( datlen <= (size_t)(p-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( p, pmodel::PMProfileModel::dataversion )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Inappropriate profile version number." );


    //read description line
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( pdata, patstrDESC )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    descp = p + strlen( patstrDESC );
    for( ; *descp == ' ' || *descp == '\t'; descp++ );
    if( datlen <= (size_t)(descp-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );
    desc->assign(descp, datlen - (size_t)(descp-pdata));


    //read filename
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = ( char* )strstr( pdata, patstrFILE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    filep = p + strlen( patstrFILE );
    for( ; *filep == ' ' || *filep == '\t' ; filep++ );
    if( datlen <= (size_t)(filep-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );
    file->assign(filep, datlen - (size_t)(filep-pdata));


    //read command line without analysis
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    //read profile length
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( pdata, patstrLEN )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrLEN );

    if( datlen <= (size_t)(p-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), prolen, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( *prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length." );

    if( MAXCOLUMNS < *prolen )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Too large profile length." );


    //read number of sequences w/o analysis
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    //read effective number of sequences
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( pdata, patstrEFFNOS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrEFFNOS );

    if( datlen <= (size_t)(p-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_float( p, datlen - (size_t)(p-pdata), &effnos, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( effnos <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid effective number of observations." );


    //read scale factor
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( pdata, patstrSCALE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p += strlen( patstrSCALE );

    if( datlen <= (size_t)(p-pdata))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), scale, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( *scale < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid scale factor." );


    //read residue letters; omit checking the order
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));


    //{{read background probabilities
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    if(( p = strstr( pdata, patstrNULL )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No background probabilities." );

    p += strlen( patstrNULL );

    consv = 0.0f;
    for( r = 0; r < noress; r++ ) {
        if( datlen <= (size_t)(p-pdata))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete background probability data." );

        if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || *scale < intval )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Invalid background probability value." );

        value = (float)intval / (float)*scale;
        consv += value;

        //NOTE: SET profile-specific data
        pm2dvec.SetField( pps2DBkgPrbs + r, 0, (FPTYPE)value );
    }

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + 
        "Invalid profile data: Invalid background probabilities." );

    //NOTE: SET profile-specific data
    pm2dvec.SetField( pps2DENO, 0, (FPTYPE)effnos );
    pm2dvec.SetField( pps2DLen, 0, (INTYPE)*prolen );
    //}}


    //read posterior (generalized target) probabilities, but do not analyze them
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    return true;
}

// -------------------------------------------------------------------------
// TextReadProfileDataBufferless: read profile data from file or stream 
// without using an additional buffer;
// distn, distance in positions to profile to write to buffer;
// profnr, profile serial number;
// prolen, profile length;
// scale, scale factor for profile data;
// pmorgaddr, 2D vector of original addresses for a batch of profile models,
//          used to calculate a shift for CV and SS data from the beginning;
// pmdata is a 2D vector representing a batch of profile models, each of 
// which includes: 
// seqn, amino acid sequence the profile has been coonstructed for; 
//       required for host;
// prodata, position-specific profile model data required for device and 
//          host;
// ptrdata, position-specific profile model transition data required for 
//          device and host;
// NOTE: on return, the individual vectors of pmdata point to the end of 
// data;
//
template<typename T>
void TextReadProfileDataBufferless(T* fp, 
    unsigned int distn,
    int profnr,
    int prolen,
    int scale, 
    char** pmorgaddr, 
    char** pmdata )
{
    MYMSG( "TextReadProfileDataBufferless", 5 );
    const mystring preamb = "TextReadProfileDataBufferless: ";

    if( !fp )
        throw MYRUNTIME_ERROR( preamb + "Null file descriptor.");
    if( !pmorgaddr || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");

    //const float accuracy = 1.0e-3f;
    const int maxnocls = 1000;
    const int noress = NUMAA;
    char msgbuf[BUF_MAX];
    char* pdata = NULL;
    const char* p;
    size_t rbts, datlen = 0;
    bool lineread;
    int emsg;
    int intval;
    char res, sss;
    float value/*, consv*/;
    size_t noppps;
    float HDP1prb;
    int HDP1ind;
    bool HDPctset = false;//HDP ctx data's present
    bool CtxVecset = false;//CV's present
    extspsl::Pslvector ctxvec;//context vector
    float lpprb = 0.0f, norm2 = 0.0f;
    int vsize = 0;//cv size
    bool SSSset = false;//SS information's present
    bool SSSP3 = false;//SS information format
    float sssprob[SS_NSTATES];
    int m, r, tt, t, t0 = 0, tn = 0, n = 0;

    PM2DVector pm2dvec( pmdata );

    //{{beginning transition probabilities
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );

    p = pdata;

    for( t = tt = 0; t < P_NSTATES; t++ )
    {
        if( datlen <= (size_t)(p-pdata))
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format: Incomplete transition probability data." );

        if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 )
            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

        p += rbts;

        value = (float)intval / (float)scale;

        //NOTE: SET the beginning transitions of profile model and 
        // MOVE the pointer to the next position
//         tt = t;
//         if( P_DI < t )
//             tt -= 2;
//         else if( P_ID < t && t < P_DI )
//             tt -= 1;
//         else if( P_ID <= t )
//             continue;
        if( P_ID == t || t == P_DI )
            continue;
        value = value? (P_MD<t? logf(value)*MAIN_TRNPRB_LOGFCT: logf(value)): -32768.0f;
        pm2dvec.SetFieldNext( ptr2DTrnPrbs + tt, 0, (FPTYPE)value );
//         pm2dvec.SetFieldNext( ptr2DTrnPrbsExp + tt, 0, (FPTYPE)((value<SLC_LOG_SP_MIN)? 0.0f: expf(value)));
        tt++;
    }

    //NOTE: SET profile-specific data
    pm2dvec.SetFieldNext( pps2DDist, 0, (LNTYPE)distn );
    //}}


    //beginning MID state expected observations w/o analysis
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format." );


    lineread = false;

    for( m = 0; m < prolen; m++ )
    {
        //target probabilities
        if( !lineread ) {
            if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || datlen < 1 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }
        }

        lineread = false;

        if(( emsg = read_integer( p = pdata, datlen, &n, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( n != m + 1 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: Wrong numbering.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p += rbts;

        if( datlen <= (size_t)(p-pdata)) {
            sprintf( msgbuf, "Wrong profile format at pos %d: No residue.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_symbol( p, datlen - (size_t)(p-pdata), &res, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        p += rbts;

        //NOTE: SET amino acid and MOVE to the next position
        //NOTE:symbol is not hashed!
        //res = HashAlphSymbol( res );
        pm2dvec.SetFieldNext( pmv2Daa, 0, (CHTYPE)res );

        for( r = 0; r < noress; r++ )
        {
            if( datlen <= (size_t)(p-pdata)) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete target probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            value = (float)intval / (float)scale;

            //NOTE: SET profile model position-specific data and 
            // MOVE the pointer to the next position; 
            //target frequencies are always present and occupy the same address in the vector
            pm2dvec.SetFieldNext( pmv2DTrgFrqs + r, 0, (FPTYPE)value );
        }


        //observed frequencies w/o analysis
        //NOTE: condensed format does not include this section
//         if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
//             sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
//         }
// 
//         if( feof( fp ) || datlen < 1) {
//             sprintf( msgbuf, "Wrong profile format at pos %d.", m );
//             throw MYRUNTIME_ERROR( preamb + msgbuf );
//         }


        //transition probabilities
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || datlen < 1) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p = pdata;

        for( t = tt = 0; t < P_NSTATES; t++ )
        {
            if( datlen <= (size_t)(p-pdata)) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                 "Incomplete transition probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            p += rbts;

            value = (float)intval / (float)scale;

            //NOTE: SET the transitions of profile model at a given position and
            // MOVE to the next
//             tt = t;
//             if( P_DI < t )
//                 tt -= 2;
//             else if( P_ID < t && t < P_DI )
//                 tt -= 1;
//             else if( P_ID <= t )
//                 continue;
            if( P_ID == t || t == P_DI )
                continue;
            value = value? (P_MD<t? logf(value)*MAIN_TRNPRB_LOGFCT: logf(value)): -32768.0f;
            pm2dvec.SetFieldNext( ptr2DTrnPrbs + tt, 0, (FPTYPE)value );
//             pm2dvec.SetFieldNext( ptr2DTrnPrbsExp + tt, 0, (FPTYPE)((value<SLC_LOG_SP_MIN)? 0.0f: expf(value)));
            tt++;
        }


        //MID state observations w/o analysis
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || datlen < 1) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }


        //{{CLUSTERS: HDP1
        //bck posterior probability, no. posterior probability values
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( feof( fp ) || datlen < 1) {
            sprintf( msgbuf, "Wrong profile format at pos %d.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        p = pdata;

        if( datlen <= (size_t)(p-pdata)) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative HDP1 background posterior probability.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        //bppprob = (float)intval / (float)scale;

        p += rbts;

        if( datlen <= (size_t)(p-pdata)) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "No number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
        }

        if( intval < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Negative number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }
        if( maxnocls < intval ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Too large number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        noppps = (size_t)intval;

        if( noppps < 1 ) {
            HDP1prb = 0.0f;
            HDP1ind = -1;
        }
        else {
            //posterior predictive probabilities
            if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || datlen < 1) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p = pdata;

            //read only the first probability value
            if( datlen <= (size_t)(p-pdata)) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Incomplete HDP1 posterior probability data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( intval < 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Negative HDP1 posterior probability value.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p += rbts;

            HDP1prb = (float)intval / (float)scale;

            //indices of posteriors
            if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || datlen < 1) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                  "No HDP1 cluster indices.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p = pdata;

            //read only the index of the first cluster
            if( datlen <= (size_t)(p-pdata)) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Incomplete HDP1 cluster index data.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( intval != -1 && intval < 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Negative HDP1 cluster index.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            p += rbts;

            HDP1ind = intval;
        }//if( noppps )
        //}}HDP1


        //{{HDP ctx
        if( m < 1 )
            //NOTE: do not expect currently HDP ctx information!
            ;//HDPctset = true;//probe for HDP ctx information
        if( HDPctset ) {
            //bck posterior probability, no. posterior probability values
            if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
            }

            if( feof( fp ) || datlen < 1) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( pdata, patstrCT )) == NULL ) {
                if( m < 1 )
                    HDPctset = false;//no HDP ctx information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No HDP ctx probabilities.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCT );
                //do not read and check data futher in this line dedicated to HDP ctx information
            }//else
        }//if( HDPctset )
        //}}HDP ctx


        //{{Context vector
        if( m < 1 )
            CtxVecset = true;//probe for context vector data
        //Initialize even in the case of no use of CV scoring;
        //init each entry to ensure CVS2S score==0 when CVS information is not to be in use
        ctxvec.Reserve(pmv2DNoCVEls);
        ctxvec.AssignTo(0.17302947f);
        lpprb = 0.0f;
        norm2 = 0.0f;
        if( CtxVecset ) {
            if( !lineread )
                if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || datlen < 1) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( pdata, patstrCV )) == NULL ) {
                if( m < 1 )
                    CtxVecset = false;//no context vector data
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrCV );

                if( datlen <= (size_t)(p-pdata)) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No prior for context vector.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                lpprb = (float)intval / (float)scale;

                p += rbts;

                if( datlen <= (size_t)(p-pdata)) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                norm2 = (float)intval / (float)scale;

                p += rbts;

                if( datlen <= (size_t)(p-pdata)) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "No context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                if(( emsg = read_integer( p, datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( 1000 < intval ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Too large context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( !vsize && intval != pmv2DNoCVEls ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Only context vectors of size %d are supported.", m, pmv2DNoCVEls );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                if( vsize && vsize != intval && vsize != pmv2DNoCVEls ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Inconsistent context vector size.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }

                vsize = intval;

                p += rbts;

                if( vsize ) {
                    ctxvec.Allocate( vsize );
                    ctxvec.Clear();

                    //vector elements
                    for( t = 0; t < vsize; t++ )
                    {
                        if( datlen <= (size_t)(p-pdata)) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "Incomplete context vector data.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }

                        if(( emsg = read_integer( p, 
                                    datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                            sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                        }

                        p += rbts;

                        ctxvec.Push( (float)intval / (float)scale );
                    }
                }//if( vsize )
            }//else
        }//if( CtxVecset )
        //}}Context vector


        //{{SS state
        if( m < 1 ) {
            SSSset = true;//probe SS information existance
            SSSP3 = true;
            t0 = 0; tn = SS_NSTATES-1;
        }
        //Initialize even in the case of no use of SS information
        memset( sssprob, 0, SS_NSTATES*sizeof(float));
        sss = 0;
        if( SSSset ) {
            if( !lineread )
                if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

            if( feof( fp ) || datlen < 1) {
                sprintf( msgbuf, "Wrong profile format at pos %d.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            lineread = true;

            if(( p = strstr( pdata, patstrSS )) == NULL ) {
                if( m < 1 )
                    SSSset = false;//no SS information
                else {
                    sprintf( msgbuf, "Wrong profile format at pos %d: No SS state.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
            }
            else {
                lineread = false;
                p += strlen( patstrSS );

                if(( emsg = read_symbol( p, datlen - (size_t)(p-pdata), &sss, &rbts )) != 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                }

                p += rbts;

                sss = HashSSState( sss );

                if( !SSSP3 ) {
                    t0 = sss; tn = sss;
                }

                for( t = t0; t <= tn; t++ )
                {
                    if( datlen <= (size_t)(p-pdata)) {
                        if( m < 1 && t == t0+1 ) {
                            SSSP3 = false;
                            sssprob[(int)sss] = (float)intval / (float)scale;
                            break;
                        }
                        else {
                            sprintf( msgbuf, "Wrong profile format at pos %d: "
                                             "No SS state probability.", m );
                            throw MYRUNTIME_ERROR( preamb + msgbuf );
                        }
                    }

                    if(( emsg = read_integer( p, 
                                datlen - (size_t)(p-pdata), &intval, &rbts )) != 0 ) {
                        sprintf( msgbuf, "Wrong profile format at pos %d: ", m );
                        throw MYRUNTIME_ERROR( preamb + msgbuf + TranslateReadError( emsg ));
                    }

                    p += rbts;
                    for( ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n'; p++ );

                    sssprob[t] = (float)intval / (float)scale;
                }
            }
        }//if( SSSset )
        //}}SS state


        //{{NOTE: SET profile model position-specific data and 
        // MOVE all pointers to the next position; 
        //target frequencies have been set already
        if( CtxVecset ) {
            //first, set the address
            pm2dvec.SetFieldNext( pmv2DAddrCV, 0, 
                (LNTYPE)(pm2dvec.GetFieldVector(pmv2DCVentrs)-pmorgaddr[pmv2DCVentrs]) / (LNTYPE)SZFPTYPE
            );
        }
        else
            pm2dvec.SetFieldNext( pmv2DAddrCV, 0, (LNTYPE)-1 );
        //set CV data
        for( r = 0; r < pmv2DNoCVEls; r++ )
            pm2dvec.SetFieldNext( pmv2DCVentrs + r, 0, (FPTYPE)ctxvec.GetValueAt(r));
        pm2dvec.SetFieldNext( pmv2DCVprior, 0, (FPTYPE)lpprb );
        pm2dvec.SetFieldNext( pmv2DCVnorm2, 0, (FPTYPE)norm2 );
        //
        if( SSSset ) {
            //first, set the address
            pm2dvec.SetFieldNext( pmv2DAddrSS, 0, 
                (LNTYPE)(pm2dvec.GetFieldVector(pmv2DSSsprbs)-pmorgaddr[pmv2DSSsprbs]) / (LNTYPE)SZFPTYPE
            );
        }
        else
            pm2dvec.SetFieldNext( pmv2DAddrSS, 0, (LNTYPE)-1 );
        //set SS data
        for( r = 0; r < pmv2DNoSSSps; r++ )
            pm2dvec.SetFieldNext( pmv2DSSsprbs + r, 0, (FPTYPE)sssprob[r]);
        pm2dvec.SetFieldNext( pmv2DSSstate, 0, (CHTYPE)DehashSSCodeLc(sss));
        //
        pm2dvec.SetFieldNext( pmv2DHDP1prb, 0, (FPTYPE)HDP1prb );
        pm2dvec.SetFieldNext( pmv2DHDP1ind, 0, (INTYPE)HDP1ind );
        pm2dvec.SetFieldNext( pmv2DAddrPro, 0, (INTYPE)profnr );
        //}}SET

    }//for(;m<prolen;)



    //NOTE: Initialize() for transitions



    //{{statistical parameters
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( feof( fp ) || datlen < 1)
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: No statistical parameters." );

    if(( p = strstr( pdata, patstrEXPNN )) != NULL ) {
        p += strlen( patstrEXPNN );
        //do not do further analysis
    }
    else {
        //computed
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || datlen < 1)
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No statistical parameters." );

        //omit further analysis


        //reference
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || datlen < 1)
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No reference parameters." );

        //omit further analysis


        //entropy, expected
        if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: " + TranslateReadError( emsg ));

        if( feof( fp ) || datlen < 1)
            throw MYRUNTIME_ERROR( preamb + 
            "Wrong profile format at the end: No entropy." );

        //omit further analysis
    }//else (strstr( locbuffer, patstrEXPNN )) == NULL)
    //}}statistical parameters


    //footer
    if(( emsg = skip_comments( fp, pdata, 0, &datlen )) != 0 )
        throw MYRUNTIME_ERROR( preamb + 
        "Wrong profile format at the end: " + TranslateReadError( emsg ));

    if( datlen < lenstrEND )
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: No ending." );

    for( n = 0; n < (int)lenstrEND; n++ )
        if( pdata[n] != patstrEND[n] )
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Invalid ending." );
}

// =========================================================================
// =========================================================================
// Instantiations

//NOTE: version FILE needs using pdata with allocated space which requires revision!
// template
// void TextReadProfileHeaderBufferless<FILE>(FILE* fp, 
//                     mystring* desc, mystring* file,
//                     int* prolen, int* scale, char** pmdata);
template
bool TextReadProfileHeaderBufferless<TCharStream>(TCharStream* fp, 
                    mystring* desc, mystring* file,
                    int* prolen, int* scale, char** pmdata);

//NOTE: version FILE needs using pdata with allocated space which requires revision!
// template
// void TextReadProfileDataBufferless<FILE>(FILE* fp,
//                     unsigned int distn,
//                     int profnr, int prolen, int scale,
//                     char** pmorgaddr, char** pmdata );
template
void TextReadProfileDataBufferless<TCharStream>(TCharStream* fp,
                    unsigned int distn,
                    int profnr, int prolen, int scale,
                    char** pmorgaddr, char** pmdata );

// =========================================================================
// =========================================================================
