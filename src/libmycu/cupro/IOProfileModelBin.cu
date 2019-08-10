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

// ===== INLINES ===========================================================
// ReadSetFieldsUnrollSum: read and assign values to 2D vector of profile 
// data at the given field; calculate the sum;
// NOTE: SET profile-specific data and it is IMPORTANT to begin template 
// indexing with 0 as the pointer address is updated in pm2dvec.SetFieldNext
template<int N>
inline
float ReadSetFieldsUnrollSum(
    PM2DVector& pm2dvec, TPM2DVectorFields field, const char* ptr)
{
    float value = ((float*)(ptr))[N];
    pm2dvec.SetField( field + N, 0, (FPTYPE)value );
    return value + ReadSetFieldsUnrollSum<N+1>(pm2dvec, field, ptr);
}
template<>
inline float ReadSetFieldsUnrollSum<NUMAA>(PM2DVector&, TPM2DVectorFields, const char*) {return 0.0f;}

// ReadSetFieldsNextUnroll: read and assign values to 2D vector of profile;
// move the pointers appropriately
template<int N, int LIMIT>
struct ReadSetFieldsNextUnroll {
    static inline void run(
        PM2DVector& pm2dvec, TPM2DVectorFields field, const char* ptr)
    {
        pm2dvec.SetFieldNext( field + N, 0, (FPTYPE)( ((float*)(ptr))[N] ));
        ReadSetFieldsNextUnroll<N+1,LIMIT>::run(pm2dvec, field, ptr);
    }
};
template<int LIMIT>
struct ReadSetFieldsNextUnroll<LIMIT,LIMIT> {
    static inline void run(PM2DVector&, TPM2DVectorFields, const char*) {}
};

// ReadSetTrnFieldsNextUnroll: read and set transition data and move 
// appropriately the pointers
template<int N, int F>
struct ReadSetTrnFieldsNextUnroll{
    static inline void run(
        PM2DVector& pm2dvec, TPM2DVectorFields field, const char* ptr)
    {
        pm2dvec.SetFieldNext( field + F, 0, (FPTYPE)( ((float*)(ptr))[N] ));
        ReadSetTrnFieldsNextUnroll<N+1,F+1>::run(pm2dvec, field, ptr);
    }
};
template<int F>
struct ReadSetTrnFieldsNextUnroll<P_ID,F>{
    static inline void run(PM2DVector& pm2dvec, TPM2DVectorFields field, const char* ptr) {
        ReadSetTrnFieldsNextUnroll<P_ID+1,F>::run(pm2dvec, field, ptr);
    }
};
template<int F> 
struct ReadSetTrnFieldsNextUnroll<P_DI,F>{
    static inline void run(PM2DVector& pm2dvec, TPM2DVectorFields field, const char* ptr) {
        ReadSetTrnFieldsNextUnroll<P_DI+1,F>::run(pm2dvec, field, ptr);
    }
};
template<int F> 
struct ReadSetTrnFieldsNextUnroll<P_NSTATES,F>{
    static inline void run(PM2DVector&, TPM2DVectorFields, const char*) {}
};
// =========================================================================

// -------------------------------------------------------------------------
// BinaryReadProfileHeader: read profile header data from binary file;
// desc, profile description; required for host;
// file, profile filename; required for host;
// prolen, profile length;
// scale, scale factor set to 1;
// pmdata, profile-specific data required for both host and device;
// NOTE: on return, desc, file, and pmdata do not point to the end of data;
// 
void BinaryReadProfileHeader(TCharStream* fp, 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata)
{
    MYMSG( "BinaryReadProfileHeader", 5 );
    const mystring preamb = "BinaryReadProfileHeader: ";

    if( !fp || !fp->data_ || fp->datlen_ < 1 )
        throw MYRUNTIME_ERROR( preamb + "Null data stream.");
    if( !prolen || !scale || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");

    const float accuracy = 1.0e-3f;
    const int noress = NUMAA;
    //mystring buffer;
    const char* p;
    size_t len;
    float effnos = 0.0f;
    float /*value, */consv;
    //int r;

    PM2DVector pm2dvec( pmdata );

    *scale = 1;

    p = fp->data_+fp->curpos_;

    //the first field is profile size; pass over it
    len = sizeof(int);
    p += len;
    fp->incpos(len);

    //read version number
    len = strlen(patstrDATBINVER);
    if( strncmp(p,patstrDATBINVER,len))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");

    len = strlen(pmodel::PMProfileModel::dataversion);
    if( strncmp(p,pmodel::PMProfileModel::dataversion,len))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");

    //read description line
    len = strlen(p);
    fp->incpos(len+1);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    *desc++ = p;
    p += len+1;//add 1 to pass over 0

    //read filename
    len = strlen(p);
    fp->incpos(len+1);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    *file++ = p;
    p += len+1;//add 1 to pass over 0

    //read profile length
    len = sizeof(int);
    *prolen = *(int*)p;
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    if( MAXCOLUMNS < *prolen || *prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length.");

    //read effective number of observations
    len = sizeof(float);
    effnos = *(float*)p;
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    if( effnos <= 0.0f )
        throw MYRUNTIME_ERROR( preamb + "Invalid effective number of observations.");

    //{{read background probabilities
    len = sizeof(float) * noress;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
//     consv = 0.0f;
//     for( r = 0; r < noress; r++ ) {
//         consv += ( value = *(float*)(p+r) );
//         if( value < 0.0f )
//             throw MYRUNTIME_ERROR( preamb + "Invalid probability value.");
//         pm2dvec.SetField( pps2DBkgPrbs + r, 0, (FPTYPE)value );
//     }
    //NOTE: SET profile-specific data
    consv = ReadSetFieldsUnrollSum<0>(pm2dvec, pps2DBkgPrbs, p);
    p += len;

    if( consv < 1.0f - accuracy || consv > 1.0f + accuracy )
        throw MYRUNTIME_ERROR( preamb + "Invalid background probabilities.");

    //NOTE: SET profile-specific data
    pm2dvec.SetField( pps2DENO, 0, (FPTYPE)effnos );
    pm2dvec.SetField( pps2DLen, 0, (INTYPE)*prolen );
    //}}

    //{{omit posterior (generalized target) probabilities that are 
    // present in the profile for the future
    len = sizeof(float) * noress;
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    //}}
}

// -------------------------------------------------------------------------
// TextReadProfileData: read profile data from file;
// distn, distance in positions to profile to write to buffer;
// profnr, profile serial number;
// prolen, profile length;
// scale, scale factor unused;
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
void BinaryReadProfileData(TCharStream* fp, 
    unsigned int distn,
    int profnr,
    int prolen,
    int /*scale*/, 
    char** pmorgaddr, 
    char** pmdata )
{
    MYMSG( "BinaryReadProfileData", 5 );
    const mystring preamb = "BinaryReadProfileData: ";

    if( !fp || !fp->data_ || fp->datlen_ < 1 )
        throw MYRUNTIME_ERROR( preamb + "Null data stream.");
    if( !pmorgaddr || !pmdata )
        throw MYRUNTIME_ERROR( preamb + "Null addresses for output data.");
    if( prolen < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid profile length.");

    //const float accuracy = 1.0e-3f;
    const int maxnocls = 1000;
    const int noress = NUMAA;
    char msgbuf[BUF_MAX];
    mystring sequence;
    const char* p;
    size_t len;
    char sss;
    int noppps;
    float HDP1prb;
    int HDP1ind;
    bool HDPctset = false;//HDP ctx data's present
    bool CtxVecset = false;//CV's present
    float lpprb = 0.0f, norm2 = 0.0f;
    int vsize = pmv2DNoCVEls;//cv size
    extspsl::Pslvector ctxvec;//context vector
    extspsl::Pslvector uninfctxvec(pmv2DNoCVEls);//uninformative context vector
    bool SSSset = false;//SS information's present
    extspsl::Pslvector sssvec;//SSS vector
    extspsl::Pslvector uninfsssvec(SS_NSTATES);//uninformative SSS vector
    int m, intval;

    //init each entry to ensure CVS2S score==0 when CVS information is not to be in use
    uninfctxvec.AssignTo(0.17302947f);

    //Initialize uninformative SSS information
    memset( uninfsssvec.GetVector(), 0, uninfsssvec.GetSize()*sizeof(float));

    PM2DVector pm2dvec( pmdata );

    p = fp->data_+fp->curpos_;

    //{{read the whole sequence at once
    sequence.assign(p, prolen);
    p += prolen;
    fp->incpos(prolen);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    //}}

    //{{beginning transition probabilities
    len = sizeof(float) * P_NSTATES;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
    //NOTE: SET the beginning transitions of profile model and 
    // MOVE the pointer to the next position
    ReadSetTrnFieldsNextUnroll<0,0>::run(pm2dvec, ptr2DTrnPrbs, p);
    p += len;
    //NOTE: SET profile-specific data
    pm2dvec.SetFieldNext( pps2DDist, 0, (LNTYPE)distn );
    //}}


    for( m = 0; m < prolen; m++ )
    {
        //NOTE: SET amino acid and MOVE to the next position
        //NOTE:symbol is not hashed!
        pm2dvec.SetFieldNext( pmv2Daa, 0, (CHTYPE)sequence[m]);

        //{{target probabilities
        len = sizeof(float) * noress;
        fp->incpos(len);
        if( fp->eof())
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
        //NOTE: SET profile-specific data
        ReadSetFieldsNextUnroll<0,NUMAA>::run(pm2dvec, pmv2DTrgFrqs, p);
        p += len;
        //}}

        //{{transition probabilities
        len = sizeof(float) * P_NSTATES;
        fp->incpos(len);
        if( fp->eof())
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
        //NOTE: SET the transitions of profile model at a position and 
        // MOVE the pointer to the next position
        ReadSetTrnFieldsNextUnroll<0,0>::run(pm2dvec, ptr2DTrnPrbs, p);
        p += len;
        //}}

        //{{CLUSTERS: HDP1
        //omit background probability:
        len = sizeof(float);
        p += len;
        fp->incpos(len);
        if( fp->eof())
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");

        //number of components:
        len = sizeof(int);
        noppps = *(int*)p;
        p += len;
        fp->incpos(len);
        if( fp->eof())
            throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
        if( maxnocls < noppps || noppps < 0 ) {
            sprintf( msgbuf, "Wrong profile format at pos %d: "
                             "Invalid number of HDP1 posterior probabilities.", m );
            throw MYRUNTIME_ERROR( preamb + msgbuf );
        }

        if( noppps < 1 ) {
            HDP1prb = 0.0f;
            HDP1ind = -1;
        }
        else {
            //posterior predictive probabilities
            //read the first probability value
            len = sizeof(float) * noppps;
            HDP1prb = *(float*)p;
            p += len;
            fp->incpos(len);
            if( fp->eof())
                throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
            if( HDP1prb < 0.0f ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                "Invalid HDP1 posterior probability value.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }

            //indices of posteriors
            //read the index of the first cluster
            len = sizeof(int) * noppps;
            HDP1ind = *(int*)p;
            p += len;
            fp->incpos(len);
            if( fp->eof())
                throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
            if( HDP1ind < -1 ) {
                sprintf( msgbuf, "Wrong profile format at pos %d: "
                                "Invalid HDP1 cluster index.", m );
                throw MYRUNTIME_ERROR( preamb + msgbuf );
            }
        }//if( noppps )
        //}}HDP1

        //{{HDP ctx
        if( m < 1 )
            HDPctset = true;//probe for HDP ctx information
        if( HDPctset ) {
            len = strlen(patstrCT);
            if( strncmp(p,patstrCT,len)) {
                if(m) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                        "No HDP ctx probabilities.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                else
                    HDPctset = false;//no HDP ctx information
            }
            else {
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
                intval = *(int*)p;
                if( intval < 0 ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                    "Invalid number of HDP ctx components.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                //omit information currently
                len = sizeof(int) + sizeof(float) + (sizeof(float) + sizeof(int)) * intval;
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
           }
        }//if( HDPctset )
        //}}HDP ctx

        //{{Context vector
        if( m < 1 )
            CtxVecset = true;//probe for context vector data
        //Initialize even in the case of no use of CV scoring;
        //init each entry to ensure CVS2S score==0 when CVS information is not to be in use
        ctxvec.SetVector(uninfctxvec.GetVector(), /*uninfctxvec.GetSize()*/pmv2DNoCVEls);
        lpprb = 0.0f;
        norm2 = 0.0f;
        if( CtxVecset ) {
            len = strlen(patstrCV);
            if( strncmp(p,patstrCV,len)) {
                if(m) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                        "No context vector data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                else
                    CtxVecset = false;//no context vector data
            }
            else {
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
                //CV: probability:
                len = sizeof(float);
                lpprb = *(float*)p;
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR(preamb+"Wrong profile format: No prior for context vector.");
                //CV: squared norm:
                len = sizeof(float);
                norm2 = *(float*)p;
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR(preamb+"Wrong profile format: No context vector norm.");
                if( norm2 < 0.0f ) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: "
                                     "Negative context vector norm.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                //CV: size:
                vsize = pmv2DNoCVEls;
                //NOTE: (not in profile binary format)
                //len = sizeof(int);
                //vsize = *(int*)p;
                //p += len;
                //fp->incpos(len);
                //if( fp->eof())
                //    throw MYRUNTIME_ERROR(preamb+"Wrong profile format: No context vector size.");
                //if( 1000 < vsize || vsize < 0 ) {
                //    sprintf( msgbuf, "Wrong profile format at pos %d: "
                //                     "Invalid context vector size.", m );
                //    throw MYRUNTIME_ERROR( preamb + msgbuf );
                //}

                if( vsize ) {
                    //CV: vector elements:
                    len = sizeof(float) * vsize;
                    ctxvec.SetVector((float*)p, vsize);
                    p += len;
                    fp->incpos(len);
                    if( fp->eof())
                        throw MYRUNTIME_ERROR(preamb+"Wrong profile format: No context vector data.");
                }//if( vsize )
            }//else
        }//if( CtxVecset )
        //}}Context vector

        //{{SS state
        if( m < 1 )
            SSSset = true;//probe SS information existance
        //Initialize even in the case of no use of SS information
        sssvec.SetVector(uninfsssvec.GetVector(), /*uninfsssvec.GetSize()*/SS_NSTATES);
        sss = ' ';
        if( SSSset ) {
            len = strlen(patstrSS);
            if( strncmp(p,patstrSS,len)) {
                if(m) {
                    sprintf( msgbuf, "Wrong profile format at pos %d: No SS data.", m );
                    throw MYRUNTIME_ERROR( preamb + msgbuf );
                }
                else
                    SSSset = false;//no SS information
            }
            else {
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
                //SS: state:
                len = sizeof(char);
                sss = *p;
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR(preamb+"Wrong profile format.");
                //SS: vector elements:
                len = sizeof(float) * SS_NSTATES;
                sssvec.SetVector((float*)p, SS_NSTATES);
                p += len;
                fp->incpos(len);
                if( fp->eof())
                    throw MYRUNTIME_ERROR(preamb+"Wrong profile format: No SS state probability.");
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
        ReadSetFieldsNextUnroll<0,pmv2DNoCVEls>::run(pm2dvec, pmv2DCVentrs, (char*)ctxvec.GetVector());
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
        ReadSetFieldsNextUnroll<0,pmv2DNoSSSps>::run(pm2dvec, pmv2DSSsprbs, (char*)sssvec.GetVector());
        pm2dvec.SetFieldNext( pmv2DSSstate, 0, (CHTYPE)sss);
        //
        pm2dvec.SetFieldNext( pmv2DHDP1prb, 0, (FPTYPE)HDP1prb );
        pm2dvec.SetFieldNext( pmv2DHDP1ind, 0, (INTYPE)HDP1ind );
        pm2dvec.SetFieldNext( pmv2DAddrPro, 0, (INTYPE)profnr );
        //}}SET

    }//for(;m<prolen;)



    //NOTE: Initialize() for transitions



    //footer
    len = strlen(patstrEND);
    if( strncmp(p,patstrEND,len))
        throw MYRUNTIME_ERROR( preamb + "Wrong profile format: Invalid ending.");
    p += len;
    fp->incpos(len);
    //if( fp->eof())
    //    throw MYRUNTIME_ERROR( preamb + "Wrong profile format.");
}
