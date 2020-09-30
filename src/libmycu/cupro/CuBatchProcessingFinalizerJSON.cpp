/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <math.h>
#include <cmath>

#include <utility>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/fmtdescription.h"
#include "libmycu/cucom/btckcoords.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libmycu/cuss/CuBatchSS_com.h"

#include "CuBatchProcessingFinalizer.h"


// -------------------------------------------------------------------------
// CompressResultsJSON: process results obtained; compress them for passing 
// them to the writing thread; use JSON format
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::CompressResultsJSON()
{
    MYMSG( "CuBatchProcessingFinalizer::CompressResultsJSON", 4 );
    static const unsigned int annotlen = ANNOTATIONLEN;
    static const bool printsss = MOptions::GetSSSWGT() > 0.0f;
    const unsigned int dsclen = MOptions::GetDSCLEN();
    //const int show = MOptions::GetSHOW();
    size_t szannot = 0UL;
    size_t szalns = 0UL;
    size_t szalnswodesc = 0UL;
    const char* desc;//profile description
    int written, sernr = 0;//initialize serial number here
    //
    bool qrysssinuse = 
        GetQueryFieldPos<LNTYPE>(pmv2DAddrSS, 0/*cubp_set_querposoffset_*/) != (LNTYPE)-1;

    GetSizeOfCompressedResultsJSON( &szannot, &szalns, &szalnswodesc);

    annotations_.reset();
    alignments_.reset();

    ReserveVectors( cubp_set_npros_ );

    if( szalns < szannot || 
        szalnswodesc > 2 * (cubp_set_sz_alndata_ + cubp_set_sz_alns_))
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResultsJSON: "
        "Size of compressed results is unusually large.");

    if( szannot < 1 || szalns < 1 )
        return;

    annotations_.reset((char*)std::malloc(szannot));
    alignments_.reset((char*)std::malloc(szalns));

    if( !annotations_ || !alignments_)
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResultsJSON: "
        "Not enough memory.");

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResultsJSON: "
        "Not enough memory.");

    char* annptr = annotations_.get();
    char* outptr = alignments_.get();

    for(unsigned int prondx = 0; prondx < cubp_set_npros_; prondx++)
    {
        float logeval = GetOutputAlnDataField<float>(prondx, dp2oadEvalue);
        if( cubp_set_logevthld_ < logeval )
            continue;
        unsigned int alnlen = GetOutputAlnDataField<int>(prondx, dp2oadEpA);
        if( alnlen < 1 )
            continue;
        //get the index over all (two) db pm structures
        unsigned int orgprondx = GetOutputAlnDataField<unsigned int>(prondx, dp2oadOrgProNo);
        //distance form the beginning in phase 2 (profiles passed through phase 1)
        unsigned int dbpro2dst = GetOutputAlnDataField<unsigned int>(prondx, dp2oadProNewDst);
        int dbprolen = GetDbProfileField<INTYPE>(orgprondx, pps2DLen);
        if( cubp_set_nposits_ < dbpro2dst + dbprolen)
            //the profile has not been processed due to GPU memory restrictions
            continue;
        //
        float dbproENO = GetDbProfileField<FPTYPE>(orgprondx, pps2DENO);
        //save the addresses of the annotations and the alignment records
        srtindxs_->push_back(sernr++);
        logevalues_->push_back(logeval);
        annotptrs_->push_back(annptr);
        alnptrs_->push_back(outptr);
        //
        double evalue = GetEvalue(logeval);
        float score = GetOutputAlnDataField<float>(prondx, dp2oadScore);
        //get the name and description
        GetDbProfileDesc(desc, orgprondx);
        //make an annotation
        MakeAnnotationJSON( annptr,
            desc, annotlen, 
            score, evalue);
        *annptr++ = 0;//end of record
        //
        //compress the alignment and relative information...
        written = sprintf( outptr,
                "    {\"hit_record\": {%s"
                "      \"target_description\": \"",NL);
        outptr += written;
        //put the description...
        int outpos = 0;
        if( desc )
            FormatDescriptionJSON(outptr, desc, dsclen, outpos);
        written = sprintf( outptr,"\",%s",NL);
        outptr += written;
        written = sprintf( outptr,
                "      \"query_length\": %d,%s"
                "      \"query_eno\": %.1f,%s"
                "      \"target_length\": %d,%s"
                "      \"target_eno\": %.1f,%s"
                "      \"alignment\": {%s",
                cubp_set_nqyposs_,NL,cubp_set_qyeno_,NL,
                dbprolen,NL,dbproENO,NL,NL);
        outptr += written;
        //
        FormatScoresJSON(outptr, prondx, alnlen, score, logeval, evalue);
        FormatAlignmentJSON( outptr,
            prondx, orgprondx, dbpro2dst, alnlen, dbprolen, 
            printsss, qrysssinuse );
        //if( show )
        FormatFooterJSON( outptr, prondx );
        //
        written = sprintf( outptr,"      }%s    }}",NL);//,%s",NL,NL);
        outptr += written;
        *outptr++ = 0;//end of record
    }
}

// -------------------------------------------------------------------------
// MakeAnnotationJSON: format profile description;
// NOTE: space is assumed to be pre-allocated;
// outptr, pointer to the output buffer;
// name, profile name;
// desc, profile description;
// printname, whether to include the name in the annotation;
// maxoutlen, maximum length of output description;
// width, width to wrap the profile description;
// score, alignment score;
// evalue, e-value of the alignment;
inline
void CuBatchProcessingFinalizer::MakeAnnotationJSON( 
    char*& outptr,
    const char* desc,
    const int maxoutlen,
    const float score,
    const double evalue) const
{
    int outpos = 0;
    int written;

    written = sprintf( outptr,
                "    {\"summary_entry\": {%s"
                "      \"description\": \"",NL);
    outptr += written;

    if( desc )
        FormatDescriptionJSON(outptr, desc, maxoutlen, outpos);

    written = sprintf( outptr,"\",%s",NL);
    outptr += written;

    written = sprintf( outptr,
                "      \"score\": %.0f,%s"
                "      \"evalue\": %.2g%s"
                "    }}",//,%s",
            score,NL,evalue,NL/*,NL*/);
    outptr += written;
}

// -------------------------------------------------------------------------
// outptr, pointer to the output buffer;
// prondx, profile index in the results list;
// alnlen, alignment length;
// score, alignment score;
// logeval, log e-value;
// evalue, e-value;
inline
void CuBatchProcessingFinalizer::FormatScoresJSON(
    char*& outptr,
    unsigned int prondx,
    unsigned int alnlen,
    float score,
    float logeval,
    double evalue )
{
    int written;
    int psts = (int)GetOutputAlnDataField<float>(prondx, dp2oadPstvs);
    int idts = (int)GetOutputAlnDataField<float>(prondx, dp2oadIdnts);
    int gaps = (int)GetOutputAlnDataField<float>(prondx, dp2oadNGaps);
    written = sprintf( outptr,
            "        \"score\": %.2f,%s"
            "        \"bit_score\": %.1f,%s"
            "        \"evalue\": %.2g,%s"
            "        \"pvalue\": %.2g,%s",
            score,NL,GetBitScore(logeval),NL,evalue,NL,GetPvalue(evalue),NL);
    outptr += written;
    written = sprintf( outptr,
            "        \"n_identities\": %d,%s"
            "        \"n_positives\": %d,%s"
            "        \"n_gaps\": %d,%s"
            "        \"aln_length\": %d,%s",
            idts,NL,psts,NL,gaps,NL,alnlen,NL);
    outptr += written;
}

// -------------------------------------------------------------------------
// outptr, pointer to the output buffer;
// prondx, profile index in the results list;
// orgprondx, profile index over all pm data structures;
// dbpro2dst, distance form the beginning in phase 2;
// alnlen, alignment length;
// dbprolen, length of a db profile;
// width, alignment output width;
// printsss, print SSS information;
// qrysssinuse, query contains SSS information;
inline
void CuBatchProcessingFinalizer::FormatAlignmentJSON(
    char*& outptr,
    unsigned int prondx,
    unsigned int orgprondx,
    unsigned int dbpro2dst,
    int alnlen,
    int dbprolen,
    bool printsss,
    const bool qrysssinuse )
{
    //alignment beginning coordinates:
    unsigned int alnbegcoords = GetOutputAlnDataField<unsigned int>(prondx, dp2oadBegCoords);
    unsigned int qrybeg = GetCoordY(alnbegcoords)+1;
    unsigned int trgbeg = GetCoordX(alnbegcoords)+1;
    unsigned int dbprodst = GetDbProfileField<LNTYPE>(orgprondx, pps2DDist);
    bool trgsssinuse = 
        GetDbProfileFieldPos<LNTYPE>(orgprondx, pmv2DAddrSS, dbprodst) != (LNTYPE)-1;
    //offset to the end of the db profile in phase 2:
    int dbpos2off = dbpro2dst + dbprolen-1;
    //alignment beginning position:
    int alnbeg = dbpos2off + (prondx+1) * cubp_set_nqyposs_ - alnlen;
    //
    //beginning of the alignment:
    const char* palnbeg = GetBegOfAlns() + alnbeg;
    const char* pqueryaln = GetAlnSectionAt(palnbeg, dp2oaQuery);
    const char* ptargetaln = GetAlnSectionAt(palnbeg, dp2oaTarget);
    const char* p;

    int written;
    int nbytes = alnlen;
    int ngapsqry = (int)(std::count(pqueryaln, pqueryaln+nbytes, '-'));//#gaps in query
    int ngapstrg = (int)(std::count(ptargetaln, ptargetaln+nbytes, '-'));//#gaps in target

    unsigned int qryend = qrybeg + nbytes - ngapsqry;
    unsigned int trgend = trgbeg + nbytes - ngapstrg;

    printsss = /*printsss && */qrysssinuse && trgsssinuse;

    written = sprintf( outptr,
            "        \"query_from\": %u,%s"
            "        \"query_to\": %u,%s"
            "        \"target_from\": %u,%s"
            "        \"target_to\": %u,%s",
            qrybeg,NL,nbytes<=ngapsqry? qryend: qryend-1,NL,
            trgbeg,NL,nbytes<=ngapstrg? trgend: trgend-1,NL);
    outptr += written;

    written = sprintf( outptr,
            "        \"query_secstr_available\": %d,%s"
            "        \"target_secstr_available\": %d,%s"
            "        \"query_secstr\": \"",
            printsss,NL,printsss,NL);
    outptr += written;

    if(printsss) {
        p = GetAlnSectionAt(palnbeg, dp2oaQuerySSS);
        strncpy(outptr, p, nbytes);
        outptr += nbytes;
    }
    written = sprintf( outptr,"\",%s        \"target_secstr\": \"",NL);
    outptr += written;
    if(printsss) {
        p = GetAlnSectionAt(palnbeg, dp2oaTargetSSS);
        strncpy(outptr, p, nbytes);
        outptr += nbytes;
    }
    written = sprintf( outptr,"\",%s        \"query_aln\": \"",NL);
    outptr += written;
    strncpy(outptr, pqueryaln, nbytes);
    outptr += nbytes;
    written = sprintf( outptr,"\",%s        \"target_aln\": \"",NL);
    outptr += written;
    strncpy(outptr, ptargetaln, nbytes);
    outptr += nbytes;
    written = sprintf( outptr,"\",%s        \"middle\": \"",NL);
    outptr += written;
    p = GetAlnSectionAt(palnbeg, dp2oaMiddle);
    strncpy(outptr, p, nbytes);
    outptr += nbytes;
    written = sprintf( outptr,"\",%s",NL);
    outptr += written;
}

// -------------------------------------------------------------------------
// outptr, pointer to the output buffer;
// prondx, profile index in the results list;
inline
void CuBatchProcessingFinalizer::FormatFooterJSON(
    char*& outptr,
    unsigned int prondx )
{
    //statistics:
    float lmbd_est = GetOutputAlnDataField<float>(prondx, dp2oadLmbdEst);
    float K_est = GetOutputAlnDataField<float>(prondx, dp2oadKEst);
    float minsc = GetOutputAlnDataField<float>(prondx, dp2oadMin);
    float maxsc = GetOutputAlnDataField<float>(prondx, dp2oadMax);
    float expected = GetOutputAlnDataField<float>(prondx, dp2oadE);
    float entropy = GetOutputAlnDataField<float>(prondx, dp2oadH);
    float lmbd = GetOutputAlnDataField<float>(prondx, dp2oadLmbd);
    float K = GetOutputAlnDataField<float>(prondx, dp2oadK);
    int written;
    static const char* na = "NA";
    char kbuf[BUF_MAX], lbuf[BUF_MAX];
    const char* kp, *lp;

    kp = lp = na;
    if( expected < 0.0f ) {
        if(K > 0.0f) { 
            sprintf( kbuf,"%6.4f",K);
            kp = kbuf;
        }
        if(lmbd > 0.0f) {
            sprintf( lbuf,"%6.4f",lmbd);
            lp = lbuf;
        }
    }
    written = sprintf( outptr,
            "        \"computed_K_ungapped\": \"%s\",%s"
            "        \"computed_lambda_ungapped\": \"%s\",%s",
            kp,NL,lp,NL);
    outptr += written;
    kp = lp = na;
    if( expected < 0.0f ) {
        if(K_est > 0.0f) {
            sprintf( kbuf,"%6.4f",K_est);
            kp = kbuf;
        }
        if(lmbd_est > 0.0f) {
            sprintf( lbuf,"%6.4f",lmbd_est);
            lp = lbuf;
        }
    }
    written = sprintf( outptr,
            "        \"estimated_K_gapped\": \"%s\",%s"
            "        \"estimated_lambda_gapped\": \"%s\",%s"
            "        \"entropy\": %.4f,%s"
            "        \"expected_score\": %.4f,%s"
            "        \"min_score\": %.0f,%s"
            "        \"max_score\": %.0f%s",
            kp,NL,lp,NL,
            entropy,NL,expected,NL,minsc,NL,maxsc,NL);
    outptr += written;
}



// -------------------------------------------------------------------------
// GetSizeOfCompressedResultsJSON: get total size required for annotations 
// and complete alignments; using JSON format;
// szannot, size of annotations;
// szalns, size of complete alignments (with descriptions);
// szalnswodesc, size of alignments without descriptions;
//
inline
void CuBatchProcessingFinalizer::GetSizeOfCompressedResultsJSON(
    size_t* szannot, size_t* szalns, size_t* szalnswodesc) const
{
    MYMSG( "CuBatchProcessingFinalizer::GetSizeOfCompressedResultsJSON", 5 );
    static const unsigned int sznl = (int)strlen(NL);
    static const unsigned int maxfieldlen = 40;//JSON max field length
    static const unsigned int annotlen = ANNOTATIONLEN;
    static const unsigned int annotheadlines = 2;//number of lines for score and evalue
    //static const bool printsss = MOptions::GetSSSWGT() > 0.0f;
    const unsigned int dsclen = MOptions::GetDSCLEN();
    //static const int show = MOptions::GetSHOW();
    static const unsigned int openlines = 6;//number of lines for opening a hit record
    static const unsigned int headlines = 14;//number of lines for scores, etc.
    static const unsigned int footlines = 8;//number of lines for statistical parameters
    static const unsigned int closlines = 2;//number of lines for closing a hit record
    static const unsigned int maxlinelen = 90;//maximum length of lines other than alignment lines
    static const unsigned int maxopenlinelen = 50;//maximum length of `openlines'
    static const unsigned int maxcloslinelen = 20;//maximum length of `closlines'
    int alnsize;//alignment size
    const char* desc;//profile description
    //
    bool qrysssinuse = 
        GetQueryFieldPos<LNTYPE>(pmv2DAddrSS, 0/*cubp_set_querposoffset_*/) != (LNTYPE)-1;

    if( cubp_set_bdbCpmbeg_[0] && cubp_set_bdbCpmend_[0] && !cubp_set_bdbCdesc_)
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::GetSizeOfCompressedResultsJSON: Null profile descriptions.");

#ifdef __DEBUG__
    if( !cached_queryfnames_ || !cached_querydescs_ /*||!cached_bdb1descs_*/)
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::GetSizeOfCompressedResultsJSON: Null query descriptions.");
#endif

    *szannot = 0UL;
    *szalns = 0UL;
    *szalnswodesc = 0UL;

    //on the first receive of results, they can be empty if 
    // there are no hits found
    if( !cubp_set_h_results_)
        return;

    for(unsigned int prondx = 0; prondx < cubp_set_npros_; prondx++)
    {
        float logeval = GetOutputAlnDataField<float>(prondx, dp2oadEvalue);
        if( cubp_set_logevthld_ < logeval )
            continue;
        unsigned int alnlen = GetOutputAlnDataField<unsigned int>(prondx, dp2oadEpA);
        if( alnlen < 1 )
            continue;
        //get the index over all (two) db pm structures:
        unsigned int orgprondx = GetOutputAlnDataField<unsigned int>(prondx, dp2oadOrgProNo);
        //distance form the beginning in phase 2 (profiles passed through phase 1)
        unsigned int dbpro2dst = GetOutputAlnDataField<unsigned int>(prondx, dp2oadProNewDst);
        int dbprolen = GetDbProfileField<INTYPE>(orgprondx, pps2DLen);
        if( cubp_set_nposits_ < dbpro2dst + dbprolen)
            //the profile has not been processed due to GPU memory restrictions
            continue;
        //
        //
        unsigned int varwidth = alnlen;
        //calculate the size of the alignment section...
        int alnlines = nTDP2OutputAlignment;
        if( /*printsss && */qrysssinuse ) {
            unsigned int dbprodst = GetDbProfileField<LNTYPE>(orgprondx, pps2DDist);
            bool trgsssinuse = 
                GetDbProfileFieldPos<LNTYPE>(orgprondx, pmv2DAddrSS, dbprodst) != (LNTYPE)-1;
            if( trgsssinuse )
                alnlines = nTDP2OutputAlignmentSSS;
        }
        alnsize = openlines * (maxopenlinelen+sznl);//opening
        alnsize += headlines * (maxlinelen+sznl);
        alnsize += alnlines * (varwidth + maxfieldlen + sznl + 2) + 2*(maxfieldlen + sznl + 2);//+2 (quotes)
        //if( show )
        alnsize += footlines * (maxlinelen+sznl);
        alnsize += closlines * (maxcloslinelen+sznl);
        *szalnswodesc += alnsize;
        //
        //calculate the size for description...
        GetDbProfileDesc(desc, orgprondx);
        alnlen = (unsigned int)(strlen(desc) + 2);
        alnlen = SLC_MIN(alnlen, dsclen);
        varwidth = alnlen;
        alnsize += (varwidth + maxfieldlen + sznl + 2);
        //
        *szannot += maxfieldlen + sznl/*open*/ + maxcloslinelen + sznl/*close*/ + 
                annotlen + maxfieldlen + sznl + 2/*annotation*/ + 
                annotheadlines * (maxlinelen+sznl)/*score & evalue*/;
        *szalns += alnsize;
    }

    MYMSGBEGl(5)
        char msgbuf[KBYTE];
        mystring strbuf = "CuBatchProcessingFinalizer::GetSizeOfCompressedResultsJSON: ";
        sprintf(msgbuf," szannot %zu szalns %zu (no desc. %zu)",*szannot,*szalns,*szalnswodesc);
        strbuf += msgbuf;
        MYMSG(strbuf.c_str(),5);
    MYMSGENDl
}

