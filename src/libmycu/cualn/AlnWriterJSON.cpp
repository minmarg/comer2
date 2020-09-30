/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extsp/psl.h"
#include "liblib/fmtdescription.h"

#include "AlnWriter.h"

// -------------------------------------------------------------------------
// WriteResultsJSON: write merged results to file in JSON format
//
void AlnWriter::WriteResultsJSON()
{
    MYMSG( "AlnWriter::WriteResultsJSON", 4 );
    mystring preamb = "AlnWriter::WriteResultsJSON: ";
    static const unsigned int indent = OUTPUTINDENT;
    static const unsigned int annotlen = ANNOTATIONLEN;
    const unsigned int dsclen = MOptions::GetDSCLEN();
    const unsigned int dscwidth = MOptions::GetDSCWIDTH();
    const size_t nhits = MOptions::GetNOHITS();
    const size_t nalns = MOptions::GetNOALNS();
    static const bool printname = false;
    myruntime_error mre;
    mystring filename;
    FILE* fp = NULL;
    char srchinfo[szTmpBuffer];
    char* pb = buffer_, *ptmp = srchinfo;
    int written, left;
    int size = 0;//number of characters written already in the buffer

    GetOutputFilename( filename,
        mstr_set_outdirname_, vec_qryname_[qrysernr_], qrysernr_);

    //use mode b for OS_MS_WINDOWS to not use translation
    if((fp = fopen( filename.c_str(), "wb")) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Failed to open file for writing: " + filename);

    try {
        left = szWriterBuffer;
        size += WritePrognameJSON(pb, left/*maxsize*/, dscwidth);

        left = szWriterBuffer - size;
        left = SLC_MIN(left, (int)dsclen);

        size += WriteQueryDescriptionJSON(pb, left/*maxsize*/, printname,
            vec_nqyposs_[qrysernr_], vec_qryname_[qrysernr_], vec_qrydesc_[qrysernr_],
            dscwidth );

        bool hitsfound = p_finalindxs_ != NULL && p_finalindxs_->size()>0;

        written = WriteSearchInformationJSON(ptmp, szTmpBuffer/*maxsize*/,
            mstr_set_dbname_, mstr_set_prodbsize_, mstr_set_ndbentries_,
            vec_logevthld_[qrysernr_], indent, annotlen, hitsfound );

        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        written = sprintf(srchinfo,"  \"search_summary\": [%s",NL);
        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        const char* recsep = "," NL;//record separator
        const char* finsep = NL;//separator following the last record
        const int lenrecsep = (int)strlen(recsep);
        const int lenfinsep = (int)strlen(finsep);

        if( hitsfound )
        {
            //first, buffer annotations
            for( size_t i = 0; i < p_finalindxs_->size() && i < nhits; i++ ) {
                const char* annot = 
                    ( *vec_annotptrs_[qrysernr_][allsrtvecs_[ (*p_finalindxs_)[i] ]] )
                                               [allsrtindxs_[ (*p_finalindxs_)[i] ]];
                int annotlen = (int)strlen(annot);
                bool last = !(i+1 < p_finalindxs_->size() && i+1 < nhits);
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    annot, annotlen );
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    last? finsep: recsep, last? lenfinsep: lenrecsep );
            }
        }

        written = sprintf(srchinfo,"  ],%s  \"search_hits\": [%s",NL,NL);
        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        if( hitsfound )
        {
            //now, buffer alignments
            for( size_t i = 0; i < p_finalindxs_->size() && i < nalns; i++ ) {
                const char* aln = 
                      ( *vec_alnptrs_[qrysernr_][allsrtvecs_[ (*p_finalindxs_)[i] ]] )
                                               [allsrtindxs_[ (*p_finalindxs_)[i] ]];
                int alnlen = (int)strlen(aln);
                bool last = !(i+1 < p_finalindxs_->size() && i+1 < nhits);
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    aln, alnlen );
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    last? finsep: recsep, last? lenfinsep: lenrecsep );
            }
        }

        written = sprintf(srchinfo,"  ],%s",NL);
        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        ptmp = srchinfo;

        written = WriteSummaryJSON( ptmp,
            refLambda_, refK_,
            expGappedLambda_, expGappedK_,
            vec_nqyposs_[qrysernr_], mstr_set_prodbsize_, mstr_set_ndbentries_,
            vec_sspace_[qrysernr_], vec_deltalen_[qrysernr_] );

        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        //flush the buffer to file
        if( size > 0 )
            WriteToFile( fp, buffer_, size );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( fp )
        fclose(fp);

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// WritePrognameJSON: write the programme name in JSON format;
// outptr, address of the buffer to write;
// maxsize, maximum allowed number of bytes to write;
// returns the number of bytes written;
//
int AlnWriter::WritePrognameJSON( char*& outptr, int maxsize, const int /*width*/ )
{
    //max allowed length for the section of references
    const int maxlenreferences = KBYTE;
    const int maxlenrefinf = maxlenreferences - 60;//max space for references text itself
    int written;
    int size = 0;

    written = sprintf( outptr,
            "{\"comer_search\": {%s"
            "  \"program\": \"%s\",%s"
            "  \"version\": \"%s\",%s",
            NL,PROGNAME? PROGNAME: "",
            NL,PROGVERSION? PROGVERSION: "",NL);
    outptr += written;
    size += written;

    //subtract the maximum reserved space for references
    maxsize -= size + maxlenreferences;

    if(maxsize < 1)
        return size;

    written = sprintf(outptr,"  \"references\": [%s",NL);
    outptr += written;
    size += written;

    int totlenref = 0;
    bool bwrt = true;

    for(int i = 0; !i || (bwrt && PROGREFERENCES[i]); i++)
    {
        const char* strref = PROGREFERENCES[i];
        totlenref += strref? (int)strlen(strref): 0;
        bwrt = totlenref < maxlenrefinf;
        written = sprintf( outptr,
            "    \"%s\"%s%s",
            (bwrt && strref)? strref: "",
            (bwrt && strref && PROGREFERENCES[i+1])? ",":"", NL);
        outptr += written;
        size += written;
    }

    written = sprintf(outptr,"  ],%s",NL);
    outptr += written;
    size += written;

    return size;
}

// -------------------------------------------------------------------------
// WriteQueryDescriptionJSON: write the programme name in JSON format;
// outptr, address of the buffer to write;
// maxsize, maximum allowed number of bytes to write;
// printname, include the name in the output text;
// qrylen, query length;
// name, query name;
// desc, query description;
// width, width to wrap the query description;
// returns the number of bytes written;
//
int AlnWriter::WriteQueryDescriptionJSON( 
    char*& outptr,
    int maxsize,
    bool printname,
    const int qrylen,
    const char* name,
    const char* desc,
    const int /*width*/)
{
    int written;
    int size = 0;

    if( maxsize < 60 )
        return size;

    written = sprintf(outptr,
            "  \"query\": {%s"
            "    \"length\": %d,%s",
            NL,qrylen,NL);
    outptr += written;
    size += written;

    //subtract BUF_MAX for accounting for the two fields preceding description
    maxsize -= size + BUF_MAX;

    if( maxsize < 1 )
        return size;

    char* p = outptr;
    int outpos = 0;

    written = sprintf(outptr,"    \"name\": \"");
    outptr += written;
    outpos += written;
    if( printname && name)
        FormatDescriptionJSON(outptr, name, maxsize, outpos);
    written = sprintf( outptr,"\",%s",NL);
    outptr += written;
    outpos += written;

    written = sprintf(outptr,"    \"description\": \"");
    outptr += written;
    outpos += written;
    if( desc )
        FormatDescriptionJSON(outptr, desc, maxsize, outpos);
    written = sprintf( outptr,"\"%s  },%s",NL,NL);
    outptr += written;
    outpos += written;

    written = (int)(outptr - p);
    size += written;

    return size;
}

// -------------------------------------------------------------------------
// WriteSearchInformationJSON: write database information in JSON format;
// outptr, address of the buffer to write;
// maxsize, maximum allowed number of bytes to write;
// name, database name;
// dbsize, database size;
// ndbentries, number of profiles in the database;
// logevthrld, log e-value threshold;
// indent, indentation length;
// annotlen, annotation length;
// found, whether any profiles have been found;
// returns the number of bytes written;
//
int AlnWriter::WriteSearchInformationJSON( 
    char*& outptr,
    int maxsize,
    const char* dbname,
    const size_t dbsize,
    const size_t ndbentries,
    const float logevthrld,
    const int /*indent*/,
    const int /*annotlen*/,
    const bool found )
{
    static const int envlpines = 2;//number of lines wrapping a database record
    static const int headlines = 4;//number of lines for information
    static const int maxfieldlen = 40;//JSON max field length
    static const int maxlinelen = 90;//maximum length of lines other than alignment lines
    int written, szname = 0, namelen = 0;
    int size = 0;

    maxsize -= envlpines * maxfieldlen + headlines * maxlinelen;
    if( dbname ) {
        dbname = my_basename(dbname);
        szname = namelen = (int)strlen(dbname);
    }
    if( maxsize < szname )
        szname = maxsize;

    if( szname < 1 )
        return size;

    written = sprintf(outptr,
            "  \"database\": {%s"
            "    \"name\": \"",NL);
    outptr += written;
    size += written;
    strncpy( outptr, my_basename(dbname), szname );
    outptr += szname;
    size += szname;
    if(szname < namelen && 3 < szname)
        for(int i=1; i<=3; i++) *(outptr-i) = '.';
    written = sprintf( outptr,"\",%s",NL);
    outptr += written;
    size += written;

    written = sprintf(outptr,
            "    \"number_of_profiles\": %zu,%s"
            "    \"number_of_positions\": %zu%s"
            "  },%s",
            ndbentries,NL,dbsize,NL,NL);
    outptr += written;
    size += written;

    if( found ) {
        written = sprintf( outptr,
            "  \"message\": \"Profiles found below the e-value threshold:\",%s",NL);
        outptr += written;
        size += written;
    }
    else {
        written = sprintf( outptr,
            "  \"message\": \"No profiles found below an e-value threshold of %g\",%s",
            exp(logevthrld),NL);
        outptr += written;
        size += written;
    }

    return size;
}

// -------------------------------------------------------------------------
// WriteSummaryJSON: write summary of parameters and database information in
// JSON format;
// outptr, address of the output buffer;
// refLambda, refK, and...
// expGappedLambda, expGappedK, statistical parameters;
// qrylen, query length;
// dbsize, database size;
// ndbentries, number of profiles in the database;
// sspace, search space;
// deltalen, value of length correction;
// returns the number of bytes written;
//
int AlnWriter::WriteSummaryJSON( 
    char*& outptr,
    const float refLambda, const float refK,
    const float expGappedLambda, const float expGappedK,
    const int qrylen,
    const size_t dbsize,
    const size_t ndbentries,
    const float sspace,
    const int deltalen )
{
    int written;
    int size = 0;

    int     eff_query_length = qrylen;
    size_t  eff_db_length = dbsize;

    if( deltalen < qrylen )
        eff_query_length -= deltalen;
    if( ndbentries * deltalen < dbsize )
        eff_db_length -= ndbentries * deltalen;

    written = sprintf( outptr,
            "  \"search_statistics\": {%s"
            "    \"reference_K_ungapped\": %.4f,%s"
            "    \"reference_lambda_ungapped\": %.4f,%s"
            "    \"reference_K_gapped\": %.4f,%s"
            "    \"reference_lambda_gapped\": %.4f,%s"
            "    \"query_length\": %d,%s"
            "    \"database_size\": %zu,%s"
            "    \"effective_query_length\": %d,%s"
            "    \"effective_database_size\": %zu,%s"
            "    \"effective_search_space\": %zu%s"
            "  }%s}}%s",
            NL,refK,NL,refLambda,NL,expGappedK,NL,expGappedLambda,NL,
            qrylen,NL,dbsize,NL,eff_query_length,NL,eff_db_length,NL,(size_t)sspace,NL,
            NL,NL);
    outptr += written;
    size += written;

    return size;
}
