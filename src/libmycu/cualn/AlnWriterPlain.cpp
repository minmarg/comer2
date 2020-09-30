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
// WriteResultsPlain: write merged results to file
void AlnWriter::WriteResultsPlain()
{
    MYMSG( "AlnWriter::WriteResultsPlain", 4 );
    mystring preamb = "AlnWriter::WriteResultsPlain: ";
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
        size += WritePrognamePlain(pb, left/*maxsize*/, dscwidth);

        left = szWriterBuffer - size;
        left = SLC_MIN(left, (int)dsclen);

        size += WriteQueryDescriptionPlain(pb, left/*maxsize*/, printname,
            vec_nqyposs_[qrysernr_], vec_qryname_[qrysernr_], vec_qrydesc_[qrysernr_],
            dscwidth );

        bool hitsfound = p_finalindxs_ != NULL && p_finalindxs_->size()>0;

        written = WriteSearchInformationPlain(ptmp, szTmpBuffer/*maxsize*/,
            mstr_set_dbname_, mstr_set_prodbsize_, mstr_set_ndbentries_,
            vec_logevthld_[qrysernr_], indent, annotlen, hitsfound );

        BufferData( fp,
            buffer_, szWriterBuffer, pb, size,
            srchinfo, written );

        if( hitsfound )
        {
            //first, buffer annotations
            for( size_t i = 0; i < p_finalindxs_->size() && i < nhits; i++ ) {
                const char* annot = 
                    ( *vec_annotptrs_[qrysernr_][allsrtvecs_[ (*p_finalindxs_)[i] ]] )
                                               [allsrtindxs_[ (*p_finalindxs_)[i] ]];
                int annotlen = (int)strlen(annot);
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    annot, annotlen );
            }

            written = sprintf(srchinfo,"%s%s%s",NL,NL,NL);
            BufferData( fp,
                buffer_, szWriterBuffer, pb, size,
                srchinfo, written );

            //now, buffer alignments
            for( size_t i = 0; i < p_finalindxs_->size() && i < nalns; i++ ) {
                const char* aln = 
                      ( *vec_alnptrs_[qrysernr_][allsrtvecs_[ (*p_finalindxs_)[i] ]] )
                                               [allsrtindxs_[ (*p_finalindxs_)[i] ]];
                int alnlen = (int)strlen(aln);
                BufferData( fp,
                    buffer_, szWriterBuffer, pb, size,
                    aln, alnlen );
            }
        }

        ptmp = srchinfo;

        written = WriteSummaryPlain( ptmp,
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
// WriteProgname: write the programme name;
// outptr, address of the buffer to write;
// maxsize, maximum allowed number of bytes to write;
// width, width of text to wrap;
// returns the number of bytes written;
//
int AlnWriter::WritePrognamePlain( char*& outptr, int maxsize, const int width )
{
    static const int sznl = (int)strlen(NL) * 3 + 1;
    int written;
    int size = 0;

    if( PROGNAME ) {
        written = sprintf( outptr, "%s", PROGNAME );
        outptr += written;
        size += written;
    }
    if( PROGVERSION ) {
        written = sprintf( outptr, " %s", PROGVERSION );
        outptr += written;
        size += written;
    }
    if( PROGNAME || PROGVERSION ) {
        written = sprintf( outptr, "%s%s", NL,NL);
        outptr += written;
        size += written;
    }

    maxsize -= size + sznl;

    if(maxsize < 1)
        return size;

    char* p = outptr;
    int outpos = 0;
    int linepos = 0;
    char addsep = 0;
    for(int i = 0; PROGREFERENCES[i] && outpos < maxsize; i++) {
        FormatDescription(
            outptr, PROGREFERENCES[i],
            maxsize, 0/*indent*/, width, outpos, linepos, addsep);
        outpos += PutNL(outptr);
        linepos = 0;
    }

    written = sprintf( outptr, "%s", NL);
    outptr += written;

    written = (int)(outptr - p);
    size += written;

    return size;
}

// -------------------------------------------------------------------------
// WriteProgname: write the programme name;
// outptr, address of the buffer to write;
// maxsize, maximum allowed number of bytes to write;
// printname, include the name in the output text;
// qrylen, query length;
// name, query name;
// desc, query description;
// width, width to wrap the query description;
// returns the number of bytes written;
//
int AlnWriter::WriteQueryDescriptionPlain( 
    char*& outptr,
    int maxsize,
    bool printname,
    const int qrylen,
    const char* name,
    const char* desc,
    const int width )
{
    int written;
    int size = 0;

    if( maxsize < 40 )
        return size;

    written = sprintf( outptr, "%s Query (%d positions):%s", NL, qrylen, NL);
    outptr += written;
    size += written;

    maxsize -= size;

    if( maxsize < 1 )
        return size;

    char* p = outptr;
    int outpos = 0;
    int linepos = 0;
    char addsep = 0;
    if( printname && name) {
        FormatDescription( 
            outptr, name,
            maxsize, 0/*indent*/, width, outpos, linepos, addsep);
        addsep = ' ';
    }
    if( desc )
        FormatDescription( 
            outptr, desc,
            maxsize, 0/*indent*/, width, outpos, linepos, addsep);

    written = (int)(outptr - p);
    size += written;

    return size;
}

// -------------------------------------------------------------------------
// WriteDbInformation: write database information;
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
int AlnWriter::WriteSearchInformationPlain( 
    char*& outptr,
    int maxsize,
    const char* dbname,
    const size_t dbsize,
    const size_t ndbentries,
    const float logevthrld,
    const int indent,
    const int annotlen,
    const bool found )
{
    int written, szname = 0, namelen = 0;
    int size = 0;

    maxsize -= 20 + 100 + 2 * indent + annotlen + 100;
    if( dbname ) {
        dbname = my_basename(dbname);
        szname = namelen = (int)strlen(dbname);
    }
    if( maxsize < szname )
        szname = maxsize;

    if( szname < 1 )
        return size;

    written = sprintf( outptr, "%s%s Database:%s",NL,NL,NL);
    outptr += written;
    size += written;

    strncpy( outptr, my_basename(dbname), szname );
    outptr += szname;
    size += szname;

    if(szname < namelen && 3 < szname)
        for(int i=1; i<=3; i++) *(outptr-i) = '.';

    size += PutNL(outptr);

    for(int i=0; i<indent; i++, size++ ) *outptr++ = ' ';

    written = sprintf( outptr, "%zu profiles%s",ndbentries,NL);
    outptr += written;
    size += written;

    for(int i=0; i<indent; i++, size++ ) *outptr++ = ' ';

    written = sprintf( outptr, "%zu total positions%s%s%s",dbsize,NL,NL,NL);
    outptr += written;
    size += written;

    if( found ) {
        written = sprintf( outptr," Profiles found below the e-value threshold:");
        outptr += written;
        size += written;
        for(; written < annotlen; written++, size++ ) *outptr++ = ' ';
        written = sprintf( outptr," %7s %7s%s%s","Score","E-value",NL,NL);
        outptr += written;
        size += written;
    }
    else {
        written = sprintf( outptr,
            " No profiles found below an e-value threshold of %g.%s%s%s",
            exp(logevthrld),NL,NL,NL);
        outptr += written;
        size += written;
    }

    return size;
}

// -------------------------------------------------------------------------
// WriteDbInformation: write database information;
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
int AlnWriter::WriteSummaryPlain( 
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

    written = sprintf( outptr,"%-20s  %-6s   %-6s%s","Reference values of","K","Lambda",NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"%-20s  %6.4f   %6.4f%s","Ungapped", refK, refLambda,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"%-20s  %6.4f   %6.4f%s%s","Gapped", expGappedK, expGappedLambda,NL,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"Length of query, %d%s", qrylen,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"Length of database, %zu%s", dbsize,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"Effective length of query, %d%s", eff_query_length,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"Effective length of database, %zu%s", eff_db_length,NL);
    outptr += written;
    size += written;
    written = sprintf( outptr,"Effective search space, %zu%s%s", (size_t)sspace,NL,NL);
    outptr += written;
    size += written;

    return size;
}
