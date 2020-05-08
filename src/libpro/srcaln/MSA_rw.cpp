/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <memory>

#include "extsp/psl.h"
#include "extsp/rng.h"
#include "liblib/msg.h"
#include "liblib/mysort.h"
#include "liblib/BinarySearchStructure.h"
#include "MSA.h"

// -------------------------------------------------------------------------
// SetName: set the name of the MSA
//
void MSA::SetName( const char* newname )
{
    size_t newlength = strlen( newname );
    if( !newlength )
        throw MYRUNTIME_ERROR( "MSA::SetName: No name." );

    if( name_ )
        free( name_ );

    name_ = (char*)malloc( newlength + 1 );

    if( !name_ )
        throw MYRUNTIME_ERROR( "MSA::SetName: Not enough memory." );

    strncpy( name_, newname, newlength + 1 );//include the terminal null character
}

// -------------------------------------------------------------------------
// SetDescription: set the description of the MSA
//
void MSA::SetDescription( const char* newdesc )
{
    size_t newlength = strlen( newdesc );
    if( !newlength )
        throw MYRUNTIME_ERROR( "MSA::SetDescription: No description." );

    if( description_ )
        free( description_ );

    description_ = (char*)malloc( newlength + 1 );

    if( !description_ )
        throw MYRUNTIME_ERROR( "MSA::SetDescription: Not enough memory." );

    //include the terminal null character
    strncpy( description_, newdesc, newlength + 1 );
}

// -------------------------------------------------------------------------
// AppendDescription: append text to the description of the MSA; if the 
// description is null, the text is simply copied to the description
//
void MSA::AppendDescription( const char* newdesc, size_t beg, size_t end )
{
    if( end < beg )
        return;

    size_t  from = 0;
    size_t  newlength = strlen( newdesc );

    if( !newlength )
        throw MYRUNTIME_ERROR( "MSA::AppendDescription: No description." );

    if( newlength < beg )
        return;

    if( newlength < end )
        end = newlength;

    newlength = end - beg + 1;

    if( description_ ) {
        size_t curlength = strlen( description_ );
        description_ = (char*)realloc( description_, curlength + newlength + 1 );
        from = curlength;
    } else
        description_ = (char*)malloc( newlength + 1 );

    if( !description_ )
        throw MYRUNTIME_ERROR( "MSA::AppendDescription: Not enough memory." );

    strncpy( description_ + from, newdesc + beg, newlength );
    description_[ from + newlength ] = 0;
}

// -------------------------------------------------------------------------

enum {
    eninfmtSTOCKHOLM,
    eninfmtFASTA,
    eninfmtNO
};

// -------------------------------------------------------------------------
// ReadAlignment: read multiple sequence alignment
//
void MSA::ReadAlignment( const char* filename )
{
    myruntime_error mre;
    int informat = eninfmtNO;

    try {
        //try STOCKHOLM
        ReadSTOCKHOLM1( filename, &informat );
    } catch( myexception const& ex1 ) {
        mre = ex1;
        if( informat == eninfmtSTOCKHOLM )
            throw mre;
        try {
            //try FASTA
            ReadFASTA( filename );
        } catch( myexception const& ex2 ) {
            mre = ex2;
            throw mre;
        }
    }
}

// -------------------------------------------------------------------------
// ReadFASTA: read multiple sequence alignment in FASTA format from file
// TODO: test line ends on Windows OS
// FIXME: AppendDescription, the last argument is incorrect
//
void MSA::ReadFASTA( const char* filename )
{
    const mystring preamb = "MSA::ReadFASTA: ";
    myruntime_error mre;
    FILE* fp = fopen( filename, "r" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( preamb + "Failed to open input file.");

    char        p = 0;
    char        buffer[KBYTE];
    int         toread = KBYTE;
    size_t      length = 0;     //length of sequences
    size_t      cread = 1;      //number of characters read by the last operation
    size_t      begin = 0;      //beginning of the description substring
    size_t      pos = 0;        //position counter
    BMSequence  fake( scnDefNoPoss );//fake sequence to delineate gap positions in query (the first sequence)
    BMSequence* one = NULL;
    bool        desc = false;

    clear();

    one = new BMSequence( scnDefNoPoss );
    if( one == NULL )
        throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

    try {
#if 1 
        SetName( my_basename( filename ));
#else
        // or use the posix function
        SetName( basename( filename ));
#endif
        while( !feof(fp) && cread )
        {
            cread = fread( buffer, 1, toread - 1, fp );
            buffer[cread] = 0;
            begin = 0;

            for( size_t n = 0; n < cread; n++, p = buffer[n-1] )
            {
                if(( !p || p == '\r' || p == '\n' ) && buffer[n] == '>' ) {
                    begin = n + 1;
                    desc = true;
                    pos = 0;

                    if( one->GetSize()) {
                        if( !GetDescription()) {
                            //allow null descriptions in FASTA files
                            ;//throw MYRUNTIME_ERROR( preamb + 
                             //"Wrong file format: No description.");
                        }
                        if( length ) {
                            if( one->GetSize() != length )
                                throw MYRUNTIME_ERROR( preamb + 
                                "Wrong file format: Lengths of sequences are not equal.");
                        } else
                            length = one->GetSize();

                        push( one );
                        one = new BMSequence( length );
                        if( one == NULL )
                            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
                    }
                    continue;
                }

                if( desc ) {
                    // if the end of buffer or the beginning of the next line is reached,
                    // save the description of the first sequence
                    if( n + 1 == cread || buffer[n] == '\r' || buffer[n] == '\n' )
                        if( n ) {
                            size_t  bufbeg = ( cread <= begin )? 0: begin;
                            size_t  buflen = ( (buffer[n]=='\r' || buffer[n]=='\n')? n-1: n) - 
                                              bufbeg + 1;
                            if( !length && buflen )
                                AppendDescription( buffer, bufbeg, buflen );
                            if( GetKeepSeqDescriptions() && one )
                                one->AppendDescription( buffer, bufbeg, buflen );
                        }

                    if( buffer[n] == '\r' || buffer[n] == '\n' )
                        desc = false;
                    continue;
                }

                if( !desc ) {
                    if( buffer[n]==' ' || buffer[n]=='\t' || buffer[n]=='\r' || buffer[n]=='\n' )
                        continue;

                    if( GetIgnoreGapsInQuery()) {
                        //if this is the first sequence
                        if( !length ) {
                            fake.push( HashAlphSymbol( buffer[n] ));
                        }
                        if( fake.GetResidueAt( pos ) != GAP )
                            //push alignment symbols
                            one->push( HashAlphSymbol( buffer[n] ));
                    }
                    else
                        //push alignment symbols
                        one->push( HashAlphSymbol( buffer[n] ));    
                    pos++;
                }
            }
        }

        if( !feof( fp ))
            warning(( preamb + "Some data left unprocessed in input file." ).c_str());

        if( one->GetSize()) {
            if( !GetDescription()) {
                //allow null descriptions in FASTA files
                ;//throw MYRUNTIME_ERROR( preamb + "Wrong file format: No description." );
            }
            if( length && one->GetSize() != length )
                throw MYRUNTIME_ERROR( preamb + 
                "Wrong file format: Lengths of sequences are not equal.");
            push( one );
            one = NULL;
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    fclose( fp );
    if( one )
        delete one;

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------

static const char*  patstrFMTSTO = "# STOCKHOLM";
static const char*  patstrSTOTER = "//";
static const char*  patstrSTOGF = "#=GF";
static const char*  patstrSTOGFfts[] = {"AC","DE"};
static const int    noSTOGFfts = 2;
static const char*  patstrSTOGS = "#=GS";
// static const char* patstrSTOGSDE = "DE";
static const char*  patstrSTOGC = "#=GC";
static const char*  patstrSTOseqcons = "seq_cons";
static const size_t lenstrFMTSTO = strlen( patstrFMTSTO );

// -------------------------------------------------------------------------
// ReadSTOCKHOLM1: read multiple sequence alignment in Stockholm 1.0 format
// NOTE: reading MSA in ReadSTOCKHOLM adds an additional sequence generated 
// according to the state sequence; this additional sequence affects 
// (find out how) the profile model when the MSA contains a small number of 
// sequences and if the sequence information is included in the process to 
// infer the profile model
//
void MSA::ReadSTOCKHOLM1( const char* filename, int* format )
{
    if( format )
        *format = eninfmtNO;

    const mystring preamb = "MSA::ReadSTOCKHOLM1: ";
    myruntime_error mre;
    FILE* fp = fopen( filename, "r" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( preamb + "Failed to open input file.");

    size_t          length;//, rbts;
    const size_t    defsize = KBYTE;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    bool            stoter = false;
    bool            descin = false;
    bool            circled = false;
    const char*     p, *pp;
    char            ch;
    int             emsg, n;
    size_t          sind;
    bool            statesempty;
    mystring        seqcons, states;
    mystring        line, str, *sname;
    const mystring* svn;
    SimpleVector    svhdngs( KBYTE );//headings
    SimpleVector    svnames( KBYTE );//names
    SimpleVector    svaseqs( KBYTE );//alignment sequences
    BMSequence*     svs, *seq = NULL;

    clear();
    states.reserve(TIMES4(KBYTE));

#if 1 
    SetName( my_basename( filename ));
#else
    // or use the posix function
    SetName( basename( filename ));
#endif

    //read header of file
    if(( emsg = skip_comments( fp, locbuffer, SLC_MIN( locsize, lenstrFMTSTO+1 ), &length, 0 )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

    if(( p = strstr( locbuffer, patstrFMTSTO )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong file format.");

    //STOCKHOLM format
    if( format )
        *format = eninfmtSTOCKHOLM;
    sind = 0;
    circled = false;

    try {
        while( !feof( fp )) {
            //read full line
            if(( emsg = skip_comments( fp, line, 0 )) != 0 )
                throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

            if( feof( fp ) || line.empty())
                continue;

            if( line[0] == '\n' || line[0] == '\r')
                continue;

            for( n = (int)line.length() - 1; 
                 0 <= n && ( line[n] == '\n' || line[n] == '\r' ); n-- ) line[n] = 0;

            //save family description if available 
            if(( p = (char*)strstr( line.c_str(), patstrSTOGF )) != NULL ) {
                p += strlen( patstrSTOGF );
                for( ; *p == ' ' || *p == '\t' ; p++ );
                for( n = 0; n < noSTOGFfts; n++ ) {
                    if( strncmp( p, patstrSTOGFfts[n], strlen( patstrSTOGFfts[n])) == 0 ) {
                        pp = p + strlen( patstrSTOGFfts[n]);
                        for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                        if( strlen( pp ))
                            str = mystring( pp );
                        if( !str.empty()) {
                            if( GetDescription()) 
                                AppendDescription(" ");
                            AppendDescription( str.c_str());
                            descin = true;
                        }
                        break;
                    }
                }
                continue;
            }

            //save description if not saved 
            if(( p = ( char* )strstr( line.c_str(), patstrSTOGS )) != NULL ) {
                pp = p + strlen( patstrSTOGS );
                for( p = pp; *p == ' ' || *p == '\t' ; p++ );
                for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
                str = mystring( p, n );
                for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                //if(( p = strstr( pp, patstrSTOGSDE )) != NULL ) {
                //    pp = p + strlen( patstrSTOGSDE );
                for( ; *pp && *pp != ' ' && *pp != '\t' ; pp++ );//jump over <feature>
                //
                for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                if( strlen( pp ))
                    str += ' ' + mystring( pp );
                if( !str.empty() && !descin )
                    AppendDescription( str.c_str());
                svhdngs.Push( new mystring( str ));
                //}
                descin = true;
                continue;
            }

            //if end of MSA exit the loop
            if( strlen( patstrSTOTER ) <= strlen( line.c_str()))
                if( strncmp( line.c_str(), patstrSTOTER, strlen( patstrSTOTER )) == 0 ) {
                    stoter = true;
                    if( !feof( fp ))
                        if(( emsg = skip_comments( fp, line, 0 )) != 0 )
                            throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                    break;
                }

            //read consensus sequence seq_cons
            if(( p = ( char* )strstr( line.c_str(), patstrSTOGC )) != NULL ) {
                p += strlen( patstrSTOGC );
                for( ; *p == ' ' || *p == '\t' ; p++ );
                if( strncmp( p, patstrSTOseqcons, strlen( patstrSTOseqcons )) == 0 ) {
                    p += strlen( patstrSTOseqcons );
                    for( ; *p == ' ' || *p == '\t' ; p++ );
                    for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
                    if( 0 < n )
                        seqcons += mystring( p, n );
                }
                continue;
            }

            if( line[0] == '#' )
                continue;

            //process alignment sequences
            p = line.c_str();
            for( ; *p == ' ' || *p == '\t' ; p++ );
            for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
            if( n <= 0 )
                throw MYRUNTIME_ERROR( preamb + "Wrong file format: No sequence name.");

            //check name
            sname = new mystring( p, n );
            if( sname == NULL )
                throw MYRUNTIME_ERROR( preamb + "Not enough memory.");

            if( svnames.GetSize()) {
                if(( svn = ( mystring* )svnames.GetValueAt( 0 )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Memory access error.");
                if( *svn == *sname ) {
                    sind = 0;
                    circled = true;
                }
            }

            if( circled ) {
                if( svnames.GetSize() <= sind || svaseqs.GetSize() <= sind )
                    throw MYRUNTIME_ERROR( preamb + "Invalid sequence index.");
                if(( svn = ( mystring* )svnames.GetValueAt( sind )) == NULL  || 
                   ( svs = ( BMSequence* )svaseqs.GetValueAt( sind )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Memory access error.");
                if( *svn != *sname )
                    throw MYRUNTIME_ERROR( preamb + 
                    "Wrong file format: Inconsistency between names.");
                delete sname;
                sname = NULL;
            }
            else {
                svnames.Push( sname );
                sname = NULL;
                if(( svs = new BMSequence( defsize )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
                push( svs );
                svaseqs.Push( svs );
                if(( svn = ( mystring* )svhdngs.GetValueAt( sind )) == NULL )
                    throw MYRUNTIME_ERROR( preamb + "Memory access error.");
                if( GetKeepSeqDescriptions())
                    svs->AppendDescription( svn->c_str(), 0, svn->length());
                delete svn;
                svhdngs.SetValueAt( sind, svn = NULL );
            }

            sind++;
            statesempty = (sind == 1);//states.empty();//NOTE:go on the state sequence once cycled

            //save alignment sequence
            for( p = pp; *p == ' ' || *p == '\t' ; p++ );
            for( pp = p; *pp && *pp != ' ' && *pp != '\t'; pp++ ) {
                if( 'A' <= *pp && *pp <= 'Z' )
                    ch = *pp;
                else if( 'a' <= *pp && *pp <= 'z' )
                    ch = *pp;
                else if( *pp == '-' || *pp == '.' || *pp == '_' || *pp == '~' )
                    ch = '-';
                else
                    throw MYRUNTIME_ERROR( preamb + "Unrecognized alignment symbol.");
                svs->push( HashAlphSymbol(ch));
                if( statesempty ) {
                    if(('A' <= *pp && *pp <= 'Z') || *pp == '-')
                        states.append('1');//match or delete state
                    else
                        states.append('0');
                }
            }
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    if( sname ) {
        delete sname;
        sname = NULL;
    }
    for( n = 0; n < ( int )svhdngs.GetSize(); n++ )
        if(( svn = ( mystring* )svhdngs.GetValueAt( n )) != NULL ) {
            delete svn;
            svhdngs.SetValueAt( n, NULL );
        }
    for( n = 0; n < ( int )svnames.GetSize(); n++ )
        if(( svn = ( mystring* )svnames.GetValueAt( n )) != NULL ) {
            delete svn;
            svnames.SetValueAt( n, NULL );
        }

    if( stoter && !feof( fp ))
        warning(( preamb + "Only a single MSA processed." ).c_str());
    fclose( fp );

    if( !stoter && !mre.isset())
        warning(( preamb + "Unterminated file." ).c_str());

    if( mre.isset())
        throw mre;

    if( !GetDescription())
        warning(( preamb + "No MSA description." ).c_str());

    if( 0 < GetSize()) {
        seq = GetSequenceAt(0);
        if( seq == NULL )
            throw MYRUNTIME_ERROR( preamb + "Memory access error.");
    }

    //check alignment sequence lengths
    for( n = 1; n < (int)GetSize(); n++ ) {
        svs = GetSequenceAt(n);
        if( !svs || !seq )
            throw MYRUNTIME_ERROR( preamb + "Memory access error.");
        if( svs->GetSize() != seq->GetSize())
            throw MYRUNTIME_ERROR( preamb + "Lengths of aligned sequences differ.");
    }

    //check consensus sequence length, format and save it
//     if( seq && !seqcons.empty()) {
//         if( seq->GetSize() != seqcons.length())
//             throw MYRUNTIME_ERROR( preamb + "Inconsistent length of consensus sequence.");
//         //move the first sequence to the end
//         push( seq );
//         //allocate space for consensus sequence in the first place
//         if(( svs = new BMSequence( seq->GetSize())) == NULL )
//             throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
//         SetSequenceAt( 0, svs );
//         if( GetKeepSeqDescriptions())
//             svs->AppendDescription( GetDescription(), 0, strlen( GetDescription()));
//         //translate and save consensus sequence
//         TranslateSTOConsSequence( seqcons, svs );
//     }

    //check state sequence length, format and save it;
    //NOTE: should be consistent with ssp2.pl
    if( seq ) {
        if( seq->GetSize() != states.length())
            throw MYRUNTIME_ERROR( preamb + "Inconsistent length of state sequence.");
        //move the first sequence to the end
        push( seq );
        //allocate space for consensus sequence in the first place
        if(( svs = new BMSequence( seq->GetSize())) == NULL )
            throw MYRUNTIME_ERROR( preamb + "Not enough memory.");
        SetSequenceAt( 0, svs );
        if( GetKeepSeqDescriptions())
            svs->AppendDescription( GetDescription(), 0, strlen( GetDescription()));
        //process state sequence, make it the support sequence of profile
        TranslateSTOStates( states, svs );
    }
}

// -------------------------------------------------------------------------
// TranslateSTOConsSequence: translate STOCKHOLM consensus sequence
//
//CntCmp: private function: count comparator
//
static inline
int CntCmp( const void* vector, size_t n1, size_t n2 )
{
    const size_t* vec = ( const size_t* )vector;
    if( vec == NULL )
        return 0;
    //[n2]-[n1], to sort in descending order
    return (int)(vec[n2] - vec[n1]);
}
//
void MSA::TranslateSTOConsSequence( const mystring& seqcons, BMSequence* seq )
{
    const mystring preamb = "MSA::TranslateSTOConsSequence: ";
    size_t  n, s; 
    int     slcted;
    char    ch; 
    unsigned char rr;
    char    bufstr[KBYTE];
    size_t  resdst[NUMAA];
    size_t  srtind[NUMAA];
    size_t  chind;
    extspsl::MTRng rng;
    BMSequence* svs;

    const mystring code("olach-p+sut");
    const size_t lencode = code.length();
    //mystring    codeRES[lencode];
	std::unique_ptr<mystring[]> codeRES(new mystring[lencode]);
	n = 0;
    codeRES[n++] = "ST";              //o, alcohol
    codeRES[n++] = "ILV";             //l, aliphatic
    codeRES[n++] = "FHWY";            //a, aromatic
    codeRES[n++] = "DEHKR";           //c, charged
    codeRES[n++] = "ACFGHIKLMRTVWY";  //h, hydrophobic
    codeRES[n++] = "DE";              //-, negative
    codeRES[n++] = "CDEHKNQRST";      //p, polar
    codeRES[n++] = "HKR";             //+, positive
    codeRES[n++] = "ACDGNPSTV";       //s, small
    codeRES[n++] = "AGS";             //u, tiny
    codeRES[n++] = "ACDEGHKNQRST";    //t, turnlike

    if( seq == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null argument.");

    seq->clear();

    rng.Set((unsigned int)(size_t)(&rng));// +(unsigned long)tm );

    for( n = 0; n < seqcons.length(); n++ ) {
        ch = seqcons[n];
        if( 'A' <= ch && ch <= 'Z' ) {
            seq->push( HashAlphSymbol(ch));
            continue;
        }
        if( '.' == ch ) {
            //any residue, gap
            seq->push( HashAlphSymbol('-'));
            continue;
        }
        if(( chind = code.find(ch)) == mystring::npos ) {
            sprintf( bufstr, "Unrecognized consensus symbol: %c.", ch );
            throw MYRUNTIME_ERROR( preamb + bufstr );
        }
        //find out residue distribution at the position
        memset( resdst, 0, NUMAA*sizeof(size_t));
        //the 0th position reserved for seq!
        for( s = 1; s < GetSize(); s++ ) {
            svs = GetSequenceAt(s);
            if( svs == NULL )
                throw MYRUNTIME_ERROR( preamb + "Memory access error.");
            rr = svs->GetResidueAt(n);
            if( NUMAA <= rr )
                continue;
            resdst[rr]++;
        }
        //sort by count
        HeapSortInd( srtind, resdst, NUMAA, CntCmp );
        //select a residue observed max number of times
        slcted = 0;
        for( s = 0; s < NUMAA; s++ ) {
            rr = ( unsigned char )srtind[s];
            if( NUMAA <= rr )
                throw MYRUNTIME_ERROR( preamb + "Memory access error.");
            if( resdst[rr] < 1 )
                break;
            if( codeRES[chind].find(DehashCode(rr)) != mystring::npos ) {
                seq->push( rr );
                slcted = 1;
                break;
            }
        }
        //if the consensus symbol is inconsistent with the residues observed most often,
        //randomly select a residues from the corresponding class
        if( !slcted ) {
            s = (size_t)( rng.GetSingle()*(float)(codeRES[chind].length()));
            if( codeRES[chind].length() <= s )
                s--;
            if( codeRES[chind].length() <= s )
                throw MYRUNTIME_ERROR( preamb + "Memory access error.");
            rr = HashAlphSymbol( codeRES[chind][s]);
            seq->push( rr );
            sprintf( bufstr, "Unmatched consensus symbol '%c' at position %zu.", ch, n );
            warning(( preamb + bufstr ).c_str());
        }
    }
}

// -------------------------------------------------------------------------
// TranslateSTOStates: process state sequence for making the support 
// sequence of profile
//
void MSA::TranslateSTOStates( const mystring& states, BMSequence* seq )
{
    const mystring preamb = "MSA::TranslateSTOStates: ";
    //minimum fraction of delete state symbols per position to consider a
    // position to be in the delete state
    const float frcDEL = 0.7f;
    //calculate delete states using frcDEL as a criterion for the minimum 
    // fraction of gaps per column;
    // if false, all match (or delete) states are considered match states 
    // and the delete state probability for each column is calculated while 
    // processing the profile; otherwise (true), delete states will become
    // insertions during processing of the profile
    const bool cbCALCDELSTATES = false;
    //
    size_t  n, s; 
    int     slcted;
    unsigned char rr;
    char    ch; 
    char    bufstr[KBYTE];
    size_t  resdst[NUMALPH];
    size_t  srtind[NUMALPH];
    size_t  restot;
    BMSequence* svs;

    if( seq == NULL )
        throw MYRUNTIME_ERROR( preamb + "Null argument.");

    seq->clear();

    for( n = 0; n < states.length(); n++ ) {
        ch = states[n];
        if( '1' != ch ) {
            //any residue, gap
            seq->push( HashAlphSymbol('-'));
            continue;
        }
        //find out residue distribution at the position
        memset( resdst, 0, NUMALPH*sizeof(size_t));
        restot = 0;
        //the 0th position reserved for seq!
        for( s = 1; s < GetSize(); s++ ) {
            svs = GetSequenceAt(s);
            if( svs == NULL )
                throw MYRUNTIME_ERROR( preamb + "Memory access error.");
            rr = svs->GetResidueAt(n);
            if( NUMALPH <= rr )
                continue;
            resdst[rr]++;
            restot++;
        }
        //sort by count
        HeapSortInd( srtind, resdst, NUMALPH, CntCmp );
        //if GAP is observed most often
        rr = (unsigned char)srtind[0];
        if( GAP == rr && restot && cbCALCDELSTATES ) {
            if( frcDEL < (float)resdst[rr] / (float)restot ) {
                seq->push( HashAlphSymbol('-'));
                continue;
            }
        }
        //select a residue observed max number of times
        slcted = 0;
        for( s = 0; s < NUMALPH; s++ ) {
            rr = ( unsigned char )srtind[s];
            if( GAP == rr )
                continue;
            if( NUMALPH <= rr )
                throw MYRUNTIME_ERROR( preamb + "Memory access error.");
            if( resdst[rr] < 1 )
                break;
            if( NUMAA <= rr && s+1 < NUMALPH && 
                srtind[s+1] < NUMAA && resdst[rr] == resdst[srtind[s+1]])
                //count of actual amino acid is the same
                continue;
            seq->push( rr );
            slcted = 1;
            break;
        }
        //if no residue selected
        if( !slcted ) {
            //delete state columns containing only GAP may occur
            seq->push( HashAlphSymbol('-'));
            continue;
            sprintf( bufstr, "Unmatched state at position %zu.", n );
            throw MYRUNTIME_ERROR( preamb + bufstr );
        }
    }
}





// -------------------------------------------------------------------------
// PutAlignment: write MSA to file
//
void MSA::PutAlignment( const char* filename )
{
    FILE* fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "MSA::PutAlignment: Failed to open file for writing." );

    BMSequence* one = NULL;
    size_t i, p;

    for( i = 0; i < GetSize(); i++ ){
        one = GetSequenceAt(i);

        if( i == 0 && protoprofile_ )
            one = protoprofile_;

        if( !one || !one->GetUsed())
            continue;

        if( !i )
            fprintf( fp, ">%s%s", GetDescription(), NL );
        else {
            if( GetKeepSeqDescriptions() && one->GetDescription())
                fprintf( fp, ">%s%s", one->GetDescription(), NL );
            else
                fprintf( fp, ">Sequence %zu%s", i+1, NL );
        }

        if( protoprofile_ ) {
            for( p = 0; p < protoprofile_->GetSize(); p++ ) {
                if( !protoprofile_->IsUsedAt( p ))
                    continue;
                putc( DehashCode(one->GetResidueAt(p)), fp );
            }
        }
        else {
            for( p = 0; p < one->GetSize(); p++ ) {
                putc( DehashCode(one->GetResidueAt(p)), fp );
            }
        }
        fprintf( fp, "%s", NL );
    }

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputWeightedFrequencies: output observed weighted frequencies
//
void MSA::OutputWeightedFrequencies( const char* filename )
{
    FILE*   fp = stdout;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::OutputWeightedFrequencies: Null prototype profile." );

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "MSA::OutputWeightedFrequencies: "
                               "Failed to open file for writing." );

    protoprofile_->PrintMatchWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputPSSM: output PSSM matrix
//
void MSA::OutputPSSM( const char* filename )
{
    FILE* fp = stdout;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::OutputPSSM: Null prototype profile." );

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "MSA::OutputPSSM: Failed to open file for writing." );

    protoprofile_->PrintPSSMatrix( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputSuppressedProfile: output PSSM matrix and weighted frequencies
//
void MSA::OutputSuppressedProfile( const char* filename )
{
    FILE* fp = stdout;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::OutputSuppressedProfile: Null prototype profile." );

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "MSA::OutputSuppressedProfile: Failed to open file for writing." );

    protoprofile_->PrintSuppressedPSSMandWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputProfile: output profile information in text format
//
void MSA::OutputProfile( const char* filename )
{
    FILE* fp = stdout;

    if( !protoprofile_ )
        throw MYRUNTIME_ERROR( "MSA::OutputProfile: Null prototype profile." );

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw MYRUNTIME_ERROR( "MSA::OutputProfile: Failed to open file for writing." );

    if( GetName()) {
        fprintf( fp, "Multiple alignment, %s%s", GetName(), NL );
    }
    if( GetDescription()) {
        fprintf( fp, "Description, %s%s%s", GetDescription(), NL, NL );
    }

    protoprofile_->PrintProfile( fp );

    if( fp != stdout )
        fclose( fp );
}
