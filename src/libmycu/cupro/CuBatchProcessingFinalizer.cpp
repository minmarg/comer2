/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <memory>
#include <utility>
#include <functional>
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <thread>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "liblib/fmtdescription.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/btckcoords.h"
#include "libmycu/cualn/Devices.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/CLOptions.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/PMBatchProData.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "libmycu/cuss/CuBatchSS_com.h"

#include "CuBatchProcessingFinalizer.h"

static cudaStream_t streamdummy;

// _________________________________________________________________________
// class statics
//

// _________________________________________________________________________
// Class CuBatchProcessingFinalizer
//
// Constructor
//
CuBatchProcessingFinalizer::CuBatchProcessingFinalizer(
    cudaStream_t& strcopyres,
    DeviceProperties dprop, 
    AlnWriter* writer,
    Configuration* config,
    const mystring* queryfnames,
    const mystring* querydescs,
    const mystring* bdb1fnames,
    const mystring* bdb1descs )
:   tobj_(NULL),
    req_msg_(CUBPTHREAD_MSG_UNSET),
    rsp_msg_(CUBPTHREAD_MSG_UNSET),
    //
    strcopyres_(strcopyres),
    dprop_(dprop),
    alnwriter_(writer),
    config_(config),
    //cached data:
    cached_queryfnames_(queryfnames),
    cached_querydescs_(querydescs),
    cached_bdb1fnames_(bdb1fnames),
    cached_bdb1descs_(bdb1descs),
    //data arguments:
    cubp_set_qrysernr_(-1),
    cubp_set_nqyposs_(0),
    cubp_set_qyeno_(0.0f),
    cubp_set_deltalen_(0U),
    cubp_set_sspace_(0.0f),
    cubp_set_logevthld_(0.0f),
    cubp_set_bdbC_(nullptr),
    cubp_set_ndb1pros_(0UL),
    cubp_set_ndbCpros_(0UL),
    cubp_set_querprosOmtd_(0UL),
    cubp_set_ndb1prosOmtd_(0UL),
    cubp_set_ndbCprosOmtd_(0UL),
    cubp_set_querposoffset_(0UL),
    cubp_set_bdb1posoffset_(0UL),
    cubp_set_bdbCposoffset_(0UL),
    cubp_set_nposits_(0U),
    cubp_set_npros_(0U),
    cubp_set_h_results_(NULL),
    cubp_set_sz_alndata_(0UL),
    cubp_set_sz_alns_(0UL),
    cubp_set_sz_dbalnlen2_(0),
    annotations_(nullptr),
    alignments_(nullptr)
{
    MYMSG( "CuBatchProcessingFinalizer::CuBatchProcessingFinalizer", 3 );
    memset( cubp_set_querypmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_querypmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdb1pmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdb1pmend_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdbCpmbeg_, 0, pmv2DTotFlds * sizeof(void*));
    memset( cubp_set_bdbCpmend_, 0, pmv2DTotFlds * sizeof(void*));
    tobj_ = new std::thread( &CuBatchProcessingFinalizer::Execute, this, (void*)NULL );
}

// Default constructor
//
CuBatchProcessingFinalizer::CuBatchProcessingFinalizer()
:   tobj_(NULL),
    strcopyres_(streamdummy),
    config_(NULL)
{
    throw MYRUNTIME_ERROR(
    "CuBatchProcessingFinalizer::CuBatchProcessingFinalizer: "
    "Default initialization is prohibited.");
}

// Destructor
//
CuBatchProcessingFinalizer::~CuBatchProcessingFinalizer()
{
    MYMSG( "CuBatchProcessingFinalizer::~CuBatchProcessingFinalizer", 3 );
    if( tobj_ ) {
        tobj_->join();
        delete tobj_;
        tobj_ = NULL;
    }
}

// -------------------------------------------------------------------------
// Execute: thread's starting point and execution process
//
void CuBatchProcessingFinalizer::Execute( void* )
{
    MYMSG( "CuBatchProcessingFinalizer::Execute", 3 );
    myruntime_error mre;

    try {
        MYCUDACHECK( cudaSetDevice( dprop_.devid_ ));
        MYCUDACHECKLAST;

        while(1) {
            //wait until a message arrives
            std::unique_lock<std::mutex> lck_msg(mx_dataccess_);

            cv_msg_.wait(lck_msg,
                [this]{return 
                    ((0 <= req_msg_ && req_msg_ <= cubpthreadmsgTerminate) || 
                     req_msg_ == CUBPTHREAD_MSG_ERROR
                    );}
            );

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "CuBatchProcessingFinalizer::Execute: Msg %d",req_msg_);
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //thread owns the lock after the wait;
            //read message req_msg_
            int reqmsg = req_msg_;

            //unset the message to avoid live cycle when starting over the loop
            req_msg_ = CUBPTHREAD_MSG_UNSET;

            //set response msg to error upon occurance of an exception
            rsp_msg_ = CUBPTHREAD_MSG_ERROR;
            int rspmsg = rsp_msg_;

            switch(reqmsg) {
                case cubpthreadmsgFinalize:
                        //data addresses have been written already;
                        //make sure data transfer has finished
                        MYCUDACHECK( cudaStreamSynchronize(strcopyres_));
                        MYCUDACHECKLAST;
                        ;;
                        CompressResults();
                        SortCompressedResults();
                        PassResultsToWriter();
                        ;;
                        //parent does not wait for a response nor requires data to read;
                        //unset response code
                        rspmsg = CUBPTHREAD_MSG_UNSET;//cubptrespmsgFinalizing;
                        break;
                case cubpthreadmsgTerminate:
                        rspmsg = cubptrespmsgTerminating;
                        break;
                default:
                        rspmsg = CUBPTHREAD_MSG_UNSET;
                        break;
            };

            MYMSGBEGl(3)
                char msgbuf[BUF_MAX];
                sprintf( msgbuf, "CuBatchProcessingFinalizer::Execute: Msg %d Rsp %d",reqmsg, rspmsg );
                MYMSG( msgbuf, 3 );
            MYMSGENDl

            //save response code
            rsp_msg_ = rspmsg;

//             //unlock the mutex and notify the parent using the cv
//             lck_msg.unlock();
//             cv_msg_.notify_one();

            if( reqmsg < 0 || reqmsg == cubpthreadmsgTerminate)
                //terminate execution
                break;
        }

    } catch( myruntime_error const& ex ) {
        mre = ex;
    } catch( myexception const& ex ) {
        mre = ex;
    } catch( ... ) {
        mre = myruntime_error("Unknown exception caught.");
    }

    if( mre.isset()) {
        error( mre.pretty_format().c_str());
        SetResponseError();
//         cv_msg_.notify_one();
        return;
    }
}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// PassResultsToWriter: transfer the addresses of sorted results to the 
// alignment writer;
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::PassResultsToWriter()
{
    MYMSG( "CuBatchProcessingFinalizer::PassResultsToWriter", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::PassResultsToWriter: ";

    if( alnwriter_ ) {
        alnwriter_->PushPartOfResults(
            cubp_set_qrysernr_,
            cubp_set_nqyposs_,
            cached_queryfnames_? cached_queryfnames_[cubp_set_querprosOmtd_].c_str(): NULL,
            cached_querydescs_? cached_querydescs_[cubp_set_querprosOmtd_].c_str(): NULL,
            cubp_set_deltalen_,
            cubp_set_sspace_,
            cubp_set_logevthld_,
            std::move(annotations_),
            std::move(alignments_),
            std::move(srtindxs_),
            std::move(logevalues_),
            std::move(alnptrs_),
            std::move(annotptrs_)
        );
    }
    cubp_set_bdbC_.reset();
}

// -------------------------------------------------------------------------
// SortCompressedResults: sort formatted results;
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::SortCompressedResults()
{
    MYMSG( "CuBatchProcessingFinalizer::SortCompressedResults", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::SortCompressedResults: ";

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(preamb + "Null compressed results.");

    if( srtindxs_->size() !=logevalues_->size() ||
        logevalues_->size() != alnptrs_->size() ||
        logevalues_->size() != annotptrs_->size())
        throw MYRUNTIME_ERROR(preamb + "Inconsistent result sizes.");

    std::sort(srtindxs_->begin(), srtindxs_->end(),
        [this](size_t n1, size_t n2) {
            return (*logevalues_)[n1] < (*logevalues_)[n2];
        });
}

// -------------------------------------------------------------------------
// CompressResults: process results obtained; compress them for passing 
// them to the writing thread;
// NOTE: all operations performed under lock
//
void CuBatchProcessingFinalizer::CompressResults()
{
    MYMSG( "CuBatchProcessingFinalizer::CompressResults", 4 );
    static const unsigned int indent = OUTPUTINDENT;
    static const unsigned int annotlen = ANNOTATIONLEN;
    static const bool printsss = MOptions::GetSSSWGT() > 0.0f;
    const unsigned int dsclen = MOptions::GetDSCLEN();
    const unsigned int dscwidth = MOptions::GetDSCWIDTH();
    const unsigned int alnwidth = MOptions::GetALNWIDTH();
    const unsigned int width = indent < dscwidth? dscwidth: 1;
    const int show = MOptions::GetSHOW();
    size_t szannot = 0UL;
    size_t szalns = 0UL;
    size_t szalnswodesc = 0UL;
    const mystring* name;//profile name
    const mystring* desc;//profile description
    static const bool printname = false;
    int written, sernr = 0;//initialize serial number here
    //
    bool qrysssinuse = 
        GetQueryFieldPos<LNTYPE>(pmv2DAddrSS, 0/*cubp_set_querposoffset_*/) != (LNTYPE)-1;

    GetSizeOfCompressedResults( &szannot, &szalns, &szalnswodesc);

    annotations_.reset();
    alignments_.reset();

    ReserveVectors( cubp_set_npros_ );

    if( szalns < szannot || 
        szalnswodesc > 2 * (cubp_set_sz_alndata_ + cubp_set_sz_alns_))
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResults: "
        "Size of compressed results is unusually large.");

    if( szannot < 1 || szalns < 1 )
        return;

    annotations_.reset((char*)std::malloc(szannot));
    alignments_.reset((char*)std::malloc(szalns));

    if( !annotations_ || !alignments_)
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResults: "
        "Not enough memory.");

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::CompressResults: "
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
        //save the addresses of the annotations and the alignment records
        srtindxs_->push_back(sernr++);
        logevalues_->push_back(logeval);
        annotptrs_->push_back(annptr);
        alnptrs_->push_back(outptr);
        //
        double evalue = GetEvalue(logeval);
        float score = GetOutputAlnDataField<float>(prondx, dp2oadScore);
        //get the name and description
        GetDbProfileNameDesc(name, desc, orgprondx);
        //make an annotation
        MakeAnnotation( annptr,
            name, desc, printname, annotlen, annotlen/*width*/,
            score, evalue);
        *annptr++ = 0;//end of record
        //
        //put the name and description...
        int outpos = 0;
        int linepos = 0;
        char addsep = '>';
        if( printname && name ) {
            FormatDescription( 
                outptr, name->c_str(),
                dsclen, indent, width, outpos, linepos, addsep);
            addsep = ' ';
        }
        if( desc )
            FormatDescription( 
                outptr, desc->c_str(),
                dsclen, indent, width, outpos, linepos, addsep);
        PutNL(outptr);
        //
        //compress the alignment and relative information...
        FormatScores( outptr,
            prondx, orgprondx, alnlen, score, logeval, evalue, dbprolen );
        FormatAlignment( outptr,
            prondx, orgprondx, dbpro2dst, alnlen, dbprolen, alnwidth,
            printsss, qrysssinuse );
        if( show )
            FormatFooter( outptr, prondx );
        //
        written = sprintf( outptr,"%s%s",NL,NL);
        outptr += written;
        *outptr++ = 0;//end of record
    }
}

// -------------------------------------------------------------------------
// FormatDescription: format profile description;
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
void CuBatchProcessingFinalizer::MakeAnnotation( 
    char*& outptr,
    const mystring* name,
    const mystring* desc,
    const bool printname,
    const int maxoutlen,
    const int width,
    const float score,
    const double evalue) const
{
    char* p = outptr;
    int outpos = 0;
    int linepos = 0;
    char addsep = 0;
    int written;
    if( printname && name ) {
        FormatDescription( 
            outptr, name->c_str(),
            maxoutlen, 0/*indent*/, width, outpos, linepos, addsep);
        addsep = ' ';
    }
    if( desc )
        FormatDescription( 
            outptr, desc->c_str(),
            maxoutlen, 0/*indent*/, width, outpos, linepos, addsep);

    written = (int)(outptr - p);
    if( written < maxoutlen )
        for(int i = written; i < maxoutlen; i++) *outptr++ = ' ';

    written = sprintf( outptr," %7.0f %7.1g%s",score,evalue,NL);
    outptr += written;
}

// -------------------------------------------------------------------------
// outptr, pointer to the output buffer;
// prondx, profile index in the results list;
// orgprondx, profile index over all pm data structures;
// alnlen, alignment length;
// score, alignment score;
// logeval, log e-value;
// evalue, e-value;
// dbprolen, db profile length;
inline
void CuBatchProcessingFinalizer::FormatScores(
    char*& outptr,
    unsigned int prondx,
    unsigned int orgprondx,
    unsigned int alnlen,
    float score,
    float logeval,
    double evalue,
    int dbprolen )
{
    int written;
    float dbproENO = GetDbProfileField<FPTYPE>(orgprondx, pps2DENO);
    int psts = (int)GetOutputAlnDataField<float>(prondx, dp2oadPstvs);
    int idts = (int)GetOutputAlnDataField<float>(prondx, dp2oadIdnts);
    int gaps = (int)GetOutputAlnDataField<float>(prondx, dp2oadNGaps);
    written = 
    sprintf( outptr,"  Length/ENO: Query = %d/%.1f, Sbjct = %d/%.1f%s%s",
            cubp_set_nqyposs_, cubp_set_qyeno_, dbprolen, dbproENO, NL,NL);
    outptr += written;
    written = 
    sprintf( outptr," Score = %.2f (%.1f bits),  Expect = %.2g, P-value = %.2g%s",
            score, GetBitScore(logeval), evalue, GetPvalue(evalue), NL);
    outptr += written;
    if( idts > 0 ) {
        written = 
            sprintf( outptr," Identities = %d/%d (%d%%)",idts,alnlen,idts*100/alnlen);
        outptr += written;
    }
    if( psts > 0 ) {
        if( idts ) {
            written = sprintf( outptr,",");
            outptr += written;
        }
        written = 
            sprintf( outptr," Positives = %d/%d (%d%%)",psts,alnlen,psts*100/alnlen);
        outptr += written;
    }
    if( gaps > 0 ) {
        if( idts || psts ) {
            written = sprintf( outptr,",");
            outptr += written;
        }
        written = sprintf( outptr," Gaps = %d/%d (%d%%)",gaps,alnlen,gaps*100/alnlen);
        outptr += written;
    }
    written = sprintf( outptr,"%s%s",NL,NL);
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
void CuBatchProcessingFinalizer::FormatAlignment(
    char*& outptr,
    unsigned int prondx,
    unsigned int orgprondx,
    unsigned int dbpro2dst,
    int alnlen,
    int dbprolen,
    const int width,
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
    const char* p;

    int written;
    int nbytes;
    int fgaps;

    printsss = printsss && qrysssinuse && trgsssinuse;

    for( int f = 0; f < alnlen; f += width, palnbeg += width )
    {
        nbytes = alnlen - f;
        if( width < nbytes )
            nbytes = width;
        if(printsss) {
            p = GetAlnSectionAt(palnbeg, dp2oaQuerySSS);
            written = sprintf( outptr,"%-13s","struct");
            outptr += written;
            strncpy(outptr, p, nbytes);
            outptr += nbytes;
            PutNL(outptr);
        }
        {   p = GetAlnSectionAt(palnbeg, dp2oaQuery);
            written = sprintf( outptr,"Query: %5u ", qrybeg);
            outptr += written;
            strncpy(outptr, p, nbytes);
            outptr += nbytes;
            fgaps = (int)(std::count(p, p+nbytes, '-'));
            qrybeg += nbytes - fgaps;
            written = sprintf( outptr," %-5d%s", nbytes<=fgaps? qrybeg: qrybeg-1,NL);
            outptr += written;
        }
        {   p = GetAlnSectionAt(palnbeg, dp2oaMiddle);
            written = sprintf( outptr,"%13c",' ');
            outptr += written;
            strncpy(outptr, p, nbytes);
            outptr += nbytes;
            PutNL(outptr);
        }
        {   p = GetAlnSectionAt(palnbeg, dp2oaTarget);
            written = sprintf( outptr,"Sbjct: %5u ", trgbeg);
            outptr += written;
            strncpy(outptr, p, nbytes);
            outptr += nbytes;
            fgaps = (int)(std::count(p, p+nbytes, '-'));
            trgbeg += nbytes - fgaps;
            written = sprintf( outptr," %-5d%s", nbytes<=fgaps? trgbeg: trgbeg-1,NL);
            outptr += written;
        }
        if(printsss) {
            p = GetAlnSectionAt(palnbeg, dp2oaTargetSSS);
            written = sprintf( outptr,"%-13s","struct");
            outptr += written;
            strncpy(outptr, p, nbytes);
            outptr += nbytes;
            PutNL(outptr);
        }
        PutNL(outptr);
    }
}

// -------------------------------------------------------------------------
// outptr, pointer to the output buffer;
// prondx, profile index in the results list;
inline
void CuBatchProcessingFinalizer::FormatFooter(
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
    static const char* na = "n/a";
    char kbuf[BUF_MAX], lbuf[BUF_MAX];
    const char* kp, *lp;

    written = sprintf( outptr,"%25c  %-6s   %-6s%s",' ',"K","Lambda",NL);
    outptr += written;
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
    written = sprintf( outptr,"%-25s  %6s   %6s%s","Computed  ungapped,",kp,lp,NL);
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
    written = sprintf( outptr,"%-25s  %6s   %6s%s","Estimated gapped,",kp,lp,NL);
    outptr += written;
    written = sprintf( outptr,"Entropy, %6.4f; Expected, %6.4f; Min/Max, %.0f/%-.0f%s%s",
                    entropy,expected,minsc,maxsc,NL,NL);
    outptr += written;
}



// -------------------------------------------------------------------------
// GetSizeOfCompressedResults: get total size required for annotations and 
// complete alignments;
// szannot, size of annotations;
// szalns, size of complete alignments (with descriptions);
// szalnswodesc, size of alignments without descriptions;
//
inline
void CuBatchProcessingFinalizer::GetSizeOfCompressedResults(
    size_t* szannot, size_t* szalns, size_t* szalnswodesc) const
{
    MYMSG( "CuBatchProcessingFinalizer::GetSizeOfCompressedResults", 5 );
    static const unsigned int sznl = (int)strlen(NL);
    static const unsigned int indent = OUTPUTINDENT;
    static const bool printsss = MOptions::GetSSSWGT() > 0.0f;
    const unsigned int dsclen = MOptions::GetDSCLEN();
    const unsigned int dscwidth = MOptions::GetDSCWIDTH();
    const unsigned int alnwidth = MOptions::GetALNWIDTH();
    const unsigned int width = indent < dscwidth? dscwidth-indent: 1;
    static const int show = MOptions::GetSHOW();
    static const unsigned int headlines = 3;//number of lines for scores, etc.
    static const unsigned int footlines = 4;//number of lines for statistical parameters
    static const unsigned int maxlinelen = 140;//maximum length of lines other than alignment lines
    int alnsize;//alignment size
    const mystring* name;//profile name
    const mystring* desc;//profile description
    //
    bool qrysssinuse = 
        GetQueryFieldPos<LNTYPE>(pmv2DAddrSS, 0/*cubp_set_querposoffset_*/) != (LNTYPE)-1;

    if( cubp_set_bdbCpmbeg_[0] && cubp_set_bdbCpmend_[0] && !cubp_set_bdbC_)
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::GetSizeOfCompressedResults: Null bdbC object.");

#ifdef __DEBUG__
    if( !cached_queryfnames_ || !cached_querydescs_ ||
        !cached_bdb1fnames_ || !cached_bdb1descs_ ||
        (cubp_set_bdbC_ && (!cubp_set_bdbC_->GetFnames() || !cubp_set_bdbC_->GetDescs())))
        throw MYRUNTIME_ERROR(
        "CuBatchProcessingFinalizer::GetSizeOfCompressedResults: Null data.");
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
// fprintf(stderr,"%f\n%d\n%d\n %d\n%x\n%x\n"
// "%.4f\n%.4f\n%.4f\n%x %x %x %x\n%.1f\n%.1f\n%.4f\n%.4f\n%.4f\n%.4f\n"
// " %d\n %d\n",
// GetOutputAlnDataField<float>(prondx,dp2oadScore),
// GetOutputAlnDataField<int>(prondx,dp2oadOrgProNo),
// GetOutputAlnDataField<int>(prondx,dp2oadProNewDst),
// (int)GetOutputAlnDataField<float>(prondx,dp2oadPstvs),
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadPstvs))[0],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadPstvs))[1],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadPstvs))[2],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadPstvs))[3],
// GetOutputAlnDataField<int>(prondx,dp2oadBegCoords),
// GetOutputAlnDataField<int>(prondx,dp2oadEndCoords),
// GetOutputAlnDataField<float>(prondx,dp2oadLmbdEst),
// GetOutputAlnDataField<float>(prondx,dp2oadKEst),
// GetOutputAlnDataField<float>(prondx,dp2oadEvalue),
// // (int)GetOutputAlnDataField<float>(prondx,dp2oadEpA),
// ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadEpA))[0],
// ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadEpA))[1],
// ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadEpA))[2],
// ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadEpA))[3],
// GetOutputAlnDataField<float>(prondx,dp2oadMin),
// GetOutputAlnDataField<float>(prondx,dp2oadMax),
// GetOutputAlnDataField<float>(prondx,dp2oadE),
// GetOutputAlnDataField<float>(prondx,dp2oadH),
// GetOutputAlnDataField<float>(prondx,dp2oadLmbd),
// GetOutputAlnDataField<float>(prondx,dp2oadK),
// (int)GetOutputAlnDataField<float>(prondx,dp2oadIdnts),
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadIdnts))[0],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadIdnts))[1],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadIdnts))[2],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadIdnts))[3],
// (int)GetOutputAlnDataField<float>(prondx,dp2oadNGaps)
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadNGaps))[0],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadNGaps))[1],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadNGaps))[2],
// // ((char*)((float*)cubp_set_h_results_ + nTDP2OutputAlnData*prondx+dp2oadNGaps))[3]
// );
        //
        unsigned int varwidth = alnlen < alnwidth? alnlen: alnwidth;
        //calculate the size of the alignment section...
        int alnfrags = (alnlen + alnwidth - 1)/alnwidth;
        int alnlines = nTDP2OutputAlignment;
        if( printsss && qrysssinuse ) {
            unsigned int dbprodst = GetDbProfileField<LNTYPE>(orgprondx, pps2DDist);
            bool trgsssinuse = 
                GetDbProfileFieldPos<LNTYPE>(orgprondx, pmv2DAddrSS, dbprodst) != (LNTYPE)-1;
            if( trgsssinuse )
                alnlines = nTDP2OutputAlignmentSSS;
        }
        alnsize = 3 * sznl;//alignment separator
        alnsize += headlines * (maxlinelen+sznl) + 2 * sznl;
        alnsize += alnfrags * alnlines * (varwidth + 2*indent + sznl) + (alnfrags+1) * sznl;
        if( show )
            alnsize += footlines * (maxlinelen+sznl);
        *szalnswodesc += alnsize;
        //
        //calculate the size for description...
        GetDbProfileNameDesc(name, desc, orgprondx);
// /**TEST*/fprintf(stderr,"name= %s prondx= %u cubp_set_logevthld_= %f logeval= %f cubp_set_nqyposs_= %d dbprolen= %d "
// "alnlen= %u alnsize= %d alnfrags= %d alnsize= %d\n",
// name->c_str(),prondx,cubp_set_logevthld_,logeval,cubp_set_nqyposs_,dbprolen,alnlen,alnsize,alnfrags,alnsize);
        alnlen = (unsigned int)(name->length() + desc->length() + 2);
        alnlen = SLC_MIN(alnlen, dsclen);
        varwidth = alnlen < dscwidth? alnlen: dscwidth;
        alnfrags = (alnlen + width - 1)/width;
        alnsize += alnfrags * (varwidth + sznl) + sznl;
        //
        *szannot += maxlinelen + sznl;
        *szalns += alnsize;
    }

    MYMSGBEGl(5)
        char msgbuf[KBYTE];
        mystring strbuf = "CuBatchProcessingFinalizer::GetSizeOfCompressedResults: ";
        sprintf(msgbuf," szannot %zu szalns %zu (no desc. %zu)",*szannot,*szalns,*szalnswodesc);
        strbuf += msgbuf;
        MYMSG(strbuf.c_str(),5);
    MYMSGENDl
}



// -------------------------------------------------------------------------
// PrintCompressedResults: print formatted alignments
//
void CuBatchProcessingFinalizer::PrintCompressedResults() const
{
    MYMSG( "CuBatchProcessingFinalizer::PrintCompressedResults", 4 );
    mystring preamb = "CuBatchProcessingFinalizer::PrintCompressedResults: ";

    if( !srtindxs_ || !logevalues_ || !alnptrs_ || !annotptrs_ )
        throw MYRUNTIME_ERROR(preamb + "Null compressed results.");

    if( srtindxs_->size() !=logevalues_->size() ||
        logevalues_->size() != alnptrs_->size() ||
        logevalues_->size() != annotptrs_->size())
        throw MYRUNTIME_ERROR(preamb + "Inconsistent result sizes.");

    for(size_t i = 0; i < srtindxs_->size(); i++ ) {
        fprintf(stdout,"%s",(*annotptrs_)[(*srtindxs_)[i]]);
    }

    fprintf(stdout,"%s",NL);

    for(size_t i = 0; i < srtindxs_->size(); i++ ) {
        fprintf(stdout,"%f%s",(*logevalues_)[(*srtindxs_)[i]],NL);
        fprintf(stdout,"%s",(*alnptrs_)[(*srtindxs_)[i]]);
    }
}
