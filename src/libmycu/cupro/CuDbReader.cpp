/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <sys/types.h>
#include <sys/stat.h>

#ifdef OS_MS_WINDOWS
#else
#include <unistd.h>
#endif
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liblib/mydirent.h"

#include "libpro/srcpro/PMProfileModel.h"
#include "libpro/srcpro/Db.h"
#include "libmycu/cupro/IOProfileModel.h"
#include "CuDbReader.h"

// -------------------------------------------------------------------------
// Constructor
//
CuDbReader::CuDbReader( const char* name, bool mapped )
:
    dbname_(name),
    mapped_(mapped),
    no_sequences_(0),
    db_size_(0),
    name_buffer_(NULL),
    openedFile_(Db::DbFlN)
{
    MYMSG("CuDbReader::CuDbReader",4);

    for( int n = 0; n < Db::DbFlN; n++ )
        InvalidateFileDescriptor(n);

    if( !patstrEND || lenstrEND < 1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::CompleteProfileInProfileBuffer: Invalid profile end delimiter.");

    //{{names
    if( !GetDbName())
        throw MYRUNTIME_ERROR( "CuDbReader::CuDbReader: Null database name." );

    name_buffer_ = (char*)malloc( strlen(GetDbName()) + strlen(Db::db_extensions_[Db::DbFlMAIN]) + 1 );

    if( !name_buffer_ )
        throw MYRUNTIME_ERROR( "CuDbReader::CuDbReader: Not enough memory." );

    strcpy( name_buffer_, GetDbName());
    //}}

    //{{page data
    ssize_t pgs = ObtainPageSize();

    if( pgs <= 0 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::CuDbReader: Unable to determine system page size.");

    current_strdata_.pagesize_ = pgs;
    current_strdata_.data_ = NULL;
    if( !GetMapped())
    {
        current_strdata_.data_ = (char*)malloc(pgs);
        if( !current_strdata_.data_ )
            throw MYRUNTIME_ERROR( "CuDbReader::CuDbReader: Not enough memory." );
    }
    //}}

    //{{profile buffer
    profile_buffer_.data_ = NULL;
    profile_buffer_.pagesize_ = 0;
    if( !GetMapped())
    {
        profile_buffer_.data_ = (char*)malloc(MAXPROSIZE);
        if( !profile_buffer_.data_ )
            throw MYRUNTIME_ERROR( "CuDbReader::CuDbReader: Not enough memory." );
        profile_buffer_.pagesize_ = MAXPROSIZE;
    }
    //}}
}

// Default constructor
//
CuDbReader::CuDbReader()
:
    dbname_(NULL),
    mapped_(false),
    no_sequences_(0),
    db_size_(0),
    name_buffer_(NULL),
    openedFile_(Db::DbFlN)
{
    throw MYRUNTIME_ERROR("CuDbReader::CuDbReader: "
    "Default initialization is not allowed.");
}

// -------------------------------------------------------------------------
// Destructor
//
CuDbReader::~CuDbReader()
{
    MYMSG("CuDbReader::~CuDbReader",4);
    Destroy();
}

// -------------------------------------------------------------------------
// Destroy: destroy allocated resources and close open files
//
void CuDbReader::Destroy()
{
    MYMSG("CuDbReader::Destroy",4);
    Close();
    if( name_buffer_ ) {
        free( name_buffer_ );
        name_buffer_ = NULL;
    }
    if(current_strdata_.data_) {
        if( !GetMapped())
            free(current_strdata_.data_);
        current_strdata_.data_ = NULL;
        current_strdata_.pagesize_ = 0;
        ResetCurrentPageData();
    }
    if(profile_buffer_.data_) {
        if( !GetMapped())
            free(profile_buffer_.data_);
        profile_buffer_.data_ = NULL;
        profile_buffer_.pagesize_ = 0;
        ResetProfileBufferData();
    }
}

// -------------------------------------------------------------------------
// GetDbMainName: get the main database name
//
const char* CuDbReader::GetDbMainName()
{
    strcpy( name_buffer_ + strlen(GetDbName()), Db::db_extensions_[Db::DbFlMAIN] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// GetDbMainBinName: get the name of the main binary database file
//
const char* CuDbReader::GetDbMainBinName()
{
    strcpy( name_buffer_ + strlen( GetDbName()), Db::db_extensions_[Db::DbFlMAINBin] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// GetDbProbsName: get the filename of profile pair probabilities 
//
const char* CuDbReader::GetDbProbsName()
{
    strcpy( name_buffer_ + strlen(GetDbName()), Db::db_extensions_[Db::DbFlPPRO] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// Open: open the database and initialize descriptors
//
void CuDbReader::Open()
{
    MYMSG("CuDbReader::Open",4);
    if( file_exists( GetDbMainBinName())) {
        //OpenBin();
        OpenBin_v2_2();
    }
    else {
//         throw MYRUNTIME_ERROR("Unrecognized binary format. Please convert the database using `db2bin'.");
        warning("For an n-fold read speedup, "
        "it is highly RECOMMENDED to use a database in binary format "
        "converted by `db2bin'");
        OpenText();
    }
}

// -------------------------------------------------------------------------
// OpenBin_v2_2: open the database in binary format v2.2 and initialize 
// descriptors
//
void CuDbReader::OpenBin_v2_2()
{
    MYMSG("CuDbReader::OpenBin_v2_2",4);
    myruntime_error mre;

    OpenFile(Db::DbFlMAINBin, GetDbMainBinName());

    SetOpenedFileType(Db::DbFlMAINBin);

    try {
        ResetCurrentPageData();
        ResetProfileBufferData();

        SetDbSizeInBytes( ObtainDbSizeInBytes(GetOpenedFileType()));

        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            sprintf( msgbuf, "CuDbReader::OpenBin_v2_2: file_size %zu page_size %zu",
                    GetDbSizeInBytes(), GetPageSize());
            MYMSG( msgbuf, 3 );
        MYMSGENDl

        if( GetMapped())
            MapFile(GetOpenedFileType());

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if(mre.isset()) {
        SetOpenedFileType(Db::DbFlN);
        Close();
        throw mre;
    }
}

// -------------------------------------------------------------------------
// OpenBin: open the database in binary format and initialize descriptors
//
void CuDbReader::OpenBin()
{
    MYMSG("CuDbReader::OpenBin",4);
    myruntime_error mre;

    OpenFile(Db::DbFlMAINBin, GetDbMainBinName());

    SetOpenedFileType(Db::DbFlMAINBin);

    try {
        ResetCurrentPageData();
        ResetProfileBufferData();

        SetDbSizeInBytes( ObtainDbSizeInBytes(GetOpenedFileType()));

        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            sprintf( msgbuf, "CuDbReader::OpenBin: file_size %zu page_size %zu",
                    GetDbSizeInBytes(), GetPageSize());
            MYMSG( msgbuf, 3 );
        MYMSGENDl

        if( GetMapped())
            MapFile(GetOpenedFileType());

        CacheProfile(true);//cache one page only

        ReadBinaryHeader( &profile_buffer_ );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if(mre.isset()) {
        SetOpenedFileType(Db::DbFlN);
        Close();
        throw mre;
    }
}

// -------------------------------------------------------------------------
// OpenText: open the database in text format and initialize a file 
// descriptor
//
void CuDbReader::OpenText()
{
    MYMSG("CuDbReader::OpenText",4);
    myruntime_error mre;

    OpenFile(Db::DbFlMAIN, GetDbMainName());

    SetOpenedFileType(Db::DbFlMAIN);

    try {
        ResetCurrentPageData();
        ResetProfileBufferData();

        SetDbSizeInBytes( ObtainDbSizeInBytes(GetOpenedFileType()));

        MYMSGBEGl(3)
            char msgbuf[BUF_MAX];
            sprintf( msgbuf, "CuDbReader::OpenText: file_size %zu page_size %zu",
                    GetDbSizeInBytes(), GetPageSize());
            MYMSG( msgbuf, 3 );
        MYMSGENDl

        if( GetMapped())
            MapFile(GetOpenedFileType());

        CacheProfile(true);//cache one page only

        ReadTextHeader( &profile_buffer_ );

    } catch( myexception const& ex ) {
        mre = ex;
    }

    if(mre.isset()) {
        SetOpenedFileType(Db::DbFlN);
        Close();
        throw mre;
    }
}

// -------------------------------------------------------------------------
// Close: close the database
//
void CuDbReader::Close( Db::TFile which )
{
    MYMSG("CuDbReader::Close",3);

    if( Db::DbFlN <= which || Db::DbFlN < 0 ) {
        for( int n = 0; n < Db::DbFlN; n++ ) {
            if( GetMapped())
                UnmapFile(n);
            CloseFile(n);
        }
        db_size_bytes_ = 0;
        ResetCurrentPageData();
        ResetProfileBufferData();
        openedFile_ = Db::DbFlN;
    } else {
        if( GetMapped())
            UnmapFile(which);
        CloseFile(which);
        if( which == Db::DbFlMAIN || which == Db::DbFlMAINBin ) {
            db_size_bytes_ = 0;
            ResetCurrentPageData();
            ResetProfileBufferData();
            openedFile_ = Db::DbFlN;
        }
    }
}

// =========================================================================
// ReadTextHeader: read database header from file
//
void CuDbReader::ReadTextHeader( TCharStream* fp )
{
    MYMSG("CuDbReader::ReadTextHeader",3);

    if( !fp )
        return;

    const mystring preamb = "CuDbReader::ReadTextHeader: ";
    mystring buffer;
    size_t rbts;
    const char* p;
    int emsg;
    long long int llintval;


    //read version number
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( buffer.c_str(), Db::patstrDBVER )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( Db::patstrDBVER );

    if( buffer.length() <= size_t(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( p, pmodel::PMProfileModel::GetDataVersion())) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Inconsistent database version number." );


    //read number of profiles
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( buffer.c_str(), Db::patstrNoPROFS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( Db::patstrNoPROFS );

    if( buffer.length() <= size_t(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( emsg = read_llinteger( p, buffer.length() - (size_t)(p-buffer.c_str()), &llintval, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( llintval < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles." );

    no_sequences_ = llintval;


    //read database size
    if(( emsg = skip_comments( fp, buffer )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || buffer.empty())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( buffer.c_str(), Db::patstrSIZE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( Db::patstrSIZE );

    if( buffer.length() <= size_t(p-buffer.c_str()))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( emsg = read_llinteger( p, buffer.length() - (size_t)(p-buffer.c_str()), &llintval, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( llintval < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid database size." );

    db_size_ = llintval;
}

// -------------------------------------------------------------------------
// ReadBinaryHeader: read database header from binary file
//
void CuDbReader::ReadBinaryHeader( TCharStream* fp )
{
    MYMSG("CuDbReader::ReadBinaryHeader",3);

    if( !fp || !fp->data_ || fp->datlen_ < 1 )
        return;

    const mystring preamb = "CuDbReader::ReadBinaryHeader: ";
    const char* p;
    size_t len;

    p = fp->data_+fp->curpos_;

    //read version number
    len = strlen(Db::patstrDBBINVER);
    if( strncmp(p,Db::patstrDBBINVER,len))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format.");
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format.");

    len = strlen(pmodel::PMProfileModel::GetDataVersion());
    if( strncmp(p,pmodel::PMProfileModel::GetDataVersion(),len))
        throw MYRUNTIME_ERROR( preamb + "Inconsistent database version number.");
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format.");

    //read number of profiles
    len = sizeof(size_t);
    no_sequences_ = *(size_t*)p;
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format.");
    if( no_sequences_ < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles.");

    //read database size
    len = sizeof(uint64_mt);
    db_size_ = *(uint64_mt*)p;
    p += len;
    fp->incpos(len);
    if( fp->eof())
        throw MYRUNTIME_ERROR( preamb + "Wrong database format.");
    if( db_size_ < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid database size.");
}



// =========================================================================
// NextProfileHeader: read the header of the next profile from the database; 
// return a flag indicating whether the database end has been reached; 
//
bool CuDbReader::NextProfileHeader( 
    mystring* desc, mystring* file, 
    int* prolen, int* scale, char** pmdata )
{
    MYMSG("CuDbReader::NextProfileHeader",7);

    if( Eof())
        return false;

    try {
        CacheProfile(false);//cache at least one complete profile

        if( GetOpenedFileType()==Db::DbFlMAINBin)
            BinaryReadProfileHeader( &profile_buffer_, desc, file, prolen, scale, pmdata );
        else {
            TextReadProfileHeader( &profile_buffer_, desc, file, prolen, scale, pmdata );
//             TextReadProfileHeaderBufferless( &profile_buffer_, desc, file, prolen, scale, pmdata );
        }
    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error(ex);
    }

    return true;
}

// -------------------------------------------------------------------------
// Next: read the data of the next profile from the database; 
// return a flag indicating whether the database end has been reached
//
bool CuDbReader::Next( 
    unsigned int distn,
    int profnr, int prolen, int scale, 
    char** pmorgaddr, char** pmdata )
{
    MYMSG("CuDbReader::Next",7);

    if( Eof())
        return false;

    try {
        if( GetOpenedFileType()==Db::DbFlMAINBin)
            BinaryReadProfileData(
                &profile_buffer_, distn, profnr, prolen, scale, 
                pmorgaddr, pmdata );
        else {
            TextReadProfileData(
//             TextReadProfileDataBufferless( 
                &profile_buffer_, distn, profnr, prolen, scale, 
                pmorgaddr, pmdata );
        }
    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error(ex);
    }

    return true;
}





// =========================================================================
// CacheProfileDirect: read one or more profiles (if they fit to the page 
// size) from the main database file;
// pageonly, cache only one page;
//
void CuDbReader::CacheProfileDirect(bool pageonly)
{
    MYMSG("CuDbReader::CacheProfileDirect",5);

    size_t pos = profile_buffer_.curpos_;

    while( !CompleteProfileInProfileBuffer(pos))
    {
        ReadPage(current_strdata_);
        pos = profile_buffer_.datlen_;
        MoveDataFromPageToProfileBuffer();
        if( pageonly )
            break;
    }
}

// -------------------------------------------------------------------------
// ReadPage: read at most one page from the main database file;
// character stream chstr is modified on exit and
// NOTE: the data buffer field of chstr is assumed to be preallocated to the
// page size
//
void CuDbReader::ReadPage( TCharStream& chstr )
{
    MYMSG("CuDbReader::ReadPage",6);

    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::ReadPage: Invalid file descriptor.");

    if( !chstr.data_ || chstr.pagesize_ < 1 )
        throw MYRUNTIME_ERROR("CuDbReader::ReadPage: Invalid argument.");

    if( db_size_bytes_ <= chstr.pageoff_)
        return;

    size_t nbytes = db_size_bytes_ - chstr.pageoff_;

    if( chstr.pagesize_ < nbytes )
        nbytes = chstr.pagesize_;

    ssize_t bytesread;

#ifdef OS_MS_WINDOWS
    DWORD bytesreadhlp = 0;
    if( !ReadFile(
        db_fp_[GetOpenedFileType()],
        chstr.data_,
        (DWORD)nbytes,//nNumberOfBytesToRead (DWORD)
        &bytesreadhlp,//lpNumberOfBytesRead (LPDWORD)
        NULL)//lpOverlapped
      )
        throw MYRUNTIME_ERROR(
        "CuDbReader::ReadPage: Read of the main database file failed.");

    bytesread = bytesreadhlp;
#else
    bytesread = read(db_fp_[GetOpenedFileType()], chstr.data_, nbytes );
#endif

    if(bytesread < 0)
        throw MYRUNTIME_ERROR(
        "CuDbReader::ReadPage: Failed to read the main database file.");

    chstr.datlen_ = bytesread;
    chstr.curpos_ = 0;
    if( 0 < bytesread ) {
        chstr.pagenr_++;
        chstr.pageoff_ += bytesread;
    }
}





// =========================================================================
// ReadDataMapped: get data of size count from the mapped main database 
// file starting at position filepos;
// NOTE: dst is assumed to be preallocated to contain `count' bytes
//
void CuDbReader::ReadDataMapped( size_t filepos, void* dst, size_t count )
{
    MYMSG("CuDbReader::ReadDataMapped",6);

    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataMapped: Invalid file descriptor.");

    if( dst == NULL || count < 1 )
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataMapped: Invalid arguments.");

    if( profile_buffer_.datlen_ < filepos + count )
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataMapped: Invalid file position.");

    memcpy( dst, profile_buffer_.data_ + filepos, count );
}

// -------------------------------------------------------------------------
// ReadDataDirect: read data of size count to dst from the main database 
// file starting at position filepos;
// NOTE: dst is assumed to be preallocated to contain `count' bytes
//
void CuDbReader::ReadDataDirect( size_t filepos, void* dst, size_t count )
{
    MYMSG("CuDbReader::ReadDataDirect",6);

    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataDirect: Invalid file descriptor.");

    if( dst == NULL || count < 1 )
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataDirect: Invalid arguments.");

    if( GetDbSizeInBytes() < filepos + count )
        throw MYRUNTIME_ERROR("CuDbReader::ReadDataDirect: Invalid file position.");

    ssize_t bytesread;

    SeekDbDirect( filepos );

#ifdef OS_MS_WINDOWS
    DWORD bytesreadhlp = 0;
    if( !ReadFile(
        db_fp_[GetOpenedFileType()],
        dst,
        (DWORD)count,//nNumberOfBytesToRead (DWORD)
        &bytesreadhlp,//lpNumberOfBytesRead (LPDWORD)
        NULL)//lpOverlapped
      )
        throw MYRUNTIME_ERROR(
        "CuDbReader::ReadDataDirect: Read of the main database file failed.");

    bytesread = bytesreadhlp;
#else
    bytesread = read(db_fp_[GetOpenedFileType()], dst, count );
#endif

    if(bytesread < 0)
        throw MYRUNTIME_ERROR(
        "CuDbReader::ReadDataDirect: Failed to read the main database file.");
}

// -------------------------------------------------------------------------
// SeekDbDirect: seek the main database by changing current file position 
// (pointer) for reading;
//
void CuDbReader::SeekDbDirect( size_t filepos )
{
    MYMSG("CuDbReader::SeekDbDirect",7);
#ifdef __DEBUG__
    if( GetDbSizeInBytes() < filepos )
        throw MYRUNTIME_ERROR("CuDbReader::SeekDbDirect: Invalid file position.");
    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::SeekDbDirect: Invalid file descriptor.");
#endif

#ifdef OS_MS_WINDOWS
    LARGE_INTEGER offset;
    offset.QuadPart = filepos;
    offset.LowPart = SetFilePointer(
        db_fp_[GetOpenedFileType()],
        offset.LowPart,//lDistanceToMove (LONG), lower double word
        &offset.HighPart,//lpDistanceToMoveHigh (PLONG), upper double word
        FILE_BEGIN
    );
    if( offset.LowPart == INVALID_SET_FILE_POINTER && GetLastError() != NO_ERROR) {
        offset.QuadPart = -1;
        throw MYRUNTIME_ERROR(
        "CuDbReader::SeekDbDirect: Setting file position failed.");
    }
#else
    if( lseek( db_fp_[GetOpenedFileType()], filepos, SEEK_SET ) == 
        (off_t)-1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::SeekDbDirect: Setting file position failed.");
#endif
}
