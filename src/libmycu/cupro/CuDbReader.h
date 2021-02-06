/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuDbReader_h__
#define __CuDbReader_h__

#ifdef OS_MS_WINDOWS
#   include <Windows.h>
#else
#   include <sys/types.h>
#   include <sys/stat.h>
#   include <sys/mman.h>
#   include <unistd.h>
#   include <fcntl.h>
#endif

#include <stdio.h>
#include <string.h>

#include "extsp/psl.h"
#include "liblib/mybase.h"
#include "liblib/myfiler.h"
#include "libpro/srcpro/Db.h"

extern "C" {
extern const char* patstrEND;
extern const size_t lenstrEND;
}

// _________________________________________________________________________
// Class CuDbReader
//
// Reader of a database of profiles
//
class CuDbReader
{
    enum {
        CUDBREADER_TMPBUFFSIZE = 4096,
        MAXPROSIZE = 32 * ONEM
    };

public:
    CuDbReader( const char* name, bool mapped = true );
    virtual ~CuDbReader();

    virtual void Destroy();

    virtual void Open();//open database
    virtual void Close( Db::TFile = Db::DbFlN );//close a file of the database

    bool GetMapped() const {return mapped_;}


    bool Eof() const {return GetMapped()? EofMapped(): EofDirect();}


    bool NextProfileHeader( mystring* desc, mystring* file,
            int* prolen, int* scale, char** pmdata );
    bool Next( unsigned int distn,
            int profnr, int prolen, int scale, 
            char** pmorgaddr, char** pmdata );


    const char* GetDbMainName();//get database main filename
    const char* GetDbMainBinName();//get database binary filename
    const char* GetDbProbsName();//get probabilities filename

    void GetDbPos( TCharStream* chstr ) const {
        if( GetMapped())
            GetDbPosMapped(chstr);
        else
            GetDbPosDirect(chstr);
    }
    void SetDbPos( const TCharStream* chstr ) {
        if( GetMapped())
            SetDbPosMapped(chstr);
        else
            SetDbPosDirect(chstr);
    }

    const char*     GetDbName() const { return dbname_; }
    size_t          GetNoSequences() const { return no_sequences_; }
    uint64_mt       GetDbSize() const { return db_size_; }

    size_t GetDbSizeInBytes() const {return db_size_bytes_;}

    void ReadData( size_t filepos, void* dst, size_t count ) {
        if( GetMapped())
            ReadDataMapped( filepos, dst, count );
        else
            ReadDataDirect( filepos, dst, count );
    }

    explicit CuDbReader();

protected:
    void OpenText();//open database in text format
    void OpenBin_v2_2();//open database v2.2 in binary format
    void OpenBin();//open database in binary format

    void SetOpenedFileType(Db::TFile which) {openedFile_ = which;}
    Db::TFile GetOpenedFileType() const {return openedFile_;}

    void SetDbPosToOrigin() {
        if( GetMapped())
            SetDbPosToOriginMapped();
        else
            SetDbPosToOriginDirect();
    }

    void ReadBinaryHeader( TCharStream* );
    void ReadTextHeader( TCharStream* );

    void CacheProfile(bool pageonly) {
        if( GetMapped())
            CacheProfileMapped(pageonly);
        else
            CacheProfileDirect(pageonly);
    }
    void CacheProfileMapped(bool) {};
    void CacheProfileDirect(bool pageonly);

    bool CompleteProfileInProfileBuffer(size_t endpos) {
        if( GetOpenedFileType()==Db::DbFlMAINBin)
            return CompleteProfileInProfileBufferBin(endpos);
        else
            return CompleteProfileInProfileBufferText(endpos);
    }
    bool CompleteProfileInProfileBufferBin(size_t);
    bool CompleteProfileInProfileBufferText(size_t endpos);

    void MoveDataFromPageToProfileBuffer();

    void ResetCurrentPageData() {
        current_strdata_.datlen_ = 0;
        current_strdata_.curpos_ = 0;
        current_strdata_.pagenr_ = 0;
        current_strdata_.pageoff_ = 0;
    }
    void ResetProfileBufferData() {
        profile_buffer_.datlen_ = 0;
        profile_buffer_.curpos_ = 0;
        profile_buffer_.pagenr_ = 0;
        profile_buffer_.pageoff_ = 0;
    }

    void MoveBlockOfData( 
        char* stream, size_t szstream, 
        size_t dstpos, size_t srcpos, size_t srclen );


    bool EofMapped() const {
        return
            profile_buffer_.datlen_ <= profile_buffer_.curpos_;
    }
    bool EofDirect() const {
        return
            db_size_bytes_ <= current_strdata_.pageoff_ &&
            profile_buffer_.datlen_ <= profile_buffer_.curpos_;
    }

    void GetDbPosMapped( TCharStream* ) const;
    void GetDbPosDirect( TCharStream* ) const;

    void SetDbPosMapped( const TCharStream* );
    void SetDbPosDirect( const TCharStream* );
    void SeekDbDirect( size_t );

    void SetDbPosToOriginMapped();
    void SetDbPosToOriginDirect();


    void SetDbSizeInBytes( size_t value ) {db_size_bytes_ = value;}
    size_t ObtainDbSizeInBytes(Db::TFile);

    size_t GetPageSize() const {return current_strdata_.pagesize_;}
    static ssize_t ObtainPageSize();

    TCharStream* GetProfileBuffer() {return &profile_buffer_;}
    void ResetProfileBufferPosition() {profile_buffer_.curpos_ = 0;}

    bool ValidFileDescriptor(int);
    void InvalidateFileDescriptor(int);

    void OpenFile(Db::TFile, const char*);
    void CloseFile(int);

    void MapFile(Db::TFile);
    void UnmapFile(int);

    bool ReadPage( TCharStream& );
    void ReadDataMapped( size_t filepos, void* dst, size_t count );
    void ReadDataDirect( size_t filepos, void* dst, size_t count );

private:
    const char* dbname_;//database name
    const bool mapped_;//whether a database is memory-mapped
    size_t no_sequences_;//number of profiles in the database
    uint64_mt db_size_;//database size in positions
    size_t db_size_bytes_;//database size in bytes

    char* name_buffer_;//database name buffer
    Db::TFile openedFile_;//type of the opened database main file
#ifdef OS_MS_WINDOWS
    HANDLE db_fp_[Db::DbFlN];//database file descriptors
#else
    int db_fp_[Db::DbFlN];//database file descriptors
#endif

    TCharStream current_strdata_;//current stream data
    TCharStream profile_buffer_;//stream for processing profile data
    char tmpbuff_[CUDBREADER_TMPBUFFSIZE];//buffer for temporary data
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// CompleteProfileInProfileBufferBin: check whether the complete profile is 
// present in the binary profile buffer;
// endpos, the last position up to which the check has been performed;
inline
bool CuDbReader::CompleteProfileInProfileBufferBin( size_t /*endpos*/)
{
    MYMSG("CuDbReader::CompleteProfileInProfileBufferBin",6);

    if( !profile_buffer_.data_)
        throw MYRUNTIME_ERROR(
        "CuDbReader::CompleteProfileInProfileBufferBin: Invalid profile buffer data.");

    int prosize;

    if( profile_buffer_.datlen_ < 1 || 
        profile_buffer_.datlen_ <= profile_buffer_.curpos_ + sizeof(prosize))
        return false;

    //NOTE profile size is present as the first data field in the profile
    prosize = *(int*)(profile_buffer_.data_+profile_buffer_.curpos_);
    if( prosize < 1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::CompleteProfileInProfileBufferBin: Invalid profile size.");
    return prosize <= (int)(profile_buffer_.datlen_-profile_buffer_.curpos_);
}

// CompleteProfileInProfileBufferText: check whether the complete profile is 
// present in the profile buffer;
// endpos, position down to which the profile end delimiter is searched;
inline
bool CuDbReader::CompleteProfileInProfileBufferText( size_t endpos )
{
    MYMSG("CuDbReader::CompleteProfileInProfileBufferText",6);

    if( !profile_buffer_.data_)
        throw MYRUNTIME_ERROR(
        "CuDbReader::CompleteProfileInProfileBufferText: Invalid profile buffer data.");

    if( profile_buffer_.datlen_ < 1 || profile_buffer_.datlen_ < lenstrEND)
        return false;

    for( size_t pos = profile_buffer_.datlen_-lenstrEND; 
        pos >= profile_buffer_.curpos_ && pos >= endpos; pos--)
    {
        size_t m;
        for( m = 0; m < lenstrEND && profile_buffer_.data_[pos+m] == patstrEND[m]; m++ );
        if( m == lenstrEND )
            return true;//found at pos;
        if( pos == 0 )
            break;
    }

    return false;
}

// -------------------------------------------------------------------------
// MoveDataFromPageToProfileBuffer: move data from the current page to the
// profile buffer and accordingly adjust the data within the structures
inline
void CuDbReader::MoveDataFromPageToProfileBuffer()
{
    MYMSG("CuDbReader::MoveDataFromPageToProfileBuffer",6);

    size_t leftinstr = (current_strdata_.curpos_ <= current_strdata_.datlen_)?
        current_strdata_.datlen_ - current_strdata_.curpos_: 0;

    if( leftinstr < 1 )
        return;

    if( !profile_buffer_.data_ ||
        profile_buffer_.pagesize_ < profile_buffer_.datlen_)
        throw MYRUNTIME_ERROR(
        "CuDbReader::MoveDataFromPageToProfileBuffer: Invalid profile buffer data.");

    if( profile_buffer_.pagesize_ < profile_buffer_.datlen_ + leftinstr)
    {
        if( profile_buffer_.datlen_ < profile_buffer_.curpos_)
            throw MYRUNTIME_ERROR(
            "CuDbReader::MoveDataFromPageToProfileBuffer: "
            "Invalid profile buffer's position.");

        MoveBlockOfData(
            profile_buffer_.data_, profile_buffer_.pagesize_, 
            0, profile_buffer_.curpos_, 
            profile_buffer_.datlen_ - profile_buffer_.curpos_ );
        profile_buffer_.datlen_ -= profile_buffer_.curpos_;
        profile_buffer_.curpos_ = 0;
        //
        if( profile_buffer_.pagesize_ < profile_buffer_.datlen_ + leftinstr)
            throw MYRUNTIME_ERROR(
            "CuDbReader::MoveDataFromPageToProfileBuffer: Too small profile buffer size.");
    }

    if( !current_strdata_.data_)
        throw MYRUNTIME_ERROR(
        "CuDbReader::MoveDataFromPageToProfileBuffer: Null page data.");

    memcpy(profile_buffer_.data_ + profile_buffer_.datlen_,
            current_strdata_.data_ + current_strdata_.curpos_, leftinstr);
    profile_buffer_.datlen_ += leftinstr;
    current_strdata_.datlen_ = 0;
    current_strdata_.curpos_ = 0;
}

// -------------------------------------------------------------------------
// MoveBlockOfData: move a block of data at position srcpos of the given 
// stream to the destination position;
// stream and szstream, stream and its size;
// dstpos, destination position in the stream;
// srcpos, destination position in the stream;
// srclen, size of a block to move;
inline
void CuDbReader::MoveBlockOfData( 
    char* stream, size_t szstream, 
    size_t dstpos, size_t srcpos, size_t srclen )
{
    MYMSG("CuDbReader::MoveBlockOfData",6);

    if( !stream || szstream < 1 )
        return;

    if( dstpos == srcpos || srclen < 1)
        return;

    if( szstream < srclen || szstream < dstpos + srclen )
        throw MYRUNTIME_ERROR("CuDbReader::MoveBlockOfData: Invalid parameters.");

    if( (dstpos < srcpos)? 
            dstpos + srclen < srcpos
        :   srcpos + srclen < dstpos )
    {
        memcpy(stream + dstpos, stream + srcpos, srclen);
        return;
    }

    for( size_t szmvd = 0; szmvd < srclen; ) {
        size_t szchnk = SLC_MIN((size_t)CUDBREADER_TMPBUFFSIZE, srclen - szmvd);
        size_t offset =
            (srcpos < dstpos )
            ?   srclen - szmvd - szchnk
            :   szmvd;
        memcpy(tmpbuff_, stream + srcpos + offset, szchnk);
        memcpy(stream + dstpos + offset, tmpbuff_, szchnk);
        szmvd += szchnk;
    }
}



// =========================================================================
// GetDbPos: get current position of the main database file
//
inline
void CuDbReader::GetDbPosMapped( TCharStream* chstr ) const
{
    MYMSG("CuDbReader::GetDbPosMapped",7);
#ifdef __DEBUG__
    if( !chstr )
        throw MYRUNTIME_ERROR("CuDbReader::GetDbPosMapped: Memory access error.");
#endif
    *chstr = profile_buffer_;
}
inline
void CuDbReader::GetDbPosDirect( TCharStream* chstr ) const
{
    MYMSG("CuDbReader::GetDbPosDirect",7);
#ifdef __DEBUG__
    if( !chstr )
        throw MYRUNTIME_ERROR("CuDbReader::GetDbPosDirect: Memory access error.");
#endif
    *chstr = current_strdata_;

    size_t pagesize = current_strdata_.pagesize_;
    ssize_t diff = profile_buffer_.datlen_ - profile_buffer_.curpos_;

    //if there left data unprocessed, go to the point where the next page will 
    // include the beginning of unprocessed data:
    for( ; (ssize_t)(chstr->pageoff_) > 0 && diff > 0; 
        chstr->pagenr_--, chstr->pageoff_ -= pagesize, diff -= pagesize);

    if((ssize_t)chstr->pageoff_ < 0) {
        chstr->pagenr_ = 0;
        chstr->pageoff_ = 0;
    }

    chstr->curpos_ = 0;
    if( diff < 0 )
        chstr->curpos_ = -diff;
}

// -------------------------------------------------------------------------
// SetDbPos: set current position for the main database file
inline
void CuDbReader::SetDbPosMapped( const TCharStream* chstr )
{
    MYMSG("CuDbReader::SetDbPosMapped",7);
#ifdef __DEBUG__
    if( !chstr )
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosMapped: Memory access error.");
#endif
    profile_buffer_ = *chstr;
}
inline
void CuDbReader::SetDbPosDirect( const TCharStream* chstr )
{
    MYMSG("CuDbReader::SetDbPosDirect",7);
#ifdef __DEBUG__
    if( !chstr )
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosDirect: Memory access error.");
#endif
    current_strdata_ = *chstr;

    profile_buffer_.datlen_ = 0;
    profile_buffer_.curpos_ = current_strdata_.curpos_;

    current_strdata_.datlen_ = 0;
    current_strdata_.curpos_ = 0;

    if( db_size_bytes_ < current_strdata_.pageoff_)
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosDirect: Invalid page offset.");

    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosDirect: Invalid file descriptor.");

#ifdef OS_MS_WINDOWS
    LARGE_INTEGER offset;
    offset.QuadPart = current_strdata_.pageoff_;
    offset.LowPart = SetFilePointer(
        db_fp_[GetOpenedFileType()],
        offset.LowPart,//lDistanceToMove (LONG), lower double word
        &offset.HighPart,//lpDistanceToMoveHigh (PLONG), upper double word
        FILE_BEGIN
    );
    if( offset.LowPart == INVALID_SET_FILE_POINTER && GetLastError() != NO_ERROR) {
        offset.QuadPart = -1;
        throw MYRUNTIME_ERROR(
        "CuDbReader::SetDbPosDirect: Setting file position failed.");
    }
#else
    if( lseek( db_fp_[GetOpenedFileType()], current_strdata_.pageoff_, SEEK_SET ) == 
        (off_t)-1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::SetDbPosDirect: Setting file position failed.");
#endif
}

// -------------------------------------------------------------------------
// SetDbPosToOrigin: set the current position of the main database 
// file to the origin
inline
void CuDbReader::SetDbPosToOriginMapped()
{
    MYMSG("CuDbReader::SetDbPosToOriginMapped",7);
    profile_buffer_.curpos_ = 0;
    profile_buffer_.pagenr_ = 0;
    profile_buffer_.pageoff_ = 0;
}
inline
void CuDbReader::SetDbPosToOriginDirect()
{
    MYMSG("CuDbReader::SetDbPosToOriginDirect",7);
    profile_buffer_.curpos_ = 0;
    profile_buffer_.pagenr_ = 0;
    profile_buffer_.pageoff_ = 0;

    current_strdata_.datlen_ = 0;
    current_strdata_.curpos_ = 0;
    current_strdata_.pagenr_ = 0;
    current_strdata_.pageoff_ = 0;
    //
    if( db_size_bytes_ < current_strdata_.pageoff_)
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosToOriginDirect: Invalid file size.");

    if( !ValidFileDescriptor(GetOpenedFileType()))
        throw MYRUNTIME_ERROR("CuDbReader::SetDbPosToOriginDirect: Invalid file descriptor.");

#ifdef OS_MS_WINDOWS
    LARGE_INTEGER offset;
    offset.QuadPart = current_strdata_.pageoff_;
    offset.LowPart = SetFilePointer(
        db_fp_[GetOpenedFileType()],
        offset.LowPart,//lDistanceToMove (LONG), lower double word
        &offset.HighPart,//lpDistanceToMoveHigh (PLONG), upper double word
        FILE_BEGIN
    );
    if( offset.LowPart == INVALID_SET_FILE_POINTER && GetLastError() != NO_ERROR) {
        offset.QuadPart = -1;
        throw MYRUNTIME_ERROR(
        "CuDbReader::SetDbPosToOriginDirect: Setting file position failed.");
    }
#else
    if( lseek( db_fp_[GetOpenedFileType()], current_strdata_.pageoff_, SEEK_SET ) == 
        (off_t)-1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::SetDbPosToOriginDirect: Setting file position failed.");
#endif
}

// -------------------------------------------------------------------------
// GetPageSize: get the system page size
//
inline
ssize_t CuDbReader::ObtainPageSize()
{
#ifdef OS_MS_WINDOWS
    SYSTEM_INFO systeminfo;
    GetSystemInfo(&systeminfo);
    return static_cast<ssize_t>(systeminfo.dwAllocationGranularity);
#else
    return sysconf(_SC_PAGESIZE);
#endif
}

// -------------------------------------------------------------------------
// ObtainDbSizeInBytes: get the database size in bytes
//
inline
size_t CuDbReader::ObtainDbSizeInBytes(Db::TFile which)
{
    if( !ValidFileDescriptor(which))
        throw MYRUNTIME_ERROR("CuDbReader::ObtainDbSizeInBytes: Invalid file descriptor.");
#ifdef OS_MS_WINDOWS
    LARGE_INTEGER filesize;
    if( GetFileSizeEx(db_fp_[which], &filesize) == 0 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::ObtainDbSizeInBytes: Failed to determine file size.");
    return static_cast<size_t>(filesize.QuadPart);
#else
    struct stat st;
    if( fstat(db_fp_[which], &st) < 0 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::ObtainDbSizeInBytes: Failed to determine file size.");
    return st.st_size;
#endif
}

// -------------------------------------------------------------------------
// OpenFile: open one of the database files for reading
inline
void CuDbReader::OpenFile(Db::TFile which, const char* filename )
{
    if( !filename )
        throw MYRUNTIME_ERROR(
        "CuDbReader::OpenFile: Unable to open the database.");

    if( ValidFileDescriptor(which))
        throw MYRUNTIME_ERROR(
        "CuDbReader::OpenFile: Database has been already opened.");

#ifdef OS_MS_WINDOWS
    db_fp_[which] = CreateFileA(
        filename,
        GENERIC_READ,
        FILE_SHARE_READ,
        NULL,//security descriptor: child processes cannot inherit the handle
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL);
#else
    db_fp_[which] = open( filename, O_RDONLY );
#endif

    if( !ValidFileDescriptor(which))
        throw MYRUNTIME_ERROR("CuDbReader::OpenFile: Failed to open the database.");
}

// CloseFile: close one of the database files
inline
void CuDbReader::CloseFile(int which)
{
    if( ValidFileDescriptor(which)) {
#ifdef OS_MS_WINDOWS
        CloseHandle(db_fp_[which]);
#else
        close(db_fp_[which]);
#endif
        InvalidateFileDescriptor(which);
    }
}

// -------------------------------------------------------------------------
// MapFile: memory-map the main database file for reading
inline
void CuDbReader::MapFile(Db::TFile which )
{
    if( which != Db::DbFlMAIN && which != Db::DbFlMAINBin )
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Incorrectly trying to memory-map a database file.");

    if( !ValidFileDescriptor(which))
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Invalid file descriptor.");

    size_t mapsize = GetDbSizeInBytes();

    if( mapsize < 1 )
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Invalid database file size.");

#ifdef OS_MS_WINDOWS
    profile_buffer_.hMapFile_ = CreateFileMapping(
        db_fp_[which],
        NULL,//security descriptor: cannot be inherited
        PAGE_READONLY,
        //zeros indicate that the size of file is used for...
        // the maximum size of the file mapping object:
        ((uint64_t)mapsize)>>32/*0*/,//dwMaximumSizeHigh: upper double word
        ((uint64_t)mapsize)&0xffffffff/*0*/,//dwMaximumSizeLow: lower double word
        NULL);//no name

    if( profile_buffer_.hMapFile_ == NULL)
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Failed to create a file mapping object.");

    profile_buffer_.data_ = static_cast<char*>(
        MapViewOfFile(
            profile_buffer_.hMapFile_,
            FILE_MAP_READ,
            0,//dwFileOffsetHigh
            0,//dwFileOffsetLow
            0)//dwNumberOfBytesToMap: 
        //mapping extends from the specified offset to the end of the file mapping
    );

    if(profile_buffer_.data_ == NULL)
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Failed to map a view of the file mapping object.");
#else
    profile_buffer_.data_ = (char*)
        mmap(NULL, mapsize, PROT_READ, MAP_SHARED, db_fp_[which], 0);

    if(profile_buffer_.data_ == MAP_FAILED)
        throw MYRUNTIME_ERROR(
        "CuDbReader::MapFile: Failed to memory-map an opened database file.");
#endif

    profile_buffer_.datlen_ = mapsize;
    profile_buffer_.curpos_ = 0;
    profile_buffer_.pagenr_ = 0;
    profile_buffer_.pagesize_ = mapsize;
    profile_buffer_.pageoff_ = 0;
}

// UnmapFile: unmap the main database file
inline
void CuDbReader::UnmapFile(int which)
{
    if( which != Db::DbFlMAIN && which != Db::DbFlMAINBin )
        return;

#ifdef OS_MS_WINDOWS
    if( profile_buffer_.data_ )
        UnmapViewOfFile(profile_buffer_.data_);

    if( profile_buffer_.hMapFile_ )
        CloseHandle(profile_buffer_.hMapFile_);

    profile_buffer_.hMapFile_ = NULL;
#else
    if( profile_buffer_.data_ )
        munmap(profile_buffer_.data_, profile_buffer_.pagesize_);
#endif

    profile_buffer_.data_ = NULL;
    profile_buffer_.datlen_ = 0;
    profile_buffer_.curpos_ = 0;
    profile_buffer_.pagenr_ = 0;
    profile_buffer_.pagesize_ = 0;
    profile_buffer_.pageoff_ = 0;
}

// -------------------------------------------------------------------------
// ValidFileDescriptor: get a flag of whether the file descriptor is valid
inline
bool CuDbReader::ValidFileDescriptor(int which)
{
    return (0 <= which && which < Db::DbFlN)?
#ifdef OS_MS_WINDOWS
            db_fp_[which] != INVALID_HANDLE_VALUE
#else
            db_fp_[which] >= 0
#endif
        :   false;
}
// InvalidateFileDescriptor: invalidate the given file descriptor
inline
void CuDbReader::InvalidateFileDescriptor(int which)
{
    if(0 <= which && which < Db::DbFlN)
#ifdef OS_MS_WINDOWS
        db_fp_[which] = INVALID_HANDLE_VALUE;
#else
        db_fp_[which] = -1;
#endif
}

#endif//__CuDbReader_h__
