/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "mybase.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ConfigFile.h"

const char* ConfigFile::s_equal   = "=";
const char* ConfigFile::s_bracket = "[]";
const char* ConfigFile::s_newline = NEWLINE;

//------------------------------------------------------------------------------
// Constructor
//
ConfigFile::ConfigFile( const char* section, const char* filename )
:   section_name( section ),
    filename( filename ),
    errbuffer( NULL ),
    errbuflen( 0 )
{
    errbuflen = KBYTE;
    errbuffer = ( char* )malloc( errbuflen );
    if( errbuffer == NULL )
        throw MYRUNTIME_ERROR("ConfigFile::ConfigFile: Not enough memory.");
}

//------------------------------------------------------------------------------
// Destructor
//
ConfigFile::~ConfigFile()
{
    if( errbuffer )
        free( errbuffer );
}

//------------------------------------------------------------------------------
// GetInt: get integer value of key;
//     key,  key name;
//     value, key value to return;
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetInt( const char* key, int* value )
{
    return GetCFPrivateProfileInt( key, value );
}

//------------------------------------------------------------------------------
// GetFloat: get float value of key;
//     key, key name;
//     value, key value to return;
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetFloat( const char* key, float* value )
{
    return GetCFPrivateProfileFloat( key, value );
}

//------------------------------------------------------------------------------
// GetDouble: get double value of key;
//     key, key name;
//     value, key value to return;
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetDouble( const char* key, double* value )
{
    return GetCFPrivateProfileDouble( key, value );
}

//------------------------------------------------------------------------------
// GetString: get string value of key
//     key, key name
//     value, key value to return
//     size, size of value
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetString( const char* key, char* value, unsigned size )
{
    return GetCFPrivateProfileString( key, value, size );
}

//------------------------------------------------------------------------------
// WriteInt: write integer value of key to file
//     key, key name
//     value, key value
// return true if the key has been found; if the key or section hasn't been
//     found, they are created and inserted to the file
//
bool ConfigFile::WriteInt( const char* key, int value )
{
    char        strvalue[BUF_MAX];

    sprintf( strvalue, "%d", value );
    return WriteString( key, strvalue );
}

//------------------------------------------------------------------------------
// WriteFloat: write float value of key to file;
//     key, key name;
//     value, key value;
// return true if the key has been found; if the key or section hasn't been
//     found, they are created and inserted to the file
//
bool ConfigFile::WriteFloat( const char* key, float value )
{
    char        strvalue[BUF_MAX];

    sprintf( strvalue, "%.6f", value );
    return WriteString( key, strvalue );
}

//------------------------------------------------------------------------------
// WriteDouble: write double value of key to file
//     key, key name
//     value, key value
// return true if the key has been found; if the key or section hasn't been
//     found, they are created and inserted to the file
//
bool ConfigFile::WriteDouble( const char* key, double value )
{
    char        strvalue[BUF_MAX];

    sprintf( strvalue, "%.6f", value );
    return WriteString( key, strvalue );
}

//------------------------------------------------------------------------------
// WriteString: write string value of key to file
//     key, key name
//     value, key value
// return true if the key has been found; otherwise it'll be inserted
//
bool ConfigFile::WriteString( const char* key, const char* value )
{
    return WriteCFPrivateProfileString( key, value );
}

//------------------------------------------------------------------------------
// GetCFPrivateProfileInt: open the file to search for key and its value; 
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetCFPrivateProfileInt( const char* key, int* value )
{
    char    strvalue[BUF_MAX];
    char*   p;

    if( !value )
        throw MYRUNTIME_ERROR("ConfigFile::GetCFPrivateProfileInt: Invalid argument.");

    if( !GetCFPrivateProfileString( key, strvalue, BUF_MAX ))
        return false;

    *value = strtol( strvalue, &p, 10 );
    if( errno || *p ) {
        Format_Message("ConfigFile::GetCFPrivateProfileInt: Invalid integer value of key ", key );
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    return true;
}

//------------------------------------------------------------------------------
// GetCFPrivateProfileFloat: open the file for searching for key and its float 
// value; return true if the key has been found and its value has been read,
//      false otherwise
//
bool ConfigFile::GetCFPrivateProfileFloat( const char* key, float* value )
{
    char    strvalue[BUF_MAX];
    char*   p;

    if( !value )
        throw MYRUNTIME_ERROR("ConfigFile::GetCFPrivateProfileFloat: Invalid argument.");

    if( !GetCFPrivateProfileString( key, strvalue, BUF_MAX ))
        return false;

    *value = strtof( strvalue, &p );

    if( errno || *p ) {
        Format_Message("ConfigFile::GetCFPrivateProfileFloat: Invalid float value of key ", key );
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    return true;
}

//------------------------------------------------------------------------------
// GetCFPrivateProfileDouble: open the file for searching for key and its double 
// value; return true if the key has been found and its value has been read,
//      false otherwise
//
bool ConfigFile::GetCFPrivateProfileDouble( const char* key, double* value )
{
    char    strvalue[BUF_MAX];
    char*   p;

    if( !value )
        throw MYRUNTIME_ERROR("ConfigFile::GetCFPrivateProfileDouble: Invalid argument.");

    if( !GetCFPrivateProfileString( key, strvalue, BUF_MAX ))
        return false;

    *value = strtod( strvalue, &p );

    if( errno || *p ) {
        Format_Message("ConfigFile::GetCFPrivateProfileDouble: Invalid double value of key ", key );
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    return true;
}

//------------------------------------------------------------------------------
// GetCFPrivateProfileString: open the file for searching for key and its 
// value of the string type;
//     key, key name;
//     value, key value to return;
//     size, max bytes <value> can accomodate;
// return true if the key has been found and its value has been read,
//     false otherwise
//
bool ConfigFile::GetCFPrivateProfileString( const char* key, char* value, unsigned size )
{
    if( !GetSectionName() || !GetFileName())
        throw MYRUNTIME_ERROR("ConfigFile::GetCFPrivateProfileString: "
                              "Null filename and/or section name.");

    FILE*       fd = NULL;
    bool        found = false;

    fd = fopen( GetFileName(), "r" );

    if( fd == NULL ) {
        Format_Message("ConfigFile::GetCFPrivateProfileString: Failed to open file, ", GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    try {
        if( FindFileSection( fd ))
            if( GetSectionKey( fd, key, value, size ))
                found =  true;

    } catch( myexception const& ex )
    {
        fclose( fd );
        throw MYRUNTIME_ERROR2( ex.what(), ex.eclass());
    }

    fclose( fd );
    return found;
}

//------------------------------------------------------------------------------
// WriteCFPrivateProfileString: open the file for writing key value;
//     key, key name;
//     value, symbolic key value;
// return true if the key has been found; if the key or section hasn't been
//     found, they are created and inserted to the file
//
bool ConfigFile::WriteCFPrivateProfileString( const char* key, const char* value )
{
    if( !GetSectionName() || !GetFileName())
        throw MYRUNTIME_ERROR("ConfigFile::WriteCFPrivateProfileString: "
                              "Null filename and/or section name.");

    FILE*       fd = NULL;
    bool        found = false;

    fd = fopen( GetFileName(), "r+w" );

    if( fd == NULL ) {
        Format_Message("ConfigFile::WriteCFPrivateProfileString: Failed to open file, ", GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    try {
        found = FindFileSection( fd, true );
        if( WriteSectionKey( fd, key, value, found ))
            found =  true;

    } catch( myexception const& ex )
    {
        fclose( fd );
        throw MYRUNTIME_ERROR2( ex.what(), ex.eclass());
    }

    fclose( fd );
    return found;
}

//------------------------------------------------------------------------------
// FindFileSection: find section in the file by name; if not found, insert it;
//     fd, pointer to file;
//     bwrite, true if section is to be appended to the file;
// return true if section has been found, false otherwise
//
bool ConfigFile::FindFileSection( FILE* fd, bool bwrite )
{
    if( !GetSectionName())
        throw MYRUNTIME_ERROR("ConfigFile::FindFileSection: Null section name.");

    if( fd == NULL )
        return false;

    const int   locsize = BUF_MAX;
    char        locbuffer[locsize] = {0};
    char*       p = locbuffer;
    const size_t seclen = strlen( GetSectionName());
    const size_t nllen = strlen( s_newline );
    size_t      len = nllen;
    size_t      nlns = 0;//number of newlines to the end of file
    bool        found = false;

    while( !feof( fd )) {
        do {
            p = fgets( locbuffer, locsize, fd );

            len = strlen( locbuffer );

            if( p == NULL )
                break;

            //if( len && strcmp( locbuffer + len - nllen, s_newline ) == 0 )
            //    break;
            if( len ) {
                size_t i = 1;
                for(; i <= nllen && i <= len && locbuffer[len-i]!='\n' && locbuffer[len-i]!='\r'; i++ )
                    break;
                if( i <= nllen && i <= len )
                    break;
            }
        } while( 1 );

        if( p == NULL )
            break;

        //omit spaces before checking for a newline
        for( p = locbuffer; *p == ' ' || *p == '\t'; p++ );

        //if( nllen <= len && strcmp( p, s_newline ) == 0 )
        //    nlns += nllen + ( p - locbuffer ) - ( p != locbuffer );
        //else
        //    nlns = 0;
        if( len ) {
            nlns += (p - locbuffer) - (p != locbuffer);
            for(; *p=='\n' || *p=='\r'; p++, nlns++);
        }
        else
            nlns = 0;

        if( *( p = locbuffer ) != s_bracket[0] )
            continue;
        p++;

        if( strncmp( p, GetSectionName(), seclen ) == 0 && p[seclen] == s_bracket[1] ) {
            found = true;
            break;
        }
    }

    if( !bwrite || found )
        return found;

    SeekToPosition( fd, -(long)nlns, SEEK_END );

    if( fprintf( fd, "%s%c%s%c%s", 
          s_newline, s_bracket[0], GetSectionName(), s_bracket[1], s_newline ) < 0 ) {
        Format_Message( "ConfigFile::FindFileSection: Unable to write to file ", GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    return found;
}

//------------------------------------------------------------------------------
// GetSectionKey: search for key in the file by name;
//     fd, file pointer;
//     key, key name;
//     value, key value;
//     size, size of value;
//     pos, the address to contain file position at which the key has been 
//          found or the end of the section;
// return true if the key has been found, false otherwise
//
bool ConfigFile::GetSectionKey( 
    FILE* fd, const char* key, char* value, size_t size, long* pos, long* lastpos )
{
    if( fd == NULL || key == NULL )
        return false;

    const int   locsize = BUF_MAX;
    char        locbuffer[locsize] = {0};
    char*       p = locbuffer;
    char*       lastline = locbuffer;
    size_t      keylen = 0;
    size_t      equlen = strlen(s_equal);
    size_t      nllen = strlen(s_newline);
    size_t      len = 1;
    size_t      plen = 0;
    bool        found = false;
    bool        nl_assigned = false;


    keylen = (int)strlen(key);

    while( !feof( fd )) {

        if( pos && !nl_assigned )
            *pos = TellPosition( fd );

        lastline = NULL;

        do {
            p = fgets( locbuffer, locsize, fd );

            len = strlen( locbuffer );

            if( p == NULL )
                //if the end of file
                break;

            if( !lastline && len && *locbuffer == s_bracket[0] )
                //if new section begins
                return found;

            if(!lastline ) {
                //omit spaces
                for( lastline = locbuffer; *lastline == ' ' || *lastline == '\t'; lastline++ );
            }

            if( len ) {
                size_t i = 1;
                for(; i <= nllen && i <= len && locbuffer[len-i]!='\n' && locbuffer[len-i]!='\r'; i++)
                    break;
                if( i <= nllen && i <= len )
                    break;
            }
        //} while( len && strcmp( locbuffer + len - nllen, s_newline ));
        } while(1);

        //if( lastline && strcmp( lastline, s_newline ) == 0 )
        if( lastline && (*lastline=='\n' || *lastline=='\r'))
                    nl_assigned = true;
            else    nl_assigned = false;

        if( lastpos && !nl_assigned )
            *lastpos = TellPosition( fd );

        if( p == NULL )
            break;

        //allow spaces in front of key
        for( p = locbuffer; *p == ' ' || *p == '\t'; p++ );

        if( strncmp( key, p, keylen ))
            continue;

        //omit spaces before the equality sign
        for( p += keylen; *p == ' ' || *p == '\t'; p++ );
        if( strncmp( p, s_equal, equlen ))
            continue;

        //omit spaces after equality
        for( p += equlen; *p == ' ' || *p == '\t'; p++ );

        plen = strlen( p );
        if( size < plen )
            plen = size;

        if( value && plen ) {
            strncpy( value, p, plen - 1 );
            value[plen-1] = 0;
            ProcessValue( value );

            found = true;
            break;
        }
    }

    return found;
}

//------------------------------------------------------------------------------
// WriteSectionKey: try to write key in the file; the key is searched at first;
//     fd, pointer to file;
//     key, key name;
//     value, key value;
// return true if the key has been found, false otherwise
//
bool ConfigFile::WriteSectionKey( FILE* fd, const char* key, const char* value, bool section_found )
{
    if( fd == NULL || key == NULL || value == NULL )
        return false;

    char    strvalue[BUF_MAX];
    long    file_pos = 0;
    long    last_pos = 0;
    bool    found = false;

    found = GetSectionKey( fd, key, strvalue, BUF_MAX, &file_pos, &last_pos );
    //the position after which the key is to be inserted is found
    InsertSectionKeyAfter( fd, key, value, file_pos, last_pos, section_found );

    return found;
}

//------------------------------------------------------------------------------
// InsertSectionKeyAfter: try to insert key with its value at the given 
//     position of the file;
//     fd, file pointer;
//     key, key name;
//     value, key value;
//     file_pos, file position to insert the key at;
//
void ConfigFile::InsertSectionKeyAfter(
        FILE* fd,
        const char* key,    const char* value,
        long file_pos,      long last_pos,
        bool section_found )
{
    if( fd == NULL || key == NULL || value == NULL )
        return;

    SeekToPosition( fd, last_pos, SEEK_SET );

    const int   size = KBYTE;
    char*       inner = ( char* )malloc( size );
    char*       p = inner;
    int         cursize = size;
    int         prvsize = 0;
    int         amount = 0;
    size_t      lastread = 0;

    if( inner == NULL )
        throw MYRUNTIME_ERROR("ConfigFile::InsertSectionKeyAfter: Not enough memory.");

    *inner = 0;

    while( !feof( fd ) && section_found ) {
        lastread = fread( p, 1, size, fd );
        amount += (int)lastread;

        if( ferror( fd )) {
            free( inner );
            Format_Message("ConfigFile::InsertSectionKeyAfter: Read error in file ", 
                            GetFileName());
            throw MYRUNTIME_ERROR( GetErrBuffer());
        }

        if( feof( fd ))
            break;

        inner = ( char* )realloc( inner, cursize += size );

        if( inner == NULL )
            throw MYRUNTIME_ERROR("ConfigFile::InsertSectionKeyAfter: Not enough memory.");

        //each time the inner address changes!!
        p = inner + amount;
    }

    SeekToPosition( fd, file_pos, SEEK_SET );

    if( fprintf( fd, "%s %s %s%s", key, s_equal, value, s_newline ) < 0 ) {
        free( inner );
        Format_Message("ConfigFile::InsertSectionKeyAfter: Unable to insert key in file ", 
                        GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    for( p = inner, prvsize = size, cursize = 0; cursize < amount; cursize += size ) {
        if( amount < cursize + size )
            prvsize = amount - cursize;

        fwrite( p, 1, prvsize, fd );

        if( ferror( fd )) {
            free( inner );
            Format_Message("ConfigFile::InsertSectionKeyAfter: Write error in file ", 
                           GetFileName());
            throw MYRUNTIME_ERROR( GetErrBuffer());
        }

        p += prvsize;
    }

    free( inner );

    if( section_found )
        FlushToEOF( fd );
}

//------------------------------------------------------------------------------
// FlushToEOF: flush file from the current position to the end of file
//
void ConfigFile::FlushToEOF( FILE* fd )
{
    if( !fd )
        return;

    char*   temp = NULL;
    long    curpos = TellPosition( fd );    SeekToPosition( fd, 0, SEEK_END );
    long    endpos = TellPosition( fd );
    size_t  nllen = strlen( s_newline );

    long    diff = endpos - curpos;

    if( 0 < diff && diff < (long)nllen )
        diff = (long)nllen;

    if( 0 < diff ) {
        SeekToPosition( fd, curpos, SEEK_SET );
        temp = ( char* )malloc( diff + 1 ); //one for null character
        if( !temp )
            throw MYRUNTIME_ERROR("ConfigFile::FlushToEOF: Not enough memory.");
        memset( temp, ' ', diff );
        strcpy( temp + diff - nllen, s_newline );
        fwrite( temp, 1, diff, fd );
        free( temp );
    }
}

//------------------------------------------------------------------------------
// TellPosition: tell the current position of the file
//
long ConfigFile::TellPosition( FILE* fd )
{
    long    pos = 0;

    if( !fd )
        return pos;

    pos = ftell( fd );

    if( pos < 0 ) {
        Format_Message("ConfigFile::TellPosition: Unable to determine position in file ", 
                       GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }

    return pos;
}

//------------------------------------------------------------------------------
// SeekToPosition: seek to position in the file
//
void ConfigFile::SeekToPosition( FILE* fd, long pos, int whence )
{
    if( fseek( fd, pos, whence ) < 0 ) {
        Format_Message("ConfigFile::SeekToPosition: Unable to seek in file ", 
                       GetFileName());
        throw MYRUNTIME_ERROR( GetErrBuffer());
    }
}

//------------------------------------------------------------------------------
// ProcessValue: process key value by removing leading and trailing spaces;
// comments at the end of the value are flushed;
// NOTE: argument must be a null-terminated string
//
void ConfigFile::ProcessValue( char* value )
{
    if( !value )
        return;

    size_t  len = strlen(value);
    bool    text = false;

    for( int n = (int)len - 1; 0 <= n; n-- ) {
        if( value[n] == '#' ) {
            value[n] = 0;
            text = false;
            continue;
        }
        if( !text && ( value[n] == ' ' || value[n] == '\t' )) {
            value[n] = 0;
            continue;
        }
        text = true;
    }
}

//------------------------------------------------------------------------------
// Format_Message: append the filename to the end of message
//
void ConfigFile::Format_Message( const char* msg, const char* cval )
{
    int     length = GetErrBufLength();
    char*   p = GetErrBuffer();

    if( length <= 0 )
        throw MYRUNTIME_ERROR("ConfigFile::Format_Message: Unallocated message buffer.");

    *p = 0;

    if( !msg || !cval )
        return;

    int prelen = (int)strlen(msg);
    int fnmlen = (int)strlen(cval);

    if( length <= prelen )
        prelen  = length - 1;

    memcpy( p, msg, prelen );
    length -= prelen;
    p += prelen;

    if( length < 2 ) {
        *p = 0;
        return;
    }

    if( fnmlen + 1 < length )
        length = fnmlen + 1;

    memcpy( p, cval, length - 1 );
    p[length - 1] = 0;
}

