/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "myexception.h"
#include "mystring.h"


const size_t mystring::npos = static_cast<size_t>( -1 );

// -------------------------------------------------------------------------
// CLASS mystring
//
// Default constructor
//
mystring::mystring()
:   _string( NULL ),
    _length( 0 ),
    _capacity( 0 )
{
}

// Copy constructor
//
mystring::mystring( const mystring& str )
:   _string( NULL ),
    _length( 0 ),
    _capacity( 0 )
{
    assign( str.c_str(), str.length());
}

// Construction from C string
//     s, character C string
//     noc, number of characters to use in string construction
//
mystring::mystring( const char* s, size_t noc )
:   _string( NULL ),
    _length( 0 ),
    _capacity( 0 )
{
    if( noc == npos )
        assign( s, s? strlen( s ): 0 );
    else
        assign( s, noc );
}

// Construction from a character
//
mystring::mystring( char c )
:   _string( NULL ),
    _length( 0 ),
    _capacity( 0 )
{
    assign( &c, 1 );
}

// Destructor
//
mystring::~mystring()
{
    destroy();
}

// -------------------------------------------------------------------------
// clear: clears object's data but the space allocated remains
//
void mystring::clear()
{
    if( _string )
        memset( _string, 0, _capacity );
    _length = 0;
}

// -------------------------------------------------------------------------
// destroy: destroys all object's data
//
void mystring::destroy()
{
    if( _string )
        free( _string );
    _string = NULL;
    _length = 0;
    _capacity = 0;
}

// -------------------------------------------------------------------------
// reallocate: reallocates space for string so that new capacity is greater
//     than was before
//
void mystring::reallocate( size_t newcap )
{
    if( newcap <= capacity())
        //nothing doing, string has enough space
        return;

    if( _string == NULL )
        _string = ( char* )malloc( newcap );
    else
        _string = ( char* )realloc( _string, newcap );

    if( _string == NULL ) {
        _capacity = 0;
        _length = 0;
        throw myruntime_error( __EXCPOINT__ );
        return;
    }

    memset( _string + _capacity, 0, newcap - _capacity );
    _capacity = newcap;
}

// -------------------------------------------------------------------------
// Case conversion
//
void mystring::upper()
{
    for( size_t n = 0; n < _length; n++ )
        if( 0 < _string[n])
            _string[n] = ( char )toupper( _string[n]);
}

void mystring::lower()
{
    for( size_t n = 0; n < _length; n++ )
        if( 0 < _string[n])
            _string[n] = ( char )tolower( _string[n]);
}

// -------------------------------------------------------------------------
// Assignment methods
//
mystring& mystring::assign( const mystring& str )
{
    assign( str.c_str(), str.length());
    return *this;
}

// assign C string
//
mystring& mystring::assign( const char* s )
{
    assign( s, s? strlen( s ): 0 );
    return *this;
}

// assign a character
//
mystring& mystring::assign( char c )
{
    assign( &c, 1 );
    return *this;
}

// helper assign method
//
void mystring::assign( const char* buf, size_t len )
{
    if( _string == buf )
        return;

    if( buf == NULL ) {
        clear();
        return;
    }

    reallocate( len + 1 );

    if( len )
        memcpy( _string, buf, len );

    _length = len;
    _string[_length] = 0;
}

// -------------------------------------------------------------------------
// Append methods
//
// append substring
//
mystring& mystring::append( const mystring& str, size_t pos, size_t len )
{
    if( pos < str.length()) {
        if( str.length() + 1 <= pos + len )
            len = str.length() - pos;
        append( str.c_str() + pos, len );
    }
    return *this;
}

// append C substring
//
mystring& mystring::append( const char* s, size_t pos, size_t len )
{
    size_t  slength = 0;

    if( s )
        slength = strlen( s );

    if( pos < slength ) {
        if( slength + 1 <= pos + len )
            len = slength - pos;
        append( s + pos, len );
    }
    return *this;
}

// append mystring
//
mystring& mystring::append( const mystring& str )
{
    append( str.c_str(), str.length());
    return *this;
}

// append C string
//
mystring& mystring::append( const char* s )
{
    append( s, s? strlen( s ): 0 );
    return *this;
}

// append a character
//
mystring& mystring::append( char c )
{
    append( &c, 1 );
    return *this;
}

// helper append method
//
void mystring::append( const char* buf, size_t len )
{
    if( _string == buf )
        return;

    if( buf == NULL || len == 0 ) {
        return;
    }

//    reallocate( _length + _length + len + 1 );
    reallocate( _length + len + len + 1 );

    memcpy( _string + _length, buf, len );

    _length += len;
    _string[_length] = 0;
}

// -------------------------------------------------------------------------
// Insert methods
//
mystring& mystring::insert( size_t pos, const mystring& str )
{
    insert( pos, str.c_str(), str.length());
    return *this;
}

// insert C string
//
mystring& mystring::insert( size_t pos, const char* s )
{
    insert( pos, s, s? strlen( s ): 0 );
    return *this;
}

// insert a character
//
mystring& mystring::insert( size_t pos, char c )
{
    insert( pos, &c, 1 );
    return *this;
}

// helper insert method
//     pos, position to insert character string at
//     buf, character string
//     len, length of character string
//
void mystring::insert( size_t pos, const char* buf, size_t len )
{
    if( _string == buf )
        return;

    if( buf == NULL || len == 0 || length() <= pos ) {
        return;
    }

    reallocate( _length + len + 1 );

    for( size_t n = _length + len - 1; n >= pos + len; n-- )
        _string[n] = _string[n-len];

    for( size_t m = pos; m < pos + len; m++ )
        _string[m] = *buf++;

    _length += len;
    _string[_length] = 0;
}

// -------------------------------------------------------------------------
// helper character string find method
//
size_t mystring::find( const char* buf, size_t len ) const
{
    if( pstring() == buf )
        return 0;

    if( buf == NULL || len == 0 )
        return npos;

    if( pstring() == NULL || length() < len )
        return npos;

    const char* p = pstring();
    const char* s = buf;
    size_t      m = 0;

    for( size_t n = 0; n + len <= length(); n++ ) {
        for( m = 0; m < len && p[n+m] == s[m]; m++ );
        if( m == len )
            return n;
    }

    return npos;
}

// -------------------------------------------------------------------------
// helper method for reverse find of character string
//
size_t mystring::rfind( const char* buf, size_t len ) const
{
    if( pstring() == buf )
        return 0;

    if( buf == NULL || len == 0 )
        return npos;

    if( pstring() == NULL || length() < len )
        return npos;

    const char* p = pstring();
    const char* s = buf;
    size_t      m = 0;

    for( size_t n = length() - len; ; n-- ) {
        for( m = 0; m < len && p[n+m] == s[m]; m++ );
        if( m == len )
            return n;
        if( n == 0 )
            break;
    }

    return npos;
}

// -------------------------------------------------------------------------
// substr: returns string containing substring of this string
//     pos, position of substring
//     noc, number of characters in the substring
//
mystring mystring::substr( size_t pos, size_t noc ) const
{
    if( pstring() == NULL || length() <= pos || noc == 0 )
        return mystring();

    if( noc != npos )
        if( length() < pos + noc )
            noc = length() - pos;

    mystring    str( pstring() + pos, noc );
    return str;
}

// -------------------------------------------------------------------------
// erase: modifies string by erasing a segment from the position specified
//     pos, position to begin erase from
//     noc, number of characters to erase
//
mystring& mystring::erase( size_t pos, size_t noc )
{
    if( pstring() == NULL || length() <= pos || noc == 0 )
        return *this;

    if( noc != npos )
        if( length() < pos + noc )
            noc = npos;

    if( noc == npos ) {
        _length = pos;
        _string[_length] = 0;
        return *this;
    }

    for( size_t n = pos; n < _length - noc; n++ )
        _string[n] = _string[n+noc];

    _length -= noc;
    _string[_length] = 0;
    return *this;
}

// -------------------------------------------------------------------------
// ncompare: compares just fraction of strings determined by their minimum
//     length
//
int mystring::ncompare( const mystring& str ) const
{
    return ncompare( str.c_str());
}
// ncompare: ncompare to C string
//
int mystring::ncompare( const char* s ) const
{
    if( s == NULL )
        return (int)length();

    const size_t    tsize = length();
    const size_t    osize = strlen( s );
    size_t          minlen = tsize;

    if( osize < tsize )
        minlen = osize;

    int result = memcmp( c_str(), s, minlen );

    return result;
}

// compare: compares this string to aonther
//
int mystring::compare( const mystring& str ) const
{
    return compare( str.c_str());
}

// compare: compares this string to C string
//
int mystring::compare( const char* s ) const
{
    if( s == NULL )
        return (int)length();

    const size_t    tsize = length();
    const size_t    osize = strlen( s );
    size_t          minlen = tsize;

    if( osize < tsize )
        minlen = osize;

    int result = memcmp( c_str(), s, minlen );

    if( result == 0 )
        result = (int)(tsize - osize);

    return result;
}
