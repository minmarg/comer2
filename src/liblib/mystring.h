/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __mystring_h__
#define __mystring_h__

#include "debug.h"
#include <string.h>


// _________________________________________________________________________
// CLASS mystring
// text manipulation 
//
class mystring
{
public:
    mystring();
    mystring( const mystring& );
    mystring( const char*, size_t = npos );
    mystring( char );
    ~mystring();

    void            reserve( size_t size )  { reallocate( size ); }

    void move(mystring& str) {
        destroy();
        _string = str._string;
        _length = str._length;
        _capacity = str._capacity;
        str._string = NULL;
        str._length = 0;
        str._capacity = 0;
    }

    mystring&       operator=( const mystring& str )    { return assign( str ); }
    mystring&       operator=( const char* s )          { return assign( s ); }
    mystring&       operator=( char c )                 { return assign( c ); }

    char            operator[]( size_t pos ) const;
    char&           operator[]( size_t pos );

    const char*     c_str() const       { return _string? _string: ""; }
    size_t          length() const      { return _length; }
    bool            empty() const       { return length() == 0; }

    void            upper();
    void            lower();

    mystring&       assign( const mystring& str );
    mystring&       assign( const char* s );
    mystring&       assign( char c );
    void            assign( const char* buf, size_t len );

    mystring&       append( const mystring& str, size_t pos, size_t len );
    mystring&       append( const char* s, size_t pos, size_t len );

    mystring&       append( const mystring& str );
    mystring&       append( const char* s );
    mystring&       append( char c );

    mystring&       insert( size_t pos, const mystring& str );
    mystring&       insert( size_t pos, const char* s );
    mystring&       insert( size_t pos, char c );

    mystring&       operator+=( const mystring& str )    { return append( str ); }
    mystring&       operator+=( const char* s )          { return append( s ); }
    mystring&       operator+=( char c )                 { return append( c ); }

    size_t          find( const mystring& str ) const   { return find( str.c_str(), str.length()); }
    size_t          find( const char* s ) const         { return find( s, s? strlen( s ): 0 ); }
    size_t          find( char c ) const                { return find( &c, 1 ); }

    size_t          rfind( const mystring& str ) const  { return rfind( str.c_str(), str.length()); }
    size_t          rfind( const char* s ) const        { return rfind( s, s? strlen( s ): 0 ); }
    size_t          rfind( char c ) const               { return rfind( &c, 1 ); }

    mystring        substr( size_t pos, size_t noc = npos ) const;
    mystring&       erase( size_t pos = 0, size_t noc = npos );

    int             ncompare( const mystring& str ) const;
    int             ncompare( const char* s ) const;

    int             compare( const mystring& str ) const;
    int             compare( const char* s ) const;

    static const size_t npos;       //value returned by the methods when they fail

protected:
    void            clear();
    void            destroy();
    void            reallocate( size_t );

    void            append( const char* buf, size_t len );
    void            insert( size_t pos, const char* buf, size_t len );

    size_t          find(  const char* buf, size_t len ) const;
    size_t          rfind( const char* buf, size_t len ) const;

    const char*     pstring() const     { return _string; }
    size_t          capacity() const    { return _capacity; }

private:
    char*   _string;
    size_t  _length;
    size_t  _capacity;
};



// INLINES
// -------------------------------------------------------------------------
// operator[]
//
inline
char mystring::operator[]( size_t pos ) const
{
    if( _length <= pos )
        return 0;
    return _string[pos];
}

static char _mystring_unused_char = 0;

inline
char& mystring::operator[]( size_t pos )
{
    if( _length <= pos )
        return _mystring_unused_char;
    return _string[pos];
}

// -------------------------------------------------------------------------
// operator +
//
inline
mystring operator+( const mystring& left, const mystring& right )
{
    mystring    str( left );
    str.append( right );
    return str;
}

inline
mystring operator+( const char* s, const mystring& right )
{
    mystring    str( s );
    str.append( right );
    return str;
}

inline
mystring operator+( char c, const mystring& right )
{
    mystring    str( c );
    str.append( right );
    return str;
}

inline
mystring operator+( const mystring& left, const char* s )
{
    mystring    str( left );
    str.append( s );
    return str;
}

inline
mystring operator+( const mystring& left, char c )
{
    mystring    str( left );
    str.append( c );
    return str;
}

// -------------------------------------------------------------------------
// operator ==
//
inline
bool operator==( const mystring& left, const mystring& right ) { return left.compare( right ) == 0; }

inline
bool operator==( const char* ls, const mystring& right ) { return right.compare( ls ) == 0; }

inline
bool operator==( const mystring& left, const char* rs ) { return left.compare( rs ) == 0; }

// -------------------------------------------------------------------------
// operator !=
//
inline
bool operator!=( const mystring& left, const mystring& right ) { return left.compare( right ) != 0; }

inline
bool operator!=( const char* ls, const mystring& right ) { return right.compare( ls ) != 0; }

inline
bool operator!=( const mystring& left, const char* rs ) { return left.compare( rs ) != 0; }

// -------------------------------------------------------------------------
// operator <
//
inline
bool operator<( const mystring& left, const mystring& right ) { return left.compare( right ) < 0; }

inline
bool operator<( const char* ls, const mystring& right ) { return right.compare( ls ) > 0; }

inline
bool operator<( const mystring& left, const char* rs ) { return left.compare( rs ) < 0; }

// -------------------------------------------------------------------------
// operator >
//
inline
bool operator>( const mystring& left, const mystring& right ) { return left.compare( right ) > 0; }

inline
bool operator>( const char* ls, const mystring& right ) { return right.compare( ls ) < 0; }

inline
bool operator>( const mystring& left, const char* rs ) { return left.compare( rs ) > 0; }

// -------------------------------------------------------------------------
// operator <=
//
inline
bool operator<=( const mystring& left, const mystring& right ) { return left.compare( right ) <= 0; }

inline
bool operator<=( const char* ls, const mystring& right ) { return right.compare( ls ) >= 0; }

inline
bool operator<=( const mystring& left, const char* rs ) { return left.compare( rs ) <= 0; }

// -------------------------------------------------------------------------
// operator >=
//
inline
bool operator>=( const mystring& left, const mystring& right ) { return left.compare( right ) >= 0; }

inline
bool operator>=( const char* ls, const mystring& right ) { return right.compare( ls ) <= 0; }

inline
bool operator>=( const mystring& left, const char* rs ) { return left.compare( rs ) >= 0; }



#endif//__mystring_h__
