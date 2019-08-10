/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __ConfigFile__
#define __ConfigFile__

#include <stdio.h>

// _________________________________________________________________________
// Class ConfigFile
//
// processing configuration (options) files
//
class ConfigFile {

public:
    ConfigFile( const char* section, const char* filename );
    ~ConfigFile();

    bool    GetInt( const char* key, int* value );          //read value of key
    bool    GetFloat( const char* key, float* value );      //read float value of key
    bool    GetDouble( const char* key, double* value );    //read double value of key
                                                            //return key value in the section
    bool    GetString( const char* key, char* value, unsigned size );


    bool    WriteInt( const char* key, int value );         //write key value when key is integer
    bool    WriteFloat( const char* key, float value );     //write float value of key
    bool    WriteDouble( const char* key, double value );   //write double value of key
                                                            //write key value when key is string
    bool    WriteString( const char* key, const char* value );

    const char* GetSectionName() const      { return section_name;  }
    const char* GetFileName() const         { return filename;      }

protected:
                                                            //read key from file given section name
    bool    GetCFPrivateProfileInt( const char* key, int* value );
                                                            //read key of float type from file given section name
    bool    GetCFPrivateProfileFloat( const char* key, float* value );
                                                            //read key of double type from file given section name
    bool    GetCFPrivateProfileDouble( const char* key, double* value );

                                                            //read string value of key given section name
    bool    GetCFPrivateProfileString( const char* key, char* value, unsigned size );
                                                            //write string value of key
    bool    WriteCFPrivateProfileString( const char* key, const char* value );

protected:
                                                            //find section in file by name
    bool    FindFileSection( FILE* fd, bool bwrite = false );
                                                            //find key in file by name
    bool    GetSectionKey( FILE*, const char* key, char* value, size_t size, long* pos = NULL, long* lastpos = NULL );
                                                            //write key with its value in the file
    bool    WriteSectionKey( FILE*, const char* key, const char* value, bool );
                                                            //insert key in file after the position specified
    void    InsertSectionKeyAfter( FILE*, const char* key, const char* value, long pos, long lastpos, bool );
                                                            //helper method to format error messages
    void    Format_Message( const char* message, const char* );

    void    FlushToEOF( FILE* );                            //flush file
    long    TellPosition( FILE* );                          //tell position
    void    SeekToPosition( FILE*, long pos, int whence );  //seek to position


    void    ProcessValue( char* value );                    //process key value by removing leading and trailing spaces

    char*       GetErrBuffer()          { return errbuffer; }
    const char* GetErrBuffer() const    { return errbuffer; }
    int         GetErrBufLength() const { return errbuflen; }

private:
    const char* section_name;       //section name in file
    const char* filename;           //configuration file name

    char*       errbuffer;          //buffer for error message
    int         errbuflen;          //length of buffer

private:
    static const char*  s_equal;            //sign of equality
    static const char*  s_bracket;          //sign of bracket to indicate the beginning of section name
    static const char*  s_newline;          //indication of newline
};

#endif//__ConfigFile__
