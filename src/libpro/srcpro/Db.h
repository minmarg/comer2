/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __Db_h__
#define __Db_h__

#include "liblib/mybase.h"

#include <stdio.h>

#include <fstream>

#include "liblib/BinarySearchStructure.h"
#include "libpro/srcpro/SEGProfile.h"
#include "DbProProbs.h"

namespace pmodel {
class PMTransModel;
class PMProfileModel;
}

// _________________________________________________________________________
// Class Db
//
// Database of profiles
//
class Db
{
public:
    //enumerations for signature...
    enum {  BEGTEXT,    //signature text
            DATAVER,    //version number
            END
    };
    //enumeration for database extensions...
    enum TFile{
            DbFlMAIN,//main profile database file
            DbFlPPRO,//file of profile pair probabilities
            DbFlMAINBin,//main profile database file in binary format
            DbFlN
    };

    //Structure type for keeping information about profile file
    struct DbP1 {
        DbP1( mystring& name, size_t len ): fname_(name), prolen_(len) {}
        DbP1( const char* name, size_t len ): fname_(name), prolen_(len) {}
        mystring fname_;//filename of profile
        size_t  prolen_;//profile length
    };

public:
    //typedefs...
    typedef void ( Db::*PMETHOD )( const char*, void* );
    //
    Db( const char* name/*, int = FrequencyStore::TFreqSimpleVector */);
    Db( const char* output, const char* directory/*, int = FrequencyStore::TFreqSimpleVector */);
    Db( const char* output, char* arguments[], int no_args/*, int = FrequencyStore::TFreqSimpleVector */);
    virtual ~Db();

//     void                    ReadInFrequencies();        //read frequencies

    virtual void    Open();//open database
    virtual void    Close( TFile = DbFlN );//close a file of the database
                    //read one profile information...
    bool            Next( pmodel::PMProfileModel&, pmodel::PMTransModel&/*, int goc, int gec, bool */);
    bool            Eof() const { return GetNextSym() == EOF; }

    void            MakeBin();//make binary database
    void            Make();//make the database: the main function

    const char*     GetDbMainName();//get database main filename
    const char*     GetDbMainBinName();//get database binary main filename
    const char*     GetDbProbsName();//get probabilities filename

    void        GetDbPos( fpos_t* ) const;
    void        SetDbPos( const fpos_t* );

    const char*     GetDbName() const           { return dbname_; }
//     size_t          GetNoVectors() const        { return no_vectors_; }
    size_t          GetNoSequences() const      { return no_sequences_; }
    uint64_mt       GetDbSize() const           { return db_size_; }

    bool            GetUsingSeg() const         { return usingseg_; }
    size_t          GetSegWinLength() const     { return segwinlen_; }
    float           GetSegLowEntropy() const    { return seglowentropy_; }
    float           GetSegHighEntropy() const   { return seghighentropy_; }
    float           GetSegDistance() const      { return segdistance_; }
    void            SetSegParameters( size_t winlen, float lowent, float highent, float distance )  {
                        SetUsingSeg( true );
                        SetSegWinLength( winlen );
                        SetSegLowEntropy( lowent );
                        SetSegHighEntropy( highent );
                        SetSegDistance( distance );
                    }

//     const FrequencyStore*   GetStore() const { return store; }

    const DbProProbs* GetProProbs() const { return &proprobs_; }


//     static mystring                 GetDistributionText( int type );
//     static TFVectorProbabilities    GetDistributionType( const mystring& distrstr );
//     static TFVectorProbabilities    GetDistributionType()           { return FrequencyStore::GetDistributionType(); }
//     static void SetDistributionType( TFVectorProbabilities type )   { FrequencyStore::SetDistributionType( type ); }

protected:
    explicit Db( /*int = FrequencyStore::TFreqSimpleVector */);

    void    Init();//initialization

    void    ProcessInput( PMETHOD, void* );//process profiles given in directory
    bool    PreprocessNeeded() const;
    void    PreprocessProfile( pmodel::PMProfileModel&, pmodel::PMTransModel& );
    void    ProcessFile( const char* filename, void* );//process profile given in file
    void    WriteFile( const char* filename, void* );//append profile to the database
//     void    WriteFrequencies( FILE* fd );//write frequency vectors to file descriptor


    void        SetUsingSeg( bool value )           { usingseg_ = value; }
    void        SetSegWinLength( size_t value )     { segwinlen_ = value; }
    void        SetSegLowEntropy( float value )     { seglowentropy_ = value; }
    void        SetSegHighEntropy( float value )    { seghighentropy_ = value; }
    void        SetSegDistance( float value )       { segdistance_ = value; }


    FILE*       GetDbDesc() const { return db_fp_[DbFlMAIN]; }
    void        SetDbDesc( FILE* fp ) { 
        if( db_fp_[DbFlMAIN])
            fclose( db_fp_[DbFlMAIN]); 
        db_fp_[DbFlMAIN] = fp; 
    }

    const char* GetDataDir() const { return data_dir_; }
    const char* const* GetProfiles() const { return profiles_; }
    const char* GetProfile( int ) const; 
    int         GetNoProfs() const { return no_profs_; }

    float   GetProbNormTerm() const { return probnorm_; }
    void    IncProbNormTerm( float prob ) { probnorm_ += prob; }
    void    ResetProbNormTerm() { probnorm_ = 0.0f; }

//     void    SetNoVectors( size_t value ) { no_vectors_ = value; }
//     void    IncNoVectors() { no_vectors_++; }
    void    SetNoSequences( size_t value ) { no_sequences_ = value; }
    void    IncNoSequences() { no_sequences_++; }
    void    SetDbSize( uint64_mt value ) { db_size_ = value; }
    void    IncDbSize( size_t amount ) { db_size_ += amount; }

//     int     GetFreqStrType() const { return fstrtype_; }

//     bool        AreVectorsConsistent( float reps ) const;//verify the consistency of the vectors

    void        WriteBinaryHeader( std::ofstream& );
    void        WriteTextHeader( FILE* );
    void        ReadTextHeader( FILE* );

    void        PutHeader( FILE* );//put database header
    void        GetHeader( FILE* );//get database header

    void        PutSignature( FILE* );//put database signature
    void        GetSignature( FILE* );//get database signature

    int     GetNextSym() const { return nextsym_; }
    void    SetNextSym( int value ) { nextsym_ = value; }

private:
    static int ProlenComparer( const void* key1, const void* key2, void* pars );
    void    NewSortedFiles();
    void    DestroySortedFiles();

private:
    const char*         dbname_;            //database name
    const char*         data_dir_;          //directory of profiles
    const char* const*  profiles_;          //names of profiles used to make up the database
    const int           no_profs_;          //number of profiles given as argument

    BinarySearchStructure*  srtfiles_;      //profile filenames sorted by profile length

    float               probnorm_;          //normalization constant for probabilities
//     size_t              no_vectors_;        //overall number of frequency vectors within the database
    size_t              no_sequences_;      //number of profiles in the database
    uint64_mt           db_size_;           //size of database

    int                 nextsym_;           //next symbol from file
    FILE*               db_fp_[DbFlN];      //file descriptors of the database
//     Serializer          serializer;         //object to serialize data to file

    char*               name_buffer_;       //auxiliary buffer


    bool                usingseg_;          //seg in use
    size_t              segwinlen_;         //seg window length
    float               seglowentropy_;     //seg low entropy threshold
    float               seghighentropy_;    //seg high entropy threshold
    float               segdistance_;       //seg distance to consider vectors equivalent


//     FrequencyStore*     store_;             //store of frequency vectors
//     const int           fstrtype_;          //frequency structure type

    DbProProbs          proprobs_;          //profile pair probabilities

public:
    static const char*  db_signature_[];//database signature
    static const char*  db_extensions_[];//extensions of the database files
    //
    static const char*  patstrDBBINVER;
    static const char*  patstrDBVER;
    static const char*  patstrNoPROFS;
    static const char*  patstrSIZE;
    // static const char* patstrDISTR = "DISTR:";
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* Db::GetProfile( int n ) const
{
#ifdef __DEBUG__
    if( n < 0 || no_profs_ <= n )
        throw MYRUNTIME_ERROR( "Db::GetProfile: Memory access error." );
#endif
    return profiles_[ n ];
}

// -------------------------------------------------------------------------
// GetDbPos: get current position of the main database file
//
inline
void Db::GetDbPos( fpos_t* pos ) const
{
#ifdef __DEBUG__
    if( !GetDbDesc())
        throw MYRUNTIME_ERROR( "Db::GetDbPos: Null file descriptor." );
    if( !pos )
        throw MYRUNTIME_ERROR( "Db::GetDbPos: Memory access error." );
#endif
    if( fgetpos( GetDbDesc(), pos ))
        warning( "Db::GetDbPos: Getting file position failed.", true, 0 );
}

// SetDbPos: set current position for the main database file
inline
void Db::SetDbPos( const fpos_t* pos )
{
#ifdef __DEBUG__
    if( !GetDbDesc())
        throw MYRUNTIME_ERROR( "Db::SetDbPos: Null file descriptor." );
    if( !pos )
        throw MYRUNTIME_ERROR( "Db::SetDbPos: Memory access error." );
#endif
    if( fsetpos( GetDbDesc(), pos ))
        warning( "Db::SetDbPos: Setting file position failed.", true, 0 );

    int nsym = 0;
    nsym = fgetc( GetDbDesc());
    ungetc( nsym, GetDbDesc());
    SetNextSym( nsym );
}

// // -------------------------------------------------------------------------
// // AreVectorsConsistent: verify whether the number of vectors is valid
// //
// inline
// bool Db::AreVectorsConsistent( float reps ) const
// {
//     if( !GetStore())
//         return false;
//     return GetStore()->IsConsistent( reps, no_vectors_ );
// }

// // -------------------------------------------------------------------------
// // GetDistributionText: get distribution description given type
// //
// inline
// mystring Db::GetDistributionText( int distr )
// {
//     return( FrequencyStore::GetDistributionText( distr ));
// }

// // GetDistributionType: get distribution type given description
// //
// inline
// TFVectorProbabilities Db::GetDistributionType( const mystring& distrstr )
// {
//     return( FrequencyStore::GetDistributionType( distrstr ));
// }

// -------------------------------------------------------------------------
// Db::ProlenComparer: compare the lengths of two profiles encoded by the 
// pointers
inline
int Db::ProlenComparer( const void* key1, const void* key2, void* )
{
    const DbP1* p11 = (const DbP1*)key1;
    const DbP1* p12 = (const DbP1*)key2;
#ifdef __DEBUG__
    if( p11 == NULL || p12 == NULL )
        throw MYRUNTIME_ERROR( "Db::ProlenComparer: Null arguments." );
#endif
    //sorting in reverse order
    return (int)((ssize_t)p12->prolen_ - (ssize_t)p11->prolen_);
}
// -------------------------------------------------------------------------
// Db::NewSortedFiles: create a new object for sorted files
inline
void Db::NewSortedFiles()
{
    DestroySortedFiles();
    const size_t initsize = 8192;
    srtfiles_ = new BinarySearchStructure( &ProlenComparer, initsize, true/*keep*/, this );
    if( srtfiles_ == NULL )
        throw MYRUNTIME_ERROR( "Db::NewSortedFiles: Not enough memory." );
}
// -------------------------------------------------------------------------
// Db::DestroySortedFiles: destroy the object of sorted files
inline
void Db::DestroySortedFiles()
{
    if( srtfiles_ == NULL )
        return;
    for( size_t n = 0; n < srtfiles_->GetSize(); n++ ) {
        const DbP1* p1 = (const DbP1*)srtfiles_->GetValueAt(n);
        if( p1 )
            delete p1;
    }
    delete srtfiles_;
    srtfiles_ = NULL;
}

#endif//__Db_h__
