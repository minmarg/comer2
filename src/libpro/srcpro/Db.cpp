/***************************************************************************
 *   Copyright (C) 2013-2019 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>

#include "liblib/mydirent.h"

#include "liblib/BinarySearchStructure.h"
#include "libpro/srcpro/PMTransModel.h"
#include "libpro/srcpro/PMProfileModel.h"
#include "libpro/srcpro/SEGProfile.h"
#include "DbProProbs.h"
#include "Db.h"

const char* Db::db_signature_[] = {
    "PROFILE DATABASE VERSION ",
    pmodel::PMProfileModel::GetDataVersion(),//dbversion,
    NULL
};

//extensions of the database files
const char* Db::db_extensions_[] = {
    ".prd",
    ".prb",
    ".bin",
    NULL
};

// -------------------------------------------------------------------------
// Constructor
//
Db::Db( const char* name/*, int fstype */)
:
    dbname_( name ),
    data_dir_( NULL ),
    profiles_( NULL ),
    no_profs_( 0 ),

    srtfiles_( NULL ),

    probnorm_( 0.0f ),
//     no_vectors_( 0 ),
    no_sequences_( 0 ),
    db_size_( 0 ),
    nextsym_( 0 ),
    name_buffer_( NULL ),

    usingseg_( false ),
    segwinlen_( 0 ),
    seglowentropy_( -1.0f ),
    seghighentropy_( -1.0f ),
    segdistance_( 0.0f )

//     store_( NULL ),
//     fstrtype_( fstype )
{
    MYMSG("Db::Db",4);
    Init();
}

// Constructor
//
Db::Db( const char* output, const char* directory/*, int fstype */)
:
    dbname_( output ),
    data_dir_( directory ),
    profiles_( NULL ),
    no_profs_( 0 ),

    srtfiles_( NULL ),

    probnorm_( 0.0f ),
//     no_vectors_( 0 ),
    no_sequences_( 0 ),
    db_size_( 0 ),
    nextsym_( 0 ),
    name_buffer_( NULL ),

    usingseg_( false ),
    segwinlen_( 0 ),
    seglowentropy_( -1.0f ),
    seghighentropy_( -1.0f ),
    segdistance_( 0.0f )

//     store_( NULL ),
//     fstrtype_( fstype )
{
    MYMSG("Db::Db",4);
    Init();
}

// Constructor
//
Db::Db( const char* output, char* arguments[], int no_args/*, int fstype */)
:
    dbname_( output ),
    data_dir_( NULL ),
    profiles_( arguments ),
    no_profs_( no_args ),

    srtfiles_( NULL ),

    probnorm_( 0.0f ),
//     no_vectors_( 0 ),
    no_sequences_( 0 ),
    db_size_( 0 ),
    nextsym_( 0 ),
    name_buffer_( NULL ),

    usingseg_( false ),
    segwinlen_( 0 ),
    seglowentropy_( -1.0f ),
    seghighentropy_( -1.0f ),
    segdistance_( 0.0f )

//     store_( NULL ),
//     fstrtype_( fstype )
{
    MYMSG("Db::Db",4);
    Init();
}

// Default constructor
//
Db::Db(/* int fstype */)
:
    dbname_( NULL ),
    data_dir_( NULL ),
    profiles_( NULL ),
    no_profs_( 0 ),
    srtfiles_( NULL ),
    probnorm_( 0.0f ),
//     no_vectors_( 0 ),
    no_sequences_( 0 ),
    db_size_( 0 ),
    nextsym_( 0 ),
    name_buffer_( NULL ),

    usingseg_( false ),
    segwinlen_( 0 ),
    seglowentropy_( -1.0f ),
    seghighentropy_( -1.0f ),
    segdistance_( 0.0f )

//     store_( NULL ),
//     fstrtype_( fstype )
{
    throw MYRUNTIME_ERROR( "Db::Db: Default initialization is not allowed." );
}

// -------------------------------------------------------------------------
// Destructor
//
Db::~Db()
{
    MYMSG("Db::~Db",4);
    Close();
    if( name_buffer_ )
        free( name_buffer_ );
//     if( store_ )
//         delete store_;
    DestroySortedFiles();
}

// -------------------------------------------------------------------------
// Init: initialize members
//
void Db::Init()
{
    MYMSG("Db::Init",3);

    for( int n = 0; n < DbFlN; n++ )
        db_fp_[n] = NULL;

    if( !GetDbName())
        return;

    name_buffer_ = (char*)malloc( strlen(GetDbName()) + strlen(db_extensions_[DbFlMAIN]) + 1 );

    if( !name_buffer_ )
        throw MYRUNTIME_ERROR( "Db::Init: Not enough memory." );

    strcpy( name_buffer_, GetDbName());

//     store_ = new FrequencyStore( GetFreqStrType());

//     if( !store_ )
//         throw MYRUNTIME_ERROR( "Db::Init: Not enough memory." );
}

// -------------------------------------------------------------------------
// GetDbMainName: get the main database name
//
const char* Db::GetDbMainName()
{
    strcpy( name_buffer_ + strlen(GetDbName()), db_extensions_[DbFlMAIN] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// GetDbMainBinName: get the name of the main binary database file
//
const char* Db::GetDbMainBinName()
{
    strcpy( name_buffer_ + strlen( GetDbName()), db_extensions_[DbFlMAINBin] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// GetDbProbsName: get the filename of profile pair probabilities 
//
const char* Db::GetDbProbsName()
{
    strcpy( name_buffer_ + strlen(GetDbName()), db_extensions_[DbFlPPRO] );
    return  name_buffer_;
}

// -------------------------------------------------------------------------
// Open: open the database and initialize a file descriptor
//
void Db::Open()
{
    MYMSG("Db::Open",3);

    if( !GetDbMainName())
        throw MYRUNTIME_ERROR( "Db::Open: Unable to open database." );

    if( db_fp_[DbFlMAIN] != NULL )
        throw MYRUNTIME_ERROR( "Db::Open: Database has been already opened." );

//     db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "rb" );//obsolete
    db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "r" );

    if( db_fp_[DbFlMAIN] == NULL )
        throw MYRUNTIME_ERROR( "Db::Open: Failed to open database." );

    try {
        SetNextSym( 0 );
//         GetHeader( db_fp_[DbFlMAIN] );//obsolete
        ReadTextHeader( db_fp_[DbFlMAIN] );

        //read profile pair probabilities
        proprobs_.ReadProbs( GetDbProbsName());

    } catch( myexception const& ex )
    {
        Close();
        throw ex;
    }
}

// -------------------------------------------------------------------------
// Close: close the database
//
void Db::Close( TFile which )
{
    MYMSG("Db::Close",3);

    if( DbFlN <= which || DbFlN < 0 ) {
        for( int n = 0; n < DbFlN; n++ )
            if( db_fp_[n] != NULL ) {
                fclose( db_fp_[n] );
                db_fp_[n] = NULL;
            }
    } else
        if( db_fp_[which] != NULL ) {
            fclose( db_fp_[which] );
            db_fp_[which] = NULL;
        }
}

// -------------------------------------------------------------------------
// Next: read next profile from the database; 
// return a flag indicating whether the database end has been reached
//
bool Db::Next( pmodel::PMProfileModel& pm, pmodel::PMTransModel& tm/*,
                    int gapopencost, int gapextncost, bool fixedcosts */)
{
    MYMSG("Db::Next",3);

    int nsym = 0;

    if( GetNextSym() == EOF )
        return false;

    if( GetDbDesc() == NULL )
        throw MYRUNTIME_ERROR( "Db::Next: Db is not opened." );

    try {
//         serializer.DeserializeProfile( pm, tm, GetDbDesc());//obsolete
        TextReadProfile( GetDbDesc(), pm, tm );

    } catch( myexception const& ex )
    {
        Close();
        throw ex;
    }

//     pm.CheckForAllZeros();    //necessary to ensure valid scores

    tm.Prepare();

    if( GetNextSym() != EOF ) {
        nsym = fgetc( GetDbDesc());
        ungetc( nsym, GetDbDesc());
        SetNextSym( nsym );
    }

    return true;
}

// -------------------------------------------------------------------------
// MakeBin: make binary representation of profile database;
//
void Db::MakeBin()
{
    MYMSG("Db::MakeBin",3);
    const mystring preamb = "Db::MakeBin: ";
    myruntime_error mre;
    pmodel::PMProfileModel pm;
    pmodel::PMTransModel tm;
    const char* biname = GetDbMainBinName();
    std::ofstream fp( biname,std::ios::binary|std::ios::out);

    if(fp.bad() || fp.fail())
        throw MYRUNTIME_ERROR(preamb+"Failed to open file for writing: "+biname);

    Open();

    try {
        WriteBinaryHeader(fp);
        while( Next(pm,tm)) {
            BinaryWriteProfile(fp,pm,tm);
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    Close();

    if( mre.isset())
        throw mre;
}

// -------------------------------------------------------------------------
// Make: make a database of profiles;
// the main function to create a database
//
void Db::Make()
{
    MYMSG("Db::Make",3);

    const mystring preamb = "Db::Make: ";

    if( !GetDbMainName())
        throw MYRUNTIME_ERROR( preamb + "Unable to make database." );

    if( db_fp_[DbFlMAIN] != NULL /*|| db_fp_[DbFlFREQ] != NULL */)
        throw MYRUNTIME_ERROR( preamb + "Null file descriptors." );

    message( "Preprocessing..." );

    try {
        NewSortedFiles();

        ResetProbNormTerm();

        //clear marginal counts before updating them
        proprobs_.ClearMargs();

        //1st PASS: find out attributes of a database to be created
        ProcessInput( &Db::ProcessFile, NULL );

        //calculate probabilities after marginal counts have been updated
        proprobs_.CalcProbs();

    //     db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "wb" );//obsolete
        db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "w" );

        if( db_fp_[DbFlMAIN] == NULL )
            throw MYRUNTIME_ERROR( preamb + "Failed to open file for writing.");

    //     db_fp_[DbFlFREQ] = fopen( GetDbFreqName(), "wb" );//obsolete
//         db_fp_[DbFlFREQ] = fopen( GetDbFreqName(), "w" );

//         if( db_fp_[DbFlFREQ] == NULL )
//             throw MYRUNTIME_ERROR( preamb + "Failed to open file for writing.");

        message( "Writing..." );
//         PutHeader( db_fp_[DbFlMAIN] );//obsolete
        WriteTextHeader( db_fp_[DbFlMAIN] );

        //{{2nd PASS: write database attributes and profiles
        if( srtfiles_ ) {
            for( size_t n = 0; n < srtfiles_->GetSize(); n++ ) {
                const DbP1* p1 = (const DbP1*)srtfiles_->GetValueAt(n);
                if( p1 )
                    WriteFile( p1->fname_.c_str(), db_fp_ );
            }
        }
        else
            ProcessInput( &Db::WriteFile, db_fp_ );
        //}}

//         message( "Writing frequencies..." );
//         WriteFrequencies( db_fp_[DbFlFREQ] );

        message( "Writing probabilities..." );
        proprobs_.WriteProbs( GetDbProbsName());

    } catch( myexception const& ex )
    {
        Close();
        throw ex;
    }

    message( "Finished." );
    Close();
}

// -------------------------------------------------------------------------
// ProcessInput: process profiles given explicitly or found in directory
//
void Db::ProcessInput( PMETHOD processing_method, void* args )
{
    MYMSG("Db::ProcessInput",3);

    if( GetDataDir()) {
        struct stat     info;
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetDataDir());
        mystring        dirname = GetDataDir() + mystring( DIRSEP );
        mystring        fullname;

        if( !direct )
            throw MYRUNTIME_ERROR( 
                mystring("Db::ProcessInput: Unable to process directory ") + GetDataDir());

        while(( entry = readdir( direct )))
        {
            fullname = dirname + entry->d_name;

            if( stat( fullname.c_str(), &info ) == -1 ) {
                error( mystring(mystring("Cannot stat '")+mystring(entry->d_name)+"'; skipping.").c_str());
                continue;
            }
            if(( S_IFMT & info.st_mode ) != S_IFREG )
//             if( !S_ISREG( info.st_mode ))
                continue;

            (this->*processing_method)( fullname.c_str()/*entry->d_name*/, args );
        }
        closedir( direct );

    } else
        if( GetProfiles()) {
            bool option = false;
            bool passed = false;//options passed

            for( int n = 0; n < GetNoProfs(); n++ ) {
                if( GetProfile(n) == NULL )
                    continue;
                if( !passed &&
                    GetProfile(n)[0] == '-' && GetProfile(n)[1] == '-' && GetProfile(n)[2] == 0 ) {
                    passed = true;
                    continue;
                }
                if( !passed && option ) {
                    option = false;
                    continue;
                }
                if( !passed && GetProfile(n)[0] == '-' ) {
                    if( GetProfile(n)[1] != 'v' && GetProfile(n)[2] == 0 )
                        option = true;//assume it is an option
                    continue;
                }
                (this->*processing_method)( GetProfile(n), args );
            }
        }
}

// -------------------------------------------------------------------------
// PreprocessNeeded: verify whether preprocessing applies
//
inline
bool Db::PreprocessNeeded() const
{
    return GetUsingSeg();
}
// PreprocessProfile: apply SEG algorithm to the profile if SEG is in use
//
void Db::PreprocessProfile( pmodel::PMProfileModel& pm, pmodel::PMTransModel& tm )
{
    MYMSG("Db::PreprocessProfile",3);

    if( GetUsingSeg()) {
        //SEG logic
        SEGProfile  segpro(
                pm,
                GetSegWinLength(),
                GetSegLowEntropy(),
                GetSegHighEntropy()
        );
        segpro.SetDistance( GetSegDistance());
        segpro.Run();
        segpro.MaskSeggedPositions( pm, tm );
    }
}

// -------------------------------------------------------------------------
// ProcessFile: read data from file and process it
//
void Db::ProcessFile( const char* filename, void* )
{
    MYMSG("Db::ProcessFile",3);

    pmodel::PMProfileModel  pm;
    pmodel::PMTransModel    tm;
    FILE* fp = NULL;
//     float prob = 0.0f;

    try {
        if(( fp = fopen( filename, "r" )) ==  NULL )
            throw MYRUNTIME_ERROR("Db::ProcessFile: Failed to open file.");

//         serializer.DeserializeProfile( pm, tm, filename );//obsolete
//         serializer.ReadProfile( filename, pm, tm );
        TextReadProfile( fp, pm, tm );
        PreprocessProfile( pm, tm );
        IncNoSequences();
        IncDbSize( pm.GetSize());

        //update marginal counts
        proprobs_.UpdateMargs( pm.GetEffNoSequences(), pm.GetSize());

        //sort filenames on the fly by profile length
        if( srtfiles_ )
            srtfiles_->Push( new DbP1( filename, pm.GetSize()));

//         //store all different observed weighted frequency vectors
//         pm.CheckForAllZeros();//necessary to ensure valid scores
// 
//         for( int n = 0; n < pm.GetSize(); n++ )
//         {
//             //omit columns representing Xs
//             if( pm.GetResidueAt(n) == X )
//                 continue;
// 
//             FrequencyVector vect(   *pm.GetObsFreqsAt(n), pm.GetFrequencyWeightAt(n), 
//                                      pm.GetMIDExpNoObservationsAt( n, PS_M ),
//                                     *pm.GetVectorAt(n), pm.GetInformationAt(n));
//             IncNoVectors();
//             FrequencyVector quest_vect = ( const char* )store_->Store( vect );
// 
//             if( quest_vect.GetVector() == vect.GetVector()) {
//                 //vector has been stored;
//                 //compute probability of this frequency vector
//                 prob = vect.ComputeProbability();
//                 //adjust probability-normalizing term
//                 IncProbNormTerm( prob );
//             }
//             else {
//                 //if the vector has not been stored, update the probability of the vector...
//                 quest_vect.UpdateProbability();
//                 //and destroy the duplicate
//                 vect.Destroy();
//             }
//         }

    } catch( myexception const& ex )
    {
        error( ex.what(), false );
        error( mystring(mystring("Db::ProcessFile: '")+filename+mystring("' skipped.")).c_str());
    }

    if( fp )
        fclose(fp);
}

// -------------------------------------------------------------------------
// WriteFile: append the file to the database
//
void Db::WriteFile( const char* filename, void* args )
{
    MYMSG("Db::WriteFile",3);

    FILE** FP = ( FILE** )args;

    if( !FP || !FP[DbFlMAIN] )
        return;

    const mystring preamb = "Db::WriteFile: ";
    pmodel::PMProfileModel  pm;
    pmodel::PMTransModel    tm;
    FILE* fp = NULL;
    mystring buffer;
    int emsg;

    try {
        if(( fp = fopen( filename, "r" )) ==  NULL )
            throw MYRUNTIME_ERROR( preamb + "Failed to open file.");

        if( PreprocessNeeded()) {
//             serializer.DeserializeProfile( pm, tm, filename );//obsolete
//             serializer.ReadProfile( filename, pm, tm );
            TextReadProfile( fp, pm, tm );
            PreprocessProfile( pm, tm );
//             serializer.SerializeProfile( pm, tm, FP[DbFlMAIN] );//obsolete
            //TextWriteProfile( FP[DbFlMAIN], pm, tm );
            TextWriteProfileCondensed( FP[DbFlMAIN], pm, tm );
        }
        else {
            //otherwise, copy data...
            while( !feof(fp)) {
                //read full line 
                if(( emsg = skip_comments( fp, buffer )) != 0 )
                    throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
                //write the line 
                if( fputs( buffer.c_str(), FP[DbFlMAIN]) < 0 )
                    throw MYRUNTIME_ERROR( preamb + "Failed to write data to file.");
            }
        }
    } catch( myexception const& ex )
    {
        error( ex.what(), false );
        error( mystring(mystring("Db::WriteFile: '")+filename+mystring("' skipped.")).c_str());
    }

    if( fp )
        fclose(fp);
}

// // -------------------------------------------------------------------------
// // WriteFrequencies: write frequency vectors to the database
// //
// void Db::WriteFrequencies( FILE* fd )
// {
//     if( !fd || !store_ )
//         return;
// 
//     const SimpleVector* frequencies = store_->GetFrequencies();
//     size_t              size = 0;
//     size_t              n;
// 
//     if( frequencies )
//         size = frequencies->GetSize();
// 
//     //write the total number of frequency vectors and the number of distinct
//     //frequency vectors, respectively
// //     Serializer::Write( fd, ( char* )&no_vectors_, sizeof( no_vectors_ ), 1 );//obsolete
// //     Serializer::Write( fd, ( char* )&size, sizeof( size ), 1 );//obsolete
// 
//     fprintf( fd, "%ld %ld%s", GetNoVectors(), size, NL );
// 
//     for( n = 0; n < size; n++ ) {
//         FrequencyVector vector(( const char* )frequencies->GetValueAt( n ));
//         //normalize probability before writing frequency vector
//         vector.NormalizeProbability( GetProbNormTerm());
// //         serializer.SerializeFrequencies( vector, fd );//obsolete
//         serializer.WriteVector( fd, vector );
//     }
// 
//     char    strbuf[UCHAR_MAX];
// 
//     sprintf( strbuf, "distinct vectors,     %9d", size );                       message( strbuf, false );
//     sprintf( strbuf, "1st-level collisions, %9d", store_->GetCollisions1());    message( strbuf, false );
//     sprintf( strbuf, "2nd-level collisions, %9d", store_->GetCollisions2());    message( strbuf );
// }

// // -------------------------------------------------------------------------
// // ReadInFrequencies: read frequency vectors
// //
// void Db::ReadInFrequencies()
// {
//     if( !store_ )
//         return;
// 
//     const mystring  preamb = "Db::ReadInFrequencies: ";
//     myruntime_error mre;
//     size_t          length, rbts;
//     const size_t    locsize = KBYTE;
//     char            locbuffer[locsize+1] = {0};
//     char*           p;
//     size_t          n;
//     int             emsg;
//     int             intval;
//     long long int   llintval;
//     size_t      size = 0;       //number of distinct vectors
//     float       reps = 0.0f;     //number of occurences in the database
// 
//     if( !GetDbFreqName())
//         throw MYRUNTIME_ERROR( preamb + "Unable to open database file of frequency vectors." );
// 
//     if( db_fp_[DbFlFREQ] != NULL )
//         throw MYRUNTIME_ERROR( preamb + "Frequency vectors file has been already opened." );
// 
// //     db_fp_[FREQ] = fopen( GetDbFreqName(), "rb" );//obsolete
//     db_fp_[DbFlFREQ] = fopen( GetDbFreqName(), "r" );
// 
//     if( db_fp_[DbFlFREQ] == NULL )
//         throw MYRUNTIME_ERROR( preamb + "Null descriptor of frequency vectors file." );
// 
//     try {
// //         Serializer::Read( db_fp_[DbFlFREQ], ( char* )&no_vectors_, sizeof( no_vectors_ ), 1 );//obsolete
// //         Serializer::Read( db_fp_[DbFlFREQ], ( char* )&size, sizeof( size ), 1 );//obsolete
// 
//         if(( emsg = skip_comments( db_fp_[DbFlFREQ], locbuffer, locsize, &length )) != 0 )
//             throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
// 
//         if( feof( db_fp_[DbFlFREQ] ) || !length )
//             throw MYRUNTIME_ERROR( preamb + "Wrong format of database file of frequency vectors." );
// 
//         p = locbuffer;
// 
//         //number of vectors
//         if( length <= size_t( p - locbuffer ))
//             throw MYRUNTIME_ERROR( preamb + "Wrong format of frequency vectors." );
// 
//         if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
//             throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
// 
//         if( llintval < 1 )
//             throw MYRUNTIME_ERROR( preamb + "Wrong format of frequency vectors" );
// 
//         SetNoVectors( llintval );
//         p += rbts;
// 
//         //number of distinct vectors
//         if( length <= size_t( p - locbuffer ))
//             throw MYRUNTIME_ERROR( preamb + "Wrong format of frequency vectors." );
// 
//         if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
//             throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
// 
//         if( llintval < 1 )
//             throw MYRUNTIME_ERROR( preamb + "Wrong format of frequency vectors" );
// 
//         size = llintval;
//         p += rbts;
// 
// 
//         //allocate memory at once to avoid sequential reallocations
//         store_->SetNoFrequencyVectors( size );
// 
//         for( n = 0; n < size; n++ ) {
//             FrequencyVector vector;
//             //function is smart to allocate required memory
// //             serializer.DeserializeFrequencies( vector, db_fp_[DbFlFREQ] );//obsolete
//             serializer.ReadVector( db_fp_[DbFlFREQ], vector );
// 
//             if( store_->Store( vector ) != vector.GetVector())
//                 throw MYRUNTIME_ERROR( preamb + "Failed to store frequency vector." );
// 
//             reps += vector.GetProbability() + 1.0f;
//             vector.FinalProbability( GetNoVectors());
//         }
// 
//         if( ! AreVectorsConsistent( reps ))
//             throw MYRUNTIME_ERROR( preamb + "Frequency vectors corrupted.");
// 
//     } catch( myexception const& ex )
//     {
//         mre = ex;
//     }
// 
//     Close( DbFlFREQ );
// 
//     if( mre.isset())
//         throw mre;
// }



// -------------------------------------------------------------------------
// PutHeader: write header information to the database
//
void Db::PutHeader( FILE* fp )
{
    MYMSG("Db::PutHeader",3);

//     TFVectorProbabilities   distrib = GetDistributionType();
    PutSignature( fp );
//     Serializer::Write( fp, ( char* )&no_sequences_, sizeof( no_sequences_ ), 1 );
//     Serializer::Write( fp, ( char* )&db_size_, sizeof( db_size_ ), 1 );
//     Serializer::Write( fp, ( char* )&distrib, sizeof( distrib ), 1 );

    size_t  ret;
    ret = fwrite((char*)&no_sequences_, sizeof(no_sequences_), 1, fp );
    if( ret != 1 )
        throw MYRUNTIME_ERROR( "Db::PutHeader: Write failed." );

    ret = fwrite((char*)&db_size_, sizeof(db_size_), 1, fp );
    if( ret != 1 )
        throw MYRUNTIME_ERROR( "Db::PutHeader: Write failed." );
}

// -------------------------------------------------------------------------
// GetHeader: read header information from the database
//
void Db::GetHeader( FILE* fp )
{
    MYMSG("Db::GetHeader",3);

//     TFVectorProbabilities   distrib;

    GetSignature( fp );
//     Serializer::Read( fp, ( char* )&no_sequences_, sizeof( no_sequences_ ), 1 );
//     Serializer::Read( fp, ( char* )&db_size_, sizeof( db_size_ ), 1 );
//     Serializer::Read( fp, ( char* )&distrib, sizeof( distrib ), 1 );

    size_t  ret;
    ret = fread((char*)&no_sequences_, sizeof(no_sequences_), 1, fp );
    if( ret != 1 )
        throw MYRUNTIME_ERROR( "Db::GetHeader: Read failed." );

    ret = fread((char*)&db_size_, sizeof(db_size_), 1, fp );
    if( ret != 1 )
        throw MYRUNTIME_ERROR( "Db::GetHeader: Read failed." );

//     switch( distrib ) {
//         case DISCRETE:
//         case PROVECTOR:
//         case MULTINOMIAL:
//             break;
// 
//         default:
//             throw MYRUNTIME_ERROR( "Db::GetHeader: Unknown vector distribution type." );
//     }
// 
//     SetDistributionType( distrib );
}

// -------------------------------------------------------------------------
// PutSignature: write database signature
//
void Db::PutSignature( FILE* fp )
{
    MYMSG("Db::PutSignature",3);

    size_t ret, len;
    for( int n = 0; db_signature_[n]; n++ ) {
//         Serializer::Write( fp, ( char* )db_signature_[n], 1, strlen( db_signature_[n] ));
        len = strlen(db_signature_[n]);
        ret = fwrite((char*)db_signature_[n], 1, len, fp );
        if( ret != len )
            throw MYRUNTIME_ERROR( "Db::PutSignature: Write failed." );
    }
}

// -------------------------------------------------------------------------
// GetSignature: read database signature and check it
//
void Db::GetSignature( FILE* fp )
{
    MYMSG("Db::GetSignature",3);

    char locsig[BUF_MAX];
    size_t ret, len;
    //
    for( int n = 0; db_signature_[n]; n++ ) {
//         Serializer::Read( fp, locsig, 1, strlen( db_signature_[n] ));
        len = strlen(db_signature_[n]);
        ret = fread((char*)locsig, 1, len, fp );
        if( ret != len )
            throw MYRUNTIME_ERROR( "Db::GetSignature: Read failed." );

        if( strncmp( db_signature_[n], locsig, len )) {
            if( n == DATAVER )
                warning( "Db::GetSignature: Inconsistent database version number." );
            else
                throw MYRUNTIME_ERROR( "Db::GetSignature: Invalid database format." );
        }
    }
}



// =========================================================================

const char* Db::patstrDBBINVER = "COMER profile bin database v";
const char* Db::patstrDBVER = "COMER profile database v";
const char* Db::patstrNoPROFS = "PROFS:";
const char* Db::patstrSIZE = "SIZE:";

// =========================================================================
// WriteTextHeader: write database header to file descriptor
//
void Db::WriteBinaryHeader( std::ofstream& fp )
{
    MYMSG("Db::WriteBinaryHeader",3);
    fp.write(patstrDBBINVER,strlen(patstrDBBINVER));
    fp.write(pmodel::PMProfileModel::GetDataVersion(),strlen(pmodel::PMProfileModel::GetDataVersion()));
    fp.write(reinterpret_cast<const char*>(&no_sequences_),sizeof(no_sequences_));
    fp.write(reinterpret_cast<const char*>(&db_size_),sizeof(db_size_));
}

// =========================================================================
// WriteTextHeader: write database header to file descriptor
//
void Db::WriteTextHeader( FILE* fp )
{
    MYMSG("Db::WriteTextHeader",3);

    if( !fp )
        return;

//     TFVectorProbabilities   distrib = GetDistributionType();
//     mystring                distrstr = GetDistributionText( distrib );

//     if( distrstr.empty())
//         throw MYRUNTIME_ERROR( "Db::WriteTextHeader: Unknown vector distribution type." );

    fprintf( fp, "%s%s%s", patstrDBVER, pmodel::PMProfileModel::GetDataVersion(), NL );
    fprintf( fp, "%-9s%zu%s", patstrNoPROFS, GetNoSequences(), NL );
    fprintf( fp, "%-9s%llu%s", patstrSIZE, GetDbSize(), NL );
//     fprintf( fp, "%-9s%s%s", patstrDISTR, distrstr.c_str(), NL );
}

// -------------------------------------------------------------------------
// ReadTextHeader: read database header from file descriptor
// TODO: reimplement local variables to use heap: locbuffer
//
void Db::ReadTextHeader( FILE* fp )
{
    MYMSG("Db::ReadTextHeader",3);

    if( !fp )
        return;

    const mystring  preamb = "Db::ReadTextHeader: ";
    size_t          length, rbts;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    int             emsg;
//     int             intval;
    long long int   llintval;
//     mystring        distrstr;
//     TFVectorProbabilities   distrib;


    //read version number
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( locbuffer, patstrDBVER )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( patstrDBVER );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( p, pmodel::PMProfileModel::GetDataVersion())) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Inconsistent database version number." );


    //read number of profiles
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( locbuffer, patstrNoPROFS )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( patstrNoPROFS );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( llintval < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid number of profiles." );

    SetNoSequences( llintval );


    //read database size
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( p = strstr( locbuffer, patstrSIZE )) == NULL )
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    p += strlen( patstrSIZE );

    if( length <= size_t( p - locbuffer ))
        throw MYRUNTIME_ERROR( preamb + "Wrong database format." );

    if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
        throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));

    if( llintval < 1 )
        throw MYRUNTIME_ERROR( preamb + "Invalid database size." );

    SetDbSize( llintval );


//     //read distribution description
//     if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
//         throw MYRUNTIME_ERROR( preamb + TranslateReadError( emsg ));
// 
//     if( feof( fp ) || !length )
//         throw MYRUNTIME_ERROR( preamb + "Wrong database format." );
// 
//     if(( p = strstr( locbuffer, patstrDISTR )) == NULL )
//         throw MYRUNTIME_ERROR( preamb + "Wrong database format." );
// 
//     p += strlen( patstrDISTR );
//     for( ; *p == ' ' || *p == '\t' ; p++ );
// 
//     if( length <= size_t( p - locbuffer ))
//         throw MYRUNTIME_ERROR( preamb + "Wrong database format." );
// 
//     distrstr = p;
//     distrib = GetDistributionType( distrstr );
// 
//     if( DTypeUNKNOWN <= distrib )
//         throw MYRUNTIME_ERROR( preamb + "Unrecognized distribution for frequency vectors." );
// 
//     SetDistributionType( distrib );
}
