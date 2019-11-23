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
#include "libmycu/cupro/PM2DVectorFields.h"
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

    nmaxopenfiles_(cMaxNOpenedFiles),
    nopenfiles_(0),
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

    nmaxopenfiles_(cMaxNOpenedFiles),
    nopenfiles_(0),
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

    nmaxopenfiles_(cMaxNOpenedFiles),
    nopenfiles_(0),
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
    nmaxopenfiles_(cMaxNOpenedFiles),
    nopenfiles_(0),
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
    ClearFPVector();
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

    db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "r" );

    if( db_fp_[DbFlMAIN] == NULL )
        throw MYRUNTIME_ERROR( "Db::Open: Failed to open database." );

    try {
        SetNextSym( 0 );
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
// MakeBin_v2_2: make binary representation of profile database; v2.2
//
void Db::MakeBin_v2_2()
{
    MYMSG("Db::MakeBin_v2_2",3);
    const mystring preamb = "Db::MakeBin_v2_2: ";
    myruntime_error mre;
    pmodel::PMProfileModel pm;
    pmodel::PMTransModel tm;
    const char* biname = GetDbMainBinName();
    std::ofstream fp( biname,std::ios::binary|std::ios::out);

    if(fp.bad() || fp.fail())
        throw MYRUNTIME_ERROR(preamb+"Failed to open file for writing: "+biname);

    size_t addrfields[pmv2DTotFlds];
    size_t addrdesc_end_addrs = 0;
    size_t addrdesc = 0;
    std::streamoff fppos = -1;

    memset( addrfields, 0, pmv2DTotFlds * sizeof(addrfields[0]));

    Open();

    MYMSG("Making binary database...",1);

    try {
        WriteBinaryHeader(fp);

        if((fppos = fp.tellp()) == (std::streamoff)-1 || fp.fail())
            throw MYRUNTIME_ERROR(preamb+"Failed to determine file position.");

        pmodel::CalculateAndWriteAddressTable_v2_2(
            fp,
            (size_t)fppos/*szpreamble*/,
            addrfields,
            addrdesc_end_addrs, addrdesc,
            no_sequences_/*nprofiles*/, db_size_/*totalposs*/);

        while( Next(pm,tm)) {
            BinaryWriteProfile_v2_2(
                fp,
                addrfields,
                addrdesc_end_addrs, addrdesc,
                pm,tm);
        }
    } catch( myexception const& ex ) {
        mre = ex;
    }

    Close();

    if( mre.isset())
        throw mre;

    MYMSG("Done.",1);
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
        nopenfiles_ = 0;

        ClearFPVector();
        fps_.resize(GetMaxNumberOfOpenFiles(), NULL);

        NewSortedFiles();

        ResetProbNormTerm();

        //clear marginal counts before updating them
        proprobs_.ClearMargs();

        //1st PASS: find out attributes of a database to be created
        ProcessInput( &Db::ProcessFile, NULL );

        //calculate probabilities after marginal counts have been updated
        proprobs_.CalcProbs();

        db_fp_[DbFlMAIN] = fopen( GetDbMainName(), "w" );

        if( db_fp_[DbFlMAIN] == NULL )
            throw MYRUNTIME_ERROR( preamb + "Failed to open file for writing.");

        message( "Writing..." );
        WriteTextHeader( db_fp_[DbFlMAIN] );

        //{{2nd PASS: write database attributes and profiles
        for( size_t n = 0; n < srtfiles_->GetSize(); n++ ) {
            const DbP1* p1 = (const DbP1*)srtfiles_->GetValueAt(n);
            if( p1 )
                WriteFile( *p1, db_fp_ );
        }
        //}}

        message( "Writing probabilities..." );
        proprobs_.WriteProbs( GetDbProbsName());

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error(ex);
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
// ProcessFile: read data from file and process it;
// sernr, serial number of the given file
//
void Db::ProcessFile( const char* filename, void* )
{
    MYMSG("Db::ProcessFile",3);

    pmodel::PMProfileModel pm;
    pmodel::PMTransModel tm;
    bool leaveopened = false;

    mystring fname = filename;
    if( fname.rfind(Db::db_extensions_[DbFlPPRO]) != mystring::npos ||
        fname.rfind(Db::db_extensions_[DbFlMAINBin]) != mystring::npos )
        return;

    if( srtfiles_ == NULL )
        throw MYRUNTIME_ERROR("Db::ProcessFile: Null file list.");

    FILE* fp = NULL;

    try {
        size_t nseqs = 0;
        uint64_mt dbsize = 0;
        long fpos = 0L;

        if(( fp = fopen( filename, "r" )) ==  NULL )
            throw MYRUNTIME_ERROR("Db::ProcessFile: Failed to open file.");

        leaveopened = (nopenfiles_ < GetMaxNumberOfOpenFiles());
    
        try {
            ReadTextHeader( fp, &nseqs, &dbsize );
        } catch(myexception const& )
        {
            nseqs = 0;
            dbsize = 0;
            if( fseek( fp, 0L, SEEK_SET ) != 0 )
                throw MYRUNTIME_ERROR("Db::ProcessFile: Failed to reset file position.");
        }

        for( ;!feof(fp); )
        {
            fpos = ftell(fp);
            TextReadProfile( fp, pm, tm );
            //PreprocessProfile( pm, tm );
            IncNoSequences();
            IncDbSize( pm.GetSize());

            //update marginal counts
            proprobs_.UpdateMargs( pm.GetEffNoSequences(), pm.GetSize());

            //sort filenames on the fly by profile length
            srtfiles_->Push( new DbP1(filename, pm.GetSize(), leaveopened? fp: NULL, fpos));

            int nsym = 0;
            nsym = fgetc(fp);
            ungetc( nsym, fp);
            if( nsym == EOF )
                break;
        }

        if( leaveopened ) {
            if( nopenfiles_ < (int)fps_.size())
                fps_[nopenfiles_] = fp;
            nopenfiles_++;
        }

    } catch( myexception const& ex )
    {
        leaveopened = false;
        error( ex.what(), false );
        error( mystring(mystring("Db::ProcessFile: ")+filename+mystring(": skipped.")).c_str());
    }

    if( fp && !leaveopened )
        fclose(fp);
}

// -------------------------------------------------------------------------
// WriteFile: append the file to the database
//
void Db::WriteFile( const DbP1& p1, void* args )
{
    MYMSG("Db::WriteFile",3);

    FILE** FP = ( FILE** )args;

    if( !FP || !FP[DbFlMAIN] )
        return;

    const mystring preamb = "Db::WriteFile: ";
    pmodel::PMProfileModel pm;
    pmodel::PMTransModel tm;
    mystring buffer;
    FILE* fp = NULL;
    bool openfile = false;
    int emsg;

    try {
        if( p1.fp_ ) {
            fp = p1.fp_;
            openfile = true;
        }
        else if(( fp = fopen( p1.fname_.c_str(), "r" )) ==  NULL )
            throw MYRUNTIME_ERROR( preamb + "Failed to open file.");

        if( fseek( fp, p1.fpos_, SEEK_SET ) != 0 )
            throw MYRUNTIME_ERROR( preamb + "Failed to set a file position.");

        if( PreprocessNeeded()) {
            TextReadProfile( fp, pm, tm );
            PreprocessProfile( pm, tm );
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
                if( strncmp(buffer.c_str(), patstrEND, lenstrEND) == 0 )
                    break;
            }
        }
    } catch( myexception const& ex )
    {
        error( ex.what(), false );
        error( mystring(preamb + "'" + p1.fname_ + mystring("' skipped.")).c_str());
    }

    if( fp && !openfile )
        fclose(fp);
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
    fp.write(pmodel::PMProfileModel::GetBinaryDataVersion(),
        strlen(pmodel::PMProfileModel::GetBinaryDataVersion()));
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

    fprintf( fp, "%s%s%s", patstrDBVER, pmodel::PMProfileModel::GetDataVersion(), NL );
    fprintf( fp, "%-9s%zu%s", patstrNoPROFS, GetNoSequences(), NL );
    fprintf( fp, "%-9s%llu%s", patstrSIZE, GetDbSize(), NL );
}

// -------------------------------------------------------------------------
// ReadTextHeader: read database header from file descriptor
// TODO: reimplement local variables to use heap: locbuffer
//
void Db::ReadTextHeader( FILE* fp, size_t* nseqs, uint64_mt* dbsize )
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
    long long int   llintval;


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

    if( nseqs == NULL )
        SetNoSequences( llintval );
    else
        *nseqs = (size_t)llintval;


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

    if( dbsize == NULL )
        SetDbSize( llintval );
    else
        *dbsize = (uint64_mt)llintval;
}
