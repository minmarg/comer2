#!/usr/bin/env perl
use warnings;

##
## (C)2019-2020 Mindaugas Margelevicius
## Institute of Biotechnology, Vilnius University
##

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Spec;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(:config no_ignore_case);

my  $platform = $^O;
my  $mswin = ($platform =~ /(?:MSWin:Win32:Win64:WinXP)/i)? 1: 0;

my  $MYPROGNAME = basename($0);
my  $MYPROGDIRNAME = dirname($0);
my  $VERSION = '2.01';

my  $NGPUsVal = 1;
my  $NGPUs = "--dev-N=$NGPUsVal";
my  $DMEM = '';
my  $DEXPLVal = 50;
my  $DEXPL = "--dev-expected-length=$DEXPLVal";
my  $DP2MPVal = 10;
my  $DP2MP = "--dev-pass2memp=$DP2MPVal";
my  $NIOBUFVal = 4;
my  $NIOBUF = "--io-nbuffers=$NIOBUFVal";
my  $DIOFMP = '';
my  $DIOUNP = '';

my  $PSIPREDDIR = '/share/install/psipred';
my  $BLASTDIR = '/share/install/ncbi-blast/ncbi-blast-2.2.23+';
my  $NOPRED;

my  $usage = <<EOIN;

$MYPROGNAME $VERSION

Helper script to construct a COMER profile and search a profile 
database with the constructed profile.
(C)2019-2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$MYPROGNAME <Options>

Options:
-i <input>                  Input multiple alignment file either in 
                            FASTA or STOCKHOLM.
-d <database_name>          Database of COMER profiles (without extension).
-o <output_directory>       Directory of output files of the search.
-r <output_profile>         Filename of output profile. If not provided, the profile
                            will be saved in the directory specified by option -o.
-p <options_file>           Input file of options;
                            By default, the options file in the installation
                            directory is used.

COMER options (type comer -h for their complete description):
--dev-N=(<number>|,<id_list>) Maximum number of GPUs to use. This can be
                            specified by a number or given by a comma-separated
                            list of GPU identifiers, which should start with a comma.
                        Default=$NGPUsVal (most powerful GPU)
--dev-mem=<megabytes>       Maximum amount of GPU memory (MB) that can be used.
                        Default=[all memory of GPU(s)]
--dev-expected-length=<length> Expected length of database proteins (profiles). Its
                            values are restricted to the interval [20,200].
                        Default=$DEXPLVal
--dev-pass2memp=<percentage> GPU memory proportion dedicated to recalculation of
                            hits that passed significance threshold.
                            (Expected proportion of significant hits.)
                        Default=$DP2MPVal
--io-nbuffers=<count>       Number of buffers [1,10] used to cache data read from file.
                        Default=$NIOBUFVal
--io-filemap                Map files into memory.
--io-unpinned               Do not use pinned (page-locked) CPU memory.

Path options:
-P <path>                   Installation path of `PSIPRED'
       default=$PSIPREDDIR
-B <path>                   Installation path of `BLAST+'
       default=$BLASTDIR
--nopred                    Do not predict and include secondary structure in the 
                            profile, and ignore options -P and -B.

Other options:
-v                          Verbose mode.
-h                          This text.

EOIN

my  $INPUT;
my  $DB;
my  $OUTDIR;
my  $PROFILE;
my  $OPTIONS = '';
my  $VERBOSE = '';
my  $nargs = $#ARGV;

my  $retcode = GetOptions(
    'i=s'      => \$INPUT,
    'd=s'      => \$DB,
    'o=s'      => \$OUTDIR,
    'r=s'      => \$PROFILE,
    'p=s'      => \$OPTIONS,

    'dev-N=s'  => sub{my ($name,$value) = @_; $NGPUs = "--dev-N=$value";},
    'dev-mem=s' => sub{my ($name,$value) = @_; $DMEM = "--dev-mem=$value";},
    'dev-expected-length=s' => sub{my ($name,$value) = @_; $DEXPL = "--dev-expected-length=$value";},
    'dev-pass2memp=s' => sub{my ($name,$value) = @_; $DP2MP = "--dev-pass2memp=$value";},
    'io-nbuffers=s' => sub{my ($name,$value) = @_; $NIOBUF = "--io-nbuffers=$value";},
    'io-filemap' => sub {$DIOFMP = '--io-filemap';},
    'io-unpinned' => sub {$DIOUNP = '--io-unpinned';},

    'P=s'      => \$PSIPREDDIR,
    'B=s'      => \$BLASTDIR,
    'nopred'   => sub {$NOPRED = 1;},

    'v'        => sub {$VERBOSE = '-v';},
    'help|h'   => sub {print $usage; exit(0);}
);

exit(1) unless $retcode;
do { print $usage; exit(0); } if $nargs < 0;
do { print STDERR "ERROR: Input not specified.\n"; exit(1); } unless($INPUT);
do { print STDERR "ERROR: Database not specified.\n"; exit(1); } unless($DB);
do { print STDERR "ERROR: Output directory not specified.\n"; exit(1); } unless($OUTDIR);
do { print STDERR "ERROR: Input not found: $INPUT\n"; exit(1); } unless(-f $INPUT);
do { print STDERR "ERROR: Options file not found: $OPTIONS\n"; exit(1); } unless(!length($OPTIONS) || -f $OPTIONS);
do { print STDERR "ERROR: PSIPRED directory not found: $PSIPREDDIR\n"; exit(1); } unless($NOPRED || -d $PSIPREDDIR);
do { print STDERR "ERROR: BLAST directory not found: $BLASTDIR\n"; exit(1); } unless($NOPRED || -d $BLASTDIR);

unless($PROFILE) {
    $PROFILE = File::Spec->catfile($OUTDIR, basename($INPUT).'.pro');
}

$OPTIONS = "-p $OPTIONS" if length($OPTIONS);

my $makepro = File::Spec->catfile($MYPROGDIRNAME, 'makepro');
my $comer = File::Spec->catfile($MYPROGDIRNAME, 'comer');
my $command;

if($mswin) {
    $makepro .= ($NOPRED? '.exe': '.cmd');
    $comer .= '.exe';
} else {
    $makepro .= '.sh' unless($NOPRED);
}

do { print STDERR "ERROR: Executable not found: $makepro\n"; exit(1); } unless(-f $makepro);
do { print STDERR "ERROR: Executable not found: $comer\n"; exit(1); } unless(-f $comer);
unless(-d $OUTDIR) {
    do { print STDERR "ERROR: Failed to create directory: $OUTDIR\n"; exit(1); } unless( make_path($OUTDIR));
}

$command = "$makepro $VERBOSE -i $INPUT -o $PROFILE $OPTIONS";
$command .= " -P $PSIPREDDIR -B $BLASTDIR" unless($NOPRED);
exit(1) unless RunCommandV($command);

$command = "$comer $VERBOSE -i $PROFILE -d $DB -o $OUTDIR $OPTIONS";
$command .= " $NGPUs $DMEM $DEXPL $DP2MP $NIOBUF $DIOFMP $DIOUNP";
exit(1) unless RunCommandV($command);


## -------------------------------------------------------------------
## run system command
##
sub CheckStatus
{   
    return RunCommand();
}

sub RunCommandV
{   
    my  $cmdline = shift;
    my  $retstatus = shift;
    my  $routput = shift;##ref
    print( STDERR "CMD: $cmdline\n") if $cmdline;
    return RunCommand($cmdline, $retstatus, $routput);
}

sub RunCommand
{
    my  $cmdline = shift;
    my  $retstatus = shift;
    my  $routput = shift;##ref

    if($cmdline) {
        $$routput = `$cmdline 2>&1` if $routput;
        system( "$cmdline" ) unless $routput;
    }

    if( $? == -1 ) {
        printf( STDERR "ERROR: Failed to execute command: $!\n" );
        return 0;
    }
    if( $? & 127 ) {
        printf( STDERR "ERROR: Command terminated with signal %d (%s coredump).\n",
            ($? & 127), ($? & 128)? 'with' : 'without' );
        return 0;
    }
    else {
        if(( $? >> 8 ) != 0 ) {
            unless( $retstatus ) {
                printf( STDERR "ERROR: Command failed and exited with status %d\n", $? >> 8 );
                return 0
            }
            return( $? >> 8 );
        }
    }
    return 1;
}

##<<>>

