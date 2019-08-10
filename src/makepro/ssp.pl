#!/usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Getopt::Long;


my  $MYPROGNAME = basename( $0 );
my  $instpsipred = '/share/install/psipred';
my  $instblastplus = '/share/install/ncbi-blast/ncbi-blast-2.2.23+';
my  $instcomer = '/home/mindaugas/projects/share/comer';

my  $usage = <<EOIN;

Predict secondary structure for the first sequence from a given multiple 
alignment (either in FASTA or in STOCKHOLM) by PSIPRED 
(Jones, D.T. (1999) JMB, 292, 195-202).
2013(C)Mindaugas Margelevicius,VU IBT,Vilnius

Usage:
$MYPROGNAME <Parameters>

Parameters:

--in  <alignment>  Filename of multiple sequence alignment
--out <output>     Output file of SS predictions
--p <options>      Options file from the `comer' package (Optional)
--select           Select representative sequences from multiple alignment
                   before predicting SS
--noflush          Do not flush temporary files
--psipred <path>   Installation path of `PSIPRED'
           default=$instpsipred
--blast <path>     Installation path of `BLAST+'
           default=$instblastplus
--comer <path>     Installation path of `comer'
           default=$instcomer
--help             This text

EOIN

my  $INFILE;
my  $OUTFILE;
my  $OPTIONS;
my  $SELECTS = 0;
my  $NOFLUSH = 0;
my  $BHDNLEN = 20; ## heading substring length under blast format
my  $BHDNSEC = 40; ## heading section length under blast format

my  $result = GetOptions(
               'in=s'      => \$INFILE,
               'out=s'     => \$OUTFILE,
               'p=s'       => \$OPTIONS,
               'psipred=s' => \$instpsipred,
               'blast=s'   => \$instblastplus,
               'comer=s'   => \$instcomer,
               'select'    => sub { $SELECTS = 1; },
               'noflush'   => sub { $NOFLUSH = 1; },
               'help|h'    => sub { print $usage; exit( 0 ); }
);


do { print $usage; exit( 1 ); }  unless $result;

die "ERROR: No input filename given." unless $INFILE;
die "ERROR: No output filename given." unless $OUTFILE;
die "ERROR: Input alingment file does not exist or is invalid." unless( $INFILE && -f $INFILE );

die "ERROR: PSIPRED directory does not exist." unless( $instpsipred && -d $instpsipred );
die "ERROR: BLAST+ directory does not exist." unless( $instblastplus && -d $instblastplus );
die "ERROR: `comer' directory does not exist." unless(( $instcomer && -d $instcomer )|| !$SELECTS );

my  $select = "$instcomer/bin/select";

my  $makeblastdb = "$instblastplus/bin/makeblastdb";
my  $psiblast = "$instblastplus/bin/psiblast";

my  $chkparse = "$instpsipred/bin/chkparse";
my  $psipred = "$instpsipred/bin/psipred";
my  $psipass2 = "$instpsipred/bin/psipass2";

my  $ppweights = "$instpsipred/data/weights.dat";
my  $ppweights2 = "$instpsipred/data/weights.dat2";
my  $ppweights3 = "$instpsipred/data/weights.dat3";
my  $ppweights_p2 = "$instpsipred/data/weights_p2.dat";

printf( STDERR "\nWARNING: Assuming PSIPRED version >=3.0\n\n");

die "ERROR: `comer' executable `select' does not exist." unless(( -e $select )|| !$SELECTS );
die "ERROR: Options file does not exist." unless(( !$OPTIONS || -e $OPTIONS )|| !$SELECTS );

die "ERROR: BLAST+ executable `makeblastdb' does not exist." unless( -e $makeblastdb );
die "ERROR: BLAST+ executable `psiblast' does not exist." unless( -e $psiblast );

die "ERROR: PSIPRED executable `chkparse' does not exist." unless( -e $chkparse );
die "ERROR: PSIPRED executable `psipred' does not exist." unless( -e $psipred );
die "ERROR: PSIPRED executable `psipass2' does not exist." unless( -e $psipass2 );

die "ERROR: PSIPRED weight file does not exist." unless( -f $ppweights );
die "ERROR: PSIPRED weight file (2) does not exist." unless( -f $ppweights2 );
die "ERROR: PSIPRED weight file (3) does not exist." unless( -f $ppweights3 );
die "ERROR: PSIPRED weight file (p2) does not exist." unless( -f $ppweights_p2 );


my  $myinput = $INFILE;
my  $inbase = basename( $INFILE);
my  $outdir = dirname( $OUTFILE );
my  $infileseld = "$outdir/$inbase.sel";
my  $infileblast = "$outdir/$inbase.blast";
my  $infileblastdb = "$outdir/$inbase.blast.db";
my  $infileblastquery = "$outdir/$inbase.blast.query";
my  $infileblastchk = "$outdir/$inbase.blast.chk";
my  $infileblastout = "$outdir/$inbase.blast.out";
my  $infileppmtx = "$outdir/$inbase.pp.mtx";
my  $infileppss = "$outdir/$inbase.pp.ss";
my  $infileppss2 = "$outdir/$inbase.pp.ss2";
my  $infilepphor = "$outdir/$inbase.pp.hor";


if( $SELECTS ) {
    printf( STDERR "Selecting sequences...\n");
    exit(1) unless SelectSequences( $myinput, $OPTIONS, $infileseld );
    $myinput = $infileseld;
}

printf( STDERR "Translating input...\n");
if( IsFasta( $myinput )) {
    exit(1) unless TranslateFasta( $myinput, $infileblast, $BHDNLEN, $BHDNSEC );
}
elsif( IsStockholm( $myinput )) {
    exit(1) unless TranslateStockholm( $myinput, $infileblast, $BHDNLEN, $BHDNSEC );
}
else {
    printf( STDERR 
      "ERROR: Other input formats than FASTA and STOCKHOLM are not supported.\n");
}

printf( STDERR "Making BLAST database...\n");
exit(1) unless MakeBlastDb( $infileblast, $infileblastdb, $infileblastquery );

printf( STDERR "\nMaking BLAST profile...\n");
exit(1) unless MakeBlastProfile( $infileblast, $infileblastdb, $infileblastquery, 
            $infileblastchk, $infileblastout );

printf( STDERR "\nPredicting SS with PSIPRED...\n");
exit(1) unless PredictSS( $infileblastchk, $infileppmtx, $infileppss, $infileppss2, $infilepphor );

exit(1) unless Output( $infilepphor, $OUTFILE );

unless( $NOFLUSH ) {
    unlink <$infileblast*>;
    unlink( $infileppmtx, $infileppss, $infileppss2 );
    unlink( $infilepphor ) if $infilepphor ne $OUTFILE;
    unlink( $infileseld ) if $SELECTS;
}

printf( STDERR "\nDone.\n\n");
exit(0);

## -------------------------------------------------------------------
## Select non-redundant set of sequences from multiple alignment
##
sub SelectSequences
{
    my  $input = shift;
    my  $options = shift;
    my  $output = shift;
    my  $command;
    my  $ret = 1;

    $command = "$select -i $input -o $output ";
    $command.= "-p $options" if $options;
    $ret = 0 unless RunCommand( $command );
    return $ret;
}

## -------------------------------------------------------------------
## Check if input is in FASTA
##
sub IsFasta
{
    my  $input = shift;
    my  $ret = 0;
    return 0 unless( open( FF, $input ));
    while(<FF>) {
        next if /^#/;
        $ret = 1 if /^>/;
        last;
    }
    close( FF );
    return $ret;
}

## -------------------------------------------------------------------
## Check if input is in STOCKHOLM
##
sub IsStockholm
{
    my  $input = shift;
    my  $ret = 0;
    return 0 unless( open( FF, $input ));
    $_ = <FF>;
    $ret = 1 if /^# STOCKHOLM/;
    close( FF );
    return $ret;
}

## -------------------------------------------------------------------
## Translate alignments from FASTA to BLAST format
##
sub TranslateFasta
{
    my  $input = shift;
    my  $output = shift;
    my  $bhdnlen = shift;
    my  $bhdnsec = shift;
    my ($hdn,$seqn) = ('','');
    my ($len, $hdnlen, $counter) = (0,0,0);
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open input file: %s.\n", $input );
        return 0;
    }
    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        close(FF);
        return 0;
    }
    while(<FF>) {
        next if /^#/;
        chomp;
        if(/^>(\S*)/ || eof(FF)) {
            $seqn .= $_ if eof(FF);
            if( $hdn && $seqn ) {
                if( $len ) {
                    if( length($seqn) != $len ) {
                        printf( STDERR "ERROR: Invalid sequence length: %s.\n", $hdn );
                        $ret = 0;
                        last;
                    }
                }
                else {
                    $len = length($seqn);
                }
                printf( OF "%s %s\n", $hdn, $seqn );
            }
            $seqn = '';
            $hdn = substr( $1, 0, $bhdnlen );
            $hdn .= sprintf("_%d", $counter++ );
            $hdnlen = length( $hdn );
            if( $hdnlen < $bhdnsec ) {
                $hdn .= ' ' x($bhdnsec-$hdnlen);
            }
            next;
        }
        ##sequence
        $seqn .= $_;
    }
    close( OF );
    close( FF );
    return $ret;
}

## -------------------------------------------------------------------
## Translate alignments from STOCKHOLM to BLAST format
##
sub TranslateStockholm
{
    my  $input = shift;
    my  $output = shift;
    my  $bhdnlen = shift;
    my  $bhdnsec = shift;
    my  %aln;
    my ($hdn,$seqn) = ('','');
    my ($len, $hdnlen, $counter, $n) = (0,0,0,0);
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open input file: %s.\n", $input );
        return 0;
    }
    while(<FF>) {
        next if /^#/;
        chomp;
        if(/^(\S+)\s+(\S+)$/) {
            $hdn = $1;
            $seqn = $2;
            $seqn =~ s/[._~]/-/g;
            $seqn = uc($seqn);

            $n++;
            $aln{$hdn}{N} = $n unless exists $aln{$hdn}{N};
            $aln{$hdn}{S} .= $seqn;
        }
    }
    close( FF );
    return $ret unless $ret;

    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        return 0;
    }
    foreach $hdn( sort {$aln{$a}{N} <=> $aln{$b}{N}} keys %aln ) {
        $seqn = $aln{$hdn}{S};
        unless( $hdn && $seqn ) {
            printf( STDERR "ERROR: Invalid sequence: %s.\n", $hdn );
            $ret = 0;
            last;
        }
        if( $len ) {
            if( length($seqn) != $len ) {
                printf( STDERR "ERROR: Invalid sequence length: %s.\n", $hdn );
                $ret = 0;
                last;
            }
        }
        else {
            $len = length($seqn);
        }
        $hdn = substr( $hdn, 0, $bhdnlen );
        $hdn .= sprintf("_%d", $counter++ );
        $hdnlen = length( $hdn );
        if( $hdnlen < $bhdnsec ) {
            $hdn .= ' ' x($bhdnsec-$hdnlen);
        }
        printf( OF "%s %s\n", $hdn, $seqn );
    }
    close( OF );
    return $ret;
}

## -------------------------------------------------------------------
## Make a small blast database
##
sub MakeBlastDb
{
    my  $input = shift;
    my  $output = shift;
    my  $output2 = shift;
    my ($command, $hdn, $seqn );
    my  $ret = 1;
    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open formatted file: %s.\n", $input );
        return 0;
    }
    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        close(FF);
        return 0;
    }
    unless( open( OF2, ">$output2" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        close(OF);
        close(FF);
        return 0;
    }
    $_ = <FF>;
    if(/^(\S+)\s+(\S+)$/) {
        $hdn = $1;
        $seqn = $2;
        $seqn =~ s/-//g;
        printf( OF ">%s\n%s\n", $hdn, $seqn );
        printf( OF2 ">%s\n%s\n", $hdn, $seqn );
    }
    close( OF2 );
    close( OF );
    close( FF );
    $command = "$makeblastdb -in $output -dbtype prot";
    $ret = 0 unless RunCommand( $command );
    return $ret;
}

## -------------------------------------------------------------------
## Make BLAST profile from the given multiple alignment
##
sub MakeBlastProfile
{
    my  $inmsa = shift;
    my  $blastdb = shift;
    my  $query = shift; 
    my  $outputchk = shift;
    my  $outputblast = shift;
    my  $command;
    my  $ret = 1;

    $command = "$psiblast -db $blastdb -in_msa $inmsa -out_pssm $outputchk ".
               "-inclusion_ethresh 0.001 -num_iterations 1 -num_alignments 10 >$outputblast 2>&1";
    $ret = 0 unless RunCommand( $command );
    return $ret;
}

## -------------------------------------------------------------------
## Predict SS by using PSIPRED
##
sub PredictSS
{
    my  $chkfile = shift;##input; below -- outputs
    my  $mtxfile = shift;
    my  $ssfile = shift;
    my  $ss2file = shift;
    my  $horfile = shift;
    my  $tmpchk = "${chkfile}_tmp";
    my  $ncbicodes = 'XAXCDEFGHIKLMNPQRSTVWXYXXX';
    my  @codes = split(//, $ncbicodes );
    my  %hcds; @hcds{@codes} = 0..$#codes;
    my ($command, $n, $mm, $ln, $c, $hx );
    my  $ret = 1;

    $command = "$chkparse $chkfile >$mtxfile";
    unless( RunCommand( $command, 1 ) == 1 ) {
        ## modify chk file so that PSIPRED finds a required token
        printf( STDERR "Modifying BLAST chk file and rerunning PSIPRED...\n");
        unless( open( CF, $chkfile )) {
            printf( STDERR "ERROR: Failed to open BLAST chk file: %s.\n", $chkfile );
            return 0;
        }
        unless( open( OF, ">$tmpchk" )) {
            printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $tmpchk );
            close(CF);
            return 0;
        }
        while(<CF>) {
            if(/iupacaa/ || $mm ) {
                $n = 0;
                $mm = 1;
                $ln = $_;
                if( $ln =~ s/iupacaa/ncbistdaa/) {
                    for( $n = $+[0]; $n < length($ln); $n++ ) {
                        $c = substr($ln, $n, 1);
                        last if( $c eq '"' || $c eq "'");
                    }
                    substr($ln, $n++, 1) = "'" if( $n < length($ln));
                }
                for( ; $n < length($ln); $n++ ) {
                    $c = substr($ln, $n, 1);
                    next if $c =~ /\s/;
                    if( $c eq '"' || $c eq "'") {
                        substr($ln, $n++, 1) = "'";
                        $mm = 0;
                        last;
                    }
                    $hx = 0;
                    $hx = $hcds{$c} if exists $hcds{$c};
                    substr($ln, $n++, 1) = sprintf("%02X",$hx);
                }
                $_ = $ln;
            }
            printf( OF );
        }
        close(OF);
        close(CF);
        return 0 unless Move( $tmpchk, $chkfile );
        return 0 unless RunCommand( $command );
    }

    $command = "$psipred $mtxfile $ppweights $ppweights2 $ppweights3 >$ssfile";
    return 0 unless RunCommand( $command );

    $command = "$psipass2 $ppweights_p2 1 1.0 1.0 $ss2file $ssfile >$horfile";
    return 0 unless RunCommand( $command );

    return $ret;
}

## -------------------------------------------------------------------
## Output predictions
##
sub Output
{
    my  $horfile = shift;
    my  $output = shift;
    my (%pred, $len, $n );
    my  $wrap = 80;
    my  $ret = 1;

    unless( open( FF, $horfile )) {
        printf( STDERR "ERROR: Failed to open PSIPRED predictions file: %s.\n", $horfile );
        return 0;
    }
    while(<FF>) {
        next if /^#/;
        chomp;
        $pred{CONF} .= $1 if(/^Conf:\s+(\S+)$/);
        $pred{PRED} .= $1 if(/^Pred:\s+(\S+)$/);
        $pred{SEQN} .= $1 if(/^\s*AA:\s+(\S+)$/);
    }
    close(FF);
    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        return 0;
    }
    $len = length($pred{SEQN});
    ##previously...
    ##for( $n = 0; $n < $len; $n+=$wrap ) {
    ##    printf( OF "Conf: %s\n", substr( $pred{CONF}, $n, $wrap ));
    ##    printf( OF "SS:   %s\n", substr( $pred{PRED}, $n, $wrap ));
    ##    printf( OF "Seqn: %s\n", substr( $pred{SEQN}, $n, $wrap ));
    ##    printf( OF "\n");
    ##}
    printf( OF "## Sequence SS Probability\n##\n");
    for( $n = 0; $n < $len; $n++ ) {
        printf( OF "%3s %3s %5.1f\n", substr($pred{SEQN},$n,1), substr($pred{PRED},$n,1), (substr($pred{CONF},$n,1)/10));
    }
    close(OF);
    return $ret;
}

## -------------------------------------------------------------------
## -------------------------------------------------------------------
## Move file
##
sub Move
{
    my  $what = shift;
    my  $to = shift;
    unless( open( IF, $what )) {
        printf( STDERR "ERROR: Failed to open file to move: %s.\n", $what );
        return 0;
    }
    unless( open( OF, ">$to" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $to );
        close(IF);
        return 0;
    }
    printf( OF ) while(<IF>);
    close(OF);
    close(IF);
    unlink($what);
    return 1;
}

## -------------------------------------------------------------------
## System command
##
sub CheckStatus
{
    RunCommand();
}
sub RunCommand
{
    my  $cmdline = shift;
    my  $retstatus = shift;

##    print( "$cmdline\n" ) if $cmdline;
    system( "$cmdline" ) if $cmdline;

    if( $? == -1 ) {
        printf( STDERR "ERROR: Failed to execute command: $!\n" );
        return 0;
    }
    if( $? & 127 ) {
        printf( STDERR "%s: Command ran died with signal %d, %s coredump.\n",
                       $retstatus? 'WARNING': 'ERROR',
                       ( $? & 127 ), ( $? & 128 )? 'with' : 'without' );
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

