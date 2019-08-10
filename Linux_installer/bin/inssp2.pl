#!/usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Getopt::Long;


my  $MYPROGNAME = basename( $0 );

my  $usage = <<EOIN;

Append secondary structure prediction to `comer' profile. 
(C) 2013-2019 Mindaugas Margelevicius,VU IBT,Vilnius

Usage:
$MYPROGNAME <Parameters>

Parameters:

--in <ssfile>      Filename of secondary structure prediction
--to <profile>     Filename of `comer' profile
--help             This text

EOIN

my  $SSPFILE;
my  $PROFILE;

my  $result = GetOptions(
               'in=s'      => \$SSPFILE,
               'to=s'      => \$PROFILE,
               'help|h'    => sub { print $usage; exit( 0 ); }
);


do { print $usage; exit( 1 ); }  unless $result;

die "ERROR: No input filename given." unless $SSPFILE;
die "ERROR: No output filename given." unless $PROFILE;
die "ERROR: Input file of SS prediction does not exist or is invalid." unless( $SSPFILE && -f $SSPFILE );
die "ERROR: `comer' profile does not exist or is invalid." unless( $PROFILE && -f $PROFILE );


my  $tmpfile = "$PROFILE.tmp";
my (%SSP, %PRO );

exit(1) unless ReadSSP( $SSPFILE, \%SSP );
exit(1) unless ReadProfile( $PROFILE, \%PRO );
exit(1) unless Append( \%SSP, \%PRO );
exit(1) unless Write( \%PRO, $PROFILE );
exit(0);

## -------------------------------------------------------------------
## Read SS prediction from file
##
sub ReadSSP
{
    my  $input = shift;
    my  $rssp = shift;
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open input file: %s.\n", $input );
        return 0;
    }
    while(<FF>) {
        next if /^#/;
        if(/\s*([A-Za-z\*])\s+([A-Za-z])\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)$/) {
            $$rssp{SEQN} .= $1;
            $$rssp{SSP} .= $2;
            push( @{$$rssp{C}}, $3 );
            push( @{$$rssp{E}}, $4 );
            push( @{$$rssp{H}}, $5 );
        }
    }
    close(FF);
    return $ret;
}

## -------------------------------------------------------------------
## Read `comer' profile
##
sub ReadProfile
{
    my  $input = shift;
    my  $rpro = shift;
    my ($r, $n, $m, $t );
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open `comer' profile: %s.\n", $input );
        return 0;
    }
    $n = -1; $m = -1;
    while(<FF>) {
        if(( /^\s*K,\s+Lambda/ || /^\*$/ || eof(FF)) && !$t ) {
            $t = 1;
            ##last line exchange:
#             $$rpro{TAIL} .= $$rpro{RECD}[$n][$m] if( 0 <= $n && 0 <= $m );
#             splice( @{$$rpro{RECD}[$n]}, $m, 1 );
        }
        if(/^\s*\d+\s+([a-zA-Z])\s+-?\d+/) {
            $r = $1;
            $$rpro{SEQN} .= $r;
            $m = -1;
            $n++;
        }
        $$rpro{TAIL} .= $_ if $t;
        $$rpro{RECD}[$n][++$m] = $_ if( $r && !$t );
        $$rpro{HEAD} .= $_ unless( $r || $t );
    }
    close(FF);
    return $ret;
}

## -------------------------------------------------------------------
## Rewrite `comer' profile
##
sub Append
{
    my  $rssp = shift;
    my  $rpro = shift;
    my ($n, $r, $c, $e, $h, $p, $rr, $m );
    my  $ret = 1;
    unless( exists $$rssp{SEQN} && exists $$rssp{SSP} &&
        exists $$rssp{C} && exists $$rssp{E} && exists $$rssp{H}) {
        printf( STDERR "ERROR: No SS prediction.\n");
        return 0;
    }
    unless( length($$rssp{SEQN}) == $#{$$rssp{C}}+1 && 
            $#{$$rssp{C}} == $#{$$rssp{E}} && $#{$$rssp{C}} == $#{$$rssp{H}} &&
            length($$rssp{SEQN}) == length($$rssp{SSP})) {
        printf( STDERR "ERROR: Inconsistent SS prediction.\n");
        return 0;
    }
    unless( exists $$rpro{RECD}) {
        printf( STDERR "ERROR: No profile.\n");
        return 0;
    }
    unless( length($$rpro{SEQN}) == $#{$$rpro{RECD}}+1 ) {
        printf( STDERR "ERROR: Inconsistent profile.\n");
        return 0;
    }
    unless( length($$rssp{SEQN}) == $#{$$rpro{RECD}}+1 ) {
        printf( STDERR "ERROR: SS prediction inconsistent with profile.\n");
        return 0;
    }
    for( $n = 0; $n < length($$rssp{SEQN}); $n++ ) {
        $r = substr( $$rssp{SEQN}, $n, 1 );
        $c = $$rssp{C}[$n];
        $e = $$rssp{E}[$n];
        $h = $$rssp{H}[$n];
        $p = substr( $$rssp{SSP}, $n, 1 );
        $rr = substr( $$rpro{SEQN}, $n, 1 );
        if( $r ne '*' && uc($r) ne uc($rr)) {
            printf( STDERR "ERROR: Inconsistent sequences from SS prediction and profile.\n");
            return 0;
        }
        $m = $#{$$rpro{RECD}[$n]};
        $m++ unless( $$rpro{RECD}[$n][$m] =~ /^\s+SS:/);
        $$rpro{RECD}[$n][$m] = sprintf(" SS:%s %d %d %d\n",
              $p, int($c*10000), int($e*10000), int($h*10000));
    }
    return $ret;
}

## -------------------------------------------------------------------
## Rewrite `comer' profile
##
sub Write
{
    my  $rpro = shift;
    my  $output = shift;
    my ($rc, $rcm );
    my  $ret = 1;

    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s.\n", $output );
        return 0;
    }
    print( OF $$rpro{HEAD});
    foreach $rc( @{$$rpro{RECD}}) {
        foreach $rcm( @{$rc}) {
            print( OF $rcm );
        }
    }
    print( OF $$rpro{TAIL});
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

