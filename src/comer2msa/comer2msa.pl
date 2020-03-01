#!/usr/bin/perl -w

##
## (C)2008-2019 Mindaugas Margelevicius
## Institute of Biotechnology, Vilnius University
##

use strict;
use File::Basename;
use Getopt::Long;

my  $MYPROGNAME = basename( $0 );
my  $X = 'X';
my  $OUTFMT = "FASTA";
my  $OUTFMTCODE = 0;
my  $EVALUE = 1e-4;
my  $SBJCTs = '';
my  $WIDTH = -1;
my  $TOLXs = 0;
my  $PWALN = 0;
my  $ALTALN = 1;
my  $usage = <<EOIN;

Convert COMER output of pairwise alignments to multiple 
sequence alignment in FASTA.
(C)2008-2019 Mindaugas Margelevicius, Vilnius University

Usage:
$MYPROGNAME <Parameters>

Parameters:

-i <filename>   COMER results file of pairwise alignments.

-o <filename>   Name of output file of multiple sequence alignment.

-f [0|1]        Output format: 0, FASTA; 1, Rosetta.
        default=$OUTFMTCODE

-q <filename>   File containing query sequence in FASTA to be
                added to multiple alignment. If not specified, the
                largest aligned fraction of the query sequence will be
                added instead.
        optional

-e <e-value>    Read alignments with an expectation value less or equal
                than the specified value.
        default=$EVALUE

-l <ID_list>    Comma-separated list of the names of subjects to be 
                considered.
        By default, all subjects with an E-value less than the threshold
                are considered.

-p              Produce a number of pairwise alignments between the 
                query and subjects instead of one MSA. In this case, 
                output filenames will be constructed by inserting the 
                suffix corresponding to the subject name.

-t              Tolerate masked positions (`X' symbols) in the query
                sequence.

-a              Do not include COMER alternative alignments between the
                query and the same subjects (not in use currently).

-w <width>      Number of characters to wrap aligned sequences at.
        By default, no wrapping is used.

-h              Short description.

EOIN


my  $INPUT;
my  $OUTPUT;
my  $QUERY;
my  $Fail = 0;
my  $processed = 0;

my  $result = GetOptions(
               'i=s'      => \$INPUT,
               'o=s'      => \$OUTPUT,
               'f=i'      => \$OUTFMTCODE,
               'q=s'      => \$QUERY,
               'e=f'      => \$EVALUE,
               'l=s'      => \$SBJCTs,
               'w=i'      => \$WIDTH,
               'p'        => sub { $PWALN = 1; },
               't'        => sub { $TOLXs = 1; },
               'a'        => sub { $ALTALN = 0; },
               'help|h'   => sub { print $usage; exit( 0 ); }
);


do { print $usage; $Fail = 1; }  unless $result;
do { print STDERR "ERROR: Arguments missing.\n$usage"; $Fail = 1; } unless( $Fail || ( $INPUT && $OUTPUT ));
do { print STDERR "ERROR: File $INPUT does not exist.\n"; $Fail = 1; } unless( $Fail || -f $INPUT );
do { print STDERR "ERROR: File $QUERY does not exist.\n"; $Fail = 1; } unless( $Fail ||( !$QUERY || -f $QUERY ));

my  %HSBJCTs;

$HSBJCTs{$_} = 1 foreach split(',',$SBJCTs);

if( $OUTFMTCODE == 1 ) {
    $OUTFMT = "Rosetta";
}

## ===================================================================

if( $Fail ) {
    exit( 1 );
}

unless( AlignTopHits( $INPUT, $OUTPUT, $OUTFMT, $QUERY, $EVALUE, \%HSBJCTs, $TOLXs, $WIDTH, $PWALN, \$processed )) {
    printf( STDERR "Failed.\n" );
    exit( 1 );
}

printf( STDERR "\nTotal sequences, %d\nDone.\n", $processed );
exit( 0 );

## ===================================================================
## align top hits from blast output
##

sub AlignTopHits
{
    my  $input = shift;
    my  $output = shift;
    my  $format = shift;
    my  $queryfile = shift;
    my  $evalue = shift;
    my  $rhSBJCTs = shift;
    my  $tolxs = shift;
    my  $width = shift;
    my  $pwaln = shift;
    my  $rproc = shift;
    my  %hithash;
    my  %descrip;
    my  %queryin;
    my  $alnmask;
    my  $query = 'theonly';
    my ($rec, $len ) = ( 1, 0 ); ## the first record is reserved

    if( $queryfile ) {
        return 0 unless ReadQueryFile( $queryfile, \%queryin );
        $len = length( $queryin{SEQN} );
        $hithash{$query}[0][0] = '';
    }
    return 0 unless ExtractAlignments( $input, $evalue, $rhSBJCTs, $query, \%hithash, $rec, \%descrip );
    if($pwaln) {
        for(my $r = $rec; $r <= $#{$hithash{$query}}; $r++) {
            my %lochash;
            my $locext = ($output =~ /(\.[^\.]+)$/)? $1: '';
            my $locoutp = ($output =~ /^(.+)$locext$/)? $1: '';
            my $locsbjct = ($hithash{$query}[$r][0] =~ /^([\w\.]+)/)? $1: '';
            $locoutp .= "_${locsbjct}".$locext;
            ##
            @{$lochash{$query}[$rec]} = @{$hithash{$query}[$r]};
            return 0 unless MakeMultipleAlignment( \%lochash, $rec, $tolxs );
            return 0 unless MakeTargetOriented( \%lochash, $rec, \$alnmask, $len, $tolxs );
            return 0 unless InsertQuerySequence( \%lochash, $rec, \%descrip, \%queryin, $alnmask );
            return 0 unless PrintAligned( $format, \%lochash, $locoutp, $width, $rproc );
        }
    }
    else {
        return 0 unless MakeMultipleAlignment( \%hithash, $rec, $tolxs );
        return 0 unless MakeTargetOriented( \%hithash, $rec, \$alnmask, $len, $tolxs );
        return 0 unless InsertQuerySequence( \%hithash, $rec, \%descrip, \%queryin, $alnmask );
        return 0 unless PrintAligned( $format, \%hithash, $output, $width, $rproc );
    }
    return 1;
}

## -------------------------------------------------------------------
## read query description and sequence data
##

sub ReadQueryFile
{
    my  $filename = shift;
    my  $rqueryin = shift; ## ref. to hash

    unless( open( FA, $filename )) {
        printf( STDERR "ERROR: Failed to open $filename: $!\n" );
        return 0;
    }
    while( <FA> ) {
        chomp;
        if( /^>(.+)\s*$/ ) {
            $$rqueryin{DESC} = $1;
            $$rqueryin{SEQN} = '';
            next;
        }
        s/\s//g;
        $$rqueryin{SEQN} .= uc( "$_" );
    }
    close( FA );
    return 1;
}

## -------------------------------------------------------------------
## extract top alignments from blast output
##

sub ExtractAlignments
{
    my  $input = shift;
    my  $evalue = shift;
    my  $rhSBJCTs = shift;
    my  $query = shift;
    my  $refhash = shift; ## reference to hash of hits
    my  $stindex = shift; ## start index
    my  $refdesc = shift; ## reference to hash of query descriptions

    my  $hitnum = 0;
    my  $rec = $stindex;
    my  $last;
    my  $e_val;
    my  $score;

    my ($sbjct, $sbjctname, $titletmp, $skip );

    my  $queryfasta;
    my  $sbjctfasta;

    my  $querystart;
    my  $sbjctstart;
    my  $queryend;
    my  $sbjctend;
    my  $sbjctlen;

    my ($desc ) = ( 0 );

    unless( open( IN, $input )) {
        printf( STDERR "ERROR: Failed to open $input: $!\n" );
        return 0;
    }
    while( <IN> ) {
##        next if /^$/; ## comment to make eof( IN ) to work
        chomp;
        $last = $_;

        if( $last =~ /^\s+Query\s+\(\d+/ ) {
            $$refdesc{$query} = '';
            $desc = 1;
            next;
        }
        if(($last =~ /^\s+Database:/ || $last =~ /^$/ )&& $desc ) {
            $desc = 0;
            next;
        }
        if( $desc ) {
            $last =~ s/^\s*(.+)\s*$/$1/;
            $$refdesc{$query} .= "$last ";
            next;
        }


        if( $last =~ /^>(.*)\s*$/ || 
           ($last =~ /^\s+Score\s+=/ && $queryfasta ) ||
            eof( IN ))
        {
            $titletmp = $1;
            $titletmp = '' if $last =~ /^\s+Score/;
            $hitnum++ if $last =~ /^>/;
            $skip = 0;

            if( defined( $sbjct ) && defined( $e_val )) {
                if( !$sbjct ) {
                    printf( STDERR "WARNING: No subject description (Hit no. %d). Skipped.\n", $hitnum - 1 ) if $ALTALN;
                    $skip = 1;
                }
                elsif( $e_val =~ /^\-1/ ) {
                    printf( STDERR "WARNING: No expect value (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( length( $queryfasta ) != length( $sbjctfasta ) || length( $sbjctfasta ) == 0 ) {
                    printf( STDERR "WARNING: Invalid alignment (sequence lengths, Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( $querystart == -1 || $sbjctstart == -1 ) {
                    printf( STDERR "WARNING: No alignment start positions (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( $queryend == -1 || $sbjctend == -1 ) {
                    printf( STDERR "WARNING: No alignment end positions (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }

                $sbjctname = ($sbjct =~ /^([\w\.]+)/)? $1: '';
                $skip = 1 if $evalue < $e_val || $e_val < 0.0;
                $skip = 1 if !$skip && (scalar(keys %$rhSBJCTs) && !$$rhSBJCTs{$sbjctname});

                unless( $skip ) {
                    $$refhash{$query}[$rec][0] = $sbjct;
                    $$refhash{$query}[$rec][1] = $e_val;

                    $$refhash{$query}[$rec][2] = ''; ## reserved
                    $$refhash{$query}[$rec][3] = ''; ## reserved
                    $$refhash{$query}[$rec][4] = $queryfasta;
                    $$refhash{$query}[$rec][5] = $sbjctfasta;

                    $$refhash{$query}[$rec][6] = $querystart;
                    $$refhash{$query}[$rec][7] = $queryend;
                    $$refhash{$query}[$rec][8] = $sbjctstart;
                    $$refhash{$query}[$rec][9] = $sbjctend;
                    $$refhash{$query}[$rec][10]= $sbjctlen;

                    $rec++;
                }
            }

            $sbjct = $titletmp if $last =~ /^>/ || !$ALTALN;
            $sbjctlen = -1 if $last =~ /^>/;
            $queryfasta = '';
            $sbjctfasta = '';
            $querystart = -1;
            $sbjctstart = -1;
            $queryend = -1;
            $sbjctend = -1;
            $e_val = -1;
            next if $last =~ /^>/;
        }

        if( $last =~ /^>/ ) {
            printf( STDERR "ERROR: Missed hit (Hit no. %d: %s...). Probably wrong format. Terminating.\n",
                    $hitnum, substr( $sbjct, 0, 20 ));
            close( IN );
            return 0;
        }

        next if !defined( $e_val );

        if( $last =~ /^\s+Length\/ENO:\s+Query\s+=\s+\d+\/\d+\.\d,\s+Sbjct\s+=\s+(\d+)\/\d+\.\d/ ) {
            $sbjctlen = $1;
            next;
        }

        if( $sbjctlen && $sbjctlen < 0 ) {
            $last =~ s/^\s*(.+)\s*$/$1/;
            $sbjct .= " $last";
            next;
        }

        if( $last =~ /^\s+Score\s+=\s+([0-9\.eE\-\+]+)(?:\s+\([0-9\.eE\-\+]+\s+bits\))?,
                       \s+Expect\s+=\s+([0-9\.eE\-\+]+|n\/a),\s+/x )
        {
            $score = $1;
            $e_val = $2;
            $e_val = "1$e_val" if $e_val =~ /^e/i;
            next;
        }

        if( $last =~ /^Query:\s+([0-9]+)\s+([a-zA-Z\-]+)\s+([0-9]+)\s*$/ ) {
            $querystart = $1 if $querystart == -1;
            $queryend = $3;
            $queryfasta .= uc( "$2" );
            next;
        }

        if( $last =~ /^Sbjct:\s+([0-9]+)\s+([a-zA-Z\-]+)\s+([0-9]+)\s*$/ ) {
            $sbjctstart = $1 if $sbjctstart == -1;
            $sbjctend = $3;
            $sbjctfasta .= uc( "$2" );
            next;
        }
    }

    close( IN );
    return 1;
}

## -------------------------------------------------------------------
## print multiply aligned top hits
##

sub PrintAligned
{
    my  $format = shift;##output format
    my  $rhithash = shift;##reference
    my  $filename = shift;
    my  $width = shift;
    my  $rproc = shift;

    my  $rec;
    my  $e_val;
    my  $query;
    my  $sbjct;
    my  $qufasta;
    my  $sbfasta;

    my  $querystart;
    my  $sbjctstart;
    my  $queryend;
    my  $sbjctend;

    my  $info;##info string

    unless( open( OUT, ">$filename" )) {
        printf( STDERR "ERROR: Failed to open $filename for writing.\n" );
        return 0;
    }

    if( !($format cmp "Rosetta")) {
        foreach $query( keys %{$rhithash} ) {
            for( $rec = 0; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $sbjct      = $$rhithash{$query}[$rec][0];
                $e_val      = $$rhithash{$query}[$rec][1];
                $qufasta   = \$$rhithash{$query}[$rec][4];
                $sbfasta   = \$$rhithash{$query}[$rec][5];

                $querystart = $$rhithash{$query}[$rec][6];
                $queryend   = $$rhithash{$query}[$rec][7];
                $sbjctstart = $$rhithash{$query}[$rec][8];
                $sbjctend   = $$rhithash{$query}[$rec][9];

                PrintAlignedHelperPass1(\*OUT, $format, \$info, $width, $sbjct, $sbjctstart, $sbjctend, $sbfasta, $e_val, $rec);
            }
        }
    }

    foreach $query( keys %{$rhithash} ) {
        for( $rec = 0; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
            $sbjct      = $$rhithash{$query}[$rec][0];
            $e_val      = $$rhithash{$query}[$rec][1];
            $qufasta   = \$$rhithash{$query}[$rec][4];
            $sbfasta   = \$$rhithash{$query}[$rec][5];

            $querystart = $$rhithash{$query}[$rec][6];
            $queryend   = $$rhithash{$query}[$rec][7];
            $sbjctstart = $$rhithash{$query}[$rec][8];
            $sbjctend   = $$rhithash{$query}[$rec][9];

            PrintAlignedHelperPass2(\*OUT, $format, \$info, $width, $sbjct, $sbjctstart, $sbjctend, $sbfasta, $e_val, $rec);
            $$rproc++;
        }
    }

    close( OUT );
    return 1;
}

sub PrintAlignedHelperPass1
{
    my  $reffile = shift;##file reference
    my  $format = shift;
    my  $rinfo = shift;
    my  $width = shift;
    my  $sbjct = shift;
    my  $sbjctstart = shift;
    my  $sbjctend = shift;
    my  $sbfasta = shift;
    my  $e_val = shift;
    my  $rec = shift;

    if( !($format cmp "Rosetta")) {
        $sbjct =~ s/^(\S+)\s*(.*)$/$1 ($sbjctstart-$sbjctend) $2  Expect=$e_val/ if $rec;
        printf( $reffile "## %s\n", $sbjct );
        $$rinfo .= " $1" if $sbjct =~ /^(\S+)/;
    }
}

sub PrintAlignedHelperPass2
{
    my  $reffile = shift;##file reference
    my  $format = shift;
    my  $rinfo = shift;
    my  $width = shift;
    my  $sbjct = shift;
    my  $sbjctstart = shift;
    my  $sbjctend = shift;
    my  $sbfasta = shift;
    my  $e_val = shift;
    my  $rec = shift;

    if( !($format cmp "FASTA")) {
        $sbjct =~ s/^(\S+)\s*(.*)$/$1 ($sbjctstart-$sbjctend) $2  Expect=$e_val/ if $rec;
        printf( $reffile ">%s\n", $sbjct );
    } elsif( !($format cmp "Rosetta")) {
        printf( $reffile "##%s\nscores from program: 0\n", $$rinfo) unless $rec;
        printf( $reffile "0 ");
    }
    WrapFasta( $reffile, $sbfasta, $width );
}

## -------------------------------------------------------------------
## make multiple alignment target-oriented
##

sub MakeTargetOriented
{
    my  $rhithash = shift;  ## reference
    my  $stindex = shift;   ## start index
    my  $ralnmask = shift;  ## ref. to alignment mask of query
    my  $length = shift;    ## length of mask
    my  $tolxs = shift;     ## tolerate Xs in query sequencies
    my  $checkonly = shift; ## flag of iindication to check only
    my  $query;
    my ($qufasta, $sbfasta, $qufasta1st );
    my ($querystart, $queryend );
    my ($qustart, $quend );
    my ($qures, $qu1st, $l, $i, $rec );

    unless( $ralnmask && ref( $ralnmask )) {
        printf( STDERR "ERROR: Reference expected.\n" );
        return 0;
    }
    $$ralnmask = '';
    if( $length && 0 < $length ) {
        $$ralnmask .= '-' for( 0 .. $length-1 );
    }

    foreach $query( keys %{$rhithash} ) {
        for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
            unless( $qufasta1st ) {
                $qufasta1st = \$$rhithash{$query}[$stindex][4];
                $querystart = \$$rhithash{$query}[$stindex][6];
                $queryend   = \$$rhithash{$query}[$stindex][7];
            }
            $qufasta = \$$rhithash{$query}[$rec][4];
            $sbfasta = \$$rhithash{$query}[$rec][5];
            $qustart =  $$rhithash{$query}[$rec][6];
            $quend   =  $$rhithash{$query}[$rec][7];

            if( $$ralnmask ) {
                for( $i = $qustart; $i <= $quend; $i++ ) {
                    substr( $$ralnmask, $i-1, 1 ) = $X if 0 < $i && $i <= length( $$ralnmask );
                }
            }
            next if $rec <= $stindex;

            if( length( $$qufasta ) != length( $$sbfasta )) {
                printf( STDERR "ERROR: Lengths of query and subject are not equal.\n" );
                return 0;
            }
            if( length( $$qufasta ) != length( $$qufasta1st )) {
                printf( STDERR "ERROR: Lengths of consecutive query sequences are not equal.\n" );
                return 0;
            }

            for( $l = 0; $l < length( $$qufasta ); $l++ ) {
                $qures = substr( $$qufasta, $l, 1 );
                $qu1st = substr( $$qufasta1st, $l, 1 );
                next if( $qures eq '-' && $qu1st eq '-' );
                if( $qures ne '-' && $qu1st ne '-' && $qures ne $qu1st && 
                  ( $tolxs?( $qures ne $X && $qu1st ne $X ): 1 )) {
                    printf( STDERR "ERROR: Inconsistent query sequences in alignment.\n" );
                    return 0;
                }
                next if $checkonly; ## do not change alignments below if check only
                if( $qures ne '-' ) {
                    substr( $$qufasta1st, $l, 1 ) = $qures;
                    $$querystart = $qustart if $qustart < $$querystart;
                    $$queryend = $quend if $$queryend < $quend;
                }
            }
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## align top hits from coma output in fasta
##

sub MakeMultipleAlignment
{
    my  $rhithash = shift; ## reference
    my  $stindex = shift;  ## start index
    my  $tolxs = shift;    ## tolerate Xs in query sequences
    my  $query;
    my (%minstart, %maxend, %posits );
    my ($qufasta, $sbfasta );
    my ($querystart, $queryend );
    my ($pos, $qgap, $allmatch );
    my ($qsym, $sym, $c, $l, $rec );

    foreach $query( keys %{$rhithash} ) {
        $minstart{$query} = 99999; ##large enough
        $maxend{$query} = 0;
        for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
            $qufasta    = \$$rhithash{$query}[$rec][4];
            $querystart =  $$rhithash{$query}[$rec][6];
            $queryend   =  $$rhithash{$query}[$rec][7];

            $posits{$query}[$rec] = $querystart;

            $minstart{$query} = $querystart if( $querystart < $minstart{$query} );
            $maxend{$query} = $queryend if( $maxend{$query} < $queryend );
        }
    }

    foreach $query( keys %{$rhithash} ) {
        for( $c = $minstart{$query}, $l = 0; $c <= $maxend{$query}; $l++ ) {
            undef $qsym;
            $qgap = 0;
            $allmatch = 1;
            for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $qufasta    = \$$rhithash{$query}[$rec][4];
                $querystart =  $$rhithash{$query}[$rec][6];
                $queryend   =  $$rhithash{$query}[$rec][7];
                $pos = $posits{$query}[$rec];

                next if( $pos != $c );
                unless( $qsym ) {
                    $qsym = substr( $$qufasta, $l, 1 );
                    $qgap = 1 if $qsym eq '-';
                }
                $sym = substr( $$qufasta, $l, 1 );
                ## added protection from livelock
                if( !$tolxs && $sym ne '-' && $qsym ne '-' && $sym ne $qsym ) {
                    printf( STDERR "ERROR: Inconsistent query sequences in alignment.\n" );
                    return 0;
                }
                do { $allmatch = 0; last; } if(( $sym ne $qsym ) && !( $tolxs?(( $qsym eq $X && $sym ne '-')||( $sym eq $X && $qsym ne '-')): 0 ));
            }

            for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $qufasta = \$$rhithash{$query}[$rec][4];
                $sbfasta = \$$rhithash{$query}[$rec][5];

                $querystart =  $$rhithash{$query}[$rec][6];
                $queryend   =  $$rhithash{$query}[$rec][7];

                $pos = $posits{$query}[$rec];
                if( $pos != $c ) {
                    if( $c < $pos ) {
                        substr( $$qufasta, $l, 1 )  = '-' . substr( $$qufasta, $l, 1 );
                        substr( $$sbfasta, $l, 1 )  = '-' . substr( $$sbfasta, $l, 1 );
                    }
                    else {
                        substr( $$qufasta, $l, 1 ) .= '-';
                        ##substr( $$sbfasta, $l, 1 ) .= '-';
                        if( substr( $$sbfasta, $l, 1 ) ne '-' ) {
                                substr( $$sbfasta, $l, 1 )  = '-' . substr( $$sbfasta, $l, 1 );
                        } else{ substr( $$sbfasta, $l, 1 ) .= '-';
                        }
                    }
                    next;
                }

                next if $allmatch;

                $qsym = substr( $$qufasta, $l, 1 );
                $qgap = 1 if $qsym eq '-';

                if( $qsym ne '-' ) {
                    substr( $$qufasta, $l, 1 ) = '-' . substr( $$qufasta, $l, 1 );
                    substr( $$sbfasta, $l, 1 ) = '-' . substr( $$sbfasta, $l, 1 );
                }
            }
            if( $allmatch ) {
                unless( $qgap ) {
                    for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                        $queryend   =  $$rhithash{$query}[$rec][7];
                        $pos = $posits{$query}[$rec];
                        $posits{$query}[$rec]++ if( $pos == $c && $c < $queryend );
                    }
                    $c++;
                }
            }
        }
        ## add terminal gaps if needed
        for( ;; $l++ ) {
            for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $qufasta = \$$rhithash{$query}[$rec][4];
                $sbfasta = \$$rhithash{$query}[$rec][5];
                last if $l < length( $$sbfasta );
            }
            last if $#{$$rhithash{$query}} < $rec;
            for( $rec = $stindex; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $qufasta = \$$rhithash{$query}[$rec][4];
                $sbfasta = \$$rhithash{$query}[$rec][5];
                if( length( $$sbfasta ) <= $l ) {
                    substr( $$sbfasta, $l, 1 ) .= '-';
                    substr( $$qufasta, $l, 1 ) .= '-';
                }
            }
        }
    }

    return 1;
}

## -------------------------------------------------------------------
## insert aligned query sequence in the first place
##

sub InsertQuerySequence
{
    my  $rhithash = shift;
    my  $stindex = shift; ## start index
    my  $rdescrip = shift;
    my  $rqueryin = shift;
    my  $alnmask = shift; ## alignment mask of query sequence
    my  %minstart;
    my  %maxend;
    my  $error = 0;
    my ($query, $qdesc );
    my ($qufasta, $sbfasta );
    my ($orqufasta, $orsbfasta );
    my ($querystart, $queryend );
    my ($rec, $length, $alnlength );
    my ($msqur, $orqur, $qur );
    my ($c, $l, $m ) = (1,0,0);

    $stindex-- if 0 < $stindex;

    if( exists $$rqueryin{DESC} && exists $$rqueryin{SEQN} )
    {
        foreach $query( keys %{$rhithash} ) {
            $minstart{$query} = 99999; ##large enough
            $maxend{$query} = 0;
            for( $rec = $stindex+1; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                $querystart =  $$rhithash{$query}[$rec][6];
                $queryend   =  $$rhithash{$query}[$rec][7];

                $minstart{$query} = $querystart if( $querystart < $minstart{$query} );
                $maxend{$query} = $queryend if( $maxend{$query} < $queryend );
            }
        }
        foreach $query( keys %{$rhithash} ) {
            $length = length( $$rqueryin{SEQN} );
            next if $#{$$rhithash{$query}} < $stindex;
            $rec = $stindex;


            $$rhithash{$query}[$rec][0] = $$rqueryin{DESC};
            $$rhithash{$query}[$rec][1] = 0;  ## e-value

            $$rhithash{$query}[$rec][2] = ''; ## reserved
            $$rhithash{$query}[$rec][3] = ''; ## reserved
            $$rhithash{$query}[$rec][4] = $$rqueryin{SEQN};
            $$rhithash{$query}[$rec][5] = $$rqueryin{SEQN};

            $$rhithash{$query}[$rec][6] = 1;
            $$rhithash{$query}[$rec][7] = $length;
            $$rhithash{$query}[$rec][8] = 1;
            $$rhithash{$query}[$rec][9] = $length;
            $$rhithash{$query}[$rec][10]= $length;

            next if $#{$$rhithash{$query}} <= $stindex;
            $alnlength = length( $$rhithash{$query}[$rec+1][4] );

            for( $c = 1, $l = 0, $m = 0; $c <= $length || $l < $alnlength; $l++ )
            {
                $qufasta   = \$$rhithash{$query}[$stindex+1][4];
                $orqufasta = \$$rhithash{$query}[$stindex][4];
                $orsbfasta = \$$rhithash{$query}[$stindex][5];
                $qur = substr( $$qufasta, $l, 1 );
                $orqur = substr( $$orqufasta, $l, 1 );
                $msqur = '';
                $msqur = substr( $alnmask, $m, 1 ) if $alnmask && $m < length( $alnmask );

                if( $c < $minstart{$query} && $c <= $length ) {
                    for( $rec = $stindex+1; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                        $qufasta = \$$rhithash{$query}[$rec][4];
                        $sbfasta = \$$rhithash{$query}[$rec][5];
                        substr( $$qufasta, $l, 1 )  = '-' . substr( $$qufasta, $l, 1 );
                        substr( $$sbfasta, $l, 1 )  = '-' . substr( $$sbfasta, $l, 1 );
                    }
                    $alnlength++;##corrected
                    $c++; $m++;
                    next;
                }
                if( $maxend{$query} < $c && $c <= $length ) {
                    ## if query ends with gaps
                    if( $qur eq '-' ) {
                        substr( $$orqufasta, $l, 1 )  = '-' . substr( $$orqufasta, $l, 1 );
                        substr( $$orsbfasta, $l, 1 )  = '-' . substr( $$orsbfasta, $l, 1 );
                        next;
                    }
                    for( $rec = $stindex+1; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
                        $qufasta = \$$rhithash{$query}[$rec][4];
                        $sbfasta = \$$rhithash{$query}[$rec][5];
                        substr( $$qufasta, $l, 1 )  .= '-';
                        ##substr( $$sbfasta, $l, 1 )  .= '-';
                        if( substr( $$sbfasta, $l, 1 ) ne '-' ) {
                                substr( $$sbfasta, $l, 1 )  = '-' . substr( $$sbfasta, $l, 1 );
                        } else{ substr( $$sbfasta, $l, 1 ) .= '-';
                        }
                    }
                    $c++; $m++;
                    next;
                }

                if( $length < $c ) {
                    substr( $$orqufasta, $l, 1 ) .= $qur;
                    substr( $$orsbfasta, $l, 1 ) .= $qur;
                    next;
                }
                ## position's been never aligned
                if( $qur eq '-' && $msqur && $msqur eq '-' && 1 < $c ) {
                    $c++; $m++;
                    next;
                }
                if( $qur eq '-' ) {
                    substr( $$orqufasta, $l, 1 )  = '-' . substr( $$orqufasta, $l, 1 );
                    substr( $$orsbfasta, $l, 1 )  = '-' . substr( $$orsbfasta, $l, 1 );
                    next;
                }
                if( $orqur ne $qur && $orqur ne 'X' && $qur ne 'X') {
                    printf( STDERR "ERROR: Inconsistent query sequence data. Not inserted.\n" );
                    $error = 1;
                    last;
                }
                $c++; $m++;
            }
        }
        return 1 unless $error;
    }

    foreach $query( keys %{$rhithash} ) {
        next if $#{$$rhithash{$query}} <= $stindex;
        $rec = $stindex;
            
        $$rhithash{$query}[$rec][0] = $$rdescrip{$query};
        $$rhithash{$query}[$rec][1] = 0;  ## e-value

        $$rhithash{$query}[$rec][2] = ''; ## reserved
        $$rhithash{$query}[$rec][3] = ''; ## reserved
        $$rhithash{$query}[$rec][4] = $$rhithash{$query}[$rec+1][4];
        $$rhithash{$query}[$rec][5] = $$rhithash{$query}[$rec+1][4];##ok

        $$rhithash{$query}[$rec][6] = $querystart = $$rhithash{$query}[$rec+1][6];
        $$rhithash{$query}[$rec][7] = $queryend   = $$rhithash{$query}[$rec+1][7];
        $$rhithash{$query}[$rec][8] = $$rhithash{$query}[$rec+1][6];##ok
        $$rhithash{$query}[$rec][9] = $$rhithash{$query}[$rec+1][7];##ok
        $$rhithash{$query}[$rec][10]= $$rhithash{$query}[$rec+1][7];##ok

        $qdesc = \$$rhithash{$query}[$rec][0];
        $$qdesc =~ s/^(\S+)\s*(.*)$/$1 ($querystart-$queryend) $2/;
    }

    return 1;
}



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## -------------------------------------------------------------------
## wrap sequence to fragments of equal length
##

sub WrapFasta {
    my  $reffile = shift;   ## reference to file descriptor
    my  $reffasta = shift;  ## reference to sequence
    my  $width = shift;     ## width of fragment per line
    my  $padding = 0;       ## padding at the beginning of each line
    my  $line;

    $width = 999999 if $width <= 0;

    if( ref( $reffile ) ne 'GLOB' && ref( $reffile ) ne 'SCALAR' ) {
        printf( STDERR "ERROR: WrapFasta: Wrong reference.\n" );
        return 0;
    }

##    $$reffile = '' if( ref( $reffile ) eq 'SCALAR' );

    for( my $n = 0; $n < length( $$reffasta ); $n += $width ) {
        if( $n && $padding ) {
            $line = sprintf( "%${padding}s%s\n", ' ', substr( $$reffasta, $n, $width ));
        } else {
            $line = sprintf( "%s\n", substr( $$reffasta, $n, $width ));
        }
        if( ref( $reffile ) eq 'SCALAR' ) {
                 $$reffile .= $line;
        } else { printf( $reffile $line );
        }
    }
    return 1;
}

## <<>>

