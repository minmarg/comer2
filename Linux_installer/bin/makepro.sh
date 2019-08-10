#!/bin/bash

basename=$( basename $0 )
dirname=$( dirname $0 )
[[ "${dirname:0:1}" != "/" ]] && dirname="$( pwd )/$dirname"
updirname=$( dirname $dirname )


INSTPSIPRED="/share/install/psipred"
INSTBLASTP="/share/install/ncbi-blast/ncbi-blast-2.2.23+"
INSTCOMER="$updirname"
INPUT=
OUTPUT=
OPTIONS=''
NOSELECT='--select'
NOFLUSH=''
VERBOSE=''

usage="
Make \`comer' profile with inclusion of secondary structure prediction.
2015(C)Mindaugas Margelevicius,VU IBT,Vilnius

$basename <Parameters>

Parameters:

-i <input>     Input multiple alignment file either in 
               FASTA or in STOCKHOLM.
-o <output>    Output filename of profile.
-p <options>   Input file of options;
               By default, the file in installation
               directory of this package is searched.
-s             Do not preprocess input multiple alignment before 
               predicting SS.
-r             Do not flush temporary files.
-P <path>      Installation path of \`PSIPRED'
       default=$INSTPSIPRED
-B <path>      Installation path of \`BLAST+'
       default=$INSTBLASTP
-v             Enable warnings.
-h             This text.
"


while getopts "i:o:p:srP:B:vh" Option
do
    case $Option in
        i ) INPUT=${OPTARG} ;;
        o ) OUTPUT=${OPTARG} ;;
        p ) OPTIONS="-p ${OPTARG}" ;;
        s ) NOSELECT='' ;;
        r ) NOFLUSH='--noflush' ;;
        P ) INSTPSIPRED=${OPTARG} ;;
        B ) INSTBLASTP=${OPTARG} ;;
        v ) VERBOSE='-v' ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo Error: Unrecognized argument.; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

if [[ -z "$OUTPUT" ]]; then echo "ERROR: No output filename given."; exit 1; fi
if [[ -z "$INPUT" ]]; then echo "ERROR: No input filename given."; exit 1; fi
if [[ ! -f "$INPUT" ]]; then echo "ERROR: Input file not found: $INPUT"; exit 1; fi

if [[ ! -d "$INSTPSIPRED" ]]; then echo "ERROR: Directory does not exists: $INSTPSIPRED"; exit 1; fi
if [[ ! -d "$INSTBLASTP" ]]; then echo "ERROR: Directory does not exists: $INSTBLASTP"; exit 1; fi
if [[ ! -d "$INSTCOMER" ]]; then echo "ERROR: Directory does not exists: $INSTCOMER"; exit 1; fi

makepro="$INSTCOMER/bin/makepro"
ssp="$INSTCOMER/bin/ssp2.pl"
inssp="$INSTCOMER/bin/inssp2.pl"
SSFILE="$OUTPUT.ss"

if [[ ! -f "$makepro" ]]; then echo "ERROR: Executable \`makepro' from the \`comer' package not found: $makepro"; exit 1; fi
if [[ ! -f "$ssp" ]]; then echo "ERROR: A perl script from the \`comer' package not found: $ssp"; exit 1; fi
if [[ ! -f "$inssp" ]]; then echo "ERROR: A perl script from the \`comer' package not found: $inssp"; exit 1; fi

cmd="$makepro $VERBOSE -i $INPUT -o $OUTPUT $OPTIONS"
[[ -n "$VERBOSE" ]] && echo $cmd
eval $cmd || exit 1

cmd="$ssp --in $INPUT --out $SSFILE $OPTIONS $NOSELECT $NOFLUSH "
cmd+="--psipred $INSTPSIPRED --blast $INSTBLASTP --comer $INSTCOMER"
[[ -n "$VERBOSE" ]] && echo $cmd
eval $cmd || exit 1

cmd="$inssp --in $SSFILE --to $OUTPUT"
[[ -n "$VERBOSE" ]] && echo $cmd
eval $cmd || exit 1

echo
[[ -z "$NOFLUSH" ]] && rm "$SSFILE"

