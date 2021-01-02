/***************************************************************************
 *   Copyright (C) 2013-2021 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __neff_h__
#define __neff_h__

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Calculate the effective number of sequences Neff for a multiple sequence alignment in\n\
FASTA or STOCKHOLM. Neff is calculated using the formula: SUM_i w_i, where \n\
w_i = 1/k_i and k_i = SUM_j H(S(i,j)-S_thr), with S(i,j) sequence identity (normalized\n\
by the smaller sequence length) between sequences i and j, S_thr a sequence identity\n\
threshold, and H(.) the unit step function.\n\
(C)2013-2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <input> [<options>]\n\
\n\
Options:\n\
\n\
-i <input_MSA_file>         Input multiple alignment file.\n\
--idn=<list_of_thresholds>  Comma-separated list of sequence identity (integer)\n\
                            thresholds (in percentage) to calculate Neff at.\n\
                        Default=62\n\
\n\
Other options:\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
\n\
Examples:\n\
<> -i myMSA.afa\n\
<> -i myMSA.sto --idn=90,80,62\n\
\n\
";

#endif//__neff_h__
