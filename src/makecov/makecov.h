/***************************************************************************
 *   Copyright (C) 2013-2020 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __makecov_h__
#define __makecov_h__

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.02";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Calculate weighted cross-covariance matrix between the positions of multiple \n\
sequence alignment given in FASTA or STOCKHOLM.\n\
(C)2013-2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-p <options>] [<other_options>]\n\
\n\
Basic options:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file.\n\
-o <output>     [Filename]  Output file of cross-covariance matrices.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, the options file in the installation\n\
                            directory is used.\n\
\n\
Calculation options:\n\
--extents                   Process only those pairs of positions that are\n\
                            included in their respective extents.\n\
--mix                       Mix weighted positional observations with target\n\
                            frequencies.\n\
--xcorr                     Calculate cross-correlation matrices instead of\n\
                            cross-covariance matrices.\n\
--scale                     Scale (multiply) sequence weights by the number of\n\
                            contributing sequences.\n\
--mi                        Calculate and add mutual information to the end of\n\
                            calculated cross-covariance/correlation matrix\n\
                            values.\n\
\n\
Other options:\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
";

#endif//__makecov_h__
