/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __makedb_h__
#define __makedb_h__

// #define DEFAULT_DISTRIBUTION_TYPE   ( DISCRETE )

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Make profile database from input files in the given directory or listed profiles or\n\
COMER database files (with the .prd extension).\n\
(C)2013-2019 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -o <output> [<other_options>] ( -d <directory> | <file1> <file2> ... )\n\
\n\
Make options:\n\
\n\
-o <output_name>            Name of output database.\n\
-d <directory>              Directory of files to read.\n\
-n <number_of_open_files>   Maximum number of open input files.\n\
                        Default=510\n\
\n\
SEG options:\n\
-U                          Invoke low-complexity filtering of profiles.\n\
-w <window>     [Integer]   Window length.                                  ( 12)\n\
-f <low>        [Real]      Low entropy threshold.                          (2.2)\n\
-F <high>       [Real]      High entropy threshold.                         (2.5)\n\
-D <distance>   [Real]      Distance of equivalence between profile       (12.96)\n\
                            vectors.\n\
\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o my_db -d ./my_profiles\n\
<> -o my_db d_70_1_2.pro b_119_1_1.pro c_17_1_1.pro c_69_1_5.pro\n\
<> -o my_db mydb_part1.prd mydb_part2.prd mydb_part3.prd\n\
<> -o my_db *.pro\n\
\n\
";

#endif//__makedb_h__
