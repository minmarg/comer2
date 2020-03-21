/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __makepro_h__
#define __makepro_h__

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Make profile from multiple sequence alignment in FASTA or STOCKHOLM.\n\
(C)2013-2019 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-p <options>] [<other_options>]\n\
\n\
Profile construction options:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file.\n\
-o <output>     [Filename]  Output profile.\n\
-p <options>    [Filename]  Input file of options.\n\
                        By default, the options file in the installation\n\
                            directory is used.\n\
\n\
Other options:\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
";

#endif//__makepro_h__
