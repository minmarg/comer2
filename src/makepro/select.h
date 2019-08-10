/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __select_h__
#define __select_h__

// Version history:
// 2.01   Initial project

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  selinst = "\n\
<>\n\
\n\
Select non-redundant set of sequences from multiple alignment.\n\
(C)2013-2019 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-p <options>] [-t]\n\
\n\
Parameters:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file either in\n\
                            FASTA or STOCKHOLM.\n\
-o <output>     [Filename]  Output file of multiple alignment.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, the options file in the installation\n\
                            directory is used.\n\
-t                          Do not apply the selection of sequences;\n\
                            just convert the input.\n\
\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
";

#endif//__select_h__
