/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __db2bin_h__
#define __db2bin_h__

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Make binary representation of profile database.\n\
(C)2013-2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -d <database> [<other_options>]\n\
\n\
Main options:\n\
-d <database_pathname>      Name of profile database (without extension).\n\
-p <options_filename>       Input file of options.\n\
                        By default, the options file in the installation\n\
                            directory is used.\n\
\n\
Other options:\n\
-v                          Verbose mode.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -d my_db\n\
\n\
";

#endif//__db2bin_h__
