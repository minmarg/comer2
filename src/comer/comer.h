/***************************************************************************
 *   Copyright (C) 2013-2019 by Mindaugas Margelevicius                    *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __comer_h__
#define __comer_h__


// #define DEFAULT_SCORING_SCHEME          ( AbstractScoreMatrix::ProfileSpecific )
// #define DEFAULT_STATISTICAL_BEHAVIOUR   ( AbstractScoreMatrix::ComputeStatistics )
// #define DEFAULT_PRECISION               ( AbstractScoreMatrix::AutoScalling )
// #define DEFAULT_MASKING                 ( MaskToIgnore )
// #define DEFAULT_ALNALGORITHM            ( ProfileAlignment::Optimizational )

// Version history:
// 2.01   Initial project 2

static const char*  version = "2.01";
static const char*  verdate = "";

static const char*  instructs = "\n\
<>\n\
\n\
Protein remote homology search tool.\n\
(C)2013-2019 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <query> -d <database> [-o <output>] [-p <options>] [<other_options>]\n\
\n\
Alignment search options:\n\
-i <query_profile>          <> profile made by makepro.\n\
-d <database_name>          Database made by makedb from this package.\n\
-o <output_directory>       Directory of output files for each query.\n\
-p <options_file>           Input file of options;\n\
                            By default, the options file in the installation\n\
                            directory is used.\n\
\n\
Parallel run options:\n\
--dev-N=(<number>|,<id_list>) Maximum number of GPUs to use. This can be\n\
                            specified by a number or given by a comma-separated\n\
                            list of GPU identifiers, which should start with a comma.\n\
                            In the latter case, work is distributed in the\n\
                            specified order. Otherwise, more powerful GPUs are\n\
                            selected first.\n\
                            NOTE: The first symbol preceding a list is a comma.\n\
                        Default=1 (most powerful GPU)\n\
--dev-mem=<megabytes>       Maximum amount of GPU memory (MB) that can be used.\n\
                            All memory is used if a GPU has less than the specified\n\
                            amount of memory.\n\
                        Default=[all memory of GPU(s)]\n\
--dev-cachep=<percentage>   GPU memory proportion dedicated to caching the profile\n\
                            database. This does not include the size of all queries.\n\
                            The size will be calculated with respect to the first GPU\n\
                            (see options --dev-N and --dev-list) and it will be the\n\
                            same for all GPUs.\n\
                            NOTE: Large percentages (>40%) may be desirable for GPUs\n\
                            with plenty of memory (>6GB).\n\
                        Default=20\n\
--dev-expected-length=<length> Expected length of database proteins (profiles). Its\n\
                            values are restricted to the interval [20,200].\n\
                            NOTE: Increasing it reduces GPU memory requirements, but\n\
                            mispredictions may cost additional computation time.\n\
                        Default=50\n\
--dev-pass2memp=<percentage> GPU memory proportion dedicated to recalculation of\n\
                            hits that passed significance threshold.\n\
                            (Expected proportion of significant hits.)\n\
                        Default=10\n\
--io-nofilemap              Do not map files into memory.\n\
                            Mapping is advantageous when working with large files.\n\
                            This options allows to turn it off.\n\
\n\
CPU options:\n\
--cpu-freemem               Free unused memory as soon as the required size is known.\n\
\n\
Other options:\n\
--dev-list                  List all GPUs compatible and available on the system\n\
                            and exit.\n\
-v [<level_number>]         Verbose mode.\n\
-h                          This text.\n\
\n\
";

#endif//__comer_h__
