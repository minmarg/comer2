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
// 2.02   Re-engineered load distribution on GPUS leading to significantly improved performance
// 2.03   Masked X positions
// 2.04   JSON format for output introduced


static const char*  version = "2.03.03";
static const char*  verdate = "";

static const char*  instructs = "\n\
<>\n\
\n\
Protein remote homology search tool.\n\
(C)2013-2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University\n\
\n\
\n\
Usage:\n\
<> -i <query> -d <database> [-o <output>] [-p <options>] [<other_options>]\n\
\n\
Basic search options:\n\
-i <query_profiles>         <> (stacked) profile(s) made by makepro or profile\n\
                            database (use extension .prd) made by makedb+db2bin.\n\
-d <database_name_list>     Database made by makedb+db2bin from this package.\n\
                            NOTE: Multiple comma-separated database names can be\n\
                            provided.\n\
-o <output_directory>       Directory of output files for each query.\n\
-f <file_format>            Format of output files: 0, Plain; 1, JSON.\n\
                        Default=0\n\
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
                            NOTE: For a small number of queries, using a moderate\n\
                            amount of memory (~4GB) is more efficient.\n\
                        Default=[all memory of GPU(s)]\n\
--dev-ngrids=<count>        Number of grids [1,20] to launch on each GPU.\n\
 *EXPERIMENTAL option*      A grid represents a logically independent computation unit\n\
                            (like a bunch of CPU threads). The GPU memory is divided\n\
                            over all grids, and the advantage of running multiple grids\n\
                            is determined by the available number of GPU cores. Note\n\
                            that a communicating CPU thread is created for each grid.\n\
                        Default=1\n\
--dev-expected-length=<length> Expected length of database proteins (profiles). Its\n\
                            values are restricted to the interval [20,200].\n\
                            NOTE: Increasing it reduces GPU memory requirements, but\n\
                            mispredictions may cost additional computation time.\n\
                        Default=50\n\
--dev-pass2memp=<percentage> GPU memory proportion dedicated to recalculation of\n\
                            hits that passed significance threshold.\n\
                            (Expected proportion of significant hits.)\n\
                        Default=10\n\
--io-nbuffers=<count>       Number of buffers [1,10] used to cache data read from file.\n\
                            Values greater than 1 lead to increased performance at the\n\
                            expense of increased memory consumption.\n\
                        Default=4\n\
--io-filemap                Map files into memory.\n\
                            In general, the performance with or without file mapping is\n\
                            similar. In some systems, however, mapping can lead to\n\
                            increased computation time.\n\
--io-unpinned               Do not use pinned (page-locked) CPU memory.\n\
                            Pinned CPU memory provides better performance, but reduces\n\
                            system resources. If RAM memory is scarce (<2GB), using\n\
                            pinned memory may reduce overall system performance.\n\
                            By default, pinned memory is used.\n\
\n\
Other options:\n\
--dev-list                  List all GPUs compatible and available on the system\n\
                            and exit.\n\
-v [<level_number>]         Verbose mode.\n\
-h                          This text.\n\
\n\
\n\
Examples:\n\
<> -i myprofile.pro -d mydb -o my_output_directory (database used without extension)\n\
<> -i mystackedprofiles.pro -d mydb1,mydb2,mydbn -o my_out_dir (multiple stacked query profiles; multiple databases)\n\
<> -i mydb1 -d mydb2 -o my_output_directory (query database searched against mydb2)\n\
<> -i mydb1.prd -d mydb2 -o my_out_dir -p my_search_options.txt\n\
\n\
";

#endif//__comer_h__
