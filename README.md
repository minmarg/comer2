COMER2, a cross-platform software package for
protein remote homology search and alignment

(C)2013-2021 Mindaugas Margelevicius,
Institute of Biotechnology, Vilnius University

v2.3.1

# Description

   The COMER method based on sequence profile-profile comparison is a
   sensitive, specific, and accurate computational tool developed for protein
   alignment and homology search. COMER version 2.2 (COMER2) represents one
   of the fastest implementations of calculations for sensitive protein
   homology search. High COMER2 performance is achieved by harnessing the
   power of the Graphics processing unit (GPU). Hence, a GPU is expected
   to be installed on the system.

   COMER2 is licensed under GNU General Public License version 3. Please find
   the LICENSE and COPYING files included in the software package.

# Available Platforms

   The COMER2 source code should compile and run on Linux, MS Windows, and 
   macOS. Please note, however, that COMER2 was tested on and the binaries 
   are provided for the following platforms:

  *  MS Windows 64 (64 bit)
  *  Linux x64 (64 bit)

   COMER2 was also compiled with the `clang` version 6 compiler and is
   expected to install and run on macOS.

# Hardware requirements

  *  CUDA-enabled GPU(s) with compute capability 3.5 (released in 2012) or
     above
  *  2 GB of RAM or more

# Getting the COMER Software

   The package is available at:

   [http://comer2.sourceforge.net](http://comer2.sourceforge.net)

   [https://github.com/minmarg/comer2](https://github.com/minmarg/comer2)

   The Docker image is available at:

   [https://hub.docker.com/r/minmar/comer2](https://hub.docker.com/r/minmar/comer2)

   See [Running a Docker container](#running-a-docker-container) for 
   information on using the docker image.

# COMER2 Profile databases

   Up-to-date COMER2 profile databases for PDB70, SCOP70, and PFAM sequence 
   families are available for download at:
   [https://sourceforge.net/projects/comer2/files/comer2-profile-databases-2.02](https://sourceforge.net/projects/comer2/files/comer2-profile-databases-2.02)

# Structure of the Package

   The main directories are described below:

  *  build -- an empty directory to contain built files.

  *  MS_Windows_installer -- contains an installer file for MS Windows.

  *  Linux_installer -- this directory contains the necessary files to 
     install the prebuilt COMER2 software on Linux.

  *  src -- the directory of the source files.

# Installation of pre-compiled binaries

   On MS Windows, run the installer:

     MS_Windows_installer\COMER2-installer1.msi

   NOTE: system requirements for the COMER2 software installed on Windows are 
   NVIDIA driver version >=425.25 and CUDA version >=10.1.

   On Linux, run the shell script and follow the instructions:

     Linux_installer/COMER2-installer1.sh

   NOTE: system requirements for the COMER2 software installed on Linux are 
   NVIDIA driver version >=418.87 and CUDA version >=10.1.

# Installation from source code

   ## Installation on MS Windows

   ### Software requirements

   To successfully build and install the COMER2 software from the source code
   on MS Windows, these tools are required to be installed:

  *  CMake version 3.8 or greater (free software)

  *  Visual C++ compiler, e.g., Visual Studio Community (free for open 
     source projects; COMER2 is an open source project)

  *  [the NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-downloads) version 10.0 or greater 
     (free software)

   ### Installation

   Run the command (batch) file:

     BUILD_and_INSTALL_win64.cmd

   ## Installation on Linux and macOS

   ### Software requirements

   To successfully build and install the COMER2 software from the source code
   on Linux or macOS, these tools are required to be installed:

  *  CMake version 3.8 or greater

  *  GNU Make version 3.81 or greater

  *  GNU GCC compiler version 4.8 or greater, or LLVM clang compiler
     version 6 or greater (or another C++ compiler that supports C++11)

  *  [the NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-downloads) version 10.0 or greater

   ### Installation

   Run the shell script:

     BUILD_and_INSTALL_unix.sh

# Getting Started

   The COMER2 software runs on CUDA-capable GPU devices. An appropriate
   NVIDIA driver must have been installed. (It can also be installed during
   the installation of the CUDA toolkit; see "Installation from source 
   code.")

   The software package contains four main programs in the bin directory in
   the installation path:

  *  makepro and makepro.sh (makepro.cmd on MS Windows), developed for
     making COMER profiles. It is recommended to use makepro.sh for enriching
     profiles with secondary structure predictions. makepro.sh, however,
     requires the external packages PSIPRED and PSI-BLAST to be installed.
     makepro and makepro.sh (makepro.cmd) make profiles in text format.
     Profiles, therefore, can be transferred between different platforms.

  *  makedb is developed for making a COMER profile database to be searched.
     makedb makes output profile databases in text format. They are also 
     cross-platform portable.

  *  db2bin, developed for converting a COMER profile database to binary
     format. For an n-fold read speedup, it is highly RECOMMENDED to 
     convert a profile database using db2bin before conducting homology
     search with the `comer` program. Please note that the output of db2bin
     is platform-dependent, and db2bin should be invoked on every platform.

  *  comer, the main program for homology search using one or more GPUs.

   Assuming that a query profile myprofile.pro and a profile database mydb
   (or databases mydb1 and mydb2) have been obtained, the simplest way to 
   run `comer` is to type (single database):

     comer -i myprofile.pro -d mydb -o my_output_directory

   or (mutliple databases)

     comer -i myprofile.pro -d mydb1,mydb2 -o my_output_directory

   where my_output_directory is an output directory to store output
   alignments files for each query profile present in the input file
   myprofile.pro.

   `comer` allows for multiple queries in the input file. In that case,
   profiles made using makepro or makepro.sh (makepro.cmd) should be stacked
   one on top of the other. It is also possible to search all profiles in
   one profile database against the profiles of another one:

     comer -i mydb1 -d mydb2 -o my_output_directory

   or perform an all-against-all comparison:

     comer -i mydb -d mydb -o my_output_directory

   Mutually aligning two profiles requires making a database of one of the
   two profiles:

     comer -i myprofile1.pro -d myprofile2_db -o my_output_directory

   `comer` search, as well as the process of making profiles, can be
   controlled with options read from the options file options.txt in the var
   directory in the installation path:

     comer -i myprofile.pro -d mydb -o my_output_directory -p options.txt

   The user can copy the options file options.txt to a preferred location
   and modify option values.

# Input Multiple Alignment

   The program makepro accepts input multiple alignment files in FASTA or
   STOCKHOLM format.

   The FASTA format can be described as follows. The section of each 
   sequence begins with a description line, whose first character is a ">"
   delimiter. Sequence data begins on the next line and can occupy multiple
   lines. An example of a multiple alignment in FASTA is provided below:

```
>d1qhka_ d.100.1.2 (A:) N-terminal domain of RNase HI...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
>gi|28379120|ref|NP_786012.1|:(2-47) ribonuclease H (putative)...
-KYYAVRKGRQ--P---G--------IYRTW---PETQKQVSGYP---QAQYKSFTSEKDAQDFMA
>gi|84386727|ref|ZP_00989753.1|:(2-47) hypothetical ribonuclease HI...
-KYYVVWKGRT--P---G--------IFTTW---NECKSQVDGFA---GARYKSFPTLGEAESAFG
>gi|116492108|ref|YP_803843.1|:(2-47) RNase H with double-stranded...
-KFYAVKKGRK--P---G--------LYLTW---DAAKQQVDGFA---GAVYKSFLTKAEAEEWMA
>gi|6323890|ref|NP_013961.1|:(1-47) Ribonuclease H1...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
```

   The package also contains the perl script blast2fa.pl to convert
   (PSI-)BLAST output to FASTA format. Please type `blast2fa.pl -h` for more
   information.

# Running a Docker container

   The docker image encapsulates software dependencies. 
   CUDA libraries, PSIPRED and BLAST software for incorporating 
   secondary structure prediction into the profile using `makepro.sh` 
   are included in the COMER2 docker image, so that the user can easily 
   run COMER2 on GPU(s) without setting up CUDA drivers and installing 
   external software before. 

   The COMER2 image (210MB, compressed) can be pulled from the Docker 
   repository using the command:

```
docker pull minmar/comer2
```

   The user's system is supposed to be equipped with one or more GPUs 
   with an installed NVIDIA graphics driver. 
   There is no requirement for CUDA drivers to be installed. 
   Yet, 
   [nvidia-docker2](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
   should be installed to enable the container to connect to the 
   graphics driver. 
   If the user's system already has a CUDA driver installed, it should 
   be at least version 10.1 or upgraded otherwise.

   The container for constructing a COMER2 profile can be used as 
   follows.
   Let's assume that the user's directory `~/myhostdata` contains an 
   MSA `myhostmsa.afa` in aligned FASTA format. 
   Then a profile can be constructed and written in the same directory 
   (`~/myhostdata/myhostmsa.pro`) using the following command:

```
docker run --rm --name=comer2 --gpus=all --user="$(id -u):$(id -g)" -ti \
   -v ~/myhostdata:/myhome \
   minmar/comer2 \
   bin/makepro.sh -v -i /myhome/myhostmsa.afa -o /myhome/myhostmsa.pro
```

   The first line specifies the container to be removed once finished 
   execution (the image is not removed).
   The host directory `~/myhostdata` is mounted on the container's 
   directory `/myhome` for data exchange in the second line (access 
   from within the container to the host system is not possible).
   The third line specifies the image name, and the fourth line 
   corresponds to the container's command with its arguments 
   (`makepro.sh` in this case).

   COMER2 search is initiated similarly:

```
docker run --rm --name=comer2 --gpus=all --user="$(id -u):$(id -g)" -ti \
   -v ~/myhostdata:/myhome \
   -v /data/databases/profile_dbs:/mydb \
   minmar/comer2 \
   bin/comer -v -i /myhome/myhostmsa.pro -d /mydb/pdb70_210331 -o /myhome/comer_output
```

   The main difference from the preceding command is the third line, 
   which mounts the directory of profile databases 
   `/data/databases/profile_dbs` in the host system on the container's 
   directory `/mydb`.

# Final Notes

   All executables in the COMER2 software package invoked with the "-h"
   option print a list of valid command-line options.

# References

Margelevicius, M. (2016) Bayesian nonparametrics in protein remote homology
search. Bioinformatics 32(18), 2744-2752.

Margelevicius, M. (2018) A low-complexity add-on score for protein remote
homology search with COMER. Bioinformatics 34(12), 2037-2045.

Margelevicius, M. (2019) Estimating statistical significance of local protein
profile-profile alignments. BMC Bioinformatics 20, 419.

Margelevicius, M. (2020) COMER2: GPU-accelerated sensitive and specific
homology searches. Bioinformatics 36(11), 3570-3572.

# Funding

The work was supported by the European Regional Development Fund 
[grant number 01.2.2-LMT-K-718-01-0028]

---

Contact: <mindaugas.margelevicius@bti.vu.lt>
