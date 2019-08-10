@echo off

for /F %%i in ("%0") do @set basename=%%~nxi
for /F %%i in ("%0") do @set dirname=%%~dpi
for %%i in ("%dirname:~0,-1%") do set updirname=%%~dpi

set "INSTPSIPRED=C:\Users\Public\Installed_packages\psipred"
set "INSTBLASTP=C:\Users\Public\Installed_packages\ncbi-blast\ncbi-blast-2.2.23+"
set INSTCOMER=%updirname%
set INPUT=
set OUTPUT=
set OPTIONS=
set NOSELECT=--select
set NOFLUSH=
set VERBOSE=

:optionsloop
if not "%1"=="" (
    if "%1"=="-i" (set INPUT=%2& shift)
    if "%1"=="-o" (set OUTPUT=%2& shift)
    if "%1"=="-p" (set "OPTIONS=-p %2"& shift)
    if "%1"=="-s" (set NOSELECT=)
    if "%1"=="-r" (set NOFLUSH=--noflush)
    if "%1"=="-P" (set INSTPSIPRED=%2& shift)
    if "%1"=="-B" (set INSTBLASTP=%2& shift)
    if "%1"=="-v" (set VERBOSE=-v)
    if "%1"=="-h" (call :printusage & exit /b)
    shift
    goto :optionsloop
)


if "%OUTPUT%"=="" (echo ERROR: No output filename given. & exit /b)
if "%INPUT%"=="" (echo ERROR: No input filename given. & exit /b)
if not exist %INPUT% (echo ERROR: Input file not found: %INPUT% & exit /b)

if not exist %INSTPSIPRED% (echo ERROR: Directory does not exist: %INSTPSIPRED% & exit /b)
if not exist %INSTBLASTP% (echo ERROR: Directory does not exist: %INSTBLASTP% & exit /b)
if not exist "%INSTCOMER%" (echo ERROR: Directory does not exist: "%INSTCOMER%" & exit /b)

set makepro=%INSTCOMER%bin\makepro.exe
set ssp=%INSTCOMER%bin\ssp2.pl
set inssp=%INSTCOMER%bin\inssp2.pl
set SSFILE=%OUTPUT%.ss

if not exist "%makepro%" (echo ERROR: Executable `makepro' from the `comer' package not found: "%makepro%" & exit /b)
if not exist "%ssp%" (echo ERROR: A perl script from the `comer' package not found: "%ssp%" & exit /b)
if not exist "%inssp%" (echo ERROR: A perl script from the `comer' package not found: "%inssp%" & exit /b)

set cmd="%makepro%" %VERBOSE% -i %INPUT% -o %OUTPUT% %OPTIONS%
if not "%VERBOSE%"=="" (echo %cmd%)
%cmd%  ||  exit /b

set cmd=perl "%ssp%" --in %INPUT% --out %SSFILE% %OPTIONS% %NOSELECT% %NOFLUSH% ^
    --psipred %INSTPSIPRED% --blast %INSTBLASTP% --comer "%INSTCOMER%"
if not "%VERBOSE%"=="" (echo %cmd%)
%cmd%  ||  exit /b

set cmd=perl "%inssp%" --in %SSFILE% --to %OUTPUT%
if not "%VERBOSE%"=="" (echo %cmd%)
%cmd%  ||  exit /b

echo.
if "%NOFLUSH%"=="" (del %SSFILE%)


exit /b


:printusage
echo.
echo Make `comer' profile including secondary structure predictions.
echo 2013-2019^(C^)Mindaugas Margelevicius,VU IBT,Vilnius
echo.
echo %basename% ^<Parameters^>
echo.
echo Parameters:
echo.
echo -i ^<input^>     Input multiple alignment file either in
echo                FASTA or in STOCKHOLM.
echo -o ^<output^>    Output filename of profile.
echo -p ^<options^>   Input file of options;
echo                By default, the file in installation
echo                directory of this package is searched.
echo -s             Do not preprocess input multiple alignment before
echo                predicting SS.
echo -r             Do not flush temporary files.
echo -P ^<path^>      Installation path of `PSIPRED'
echo        default=%INSTPSIPRED%
echo -B ^<path^>      Installation path of `BLAST+'
echo        default=%INSTBLASTP%
echo -v             Enable warnings.
echo -h             This text.
echo.
exit /b
