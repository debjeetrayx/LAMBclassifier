@echo off
rem LAMB classifier launcher (Windows). Double-click to choose a file or folder,
rem or drop a file or folder onto this icon.
rem Written by LAMBclassifier::lamb_install_launcher(). Research use only.
setlocal
set "RSCRIPT={{RSCRIPT}}"
set "LAMB_LAUNCHER=1"

echo LAMB classifier - research use only, not a clinical diagnosis.
echo.

if not exist "%RSCRIPT%" (
    echo R was not found at:
    echo   %RSCRIPT%
    echo R may have been updated or moved. Open R and run:
    echo   LAMBclassifier::lamb_install_launcher^(^)
    echo.
    pause
    exit /b 1
)

if "%~1"=="" (
    "%RSCRIPT%" -e "invisible(suppressWarnings(suppressPackageStartupMessages(library(LAMBclassifier)))); lamb_run()"
) else (
    set "LAMB_INPUT=%~1"
    "%RSCRIPT%" -e "invisible(suppressWarnings(suppressPackageStartupMessages(library(LAMBclassifier)))); lamb_run(Sys.getenv('LAMB_INPUT'))"
)
set "STATUS=%ERRORLEVEL%"

echo.
if not "%STATUS%"=="0" (
    echo LAMB classifier stopped with an error ^(see the message above^).
    echo Nothing was changed in your data.
) else (
    echo Done.
)
pause
exit /b %STATUS%
