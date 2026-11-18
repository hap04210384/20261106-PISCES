@echo off
rem ============================================================================
rem PISCES one-command pipeline: any feature x genome presence/absence table
rem -> exact all-pairs co-occurrence counts -> Fisher exact + BH-FDR candidates
rem
rem Usage:  pisces_run.bat <input> <outdir> [minc]
rem   input   presence/absence table: roary gene_presence_absence.csv,
rem           Panaroo gene_presence_absence.Rtab, PPanGGOLiN/PanGBank HDF5,
rem           plain CSV/TSV matrix, or two-column long list
rem           (format auto-detected; see src/gbm/to_gbm1.py --help)
rem   outdir  output directory (created if missing)
rem   minc    minimum support per feature (default 5); pairs whose joint support
rem           cannot reach minc are never materialized
rem
rem Engines must be built first: src/gbm/build_gbm.bat (MSVC) or
rem   g++ -O3 -mavx2 -DNDEBUG -std=c++17 (see README, Quick start)
rem ============================================================================
setlocal enabledelayedexpansion
if "%~2"=="" (echo Usage: pisces_run.bat ^<input^> ^<outdir^> [minc] & exit /b 1)
set IN=%~1
set OUT=%~2
set MINC=%~3
if "%MINC%"=="" set MINC=5
set HERE=%~dp0
set GBM1=%OUT%\input.gbm1.txt
set PAIRS=%OUT%\pairs.tsv
set CAND=%OUT%\candidates.tsv
if not exist "%OUT%" mkdir "%OUT%"

where python >nul 2>nul
if errorlevel 1 (echo ERROR: python not found on PATH & exit /b 1)
if not exist "%HERE%src\gbm\gbm_count.exe" (echo ERROR: build the engines first ^(src\gbm\build_gbm.bat^) & exit /b 1)
if not exist "%HERE%src\gbm\gbm_stats.exe" (echo ERROR: build the engines first ^(src\gbm\build_gbm.bat^) & exit /b 1)

echo [1/3] encode input to GBM1 bitmap (deterministic, one pass)
python "%HERE%src\gbm\to_gbm1.py" "%IN%" "%GBM1%" || exit /b 1

rem GBM1 header line is: GBM1 ^<genomes N^> ^<families M^>
for /f "usebackq tokens=1,2,3" %%a in ("%GBM1%") do (set M=%%c & goto haveM)
:haveM
echo [2/3] exact all-pairs count (minc=%MINC%, AVX2 Release)
"%HERE%src\gbm\gbm_count.exe" --count "%GBM1%" "%PAIRS%" %MINC% || exit /b 1

echo [3/3] Fisher two-sided p + odds ratio + BH-FDR (screened and all-pairs)
"%HERE%src\gbm\gbm_stats.exe" --pairs "%PAIRS%" "%CAND%" %M% || exit /b 1

echo.
echo PISCES pipeline done. Candidate table: %CAND%
echo Columns: famA famB c11 c10 c01 c00 odds_ratio p_fisher padj_BH_screened padj_BH_allpairs rank
endlocal
