@echo off
rem End-to-end demo on the bundled synthetic dataset (Windows).
rem Usage: demo\run_demo.bat [exe-dir]
setlocal
set HERE=%~dp0
set ROOT=%HERE%..
set BIN=%~1
if "%BIN%"=="" set BIN=%ROOT%\dist\cpu-windows-x64
if not exist "%BIN%\gbm_count.exe" set BIN=%ROOT%\src\gbm
if not exist "%HERE%out" mkdir "%HERE%out"

echo [1/3] exact all-pairs count + hypergeometric screen (minc=3)
"%BIN%\gbm_count.exe" --count "%ROOT%\data\demo\demo_matrix.txt" "%HERE%out\demo_pairs.tsv" 3 screen || exit /b 1

echo [2/3] Fisher exact test per screened pair
"%BIN%\gbm_stats.exe" --pairs "%HERE%out\demo_pairs.tsv" "%HERE%out\demo_stats.tsv" || exit /b 1

echo [3/3] check against planted ground truth
python "%HERE%check_demo.py" "%HERE%out\demo_stats.tsv"
endlocal
