@echo off
rem E-B02/E-B03/E-B04/E-B05 real-data pipeline for K. pneumoniae pangenome (id 11756)
rem Runs: h5 -> GBM1 matrix -> bitmap indexing -> all-pairs exact count (minc=5, depletion rule) -> stats (Fisher + BH-FDR)
rem All experiment stages execute in C++ Release builds (AVX2). Python is used only for data conversion.
setlocal
set ROOT=E:\20261031-paper-Tensor core solves big problems\Manuscript_work
set H5=%ROOT%\data\raw\pangbank\klebsiella_pneumoniae_refseq.h5
set GBM1=%ROOT%\data\processed\kpn_gbm1.txt
set PAIRS=%ROOT%\data\processed\kpn_pairs.tsv
set CAND=%ROOT%\tables\tb_eb05_candidate_table.tsv
set PY=C:\Users\Lenovo\AppData\Roaming\kimi-desktop\daimon-share\daimon\runtime\python\cpython-3.12\python.exe

cd /d %ROOT%\code\src\gbm

echo [1/4] h5 -> GBM1
%PY% h5_to_gbm1.py "%H5%" "%GBM1%" || exit /b 1

echo [2/4] all-pairs exact count (minc=5, hypergeometric deviation screen z>=3.29 per amendment 4)
gbm_count.exe --count "%GBM1%" "%PAIRS%" 5 screen || exit /b 1

echo [3/4] stats (Fisher two-sided + BH-FDR, both screened-m and allpairs-m conservative)
gbm_stats.exe --pairs "%PAIRS%" "%CAND%" 54590 || exit /b 1

echo [4/4] done
dir "%CAND%"
endlocal
