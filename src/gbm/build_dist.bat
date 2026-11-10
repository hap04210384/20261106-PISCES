@echo off
rem Static-runtime CPU pack build (dist\cpu-windows-x64). DAPR.
cd /d %~dp0
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
if not exist ..\..\dist\cpu-windows-x64 mkdir ..\..\dist\cpu-windows-x64
cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:..\..\dist\cpu-windows-x64\gbm_count.exe gbm_count.cpp >build_dist.log 2>&1
if errorlevel 1 (echo COUNT-FAILED) else (echo COUNT-OK)
cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:..\..\dist\cpu-windows-x64\gbm_stats.exe gbm_stats.cpp >>build_dist.log 2>&1
if errorlevel 1 (echo STATS-FAILED) else (echo STATS-OK)
cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:..\..\dist\cpu-windows-x64\gbm_sim.exe gbm_sim.cpp >>build_dist.log 2>&1
if errorlevel 1 (echo SIM-FAILED) else (echo SIM-OK)
cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:..\..\dist\cpu-windows-x64\gbm_approx.exe gbm_approx.cpp >>build_dist.log 2>&1
if errorlevel 1 (echo APPROX-FAILED) else (echo APPROX-OK)
echo BUILD-DONE
