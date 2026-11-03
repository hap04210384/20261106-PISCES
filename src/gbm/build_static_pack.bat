@echo off
rem static-runtime rebuild of all portable exes into pack bin dir
cd /d %~dp0
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
set OUT=%~dp0..\..\..\xgpu_run\bin
if not exist "%OUT%" mkdir "%OUT%"
set NVCC="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8\bin\nvcc.exe"
set GEN=-gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90 -gencode=arch=compute_120,code=sm_120 -gencode=arch=compute_120,code=compute_120

cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:"%OUT%\gbm_count.exe" gbm_count.cpp >build_mt.log 2>&1
if errorlevel 1 (echo COUNT-FAILED & type build_mt.log) else (echo COUNT-OK)
cl /O2 /EHsc /MT /arch:AVX2 /DNDEBUG /std:c++17 /Fe:"%OUT%\gbm_stats.exe" gbm_stats.cpp >>build_mt.log 2>&1
if errorlevel 1 (echo STATS-FAILED & type build_mt.log) else (echo STATS-OK)

%NVCC% -O3 -DNDEBUG -std=c++17 %GEN% -Xcompiler "/EHsc /MT" -cudart=static -o "%OUT%\gbm_router.exe" gbm_router.cu >build_mt2.log 2>&1
if errorlevel 1 (echo ROUTER-FAILED & type build_mt2.log) else (echo ROUTER-OK)
%NVCC% -O3 -DNDEBUG -std=c++17 %GEN% -Xcompiler "/EHsc /MT" -cudart=static -o "%OUT%\gbm_pair_gpu.exe" gbm_pair_gpu.cu >>build_mt2.log 2>&1
if errorlevel 1 (echo PAIR-FAILED & type build_mt2.log) else (echo PAIR-OK)
%NVCC% -O3 -DNDEBUG -std=c++17 %GEN% -Xcompiler "/EHsc /MT" -cudart=static -o "%OUT%\gpu_probe.exe" gpu_probe.cu >>build_mt2.log 2>&1
if errorlevel 1 (echo PROBE-FAILED & type build_mt2.log) else (echo PROBE-OK)
