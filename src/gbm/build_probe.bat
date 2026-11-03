@echo off
cd /d %~dp0
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
set NVCC="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8\bin\nvcc.exe"
set GEN=-gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_89,code=sm_89 -gencode=arch=compute_90,code=sm_90 -gencode=arch=compute_120,code=sm_120 -gencode=arch=compute_120,code=compute_120
%NVCC% -O3 -DNDEBUG -std=c++17 %GEN% -Xcompiler /EHsc -o gpu_probe.exe gpu_probe.cu >build_probe.log 2>&1
if errorlevel 1 (echo PROBE-FAILED & type build_probe.log) else (echo PROBE-OK)
