@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8\bin\nvcc.exe" -O3 -DNDEBUG -std=c++17 -arch=sm_86 -Xcompiler /EHsc -o gbm_router.exe gbm_router.cu >build_router.log 2>&1
if errorlevel 1 (echo BUILD-FAILED & type build_router.log) else (echo BUILD-OK)
