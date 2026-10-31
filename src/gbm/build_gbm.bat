@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /O2 /EHsc /arch:AVX2 /DNDEBUG /std:c++17 /Fe:gbm_count.exe gbm_count.cpp >build.log 2>&1
if errorlevel 1 (echo BUILD-FAILED & type build.log) else (echo BUILD-OK)
