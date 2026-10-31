@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
cl /O2 /EHsc /arch:AVX2 /DNDEBUG /std:c++17 /Fe:gbm_approx.exe gbm_approx.cpp >build_approx.log 2>&1
if errorlevel 1 (echo BUILD-FAILED & type build_approx.log) else (echo BUILD-OK)
