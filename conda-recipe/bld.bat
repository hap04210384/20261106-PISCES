@echo off
set CXXFLAGS=%CXXFLAGS% -O3 -mavx2 -DNDEBUG

%CXX% %CXXFLAGS% /std:c++17 /Fe:%LIBRARY_BIN%\gbm_count.exe  %SRC_DIR%\src\gbm\gbm_count.cpp
%CXX% %CXXFLAGS% /std:c++17 /Fe:%LIBRARY_BIN%\gbm_stats.exe %SRC_DIR%\src\gbm\gbm_stats.cpp
%CXX% %CXXFLAGS% /std:c++17 /Fe:%LIBRARY_BIN%\gbm_sim.exe    %SRC_DIR%\src\gbm\gbm_sim.cpp
%CXX% %CXXFLAGS% /std:c++17 /Fe:%LIBRARY_BIN%\gbm_approx.exe %SRC_DIR%\src\gbm\gbm_approx.cpp

copy %SRC_DIR%\src\gbm\h5_to_gbm1.py        %SCRIPTS%\
copy %SRC_DIR%\src\gbm\pangbank_to_gbm1.py  %SCRIPTS%\
if errorlevel 1 exit 1
