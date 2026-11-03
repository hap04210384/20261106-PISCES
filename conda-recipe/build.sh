#!/bin/bash
set -euo pipefail

CXXFLAGS="${CXXFLAGS:-} -O3 -mavx2 -DNDEBUG"

${CXX} ${CXXFLAGS} -std=c++17 -o ${PREFIX}/bin/gbm_count  ${SRC_DIR}/src/gbm/gbm_count.cpp
${CXX} ${CXXFLAGS} -std=c++17 -o ${PREFIX}/bin/gbm_stats ${SRC_DIR}/src/gbm/gbm_stats.cpp
${CXX} ${CXXFLAGS} -std=c++17 -o ${PREFIX}/bin/gbm_sim    ${SRC_DIR}/src/gbm/gbm_sim.cpp
${CXX} ${CXXFLAGS} -std=c++17 -o ${PREFIX}/bin/gbm_approx ${SRC_DIR}/src/gbm/gbm_approx.cpp

# Data-preparation helpers
cp ${SRC_DIR}/src/gbm/h5_to_gbm1.py        ${PREFIX}/bin/
cp ${SRC_DIR}/src/gbm/pangbank_to_gbm1.py  ${PREFIX}/bin/
