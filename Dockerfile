# GBM reproducibility code - CPU reference build (AVX2, no CUDA required).
# Build:  docker build -t gbm .
# Run:    docker run --rm gbm            (builds engines + runs smoke test)
# Demo:   docker run --rm -it gbm bash   (engines in /opt/gbm/bin)
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        g++ python3 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/gbm
COPY src/gbm/gbm_count.cpp src/gbm/gbm_stats.cpp src/gbm/gbm_sim.cpp src/gbm/gbm_approx.cpp \
     src/gbm/h5_to_gbm1.py src/gbm/pangbank_to_gbm1.py /opt/gbm/src/
COPY tests/smoke_test.py /opt/gbm/tests/smoke_test.py

RUN mkdir -p bin \
 && g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o bin/gbm_count  src/gbm_count.cpp \
 && g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o bin/gbm_stats src/gbm_stats.cpp \
 && g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o bin/gbm_sim    src/gbm_sim.cpp \
 && g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o bin/gbm_approx src/gbm_approx.cpp

RUN python3 tests/smoke_test.py bin

CMD ["python3", "tests/smoke_test.py", "bin"]
