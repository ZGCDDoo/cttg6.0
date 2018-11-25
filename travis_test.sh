#!/bin/bash

cd build && cmake .. && make -j 4 \
&& ./IntegratorTests            \
&& ./MarkovChainTests           \
&& ./MatrixTests                \
&& ./SelfConsistencyTests       \
&& ./ObservablesTests           \
&& ./FourierTests               \
&& ./Fourier_DCATests


cd .. && mkdir build_mpi && cd build_mpi && cmake -DMPI_BUILD=ON .. && make