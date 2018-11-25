#!/bin/bash

mkdir build && cd build && cmake .. && make  \
&& ./IntegratorTests            \
&& ./MarkovChainTests           \
&& ./MatrixTests                \
&& ./SelfConsistencyTests       \
&& ./ObservablesTests           \
&& ./FourierTests               \
&& ./Fourier_DCATests


cd .. && mkdir build_mpi && cd build_mpi && cmake -DMPI_BUILD=ON .. && make