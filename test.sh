#!/bin/bash

cd build_mpi

echo "Running tests"
./test_darray
./test_pdim
./test_darraycontainer
./test_clwrapper
mpirun -np 12 ./test_darrayscatter
mpirun -np 2 ./test_darrayscatter
mpirun -np 1 ./test_darrayscatter
mpirun -np 60 ./test_darrayscatter
mpirun -np 6 ./test_darrayscatter
