#!/bin/bash

echo "Running mpi tests"
mpirun -np 12 ./test_darrayscatter
mpirun -np 2 ./test_darrayscatter
mpirun -np 1 ./test_darrayscatter
mpirun -np 60 ./test_darrayscatter
mpirun -np 6 ./test_darrayscatter
mpirun -np 8 ./test_darrayscatter
