#!/bin/bash

echo "Building serial."
mkdir -v build
cd build
cmake ../
make clean
make

echo "Running tests"
./test_darray
./test_pdim

cd ..
echo "Building parallel."
mkdir -v build_mpi
cd build_mpi
cmake ../ -DUSE_MPI=1
make clean
make
