#!/bin/bash

echo "Building serial."
mkdir -v -p build
cd build
cmake ../
# make clean
make

cd ..
echo "Building parallel."
mkdir -v -p build_mpi
cd build_mpi
cmake ../ -DUSE_MPI=1
# make clean
make
 
cd ..
echo "Building with OpenCL support"
mkdir -p build_opencl
cd build_opencl
#make clean
cmake ../ -DUSE_OPENCL=1
make

cd ..
echo "Building with OpenCL & MPI support"
mkdir -p build_opencl_mpi
cd build_opencl_mpi
make clean
cmake ../ -DUSE_OPENCL=1 -DUSE_MPI=1
make

echo "Running tests"
./test_darray
./test_pdim
./test_rgio
./test_darraycontainer
./test_clwrapper
mpirun -np 12 ./test_darrayscatter
mpirun -np 2 ./test_darrayscatter
