#!/bin/bash

echo "Building serial."
mkdir -v -p build
cd build
cmake ../
# make clean
make -j 4

../test.sh

cd ..
echo "Building parallel."
mkdir -v -p build_mpi
cd build_mpi
cmake ../ -DUSE_MPI=1
# make clean
make -j 4
 
../test.sh

cd ..
echo "Building with OpenCL support"
mkdir -p build_opencl
cd build_opencl
#make clean
cmake ../ -DUSE_OPENCL=1
make -j 4

../test.sh

cd ..
echo "Building with OpenCL & MPI support"
mkdir -p build_opencl_mpi
cd build_opencl_mpi
#make clean
cmake ../ -DUSE_OPENCL=1 -DUSE_MPI=1
make -j 4

../test.sh
