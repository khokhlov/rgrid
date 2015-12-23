#!/bin/bash

mkdir -v build
cd build
cmake ../
make clean
make

echo "Running tests"
./test_darray
./test_pdim
