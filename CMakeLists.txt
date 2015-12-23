cmake_minimum_required(VERSION 2.6)
project(RGRID)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)
add_definitions(-O2 -g -W -Wall -Wl,--no-undefined)
include_directories(./)

aux_source_directory(rgrid/src RGRID)

set(USE_MPI 1)
if(USE_MPI)

set(CMAKE_CXX_COMPILER mpicxx)
add_definitions(-DUSE_MPI)

endif(USE_MPI)

add_library(rgrid SHARED ${RGRID})

add_executable(test_darray tests/test_darray.cpp)
target_link_libraries(test_darray rgrid)

add_executable(test_pdim tests/test_pdim.cpp)
target_link_libraries(test_pdim rgrid)
