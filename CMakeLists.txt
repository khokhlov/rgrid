cmake_minimum_required(VERSION 2.6)
project(RGRID)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)
add_definitions(-g -W -Wall -Wl,--no-undefined -fpermissive)
include_directories(./)

aux_source_directory(rgrid/src RGRID)

if(USE_MPI)
set(CMAKE_CXX_COMPILER mpicxx)
add_definitions(-DUSE_MPI)
endif(USE_MPI)

if(USE_OPENCL)
add_definitions(-DUSE_OPENCL)
find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenCL_LIBRARIES})
endif(USE_OPENCL)

add_library(RGrid SHARED ${RGRID})
add_library(RGridStatic STATIC ${RGRID})

add_executable(test_darray tests/test_darray.cpp)
target_link_libraries(test_darray RGrid ${EXTRA_LIBS})

add_executable(test_darraycontainer tests/test_darraycontainer.cpp)
target_link_libraries(test_darraycontainer RGrid ${EXTRA_LIBS})

add_executable(test_darrayscatter tests/test_darrayscatter.cpp)
target_link_libraries(test_darrayscatter RGrid ${EXTRA_LIBS})

add_executable(test_pdim tests/test_pdim.cpp)
target_link_libraries(test_pdim RGrid ${EXTRA_LIBS})

add_executable(test_clwrapper tests/test_clwrapper.cpp)
target_link_libraries(test_clwrapper RGrid ${EXTRA_LIBS})

