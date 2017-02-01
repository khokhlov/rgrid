cmake_minimum_required(VERSION 3.2)
project(rgrid)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

aux_source_directory(src RGRID_SRC)
aux_source_directory(include/rgrid RGRID_INCLUDE)

# specify headers in include to make them visible in some IDEs
add_library(RGrid ${RGRID_SRC} ${RGRID_INCLUDE})

target_compile_options(RGrid PUBLIC -rdynamic -W -Wall -Wl,--no-undefined -fpermissive)

target_include_directories(RGrid PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # hardcoded for build
	$<INSTALL_INTERFACE:include>) # relative path for installed

target_include_directories(RGrid PRIVATE src) # headers inside src won't be used outside

if(USE_MPI)
	find_package(MPI REQUIRED)
	target_compile_definitions(RGrid PRIVATE USE_MPI)
	target_include_directories(RGrid PRIVATE ${MPI_CXX_INCLUDE_PATH})
	target_link_libraries(RGrid PRIVATE ${MPI_CXX_LIBRARIES})
	target_compile_options(RGrid PRIVATE ${MPI_CXX_COMPILE_FLAGS})
	set_property(TARGET RGrid APPEND PROPERTY LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
endif(USE_MPI)

if(USE_OPENCL)
	find_package(OpenCL REQUIRED)
	target_compile_definitions(RGrid PRIVATE USE_OPENCL)
	target_include_directories(RGrid PRIVATE ${OpenCL_INCLUDE_DIRS})
	target_link_libraries(RGrid PRIVATE ${OpenCL_LIBRARIES})
endif(USE_OPENCL)

# add installation for target
install(TARGETS RGrid EXPORT RGridConfig
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin # This is for Windows
	INCLUDES DESTINATION include)

# create export file which can be imported by other cmake projects
install(EXPORT RGridConfig DESTINATION lib/cmake/RGrid)

# copy public include files
install(DIRECTORY include/ DESTINATION include)

# This makes the project importable from the build directory
export(TARGETS RGrid FILE rgridConfig.cmake)

# compile fdcoef generator
add_executable(FDCGen generators/fdcgen.cpp)
target_include_directories(FDCGen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/fdc.cpp
	COMMAND FDCGen >${CMAKE_CURRENT_SOURCE_DIR}/src/fdc.cpp
	DEPENDS FDCGen)

target_sources(RGrid PRIVATE src/fdc.cpp)

# Adding tests

add_executable(test_darray tests/test_darray.cpp)
target_link_libraries(test_darray RGrid)

add_executable(test_darraycontainer tests/test_darraycontainer.cpp)
target_link_libraries(test_darraycontainer RGrid)

add_executable(test_darrayscatter tests/test_darrayscatter.cpp)
target_link_libraries(test_darrayscatter RGrid)

add_executable(test_pdim tests/test_pdim.cpp)
target_link_libraries(test_pdim RGrid)

add_executable(test_clwrapper tests/test_clwrapper.cpp)
target_link_libraries(test_clwrapper RGrid)

