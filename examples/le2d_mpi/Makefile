# Steps to build:
# 1. replace ~/code/rgrid by path to rgrid directory
# 2. replace libRGridStatic.a by path to static library
CXXFLAGS = -DUSE_MPI -g -O2 -fmax-errors=5 -fpermissive -Wall -I ~/code/rgrid/ -lm
CXX = mpicxx

le2d: le_core.h le_core.cpp main.cpp libRGridStatic.a
	$(CXX) $(CXXFLAGS) -DUSE_SOA -DUSE_FLOAT  -o le2d_soa_f main.cpp le_core.cpp libRGridStatic.a
	$(CXX) $(CXXFLAGS) -DUSE_SOA -DUSE_DOUBLE -o le2d_soa_d main.cpp le_core.cpp libRGridStatic.a
	$(CXX) $(CXXFLAGS) -DUSE_AOS -DUSE_FLOAT  -o le2d_aos_f main.cpp le_core.cpp libRGridStatic.a
	$(CXX) $(CXXFLAGS) -DUSE_AOS -DUSE_DOUBLE -o le2d_aos_d main.cpp le_core.cpp libRGridStatic.a

.PHONY: clean
clean:
	rm *.vtk le2d_*
