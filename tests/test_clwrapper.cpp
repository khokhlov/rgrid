#define CATCH_CONFIG_MAIN


#include "catch.hpp"
#include "rgrid/clwrapper.h"
#include "rgrid/darray.h"

#include <iostream>

using namespace std;

const char* kernelAdd10src = 
"__kernel void add10(int nx, int ny, int nz, int nc, __global const int* in, __global int* out) {"
"	const int i = get_global_id(0);"
"	const int j = get_global_id(1);"
"	const int k = get_global_id(2);"
"	if (i >= nx || j >= ny || k >= nz) return;"
"	for (int cn = 0; cn != nc; ++cn)"
"		out[i + j * nx + k * nx * ny + cn * nx * ny * nz] = in[i + j * nx + k * nx * ny + cn * nx * ny * nz] + 10;"
"}";

TEST_CASE(
		"CLWrapper",
		"oneDeviceDArray"
	 )
{
	const clwrapper::CLWrapper& clw = clwrapper::CLWrapper::instance();
	if (clw.getPlatformsNum() != 0) {
		if (clw.getDevicesNum() != 0) {
			rgrid::DArray<int, int> da, da2;
			da.resize(100, 200, 150);
			da.alloc(3);
			da.fill(7);
			da2.resize(100, 200, 150);
			da2.alloc(3);
			da2.fill(17);
			
			da.setCLContext(clw.getContext());
			da.setCLCQ(clw.getCommandQueue());
			da.clHtoD();
			
			cl_int err;
			cl_mem buffer1 = da.getCLBuffer();
			cl_mem buffer2 = da.getCLBuffer2();
			CHECK_CL_ERROR(err);
			cl_program program = clCreateProgramWithSource(clw.getContext(), 1, &kernelAdd10src, NULL, &err);
			CHECK_CL_ERROR(err);
			err = clBuildProgram(program, 0, NULL, "", NULL, NULL);
			// NOTE to get compiler output call here clGetProgramBuildInfo
			CHECK_CL_ERROR(err); // check clBuildProgram
			cl_kernel kernelAdd10 = clCreateKernel(program, "add10", &err);
			CHECK_CL_ERROR(err);
			
			int size[3] = {da.localSizeGhost(rgrid::X), da.localSizeGhost(rgrid::Y), da.localSizeGhost(rgrid::Z)};
			int nc = da.getNC();
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 0, sizeof(size[0]), &size[0]));
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 1, sizeof(size[1]), &size[1]));
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 2, sizeof(size[2]), &size[2]));
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 3, sizeof(nc), &nc));
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 4, sizeof(buffer1), &buffer1));
			CHECK_CL_ERROR(clSetKernelArg(kernelAdd10, 5, sizeof(buffer2), &buffer2));
			
			size_t lws[3] = {4, 4, 4}; // local work size
			size_t gws[3]; // global work size
			for (int i = 0; i != 3; ++i) {
				gws[i] = (size[i] + lws[i] - 1) / lws[i] * lws[i];
			}
			/* 
			 * NOTE lws can be NULL(it will be chosen automatically), 
			 * so you can use size[3] instead of gws[3]
			 * The only bad thing, it will be hard to use shared memory without fixed lws
			 */
			CHECK_CL_ERROR(clEnqueueNDRangeKernel(clw.getCommandQueue(), kernelAdd10, 3, NULL, gws, lws, 0, NULL, NULL));
			da.swapCLBuffers();
			
			CHECK_CL_ERROR(clReleaseKernel(kernelAdd10));
			CHECK_CL_ERROR(clReleaseProgram(program));
			
			da.clDtoH();
			
			REQUIRE(da == da2);
		}
	}
}
