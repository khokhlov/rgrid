#include "rgrid/clwrapper.h"

#ifdef USE_OPENCL

#include "rgrid/clutils.h"
#include "CL/opencl.h"

using namespace std;

// return platforms number
unsigned getPNum()
{
	unsigned pNumTest = 0;
	unsigned pNum = 0;
	while (pNumTest == pNum) {
		CHECK_CL_ERROR(clGetPlatformIDs(++pNumTest, NULL, &pNum));
	}
	return pNum;
}

// return devices number
unsigned getDNum(cl_platform_id platform, cl_device_type devType = CL_DEVICE_TYPE_ALL)
{
	unsigned dNumTest = 0;
	unsigned dNum = 0;
	while (dNumTest == dNum) {
		CHECK_CL_ERROR(clGetDeviceIDs(platform, devType, ++dNumTest, NULL, &dNum));
	}
	return dNum;
}

namespace clwrapper {

CLWrapper::CLWrapper()
{
	platformsNum = getPNum();
	platforms.resize(platformsNum);
	CHECK_CL_ERROR(clGetPlatformIDs(platformsNum, &(platforms[0]), NULL));
	
	devicesNum.resize(platformsNum);
	devices.resize(platformsNum);
	context.resize(platformsNum);
	cq.resize(platformsNum);
	for (unsigned i = 0; i != platformsNum; ++i) {
		cl_int errCode;
		// get devices for specified platform
		devicesNum[i] = getDNum(platforms[i]);
		devices[i].resize(devicesNum[i]);
		CHECK_CL_ERROR(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, devicesNum[i], &(devices[i][0]), NULL));
		
		// create context for all specified platforms
		cl_context_properties context_properties[3] = {
			CL_CONTEXT_PLATFORM, (cl_context_properties) platforms[i], 0};
		context[i] = clCreateContext(context_properties, devicesNum[i], &(devices[i][0]), NULL, NULL, &errCode);
		CHECK_CL_ERROR(errCode);
		
		// create command queues for all devices
		cq[i].resize(devicesNum[i]);
		for (unsigned j = 0; j != devicesNum[i]; ++j) {
			// TODO clCreateCommandQueue() is depreceated in OpenCL 2.0, 
			// use #define to set OpenCL version
			// and call clCreateCommandQueueWithProperties()
			cq[i][j] = clCreateCommandQueue(context[i], devices[i][j], 0, &errCode);
			CHECK_CL_ERROR(errCode);
		}
	}	
}

CLWrapper::~CLWrapper() {
	RG_ASSERT(context.size() == cq.size(), "Context and command queue sizes must be the same");
	for (unsigned i = 0; i != context.size(); ++i) {
		CHECK_CL_ERROR(clReleaseContext(context[i]));
		for (unsigned j = 0; j != cq[i].size(); ++j) {
			CHECK_CL_ERROR(clReleaseCommandQueue(cq[i][j]));
		}
	}
}

const cl_context& CLWrapper::getContext(unsigned platformNum) const
{
	RG_ASSERT(platformNum < context.size(), "Wrong platform number");
	return context[platformNum];
}

const cl_command_queue& CLWrapper::getCommandQueue(unsigned int deviceNum, unsigned int contextNum) const
{
	RG_ASSERT(contextNum < cq.size(), "Wrong context number");
	RG_ASSERT(deviceNum < cq[contextNum].size(), "Wrong device number");
	return cq[contextNum][deviceNum];
}

const cl_device_id& CLWrapper::getDevice(unsigned deviceNum, unsigned platformNum) const
{
	RG_ASSERT(platformNum < devices.size(), "Wrong context number");
	RG_ASSERT(deviceNum < devices[platformNum].size(), "Wrong device number");
	return devices[platformNum][deviceNum];
}

unsigned int CLWrapper::getDevicesNum(unsigned int platform) const
{
	RG_ASSERT(platform < devices.size(), "Wrong platform number");
	return devices[platform].size();
}

} // namespace clwrapper

#endif // USE_OPENCL
