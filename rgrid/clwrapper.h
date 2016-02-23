#ifndef CL_WRAPPER_H
#define CL_WRAPPER_H

//#define USE_OPENCL

#ifdef USE_OPENCL

#include <CL/opencl.h>
#include <vector>

namespace clwrapper {

/*
 * Make all aviable devices ready to use
 * platforms <=> contexts
 * platform contains devices
 * devices <=> command queues
 */
class CLWrapper
{
public:
	static const CLWrapper& instance()
	{
		static CLWrapper instance;
		return instance;
	}
	
	// return context for specific platform with all devices
	const cl_context& getContext (unsigned platformNum = 0) const;
	
	// return command queue for specific context (the same as platform number) and device
	const cl_command_queue& getCommandQueue (unsigned deviceNum = 0, unsigned contextNum = 0) const;
	
	// return device for specific platform
	const cl_device_id& getDevice(unsigned deviceNum = 0, unsigned platformNum = 0) const;
		
	// return number of platforms
	unsigned getPlatformsNum() const { return platformsNum; }
	
	// return number of devices for specific platform
	unsigned getDevicesNum(unsigned platform = 0) const;
	
private:
	
	// initialization of command queue and context
	CLWrapper();	
	CLWrapper(const CLWrapper&);
	CLWrapper& operator=(const CLWrapper&);
	~CLWrapper();
	
	unsigned platformsNum;
	std::vector<unsigned> devicesNum;
	std::vector<cl_platform_id> platforms;
	std::vector<std::vector<cl_command_queue> > cq;
	std::vector<std::vector<cl_device_id> > devices;
	std::vector<cl_context> context;
};

}

#endif // USE_OPENCL

#endif // CL_WRAPPER_H
