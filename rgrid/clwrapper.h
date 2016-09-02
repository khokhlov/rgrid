/**
 * \file
 * \brief Wrapper to OpenCL
 */

#ifndef CL_WRAPPER_H
#define CL_WRAPPER_H

#ifdef USE_OPENCL

#include <CL/opencl.h>
#include <vector>

namespace clwrapper {

/**
 * \brief Wrap all OpenCL stuff
 * 
 * Make all aviable OpenCL devices (GPUs) ready to use.
 * Each OpenCL platform (AMD, NVIDIA, Intel) corresponds to one OpenCL context.
 * Single platform can contain multiple devices.
 * Each device corresponds to one OpenCL command queue.
 */
class CLWrapper
{
public:
	/**
	 * \brief Get single instance of this class
	 * \details This is the only way to get CLWrapper
	 * \return instance of CLWrapper
	 */
	static const CLWrapper& instance()
	{
		static CLWrapper instance;
		return instance;
	}
	
	/**
	 * \brief Get OpenCL context
	 * \details In most cases there is only one OpenCL implementation on one machine,
	 * so there is only one platform
	 * \param[in] platformNum number of platform 
	 * \return context for specific platform with all devices
	 */
	const cl_context& getContext (unsigned platformNum = 0) const;
	
	/** 
	 * \brief Get OpenCL command queue
	 * \param[in] deviceNum number of device on specific platform
	 * \param[in] contextNum number of this specific context(platform) 
	 * \return OpenCL command queue for specific context (the same as platform number) and device
	 */
	const cl_command_queue& getCommandQueue (unsigned deviceNum = 0, unsigned contextNum = 0) const;
	
	/** 
	 * \brief Get OpenCL representation of device
	 * \param[in] deviceNum number of device on specific platform
	 * \param[in] platformNum number of this specific platform
	 * \return device for specific platform
	 */
	const cl_device_id& getDevice(unsigned deviceNum = 0, unsigned platformNum = 0) const;
		
	/** 
	 * \brief Get number of detected platforms
	 * \return number of platforms
	 */
	unsigned getPlatformsNum() const { return platformsNum; }
	
	/** 
	 * \brief Get number of detected devices on specific platform
	 * \param[in] platform
	 * \return number of devices for specific platform
	 */
	unsigned getDevicesNum(unsigned platform = 0) const;
	
private:
	
	/* 
	 * initialization of command queue and context
	 */
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
