/**
 * \file
 * \brief Useful util to work with for OpenCL
 */

#ifndef CL_UTILS_H
#define CL_UTILS_H

#include "rgrid/utils.h"
#include <iostream>

#ifdef USE_OPENCL

#include <CL/opencl.h>

namespace clwrapper
{

/// Convert OpenCL error code to readable string
const char* clErrorToStr(cl_int code);

} // namespace clwrapper

/// Check for errors in execution of OpenCL functions
#define CHECK_CL_ERROR(param) {\
	cl_int rc = param;\
	if (rc != CL_SUCCESS) {\
		RG_ASSERT(rc == CL_SUCCESS, clwrapper::clErrorToStr(rc))\
	}\
}

#endif // USE_OPENCL

#endif // CL_UTILS_H
