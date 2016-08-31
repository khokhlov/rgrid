#ifndef CL_UTILS_H
#define CL_UTILS_H

#include "rgrid/utils.h"
#include <iostream>

#ifdef USE_OPENCL

#include <CL/opencl.h>

namespace clwrapper
{

const char* clErrorToStr(cl_int code);

} // namespace clwrapper

#define CHECK_CL_ERROR(param) {\
	cl_int rc = param;\
	if (rc != CL_SUCCESS) {\
		RG_ASSERT(rc == CL_SUCCESS, clwrapper::clErrorToStr(rc))\
	}\
}

#endif // USE_OPENCL

#endif // CL_UTILS_H
