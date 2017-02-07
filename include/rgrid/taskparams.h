#ifndef TASK_PARAMS_H
#define TASK_PARAMS_H

#include "rgrid/types.h"

namespace rgrid {
	
template <typename T, typename I>
struct TaskParams {
	Dim3D<T> pml_d_max;	
	Dim3D<T> space_step;
	T time_step;
	I half_order;
};

} // namespace rgrid

#endif
