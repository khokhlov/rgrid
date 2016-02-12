/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_DARRAY_H
#define RGRID_DARRAY_H

#include <vector>
#include <iostream>

#include "rgrid/pdim.h"
#include "rgrid/utils.h"

#include "rgrid/clwrapper.h"
#include "rgrid/clutils.h"

namespace rgrid {

template <typename T, typename I>
class DArray : public PDim<I> {
	public:		
		DArray() : PDim<I>(),
#ifdef USE_OPENCL
			clContextIsSet(false), 
			clBufInitialized(false), 
			clOnDevice(false), 
#endif
			nc(1) {}
		
		void alloc() { alloc(static_cast<I>(1)); }
		void alloc(const I &nc);

		void fillGhost(const CartDir &d, const CartSide &s);
		void fillGhost(const CartDir &d) { fillGhost(d, SIDE_LEFT); fillGhost(d, SIDE_RIGHT); }
		void fillGhost() { fillGhost(X); fillGhost(Y); fillGhost(Z); }

		void fill(const T &v);

		T &operator[](const I index) { return data[index]; }

		T &val(const I i, const I j, const I k, const I cn) { return data[this->ind(i, j, k) + cn * this->localSizeGhost()]; }
		T &val(const I i, const I j, const I cn) { return val(i, j, static_cast<T>(0), cn); }
		T &val(const I i, const I cn) { return val(i, static_cast<T>(0), cn); }
		T &val(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3, const I cn) { return data[this->ind(i1, d1, i2, d2, i3, d3) + cn * this->localSizeGhost()]; }

		T &operator()(const I i, const I j, const I k, const I cn) { return val(i, j, k, cn); }
		T &operator()(const I i, const I j, const I cn) { return (*this)(i, j, static_cast<I>(0), cn); }
		T &operator()(const I i, const I cn) { return (*this)(i, static_cast<I>(0), cn); }

		const T &operator()(const I i, const I j, const I k, const I cn) const { return data[this->ind(i, j, k) + cn * this->localSizeGhost()]; }
		const T &operator()(const I i, const I j, const I cn) const { return (*this)(i, j, static_cast<I>(0), cn); }
		const T &operator()(const I i, const I cn) const { return (*this)(i, static_cast<I>(0), cn); }
#ifdef USE_OPENCL
	public:
		/* Set OpenCL context to work with buffer */
		void setCLContext(cl_context context);
		
		/* Set OpenCL command queue */
		void setCLCQ(cl_command_queue cq);
		
		/* Get OpenCL data buffer to work in kernels. 
		 * You must run setCLContext() and setCLCQ() before */
		cl_mem &getCLBuffer();
		
		/* copy data from host to device */
		void clHtoD();
		
		/* copy data from device to host */
		void clDtoH();
		
	private:
		void clInitBuffer();
		
		/* OpenCL context to work with buffer */
		cl_context clContext;
		
		/* OpenCL command queue to work with buffer */
		cl_command_queue clCQ;
		
		/* Data array on device */
		cl_mem clBuffer;
		
		bool clContextIsSet;
		bool clCQIsSet;
		bool clBufInitialized;

		/* true - data on device, false - data on host */
		bool clOnDevice;
#endif		
	public:	
		/* Data. */
		std::vector<T> data;

		/* Number of components. */
		I nc;
}; // DArray

template <typename T, typename I>
void DArray<T, I>::alloc(const I &nc)
{
	this->nc = nc;
	data.resize(this->localSizeGhost() * nc);
}

template <typename T, typename I>
void DArray<T, I>::fill(const T &v)
{
	for (I i = 0; i < this->localSizeGhost() * nc; i++) { data[i] = v; }
}

template <typename T, typename I>
void DArray<T, I>::fillGhost(const CartDir &d, const CartSide &s)
{
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	I k = 0;
	I sign;
	if (this->isOnFace(d, s)) {
		if (s == SIDE_LEFT) {
			k = 0;
			sign = 1;
		} else {
			k = this->localSize(d) - 1;
			sign = -1;
		}
	} else {
		return;
	}
	for (I cn = 0; cn < nc; cn++) {
		for (I i = 0; i < this->localSize(ort1); i++) {
			for (I j = 0; j < this->localSize(ort2); j++) {
				for (I gs = 1; gs <= this->ghost(d); gs++) {
					this->val(i, ort1, j, ort2, k - sign * gs, d, cn) = this->val(i, ort1, j, ort2, k, d, cn);
				}
			}
		}
	}
}

#ifdef USE_OPENCL

template <typename T, typename I>
void DArray<T, I>::setCLContext(cl_context context) {
	this->clContext = context;
	clContextIsSet = true;
	clBufInitialized = false;
}

template <typename T, typename I>
void DArray<T, I>::setCLCQ(cl_command_queue cq) {
	this->clCQ = cq;
	clCQIsSet = true;
	clBufInitialized = false;
}

template <typename T, typename I>
void DArray<T, I>::clInitBuffer() {
	RG_ASSERT(clContextIsSet, "Call setCLContext() first");
	RG_ASSERT(clCQIsSet, "Call setCLCQ() first");
	cl_int err_code;
	this->clBuffer = clCreateBuffer(this->clContext, CL_MEM_READ_WRITE, data.size() * sizeof(T), NULL, &err_code);
	CHECK_CL_ERROR(err_code);
	clBufInitialized = true;
}

template <typename T, typename I>
cl_mem& DArray<T, I>::getCLBuffer() {
	if (!clBufInitialized) clInitBuffer();
	return this->clBuffer;
}

template <typename T, typename I>
void DArray<T, I>::clHtoD() {
	if (!clOnDevice) {
		if (!clBufInitialized) clInitBuffer();
		RG_ASSERT(clCQIsSet, "Call setCLCQ() first");
		CHECK_CL_ERROR(clEnqueueWriteBuffer(clCQ, clBuffer, CL_TRUE, 0, data.size() * sizeof(T), &(data[0]), 0, NULL, NULL));
		clOnDevice = true;
	}
}

template <typename T, typename I>
void DArray<T, I>::clDtoH() {
	if (clOnDevice) {
		RG_ASSERT(clCQIsSet, "Call setCLCQ() first");
		CHECK_CL_ERROR(clEnqueueReadBuffer(clCQ, clBuffer, CL_TRUE, 0, data.size() * sizeof(T), &(data[0]), 0, NULL, NULL));
		clOnDevice = false;
	}
}

#endif // USE_OPENCL

}; // rgrid

#endif // RGRID_DARRAY_H
