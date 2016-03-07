/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_DARRAY_H
#define RGRID_DARRAY_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

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
			
		~DArray();
		
		void alloc() { alloc(static_cast<I>(1)); }
		void alloc(const I &nc);
		
		void fillGhost(const CartDir &d, const CartSide &s);
		void fillGhost(const CartDir &d) { fillGhost(d, SIDE_LEFT); fillGhost(d, SIDE_RIGHT); }
		void fillGhost() { fillGhost(X); fillGhost(Y); fillGhost(Z); }
		
		// copy ghost nodes from "da"
		// "d" and "s" relative to this DArray
		void copyGhost(DArray<T, I>& da, const CartDir &d, const CartSide &s);

		void fill(const T &v);

		T &operator[](const I index) { return data[index]; }
		const T &operator[](const I index) const { return data[index]; }

		T &val(const I i, const I j, const I k, const I cn) { return data[this->ind(i, j, k) + cn * this->localSizeGhost()]; }
		T &val(const I i, const I j, const I cn) { return val(i, j, static_cast<T>(0), cn); }
		T &val(const I i, const I cn) { return val(i, static_cast<T>(0), cn); }
		T &val(const I &i1, const CartDir &d1, const I &i2, const CartDir &d2, const I &i3, const CartDir &d3, const I cn) { return data[this->ind(i1, d1, i2, d2, i3, d3) + cn * this->localSizeGhost()]; }

		T &operator()(const I i, const I j, const I k, const I cn) { return val(i, j, k, cn); }
		T &operator()(const I i, const I j, const I cn) { return (*this)(i, j, static_cast<I>(0), cn); }
		T &operator()(const I i, const I cn) { return (*this)(i, static_cast<I>(0), cn); }

		const T &operator()(const I i, const I j, const I k, const I cn) const { return data[PDim<I>::ind(i, j, k) + cn * PDim<I>::localSizeGhost()]; }
		const T &operator()(const I i, const I j, const I cn) const { return (*this)(i, j, static_cast<I>(0), cn); }
		const T &operator()(const I i, const I cn) const { return (*this)(i, static_cast<I>(0), cn); }
		
		I dataSize() const { return data.size(); }
		
		// get number of components
		I getNC() const { return nc; };
		
		// write line (all nodes on X direction) with coordinates y and z into stream
		void writeLine(std::iostream& stream, I y, I z, rgio::format fmt) const;
		
		// write line (all nodes on X direction) with coordinates y and z into stream
		void readLine(std::iostream& stream, I y, I z, rgio::format fmt);
		
#ifdef USE_OPENCL
	public:
		// the same as copy ghost, but it make copy on device
		void copyGhostCL(DArray<T, I>& da, const CartDir &d, const CartSide &s);
		
		// the same as fill ghost, but working on device
		void fillGhostCL(const CartDir &d, const CartSide &s);
		
		/* Set OpenCL context to work with buffer */
		void setCLContext(cl_context context);
		
		/* Set OpenCL command queue */
		void setCLCQ(cl_command_queue cq);
		
		/* Get OpenCL data buffer to work in kernels
		 * clHtoD() and clDtoH() uses this buffer 
		 * You must run setCLContext() and setCLCQ() before */
		cl_mem &getCLBuffer();
		
		/* Get second temporary buffer with the same size as first buffer to improve performance
		 * You must run setCLContext() and setCLCQ() before */
		cl_mem &getCLBuffer2();
		
		/* swap first and second buffers */
		void swapCLBuffers() { std::swap(clBuffer, clBuffer2); }
		
		/* copy data from host to device */
		void clHtoD();
		
		/* copy data from device to host */
		void clDtoH();
		
	private:
		void clInitBuffer();
		void clReleaseBuffer();
		
		/* OpenCL context to work with buffer */
		cl_context clContext;
		
		/* OpenCL command queue to work with buffer */
		cl_command_queue clCQ;
		
		/* Data array on device */
		cl_mem clBuffer, clBuffer2;
		
		bool clContextIsSet;
		bool clCQIsSet;
		bool clBufInitialized;

		/* true - data on device, false - data on host */
		bool clOnDevice;
#endif		
	private:	
		/* Data. */
		std::vector<T> data;

		/* Number of components. */
		I nc;
}; // DArray

template <typename T, typename I>
DArray<T, I>::~DArray() {
#ifdef USE_OPENCL
	clReleaseBuffer();
#endif
}

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
void DArray<T, I>::copyGhost(DArray<T, I>& da, const CartDir &d, const CartSide &s) {
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	I k = 0, kda = 0;
	I sign;
	if (s == SIDE_LEFT) {
		k = 0;
		kda = da.localSize(d) - 1;
		sign = 1;
	} else {
		k = PDim<I>::localSize(d) - 1;
		kda = 0;
		sign = -1;
	}
	for (I cn = 0; cn < nc; cn++) 
	for (I i = 0; i < PDim<I>::localSize(ort1); i++) 
	for (I j = 0; j < PDim<I>::localSize(ort2); j++) 
	for (I gs = 0; gs < PDim<I>::ghost(d); gs++) {
		val(i, ort1, j, ort2, k - sign * (gs + 1), d, cn) = da.val(i, ort1, j, ort2, kda - sign * gs, d, cn);
	}
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

template <typename T, typename I>
void DArray<T, I>::writeLine(std::iostream& stream, I y, I z, rgio::format fmt) const {
	for (I x = 0; x != PDim<I>::localSize(X); ++x) {
		for (I cn = 0; cn != nc; ++cn) {
			if (fmt == rgio::TEXT) {
				stream << operator()(x, y, z, cn) << " ";
			} else if (fmt == rgio::BINARY) {
				T d = operator()(x, y, z, cn);
				const char *dc = (const char*)(&d);
				stream.write(dc, sizeof(T));
			}
		}
	}
}

template <typename T, typename I>
void DArray<T, I>::readLine(std::iostream& stream, I y, I z, rgio::format fmt) {
	for (I x = 0; x != PDim<I>::localSize(X); ++x) {
		for (I cn = 0; cn != nc; ++cn) {
			if (fmt == rgio::TEXT) {
				stream >> operator()(x, y, z, cn);
			} else if (fmt == rgio::BINARY) {
				T d;
				char *dc = (char*)(&d);
				stream.read(dc, sizeof(T));
				operator()(x, y, z, cn) = d;
			}
		}
	}
}

template <typename T, typename I>
bool operator==(const DArray<T, I>& lhs, const DArray<T, I>& rhs) {
	if (&lhs == &rhs) return true;
	if (static_cast<PDim<I> >(lhs) != static_cast<PDim<I> >(rhs)) return false;
	if (lhs.getNC() != rhs.getNC()) return false;
	// don't check ghost!
	for (I c = 0; c != lhs.getNC(); ++c)
	for (I k = 0; k != lhs.localSize(Z); ++k)
	for (I j = 0; j != lhs.localSize(Y); ++j)
	for (I i = 0; i != lhs.localSize(X); ++i) {
		if (lhs(i, j, k, c) != rhs(i, j, k, c)) {
			return false;
		}
	}
	return true;
}

#ifdef USE_OPENCL

template <typename T, typename I>
void DArray<T, I>::setCLContext(cl_context context) {
	this->clContext = context;
	clContextIsSet = true;
	clReleaseBuffer();
}

template <typename T, typename I>
void DArray<T, I>::setCLCQ(cl_command_queue cq) {
	this->clCQ = cq;
	clCQIsSet = true;
	clReleaseBuffer();
}

template <typename T, typename I>
void DArray<T, I>::clInitBuffer() {
	RG_ASSERT(clContextIsSet, "Call setCLContext() first");
	RG_ASSERT(clCQIsSet, "Call setCLCQ() first");
	cl_int err;
	this->clBuffer = clCreateBuffer(this->clContext, CL_MEM_READ_WRITE, data.size() * sizeof(T), NULL, &err);
	CHECK_CL_ERROR(err);
	this->clBuffer2 = clCreateBuffer(this->clContext, CL_MEM_READ_WRITE, data.size() * sizeof(T), NULL, &err);
	CHECK_CL_ERROR(err);
	clBufInitialized = true;
}

template <typename T, typename I>
void DArray<T, I>::clReleaseBuffer() {
	if (clBufInitialized) {
		clReleaseMemObject(clBuffer);
		clReleaseMemObject(clBuffer2);
		clBufInitialized = false;
	}
}

template <typename T, typename I>
cl_mem& DArray<T, I>::getCLBuffer() {
	if (!clBufInitialized) clInitBuffer();
	return this->clBuffer;
}

template <typename T, typename I>
cl_mem& DArray<T, I>::getCLBuffer2() {
	if (!clBufInitialized) clInitBuffer();
	return this->clBuffer2;
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

template <typename T, typename I>
void DArray<T, I>::copyGhostCL(DArray<T, I>& da, const CartDir &d, const CartSide &s) {
	size_t origFrom[ALL_DIRS], origTo[ALL_DIRS];
	size_t origHost[ALL_DIRS] = { 0, 0, 0 }; 
	size_t region[ALL_DIRS];
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	origFrom[ort1] = origTo[ort1] = da.ghost(ort1);
	origFrom[ort2] = origTo[ort2] = da.ghost(ort2);
	region[d] = da.ghost(d);
	region[ort1] = da.localSize(ort1);
	region[ort2] = da.localSize(ort2);
	if (s == SIDE_LEFT) {
		origTo[d] = 0;
		origFrom[d] = (PDim<I>::localSizeGhost(d) - 2 * PDim<I>::ghost(d));
	} else {
		origTo[d] = (PDim<I>::localSizeGhost(d) - PDim<I>::ghost(d));
		origFrom[d] = da.ghost(d);
	}
	I regionSize = region[X] * region[Y] * region[Z];
	T* hostMem = new T[regionSize * nc];
	// add sizeof(T) multiplier in X direction
	origFrom[X] *= sizeof(T);
	origTo[X] *= sizeof(T);
	region[X] *= sizeof(T);
	// Create Sub buffers to work with components separately
	cl_mem* subBufFrom = new cl_mem[nc];
	cl_mem* subBufTo = new cl_mem[nc];
	for (I cn = 0; cn != nc; ++cn) {
		cl_int err;
		cl_buffer_region clbr;
		clbr.size = da.localSizeGhost() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBufFrom[cn] = clCreateSubBuffer(da.clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
		clbr.size = PDim<I>::localSizeGhost() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBufTo[cn] = clCreateSubBuffer(clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
	}
	// TODO make non blocking copy
	for (I cn = 0; cn != nc; ++cn) {
		CHECK_CL_ERROR(clEnqueueReadBufferRect(da.clCQ, subBufFrom[cn], CL_TRUE, origFrom, origHost, region, da.stride(Y) * sizeof(T), da.stride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, NULL));
	}
	for (I cn = 0; cn != nc; ++cn) {
		CHECK_CL_ERROR(clEnqueueWriteBufferRect(clCQ, subBufTo[cn], CL_TRUE, origTo, origHost, region, PDim<I>::stride(Y) * sizeof(T), PDim<I>::stride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, NULL));
	}
	for (I cn = 0; cn != nc; ++cn) {
		CHECK_CL_ERROR(clReleaseMemObject(subBufFrom[cn]));
		CHECK_CL_ERROR(clReleaseMemObject(subBufTo[cn]));
	}
	delete[] subBufTo;
	delete[] subBufFrom;
	delete[] hostMem;
}

template <typename T, typename I>
void DArray<T, I>::fillGhostCL(const CartDir &d, const CartSide &s) {
	size_t origFrom[ALL_DIRS], origTo[ALL_DIRS];
	size_t region[ALL_DIRS];
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	origFrom[ort1] = origTo[ort1] = PDim<I>::ghost(ort1);
	origFrom[ort2] = origTo[ort2] = PDim<I>::ghost(ort2);
	region[d] = 1;
	region[ort1] = PDim<I>::localSize(ort1);
	region[ort2] = PDim<I>::localSize(ort2);
	if (s == SIDE_LEFT) {
		origTo[d] = 0;
		origFrom[d] = PDim<I>::ghost(d);
	} else {
		origTo[d] = (PDim<I>::localSizeGhost(d) - PDim<I>::ghost(d));
		origFrom[d] = (PDim<I>::localSizeGhost(d) - PDim<I>::ghost(d) - 1);
	}
	// add sizeof(T) multiplier in X direction
	origFrom[X] *= sizeof(T);
	origTo[X] *= sizeof(T);
	region[X] *= sizeof(T);
	// Create Sub buffers to work with components separately
	cl_mem* subBuf = new cl_mem[nc];
	for (I cn = 0; cn != nc; ++cn) {
		cl_int err;
		cl_buffer_region clbr;
		clbr.size = PDim<I>::localSizeGhost() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBuf[cn] = clCreateSubBuffer(clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
	}
	// TODO make non blocking copy
	for (I gh = 0; gh != PDim<I>::ghost(d); ++gh) {
		for (I cn = 0; cn != nc; ++cn) {
			cl_event event;
			CHECK_CL_ERROR(clEnqueueCopyBufferRect(clCQ, subBuf[cn], subBuf[cn], origFrom, origTo, region, PDim<I>::stride(Y) * sizeof(T), PDim<I>::stride(Z) * sizeof(T), PDim<I>::stride(Y) * sizeof(T), PDim<I>::stride(Z) * sizeof(T), 0, NULL, &event));
			CHECK_CL_ERROR(clWaitForEvents(1, &event));
			CHECK_CL_ERROR(clReleaseEvent(event));
		}
		origTo[d] += (d == X) ? sizeof(T) : 1;
	}
	for (I cn = 0; cn != nc; ++cn) {
		CHECK_CL_ERROR(clReleaseMemObject(subBuf[cn]));
	}
	delete[] subBuf;
}

#endif // USE_OPENCL

}; // rgrid

#endif // RGRID_DARRAY_H
