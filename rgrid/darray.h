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
#include "rgrid/rgio.h"
#include "rgrid/clwrapper.h"
#include "rgrid/clutils.h"

namespace rgrid {

template <typename T, typename I>
class DArray : public PDim<I> {
public:
	DArray() : PDim<I>() {
		init();
		alloc();
	}

public:

	virtual void resize(const T x, const T y, const T z,
	                    const T px, const T py, const T pz,
	                    const T ox, const T oy, const T oz,
	                    const T gx, const T gy, const T gz,
	                    const T cn) {
		PDim<I>::resize(x, y, z, px, py, pz, ox, oy, oz, gx, gy, gz, cn);
		alloc();
	}
	virtual void resize(const PDim<I> &p) {
		PDim<I>::resize(p);
		alloc();
	}
	virtual void resize(const T x = 1, const T y = 1, const T z = 1, const T cn = 1) {
		PDim<I>::resize(x, y, z, cn);
		alloc();
	}
	virtual void resize(const T x, const T y, const T z,
	                    const T gx, const T gy, const T gz,
	                    const T cn) {
		PDim<I>::resize(x, y, z, gx, gy, gz, cn);
		alloc();
	}
	virtual void resize_d(const CartDir d1, const T s1, const CartDir d2, const T s2, const CartDir d3, const T s3, const T cn) {
		PDim<I>::resize(d1, s1, d2, s2, d3, s3, cn);
		alloc();
	}

private:
	void init() {
#ifdef USE_OPENCL
		clContextIsSet = false;
		clBufInitialized = false;
		clOnDevice = false;
#endif
	}

	void alloc() {
		data.resize(PDim<I>::getNC() * PDim<I>::localGhostSize());
	}

public:

	~DArray();

	void fillGhost(const CartDir d, const CartSide s);
	void fillGhost(const CartDir d) {
		fillGhost(d, SIDE_LEFT);
		fillGhost(d, SIDE_RIGHT);
	}
	void fillGhost() {
		fillGhost(X);
		fillGhost(Y);
		fillGhost(Z);
	}

	// copy ghost nodes from "da"
	// "d" and "s" relative to this DArray
	void copyGhost(DArray<T, I> &da, const CartDir d, const CartSide s);

	void fill(const T &v);

	T *getDataRaw() {
		return &data.at(0);
	}

	T &operator[](const I index) {
		return data[index];
	}
	const T &operator[](const I index) const {
		return data[index];
	}

	T &val(const I i, const I j, const I k, const I cn) {
		return data.at(PDim<I>::ind(i, j, k, cn));
	}
	T &val(const I i, const I j, const I cn) {
		return val(i, j, static_cast<T>(0), cn);
	}
	T &val(const I i, const I cn) {
		return val(i, static_cast<T>(0), cn);
	}
	T &val(const I i1, const CartDir d1, const I i2, const CartDir d2, const I i3, const CartDir d3, const I cn) {
		return data[PDim<I>::ind(i1, d1, i2, d2, i3, d3, cn)];
	}

	T &operator()(const I i, const I j, const I k, const I cn) {
		return val(i, j, k, cn);
	}
	T &operator()(const I i, const I j, const I cn) {
		return val(i, j, cn);
	}
	T &operator()(const I i, const I cn) {
		return val(i, cn);
	}

	const T &operator()(const I i, const I j, const I k, const I cn) const {
		return data[PDim<I>::ind(i, j, k, cn)];
	}
	const T &operator()(const I i, const I j, const I cn) const {
		return val(i, j, cn);
	}
	const T &operator()(const I i, const I cn) const {
		return val(i, cn);
	}

	I dataSize() const {
		return data.size();
	}

	// write line (all nodes on X direction) with coordinates y and z into stream
	void readLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt);

	// write line (all nodes on X direction) with coordinates y and z into stream
	void writeLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) const;

	void loadData(std::iostream &stream);

	void saveData(std::iostream &stream, const rgio::format fmt) const;

	/*
	 * write data to buffer,
	 * start with position "start"
	 * return number of written elements
	 * TODO implement
	 */
	typename std::vector<char>::size_type
	saveData(std::vector<char> &buffer, typename std::vector<char>::size_type start, const rgio::format fmt) const;

#ifdef USE_OPENCL
public:
	// the same as copy ghost, but it make copy on device
	void copyGhostCL(DArray<T, I> &da, const CartDir &d, const CartSide &s);

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

	/* Get second temporary buffer with the same size as first buffer to improve performaPDim<I>::getNC()e
	 * You must run setCLContext() and setCLCQ() before */
	cl_mem &getCLBuffer2();

	/* swap first and second buffers */
	void swapCLBuffers() {
		std::swap(clBuffer, clBuffer2);
	}

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
}; // DArray

template <typename T, typename I>
DArray<T, I>::~DArray() {
#ifdef USE_OPENCL
	clReleaseBuffer();
#endif
}

template <typename T, typename I>
void DArray<T, I>::fill(const T &v) {
	for (typename std::vector<T>::size_type i = 0; i != data.size(); ++i) {
		data.at(i) = v;
	}
}

template <typename T, typename I>
void DArray<T, I>::copyGhost(DArray<T, I> &da, const CartDir d, const CartSide s) {
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
	for (I cn = 0; cn < PDim<I>::getNC(); cn++)
		for (I i = 0; i < PDim<I>::localSize(ort1); i++)
			for (I j = 0; j < PDim<I>::localSize(ort2); j++)
				for (I gs = 0; gs < PDim<I>::ghost(d); gs++) {
					val(i, ort1, j, ort2, k - sign * (gs + 1), d, cn) = da.val(i, ort1, j, ort2, kda - sign * gs, d, cn);
				}
}

template <typename T, typename I>
void DArray<T, I>::fillGhost(const CartDir d, const CartSide s) {
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	I k = 0;
	I sign;
	if (s == SIDE_LEFT) {
		k = 0;
		sign = 1;
	} else {
		k = this->localSize(d) - 1;
		sign = -1;
	}
	for (I cn = 0; cn < PDim<I>::getNC(); cn++) {
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
void DArray<T, I>::writeLine(std::iostream &stream, I cn, I y, I z, rgio::format fmt) const {
	for (I x = 0; x != PDim<I>::localSize(X); ++x) {
		if (fmt == rgio::TEXT) {
			stream << operator()(x, y, z, cn) << " ";
		} else if (fmt == rgio::BINARY) {
			T d = operator()(x, y, z, cn);
			const char *dc = (const char *)(&d);
			stream.write(dc, sizeof(T));
		}
	}
}

template <typename T, typename I>
void DArray<T, I>::readLine(std::iostream &stream, const I cn, const I y, const I z, rgio::format fmt) {
	for (I x = 0; x != PDim<I>::localSize(X); ++x) {
		if (fmt == rgio::TEXT) {
			stream >> operator()(x, y, z, cn);
		} else if (fmt == rgio::BINARY) {
			T d;
			char *dc = (char *)(&d);
			stream.read(dc, sizeof(T));
			operator()(x, y, z, cn) = d;
		}
	}
}

template <typename T, typename I>
void DArray<T, I>::saveData(std::iostream &stream, const rgio::format fmt) const {
	Dim3D<I> size(PDim<I>::localSize(X), PDim<I>::localSize(Y), PDim<I>::localSize(Z));
	rgio::writeHeader(stream, size, PDim<I>::getNC(), fmt);
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		for (I z = 0; z != PDim<I>::localSize(Z); ++z) {
			for (I y = 0; y != PDim<I>::localSize(Y); ++y) {
				writeLine(stream, cn, y, z, fmt);
			}
		}
	}
}

template <typename T, typename I>
void DArray<T, I>::loadData(std::iostream &stream) {
	rgio::format fmt;
	Dim3D<I> size;
	I nc;
	rgio::loadHeader(stream, size, nc, fmt);
	PDim<I>::resize(size[X], size[Y], size[Z], nc);
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		for (I z = 0; z != PDim<I>::localSize(Z); ++z) {
			for (I y = 0; y != PDim<I>::localSize(Y); ++y) {
				readLine(stream, cn, y, z, fmt);
			}
		}
	}
}

template <typename T, typename I>
bool operator==(const DArray<T, I> &lhs, const DArray<T, I> &rhs) {
	if (&lhs == &rhs) {
		return true;
	}
	if (!lhs.isEqualNoGhost(rhs)) {
		return false;
	}
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
cl_mem &DArray<T, I>::getCLBuffer() {
	if (!clBufInitialized) {
		clInitBuffer();
	}
	return this->clBuffer;
}

template <typename T, typename I>
cl_mem &DArray<T, I>::getCLBuffer2() {
	if (!clBufInitialized) {
		clInitBuffer();
	}
	return this->clBuffer2;
}

template <typename T, typename I>
void DArray<T, I>::clHtoD() {
	if (!clOnDevice) {
		if (!clBufInitialized) {
			clInitBuffer();
		}
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
void DArray<T, I>::copyGhostCL(DArray<T, I> &da, const CartDir &d, const CartSide &s) {
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
		origFrom[d] = (PDim<I>::localGhostSize(d) - 2 * PDim<I>::ghost(d));
	} else {
		origTo[d] = (PDim<I>::localGhostSize(d) - PDim<I>::ghost(d));
		origFrom[d] = da.ghost(d);
	}
	I regionSize = region[X] * region[Y] * region[Z];
	T *hostMem = new T[regionSize * PDim<I>::getNC()];
	// add sizeof(T) multiplier in X direction
	origFrom[X] *= sizeof(T);
	origTo[X] *= sizeof(T);
	region[X] *= sizeof(T);
	// Create Sub buffers to work with components separately
	cl_mem *subBufFrom = new cl_mem[PDim<I>::getNC()];
	cl_mem *subBufTo = new cl_mem[PDim<I>::getNC()];
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		cl_int err;
		cl_buffer_region clbr;
		clbr.size = da.localGhostSize() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBufFrom[cn] = clCreateSubBuffer(da.clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
		clbr.size = PDim<I>::localGhostSize() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBufTo[cn] = clCreateSubBuffer(clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
	}
	// TODO make non blocking copy
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		CHECK_CL_ERROR(clEnqueueReadBufferRect(da.clCQ, subBufFrom[cn], CL_TRUE, origFrom, origHost, region, da.localStride(Y) * sizeof(T), da.localStride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, NULL));
	}
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		CHECK_CL_ERROR(clEnqueueWriteBufferRect(clCQ, subBufTo[cn], CL_TRUE, origTo, origHost, region, PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, NULL));
	}
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
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
		origTo[d] = (PDim<I>::localGhostSize(d) - PDim<I>::ghost(d));
		origFrom[d] = (PDim<I>::localGhostSize(d) - PDim<I>::ghost(d) - 1);
	}
	// add sizeof(T) multiplier in X direction
	origFrom[X] *= sizeof(T);
	origTo[X] *= sizeof(T);
	region[X] *= sizeof(T);
	// Create Sub buffers to work with components separately
	cl_mem *subBuf = new cl_mem[PDim<I>::getNC()];
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		cl_int err;
		cl_buffer_region clbr;
		clbr.size = PDim<I>::localGhostSize() * sizeof(T);
		clbr.origin = clbr.size * cn;
		subBuf[cn] = clCreateSubBuffer(clBuffer, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &clbr, &err);
		CHECK_CL_ERROR(err);
	}
	// TODO make non blocking copy
	for (I gh = 0; gh != PDim<I>::ghost(d); ++gh) {
		for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
			cl_event event;
			CHECK_CL_ERROR(clEnqueueCopyBufferRect(clCQ, subBuf[cn], subBuf[cn], origFrom, origTo, region, PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), 0, NULL, &event));
			CHECK_CL_ERROR(clWaitForEvents(1, &event));
			CHECK_CL_ERROR(clReleaseEvent(event));
		}
		origTo[d] += (d == X) ? sizeof(T) : 1;
	}
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		CHECK_CL_ERROR(clReleaseMemObject(subBuf[cn]));
	}
	delete[] subBuf;
}

#endif // USE_OPENCL

}; // rgrid

#endif // RGRID_DARRAY_H


