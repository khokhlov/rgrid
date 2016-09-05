/**
 * \file
 * \brief Creation of rectangular structures
 */

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

/**
 * \brief Class for creation rectangular structures.
 * 
 * DArray is 3 dimensional structure and every node can have several components of same type.
 * Also DArray can represent rectangular structure (subarray) inside bigger rectangular structure.
 * So it have two sizes: global - size in bigger structure, local - size of current DArray.
 * 
 * DArray can have ghost nodes on edges for auxiliary purposes. Indexes of ghost nodes lay outside DArray.
 * 
 * \tparam T type of every grid node (i.e. double, float)
 * \tparam I type of grid indexes (i.e. int, long)
 */
template <typename T, typename I>
class DArray : public PDim<I> {
public:
	/**
	 * Call resize() if you use this constructor
	 */
	DArray() : PDim<I>() {
		init();
		alloc();
	}

public:

	/**
	 * \brief Set new sizes to DArray
	 * 
	 * The most verbose version
	 * 
	 * \param[in] x,y,z size of bigger rectangular structure
	 * \param[in] px,py,pz size of current DArray
	 * \param[in] ox,oy,oz origin of current DArray in bigger structure
	 * \param[in] gx,gy,gz number of ghost nodes on each side
	 * \param[in] cn number of components
	 */
	virtual void resize(const I x, const I y, const I z,
	                    const I px, const I py, const I pz,
	                    const I ox, const I oy, const I oz,
	                    const I gx, const I gy, const I gz,
	                    const I cn) {
		PDim<I>::resize(x, y, z, px, py, pz, ox, oy, oz, gx, gy, gz, cn);
		alloc();
	}
	/**
	 * \brief resize this PDim as the other PDim
	 */
	virtual void resize(const PDim<I> &p) {
		PDim<I>::resize(p);
		alloc();
	}
	/**
	 * \brief Set new sizes to DArray
	 * \param[in] x,y,z size of DArray
	 * \param[in] cn number of components
	 */
	virtual void resize(const I x = 1, const I y = 1, const I z = 1, const I cn = 1) {
		PDim<I>::resize(x, y, z, cn);
		alloc();
	}
	/**
	 * \brief Set new sizes to DArray
	 * \param[in] x,y,z size of DArray
	 * \param[in] gx,gy,gz number of ghost nodes on each side
	 * \param[in] cn number of components
	 */
	virtual void resize(const I x, const I y, const I z,
	                    const I gx, const I gy, const I gz,
	                    const I cn) {
		PDim<I>::resize(x, y, z, gx, gy, gz, cn);
		alloc();
	}
	/**
	 * \brief Set new sizes to DArray
	 * 
	 * \param[in] d1,d2,d3 different dimensions in any order
	 * \param[in] s1,s2,s3 size of DArray in the same order as in d1,d2,d3
	 * \param[in] cn number of components
	 */
	virtual void resize_d(const CartDir d1, const I s1, const CartDir d2, const I s2, const CartDir d3, const I s3, const I cn) {
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

	/**
	 * \brief Fill ghost nodes by values from adjacent nodes
	 * \param[in] d direction
	 * \param[in] s side
	 */
	void fillGhost(const CartDir d, const CartSide s);
	/**
	 * \brief Fill ghost nodes by values from adjacent nodes on both sides
	 * \param[in] d direction
	 */
	void fillGhost(const CartDir d) {
		fillGhost(d, SIDE_LEFT);
		fillGhost(d, SIDE_RIGHT);
	}
	/**
	 * \brief Fill all ghost nodes by values from adjacent nodes
	 */
	void fillGhost() {
		fillGhost(X);
		fillGhost(Y);
		fillGhost(Z);
	}
	/** 
	 * \brief Copy nodes from another DArray to this DArray.
	 * 
	 * Take edge nodes from another DArray and place them into this DArray
	 * 
	 * \param[in] da DArray to copy from
	 * \param[in] d direction relative to this DArray
	 * \param[in] s side relative to this DArray
	 */
	void copyGhost(DArray<T, I> &da, const CartDir d, const CartSide s);
	/**
	 * \brief Fill all DArray by one value
	 * \param[in] v value to fill
	 */
	void fill(const T &v);
	/**
	 * \brief Get raw pointer to DArray nodes
	 * \return pointer to nodes
	 */
	T *getDataRaw() {
		return &data.at(0);
	}
	/**
	 * \brief Get node component in internal linear array of nodes
	 * \param[in] index index in internal array
	 * \return node component
	 */
	T &operator[](const I index) {
		return data[index];
	}
	/**
	 * \brief Get node component in internal linear array of nodes
	 * \param[in] index index in internal array
	 * \return node component
	 */
	const T &operator[](const I index) const {
		return data[index];
	}

	/** 
	 * \brief Get node component at specified position in 3D DArray
	 * \param[in] i,j,k coordinate
	 * \param[in] cn component number
	 * \return node component
	 */
	T &val(const I i, const I j, const I k, const I cn) {
		return data.at(PDim<I>::ind(i, j, k, cn));
	}
	/** 
	 * \brief Get node component at specified position in 2D DArray
	 * \param[in] i,j coordinate
	 * \param[in] cn component number
	 * \return node component
	 */
	T &val(const I i, const I j, const I cn) {
		return val(i, j, static_cast<T>(0), cn);
	}
	/** 
	 * \brief Get node component at specified position in 1D DArray
	 * \param[in] i coordinate
	 * \param[in] cn component number
	 * \return node component
	 */
	T &val(const I i, const I cn) {
		return val(i, static_cast<T>(0), cn);
	}
	/** 
	 * \brief Get node component at specified position in 3D DArray
	 * \param[in] d1,d2,d3 different coordinates in any order
	 * \param[in] i1,i2,i3 components in coordinates order
	 * \param[in] cn component number
	 * \return node component
	 */
	T &val(const I i1, const CartDir d1, const I i2, const CartDir d2, const I i3, const CartDir d3, const I cn) {
		return data[PDim<I>::ind(i1, d1, i2, d2, i3, d3, cn)];
	}
	/**
	 * \brief Equivalent to val()
	 */
	T &operator()(const I i, const I j, const I k, const I cn) {
		return val(i, j, k, cn);
	}
	/**
	 * \brief Equivalent to val()
	 */
	T &operator()(const I i, const I j, const I cn) {
		return val(i, j, cn);
	}
	/**
	 * \brief Equivalent to val()
	 */
	T &operator()(const I i, const I cn) {
		return val(i, cn);
	}
	/**
	 * \brief Equivalent to val()
	 */
	const T &operator()(const I i, const I j, const I k, const I cn) const {
		return data[PDim<I>::ind(i, j, k, cn)];
	}
	/**
	 * \brief Equivalent to val()
	 */
	const T &operator()(const I i, const I j, const I cn) const {
		return val(i, j, cn);
	}
	/**
	 * \brief Equivalent to val()
	 */
	const T &operator()(const I i, const I cn) const {
		return val(i, cn);
	}
	/**
	 * \brief Size of entire DArray in memory
	 */
	I dataSize() const {
		return data.size();
	}

	/** 
	 * \brief Read all nodes in X direction from stream
	 * \param[in] stream
	 * \param[in] cn component number
	 * \param[in] y,z coordinates of all nodes in X direction
	 * \param[in] fmt format
	 */
	void readLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt);

	/** 
	 * \brief Write all nodes in X direction into stream
	 * \param[in] stream
	 * \param[in] cn component number
	 * \param[in] y,z coordinates of all nodes in X direction
	 * \param[in] fmt format
	 */
	void writeLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) const;

	/**
	 * \brief Load DArray from stream
	 * \param[in] stream
	 */
	void loadData(std::iostream &stream);

	/**
	 * \brief Save DArray into stream
	 * \param[in] stream
	 * \param[in] fmt format
	 */
	void saveData(std::iostream &stream, const rgio::format fmt) const;

	/**
	 * \brief Write data to buffer starting with some position
	 * \param[out] buffer
	 * \param[in] start initial position in buffer
	 * \param[in] fmt format
	 * \return number of written elements
	 * \todo not implemented
	 */
	typename std::vector<char>::size_type
	saveData(std::vector<char> &buffer, typename std::vector<char>::size_type start, const rgio::format fmt) const;

#ifdef USE_OPENCL
public:
	/** 
	 * \brief OpenCL implementation of copyGhost()
	 * 
	 * Make copy to ghost nodes of this DArray from another DArray
	 * \note Both DArrays have to be placed on device. 
	 * If event is NULL perform blocking copy, 
	 * otherwise return event which user have to wait and perform asynchronous copy
	 * \param[in] da DArray to copy from
	 * \param[in] d direction
	 * \param[in] s side
	 * \param[in,out] event if event is NULL perform blocking copy, otherwise return event which user have to wait and perform asynchronous copy
	 */
	void copyGhostCL(DArray<T, I> &da, const CartDir &d, const CartSide &s, cl_event* event = NULL);

	/** 
	 * \brief OpenCL implementation of fillGhost()
	 */
	void fillGhostCL(const CartDir &d, const CartSide &s);

	/** 
	 * \brief Set OpenCL context to work with buffer 
	 */
	void setCLContext(cl_context context);

	/** 
	 * \brief Set OpenCL command queue
	 */
	void setCLCQ(cl_command_queue cq);

	/** 
	 * \brief OpenCL data buffer to work in kernels
	 * 
	 * clHtoD() and clDtoH() uses this buffer
	 * \note You must run setCLContext() and setCLCQ() before
	 */
	cl_mem &getCLBuffer();

	/**
	 * \brief Get second temporary buffer with the same size as first buffer 
	 * Second buffer helps to work with OpenCL
	 * \note You must run setCLContext() and setCLCQ() before 
	 */
	cl_mem &getCLBuffer2();

	/**
	 * \brief Swap first and second buffers 
	 */
	void swapCLBuffers() {
		std::swap(clBuffer, clBuffer2);
	}

	/** 
	 * \brief Copy data from host to device 
	 */
	void clHtoD();

	/** 
	 * \brief Copy data from device to host 
	 */
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

/// Check is two DArrays equal
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
void DArray<T, I>::copyGhostCL(DArray<T, I> &da, const CartDir &d, const CartSide &s, cl_event* event) {
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
		origFrom[d] = (da.localGhostSize(d) - 2 * da.ghost(d));
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
	if (event == NULL) { // make blocking copy
		cl_event* read_event_list = new cl_event[PDim<I>::getNC()];
		cl_event* write_event_list = new cl_event[PDim<I>::getNC()];
		for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
			CHECK_CL_ERROR(clEnqueueReadBufferRect(da.clCQ, da.clBuffer, CL_TRUE, origFrom, origHost, region, da.localStride(Y) * sizeof(T), da.localStride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, &read_event_list[cn]));
			origFrom[X] += da.localGhostSize() * sizeof(T);
		}
		CHECK_CL_ERROR(clWaitForEvents(PDim<I>::getNC(), read_event_list));
		for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
			CHECK_CL_ERROR(clEnqueueWriteBufferRect(clCQ, clBuffer, CL_TRUE, origTo, origHost, region, PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), region[X], region[X] * region[Y], hostMem + regionSize * cn, 0, NULL, &write_event_list[cn]));
			origTo[X] += PDim<I>::localGhostSize() * sizeof(T);
		}
		CHECK_CL_ERROR(clWaitForEvents(PDim<I>::getNC(), write_event_list));
		for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
			CHECK_CL_ERROR(clReleaseEvent(read_event_list[cn]));
			CHECK_CL_ERROR(clReleaseEvent(write_event_list[cn]));
		}
		delete[] read_event_list;
		delete[] write_event_list;
		delete[] hostMem;
	} else { // make non blocking copy
		// TODO
		delete[] hostMem;
	}
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
	// TODO make non blocking copy
	for (I cn = 0; cn != PDim<I>::getNC(); ++cn) {
		for (I gh = 0; gh != PDim<I>::ghost(d); ++gh) {
			cl_event event;
			CHECK_CL_ERROR(clEnqueueCopyBufferRect(clCQ, clBuffer, clBuffer, origFrom, origTo, region, PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), PDim<I>::localStride(Y) * sizeof(T), PDim<I>::localStride(Z) * sizeof(T), 0, NULL, &event));
			CHECK_CL_ERROR(clWaitForEvents(1, &event));
			CHECK_CL_ERROR(clReleaseEvent(event));
			origTo[d] += (d == X) ? sizeof(T) : 1;
		}
		origTo[d] -= ((d == X) ? sizeof(T) : 1) * PDim<I>::ghost(d);
		origFrom[X] += PDim<I>::localGhostSize() * sizeof(T);
		origTo[X] += PDim<I>::localGhostSize() * sizeof(T);
	}
}

#endif // USE_OPENCL

}; // rgrid

#endif // RGRID_DARRAY_H


