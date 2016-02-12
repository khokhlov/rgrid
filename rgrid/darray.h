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
		T &val(const I &i1, const CartDir &d1, const I &i2, const CartDir &d2, const I &i3, const CartDir &d3, const I cn) { return data[this->ind(i1, d1, i2, d2, i3, d3) + cn * this->localSizeGhost()]; }

		T &operator()(const I i, const I j, const I k, const I cn) { return val(i, j, k, cn); }
		T &operator()(const I i, const I j, const I cn) { return (*this)(i, j, static_cast<I>(0), cn); }
		T &operator()(const I i, const I cn) { return (*this)(i, static_cast<I>(0), cn); }

		const T &operator()(const I i, const I j, const I k, const I cn) const { return data[this->ind(i, j, k) + cn * this->localSizeGhost()]; }
		const T &operator()(const I i, const I j, const I cn) const { return (*this)(i, j, static_cast<I>(0), cn); }
		const T &operator()(const I i, const I cn) const { return (*this)(i, static_cast<I>(0), cn); }
		
		// save darray to file named "name" in binary format
		void saveBinaryFile(const char* name);
		// save darray to file named "name" in text format
		void saveTextFile(const char* name);
		// load darray from file named "name" in any format
		void loadFile(const char* name);
	private:
		void writeFileHeader(std::ofstream& f);
		void loadFileHeader(std::ifstream& f);
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

template <typename T, typename I>
void DArray<T, I>::saveTextFile(const char* name)
{
	std::ofstream f(name);
	RG_ASSERT(f.is_open(), "Can't open file");
	// write header
	writeFileHeader(f);
	// write data in text format
	f << "FORMAT: text\n";
	f << "DATA START\n";
	for (I c = 0; c != this->nc; ++c)
	for (I k = 0; k != this->localSize(Z); ++k)
	for (I j = 0; j != this->localSize(Y); ++j)
	for (I i = 0; i != this->localSize(X); ++i) {
		f << (*this)(i, j, k, c) << " ";
	}
	f.close();
}

template <typename T, typename I>
void DArray<T, I>::saveBinaryFile(const char* name)
{
	std::ofstream f(name);
	RG_ASSERT(f.is_open(), "Can't open file");
	// write header
	writeFileHeader(f);
	// write data in binary format
	f << "FORMAT: binary\n";
	f << "DATA START\n";
	for (I c = 0; c != this->nc; ++c)
	for (I k = 0; k != this->localSize(Z); ++k)
	for (I j = 0; j != this->localSize(Y); ++j)
	for (I i = 0; i != this->localSize(X); ++i) {
		T d = (*this)(i, j, k, c);
		const char *dc = reinterpret_cast<const char*>(&d);
		f.write(dc, sizeof(T));
	}
	f.close();
}

template <typename T, typename I>
void DArray<T, I>::writeFileHeader(std::ofstream& f)
{
	f << "# DARRAY DATA FORMAT\n";
	f << "LOCAL SIZE: " << this->localSize(X) << " " << this->localSize(Y) << " " << this->localSize(Z) << "\n";
	f << "GLOBAL SIZE: " << this->size(X) << " " << this->size(Y) << " " << this->size(Z) << "\n";
	f << "ORIGIN: " << this->origin(X) << " " << this->origin(Y) << " " << this->origin(Z) << "\n";
	f << "GHOST SIZE: " << this->ghost(X) << " " << this->ghost(Y) << " " << this->ghost(Z) << "\n";
	f << "COMPONENTS: " << this->nc << "\n";
}

template <typename T, typename I>
void DArray<T, I>::loadFileHeader(std::ifstream& f)
{
	// darray params
	I localSize[ALL_DIRS];
	I globalSize[ALL_DIRS];
	I origin[ALL_DIRS];
	I ghostSize[ALL_DIRS];
	I components;
	
	std::string str;
	
	std::getline(f, str); 
	RG_ASSERT(0 == str.compare("# DARRAY DATA FORMAT"), "Wrong file format");
	
	std::getline(f, str, ':');	
	RG_ASSERT(0 == str.compare("LOCAL SIZE"), "Wrong entry, LOCAL SIZE expected");
	f >> localSize[X] >> localSize[Y] >> localSize[Z];
	std::getline(f, str);
	
	std::getline(f, str, ':');
	RG_ASSERT(0 == str.compare("GLOBAL SIZE"), "Wrong entry, GLOBAL SIZE expected");
	f >> globalSize[X] >> globalSize[Y] >> globalSize[Z];
	std::getline(f, str);
	
	std::getline(f, str, ':');
	RG_ASSERT(0 == str.compare("ORIGIN"), "Wrong entry, ORIGIN expected");
	f >> origin[X] >> origin[Y] >> origin[Z];
	std::getline(f, str);
	
	std::getline(f, str, ':');
	RG_ASSERT(0 == str.compare("GHOST SIZE"), "Wrong entry, GHOST SIZE expected");
	f >> ghostSize[X] >> ghostSize[Y] >> ghostSize[Z];
	std::getline(f, str);
	
	std::getline(f, str, ':');
	RG_ASSERT(0 == str.compare("COMPONENTS"), "Wrong entry, COMPONENTS expected");
	f >> components;
	std::getline(f, str);

	this->resize(
		localSize[X], localSize[Y], localSize[Z],
		globalSize[X], globalSize[Y], globalSize[Z],
		origin[X], origin[Y], origin[Z],
		ghostSize[X], ghostSize[Y], ghostSize[Z]);
	this->alloc(components);
}

template <typename T, typename I>
void DArray<T, I>::loadFile(const char* name)
{
	std::ifstream f(name);
	RG_ASSERT(f.is_open(), "Can't open file");
	
	loadFileHeader(f);
	
	std::string str;
	std::getline(f, str, ':');
	RG_ASSERT(0 == str.compare("FORMAT"), "Wrong entry, FORMAT expected");
	f >> str;
	if (0 == str.compare("text")) {
		// load data in text format
		std::getline(f, str);
		std::getline(f, str);
		RG_ASSERT(0 == str.compare("DATA START"), "Wrong entry, DATA START expected");
		for (I c = 0; c != this->nc; ++c)
		for (I k = 0; k != this->localSize(Z); ++k)
		for (I j = 0; j != this->localSize(Y); ++j)
		for (I i = 0; i != this->localSize(X); ++i) {
			f >> (*this)(i, j, k, c);
		}
	} else if (0 == str.compare("binary")) {
		// load data in binary format
		std::getline(f, str);
		std::getline(f, str);
		RG_ASSERT(0 == str.compare("DATA START"), "Wrong entry, DATA START expected");
		for (I c = 0; c != this->nc; ++c)
		for (I k = 0; k != this->localSize(Z); ++k)
		for (I j = 0; j != this->localSize(Y); ++j)
		for (I i = 0; i != this->localSize(X); ++i) {
			T d;
			char *dc = reinterpret_cast<char*>(&d);
			f.read(dc, sizeof(T));
			(*this)(i, j, k, c) = d;
		}
	} else {
		RG_ASSERT(0, "Unknown format");
	}
	
	f.close();
}

template <typename T, typename I>
bool operator==(const DArray<T, I>& lhs, const DArray<T, I>& rhs) {
	bool is1 = lhs.data == rhs.data;
	bool is2 = lhs.nc == rhs.nc;
	bool is3 = static_cast<PDim<I> >(lhs) == static_cast<PDim<I> >(rhs);
	return is1 && is2 && is3;
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
