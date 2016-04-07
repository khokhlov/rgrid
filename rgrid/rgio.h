#ifndef RG_IO_H
#define RG_IO_H

#include <iostream>
#include "rgrid/types.h"
#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"

// TODO use another place for it
namespace rgrid {
namespace rgio {

template <typename I>
void writeHeader(std::iostream& stream, const I size[ALL_DIRS], const I components, const format fmt) {
	stream << "# DARRAY DATA FORMAT" << std::endl;
	stream << "SIZE: " << size[X] << " " << size[Y] << " " << size[Z] << std::endl;
	stream << "COMPONENTS: " << components << std::endl;
	stream << "FORMAT: ";
	if (fmt == TEXT) stream << "text" << std::endl;
	else if (fmt == BINARY) stream << "binary" << std::endl;
}

}
}

#include "rgrid/darrayscatter.h"

namespace rgrid {

namespace rgio {

template <typename T, typename I> 
void saveData(std::iostream& stream, const DArray<T, I>& dArray, format fmt);
template <typename T, typename I>
void saveData(std::iostream& stream, const DArrayContainer<T, I>& dArrayContainer, format fmt);

template <typename T, typename I>
void loadData(std::iostream& stream, DArray<T, I>& dArray);

template <typename T, typename I>
void writeHeader(std::iostream& stream, const DArray<T, I>& dArray, format fmt);
template <typename T, typename I>
void writeHeader(std::iostream& stream, const DArrayContainer<T, I>& dArrayContainer, format fmt);

template <typename T, typename I>
format loadHeader(std::iostream& stream, DArray<T, I>& dArray);

/*
 * basic header loader
 * stream - IN
 * size - OUT array of sizes in header
 * components - OUT number of comonents
 */
template <typename I>
format loadHeader(std::iostream& stream, I size[ALL_DIRS], I& components);

/*
 * basic header saver
 */
template <typename I>
void writeHeader(std::iostream& stream, const I size[ALL_DIRS], const I components, const format fmt);

#ifdef USE_MPI
//template <typename T, typename I>
//void saveDataBegin(const std::string filename, DArrayScatter<T, I>& das, const format fmt, std::vector<MPI_Request>& requests);
//
//void saveDataEnd(std::vector<MPI_Request>& requests);
//
//template <typename T, typename I>
//void loadDataBegin(const std::string filename, DArrayScatter<T, I>& das, std::vector<MPI_Request>& requests);
//
//void loadDataEnd(std::vector<MPI_Request>& requests);
#endif // USE_MPI

/* Start of implementations */

template <typename T, typename I> 
void saveData(std::iostream& stream, const DArray<T, I>& dArray, format fmt) {
	writeHeader(stream, dArray, fmt);
	stream << "DATA START" << std::endl;
	for (I z = 0; z != dArray.localSize(Z); ++z) {
		for (I y = 0; y != dArray.localSize(Y); ++y) {
			dArray.writeLine(stream, y, z, fmt);
		}
	}
}

template <typename T, typename I> 
void saveData(std::iostream& stream, const DArrayContainer<T, I>& dArrayContainer, format fmt) {
	writeHeader(stream, dArrayContainer, fmt);
	stream << "DATA START" << std::endl;
	for (I z = 0; z != dArrayContainer.numNodes(Z); ++z) {
		for (I y = 0; y != dArrayContainer.numNodes(Y); ++y) {
			dArrayContainer.writeLine(stream, y, z, fmt);
		}
	}
}

#ifdef USE_MPI
template <typename T, typename I>
void saveDataBegin(const std::string filename, DArrayScatter<T, I>& das, const format fmt, std::vector<MPI_Request>& requests) {
	MPI_File fh;
	std::stringstream ss;
	I size[ALL_DIRS] = { das.numNodes(X), das.numNodes(Y), das.numNodes(Z) };
	writeHeader(ss, size, das.getNC(), fmt);
	std::string header = ss.str();
	// write header
	if (rgmpi::worldRank() == 0) {
		MPI_CHECK(MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
		MPI_CHECK(MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL));
		MPI_CHECK(MPI_File_write(fh, header.c_str(), header.size(), MPI_CHAR, MPI_STATUS_IGNORE));
		MPI_CHECK(MPI_File_close(&fh));
	}
	// write data
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
	MPI_CHECK(MPI_File_set_view(fh, header.size() * sizeof(char), rgmpi::getMPItype<T>(), das.fileViewType(), "native", MPI_INFO_NULL));
	MPI_Request req;
	requests.clear();
	requests.reserve(das.numParts());
	for (I k = 0; k != das.numParts(Z); ++k)
	for (I j = 0; j != das.numParts(Y); ++j)
	for (I i = 0; i != das.numParts(X); ++i) {
		MPI_CHECK(MPI_File_iwrite(fh, das.getDArrayBuffer(i, j, k), 1, das.getDArrayDt(i, j, k), &req));
		requests.push_back(req);
	}
	MPI_CHECK(MPI_File_close(&fh));
}

void saveDataEnd(std::vector<MPI_Request>& requests) {
	MPI_CHECK(MPI_Waitall(requests.size(), &requests.at(0), MPI_STATUSES_IGNORE));
	requests.clear();
}

template <typename T, typename I>
void loadDataBegin(const std::string filename, DArrayScatter<T, I>& das, std::vector<MPI_Request>& requests) {
	MPI_File fh;
	std::ifstream fs;
	I size[ALL_DIRS];
	I nc;
	// read header
	fs.open(filename.c_str());
	format fmt = loadHeader(fs, size, nc);
	long offset = fs.tellg();
	fs.close();
	I globalPt[ALL_DIRS] = { das.numParts(X), das.numParts(Y), das.numParts(Z) };
	I localPt[ALL_DIRS] = { das.getDAC().numParts(X), das.getDAC().numParts(Y), das.getDAC().numParts(Z) };
	I ghost[ALL_DIRS] = { das.getGhost(X), das.getGhost(Y), das.getGhost(Z) };
	das.setSizes(size, globalPt, localPt, ghost, nc);
	// read data
	MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
	MPI_CHECK(MPI_File_set_view(fh, offset * sizeof(char), rgmpi::getMPItype<T>(), das.fileViewType(), "native", MPI_INFO_NULL));
	MPI_Request req;
	requests.clear();
	requests.reserve(das.numParts);
	for (I k = 0; k != das.numParts(Z); ++k)
	for (I j = 0; j != das.numParts(Y); ++j)
	for (I i = 0; i != das.numParts(X); ++i) {
		MPI_CHECK(MPI_File_iread(fh, das.getDArrayBuffer(i, j, k), 1, das.getDArrayDt(i, j, k), &req));
		requests.push_back(req);
	}
	MPI_CHECK(MPI_File_close(&fh));
}

void loadDataEnd(std::vector<MPI_Request>& requests) {
	MPI_CHECK(MPI_Waitall(requests.size(), &requests.at(0), MPI_STATUSES_IGNORE));
	requests.clear();
}
#endif // USE_MPI

template <typename T, typename I>
void loadData(std::iostream& stream, DArray<T, I>& dArray) {
	format fmt = loadHeader(stream, dArray);
	std::string str;
	std::getline(stream, str);
	RG_ASSERT(0 == str.compare("DATA START"), "Wrong entry, DATA START expected");
	for (I z = 0; z != dArray.localSize(Z); ++z) {
		for (I y = 0; y != dArray.localSize(Y); ++y) {
			dArray.readLine(stream, y, z, fmt);
		}
	}
}

template <typename T, typename I> 
void writeHeader(std::iostream& stream, const DArray<T, I>& da, format fmt) {
	I size[ALL_DIRS] = { da.size(X), da.size(Y), da.size(Z) };
	writeHeader(stream, size, da.getNC(), fmt);
}

template <typename T, typename I> 
void writeHeader(std::iostream& stream, const DArrayContainer<T, I>& dac, format fmt) {
	writeHeader(stream, dac.getDArrayPart(0), fmt);
}

template <typename T, typename I>
format loadHeader(std::iostream& stream, DArray<T, I>& dArray) {
	I size[ALL_DIRS];
	I components;
	format ret = loadHeader(stream, size, components);
	dArray.resize(size[X], size[Y], size[Z]);
	dArray.alloc(components);
	return ret;
}

template <typename I>
format loadHeader(std::iostream& stream, I size[ALL_DIRS], I& components) {	
	std::string str;
	
	std::getline(stream, str); 
	RG_ASSERT(0 == str.compare("# DARRAY DATA FORMAT"), "Wrong file format");
	
	std::getline(stream, str, ':');	
	RG_ASSERT(0 == str.compare("SIZE"), "Wrong entry, SIZE expected");
	stream >> size[X] >> size[Y] >> size[Z];
	std::getline(stream, str);
	
	std::getline(stream, str, ':');
	RG_ASSERT(0 == str.compare("COMPONENTS"), "Wrong entry, COMPONENTS expected");
	stream >> components;
	std::getline(stream, str);
	
	format ret = TEXT;
	std::getline(stream, str, ':');
	RG_ASSERT(0 == str.compare("FORMAT"), "Wrong entry, FORMAT expected");
	stream >> str;
	if (0 == str.compare("text")) ret = TEXT;
	else if (0 == str.compare("binary")) ret = BINARY;
	else {
		RG_ASSERT(0, "Wrong format");
	}
	std::getline(stream, str);
	return ret;
}



} // namespace rgio

} // namespace rgrid

#endif // RG_IO_H
