#ifndef RG_IO_H
#define RG_IO_H

#include <iostream>

#include "rgrid/types.h"
#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"

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



template <typename T, typename I> 
void saveData(std::iostream& stream, const DArray<T, I>& dArray, rgio::format fmt) {
	writeHeader(stream, dArray, fmt);
	stream << "DATA START" << std::endl;
	for (I z = 0; z != dArray.localSize(Z); ++z) {
		for (I y = 0; y != dArray.localSize(Y); ++y) {
			dArray.writeLine(stream, y, z, fmt);
		}
	}
}

template <typename T, typename I> 
void saveData(std::iostream& stream, const DArrayContainer<T, I>& dArrayContainer, rgio::format fmt) {
	writeHeader(stream, dArrayContainer, fmt);
	stream << "DATA START" << std::endl;
	for (I z = 0; z != dArrayContainer.numNodes(Z); ++z) {
		for (I y = 0; y != dArrayContainer.numNodes(Y); ++y) {
			dArrayContainer.writeLine(stream, y, z, fmt);
		}
	}
}

template <typename T, typename I> 
void loadData(std::iostream& stream, DArray<T, I>& dArray) {
	rgio::format fmt = loadHeader(stream, dArray);
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
void writeHeader(std::iostream& stream, const DArray<T, I>& da, rgio::format fmt) {
	stream << "# DARRAY DATA FORMAT" << std::endl;
	stream << "SIZE: " << da.size(X) << " " << da.size(Y) << " " << da.size(Z) << std::endl;
	stream << "COMPONENTS: " << da.getNC() << std::endl;
	stream << "FORMAT: ";
	if (fmt == rgio::TEXT) stream << "text" << std::endl;
	else if (fmt == rgio::BINARY) stream << "binary" << std::endl;
}

template <typename T, typename I> 
void writeHeader(std::iostream& stream, const DArrayContainer<T, I>& dac, rgio::format fmt) {
	writeHeader(stream, dac.getDArrayPart(0), fmt);
}

template <typename T, typename I>
rgio::format loadHeader(std::iostream& stream, DArray<T, I>& dArray)
{
	// darray params
	I size[ALL_DIRS];
	I components;
	
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

	dArray.resize(size[X], size[Y], size[Z]);
	dArray.alloc(components);
	
	rgio::format ret = rgio::TEXT;
	std::getline(stream, str, ':');
	RG_ASSERT(0 == str.compare("FORMAT"), "Wrong entry, FORMAT expected");
	stream >> str;
	if (0 == str.compare("text")) ret = rgio::TEXT;
	else if (0 == str.compare("binary")) ret = rgio::BINARY;
	else {
		RG_ASSERT(0, "Wrong format");
	}
	std::getline(stream, str);
	return ret;
}

} // namespace rgio

} // namespace rgrid

#endif // RG_IO_H
