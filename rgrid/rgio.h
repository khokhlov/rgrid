/**
 * \file
 * \brief Input/output operations
 */

#ifndef RG_IO_H
#define RG_IO_H

#include <iostream>
#include "rgrid/types.h"
#include <cstdio>
#include "rgrid/rgmpi.h"

namespace rgrid {

namespace rgio {

/**
 * \brief Basic header loader (own DArray format)
 * \param[in] stream
 * \param[out] size array of sizes in header
 * \param[out] components number of comonents
 * \param[out] fmt format
 */
template <typename I>
void loadHeader(std::iostream& stream, Dim3D<I>& size, I& components, format& fmt) {
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
	
	std::getline(stream, str, ':');
	RG_ASSERT(0 == str.compare("FORMAT"), "Wrong entry, FORMAT expected");
	stream >> str;
	if (0 == str.compare("text")) fmt = TEXT;
	else if (0 == str.compare("binary")) fmt = BINARY;
	else {
		RG_ASSERT(0, "Wrong format");
	}
	std::getline(stream, str);

	std::getline(stream, str);
	RG_ASSERT(0 == str.compare("DATA START"), "Wrong entry, DATA START expected");
}

#ifdef USE_MPI

/**
 * \brief read one line from file opened by single process
 * \param[in] fh file handler
 * \param[out] str string
 * \param[out] delim symbol to stop reading
 */
long getLineMPI(MPI_File fh, std::string& str, char delim = '\n');

/**
 * \brief Load header using MPI by single process
 * \return offset of data in file
 * \param[in] filename
 * \param[out] size array of sizes in header
 * \param[out] components number of comonents
 * \param[out] fmt format
 */
template <typename I>
long loadHeaderMPI(const std::string& filename, Dim3D<I>& size, I& components, format& fmt) {
	MPI_File fh;
	MPI_CHECK(MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
	MPI_CHECK(MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL));
	
	long offset = 0;
	std::string str;
	
	offset += getLineMPI(fh, str);
	RG_ASSERT(0 == str.compare("# DARRAY DATA FORMAT"), "Wrong file format");
	
	offset += getLineMPI(fh, str, ':');
	RG_ASSERT(0 == str.compare("SIZE"), "Wrong entry, SIZE expected");
	
	offset += getLineMPI(fh, str);
	sscanf(str.c_str(), "%d %d %d", &size[X], &size[Y], &size[Z]);
	
	offset += getLineMPI(fh, str, ':');	
	RG_ASSERT(0 == str.compare("COMPONENTS"), "Wrong entry, COMPONENTS expected");
	offset += getLineMPI(fh, str);
	sscanf(str.c_str(), "%d", &components);
	
	offset += getLineMPI(fh, str, ':');	
	RG_ASSERT(0 == str.compare("FORMAT"), "Wrong entry, FORMAT expected");

	offset += getLineMPI(fh, str, ' ');
	offset += getLineMPI(fh, str);
	if (0 == str.compare("text")) fmt = TEXT;
	else if (0 == str.compare("binary")) fmt = BINARY;
	else {
		RG_ASSERT(0, "Wrong format");
	}

	offset += getLineMPI(fh, str);
	RG_ASSERT(0 == str.compare("DATA START"), "Wrong entry, DATA START expected");
	
	MPI_CHECK(MPI_File_close(&fh));
	return offset;
}

#endif

/**
 * \brief Basic header saver (own DArray format)
 * \param[out] stream
 * \param[in] size
 * \param[in] components
 * \param[in] fmt format
 */
template <typename I>
void writeHeader(std::iostream& stream, const Dim3D<I>& size, const I components, const format fmt) {
	stream << "# DARRAY DATA FORMAT" << std::endl;
	stream << "SIZE: " << size[X] << " " << size[Y] << " " << size[Z] << std::endl;
	stream << "COMPONENTS: " << components << std::endl;
	stream << "FORMAT: ";
	if (fmt == TEXT) stream << "text" << std::endl;
	else if (fmt == BINARY) stream << "binary" << std::endl;
	stream << "DATA START" << std::endl;
}

} // namespace rgio

} // namespace rgrid

#endif // RG_IO_H
