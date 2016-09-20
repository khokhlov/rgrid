#ifndef VTK_SAVER_H
#define VTK_SAVER_H

#include "rgrid/darrayscatter.h"

#include "rgrid/utils.h"
#include "rgrid/types.h"

#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <fstream>

namespace rgrid {


template <typename T> std::string getVTKType();

#ifdef USE_MPI
template <typename T, typename I>
class VTKSaver {
public:
	VTKSaver() {}
	
	~VTKSaver()
	{
		typename std::list<Entry>::iterator it;
		for (it = entries.begin(); it != entries.end(); ++it)	{
			delete it->header;
		}
	}
	
	void setHeaderStructPoints(
		const std::string& title,
		const Dim3D<I>& dims,
		const Dim3D<T>& origin,
		const Dim3D<T>& spacing)
	{
		header << "# vtk DataFile Version 3.0" << std::endl;
		header << title << std::endl;
		header << "BINARY" << std::endl;
		header << "DATASET STRUCTURED_POINTS" << std::endl;
		header << "DIMENSIONS " << dims.x << " " << dims.y << " " << dims.z << std::endl;
		
		header << "ORIGIN " << origin.x << " " << origin.y << " " << origin.z << std::endl;
		header << "SPACING " << spacing.x << " " << spacing.y << " " << spacing.z << std::endl;
	}
	
	void setHeaderStructGrid(
		const std::string& title, 
		const Dim3D<I>& dims,
		const std::vector<Dim3D<T> >& points) 
	{
		RG_ASSERT(points.size() == dims.x * dims.y * dims.z, "Wrong sizes in VTK header");
		header << "# vtk DataFile Version 3.0" << std::endl;
		header << title << std::endl;
		header << "BINARY" << std::endl;
		header << "DATASET STRUCTURED_GRID" << std::endl;	
		header << "DIMENSIONS " << dims.x << " " << dims.y << " " << dims.z << std::endl;	
		
		header << "POINTS " << points.size() << " " << getVTKType<T>() << std::endl;
		
		for (typename std::vector<Dim3D<T> >::size_type i = 0; i != points.size(); ++i) {
			for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
				T val = points.at(i)[d];
				char* valPtr = (char*)(&val);
				for (int j = 0; j != sizeof(T); ++j) {
					header << valPtr[j]; // as big endian
				}
			}
		}
		header << std::endl;
	}
	
	void setHeaderRectlinGrid(
		const std::string& title,
		const Dim3D<I>& dims,
		const Dim3D<std::vector<T> >& points)
	{
		RG_ASSERT(dims.x == points.x.size(), "Wrong size X in VTK header");
		RG_ASSERT(dims.y == points.y.size(), "Wrong size Y in VTK header");
		RG_ASSERT(dims.z == points.z.size(), "Wrong size Z in VTK header");
		
		header << "# vtk DataFile Version 3.0" << std::endl;
		header << title << std::endl;
		header << "BINARY" << std::endl;
		header << "DATASET RECTILINEAR_GRID" << std::endl;
		header << "DIMENSIONS " << dims.x << " " << dims.y << " " << dims.z << std::endl;
		
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			char c = d + 'X';
			header << c << "_COORDINATES " << points[d].size() << " " << getVTKType<T>() << std::endl;

			for (typename std::vector<T>::size_type i = 0; i != points[d].size(); ++i) {
				T val = points[d].at(i);
				char* valPtr = (char*)(&val);
				for (int j = 0; j != sizeof(T); ++j) {
					header << val[j]; // as big endian
				}
			}
			header << std::endl;
		}
	}
	
	enum DataPlacement {
		POINT_DATA,
		CELL_DATA
	};
	
	void appendData(const std::string& name, DArrayScatter<T, I>& das, const DataPlacement dp) {		
		Entry e;
		e.header = new std::stringstream;
		*(e.header) << (dp == POINT_DATA ? "POINT_DATA " : "CELL_DATA ") << das.numNodes() << std::endl;
		*(e.header) << "SCALARS " << name << " float 1" << std::endl;
		*(e.header) << "LOOKUP_TABLE default" << std::endl;
		e.das = &das;
		entries.push_back(e);
	}
	
	void save(std::string filename) {
#ifdef USE_MPI
		if (rgmpi::worldRank() == 0) {
#endif
			// save main header
			std::ofstream ofs;
			ofs.open(filename.c_str());
			ofs << header.rdbuf();
			header.rdbuf()->pubseekpos(0);
			ofs.close();
#ifdef USE_MPI
		}
#endif
		typename std::list<Entry>::iterator it;
		for (it = entries.begin(); it != entries.end(); ++it) {
#ifdef USE_MPI
			if (rgmpi::worldRank() == 0) {
#endif
				std::ofstream ofs;
				ofs.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
				ofs << it->header->rdbuf();
				it->header->rdbuf()->pubseekpos(0);
				ofs.close();
#ifdef USE_MPI
			}
#endif
			it->das->appendData(filename);
		}
	}
		
private:
	
	struct Entry {
		DArrayScatter<T, I>* das;
		std::stringstream* header;
	};
	
	std::stringstream header;
	std::list<Entry> entries;
	
	VTKSaver(VTKSaver const &);
	VTKSaver &operator=(VTKSaver const &);
};
#endif

} // namespace rgrid

#endif // VTK_SAVER_H
