#ifndef DARRAY_CONTAINER_H
#define DARRAY_CONTAINER_H

#include "rgrid/darray.h"

#include <vector>

namespace rgrid {

/*
 * DArrayContainer splits DArray to work with small parts
 */
template <typename T, typename I>
class DArrayContainer {
public:
	DArrayContainer() {}
	DArrayContainer(const DArray<T, I>& da, const I numParts) { setDArray(da, numParts); }
	DArrayContainer(const DArray<T, I>& da, const I px, const I py, const I pz) { setDArray(da, px, py, pz); };
	
	~DArrayContainer() {}
	
	void setDArray(const DArray<T, I>& da, const I numParts); 
	void setDArray(const DArray<T, I>& da, const I px, const I py, const I pz);
	
	// get entire DArray
	void getDArray(DArray<T, I>&) const;
	// get specific part of DArray
	DArray<T, I>& getDArrayPart(const I partNum) {
		RG_ASSERT(dArray.size() > partNum, "Out of range");
		return dArray[partNum];
	}
	DArray<T, I>& getDArrayPart(const I& i, const I& j, const I& k) {
		I ind = k * parts[X] * parts[Y] + j * parts[X] + i;
		return dArray[ind];
	}
	
	DArray<T, I>& getDArrayPart(const I px, I py, I pz) {
		return getDArrayPart(pz * parts[X] * parts[Y] + py * parts[X] + px);
	}
	
	DArray<T, I>& operator()(const I px, I py, I pz) {
		return getDArrayPart(px, py, pz);
	}
	
	I numParts() const { return dArray.size(); }
	// synchronize ghost nodes
	void synchronize();
#ifdef USE_OPENCL
	// synchronize ghost nodes between all devices without copy to host
	void clDeviceSynchronize();
#endif // USE_OPENCL
	
private:
	DArrayContainer(const DArrayContainer<T, I>& rhs);
	DArrayContainer& operator=(const DArrayContainer<T, I>& rhs);
	
	// local size of entire darray
	I ls[ALL_DIRS];
	
	// number of darrays in this direction
	I parts[ALL_DIRS];
	
	// if part index in direction dir < iml[dir] then number of nodes in part = ml[dir] + 1
	// if index >= iml[dir] then number of nodes = ml[dir]
	I ml[ALL_DIRS];
	I iml[ALL_DIRS];
	
	std::vector<DArray<T, I> > dArray;
};

} // namespace rgrid

#endif // DARRAY_CONTAINER_H