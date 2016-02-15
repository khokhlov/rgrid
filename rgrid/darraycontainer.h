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
	DArray<T, I> getDArray() const;
	// get specific part of DArray
	DArray<T, I>& getDArrayPart(const I partNum) {
		RG_ASSERT(dArray.size() > partNum, "Out of range");
		return dArray[partNum];
	}
	
	DArray<T, I>& getDArrayPart(const I px, I py, I pz) {
		return getDArrayPart(pz * parts[X] * parts[Y] + py * parts[X] + px);
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
	
	I parts[ALL_DIRS];
	std::vector<DArray<T, I> > dArray;
};

} // namespace rgrid

#endif // DARRAY_CONTAINER_H