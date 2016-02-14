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
	DArrayContainer();
	DArrayContainer(const DArray<T, I>& da, const I numParts);
	DArrayContainer(const DArray<T, I>& da, const I px, const I py, const I pz);
	
	~DArrayContainer();
	
	void setDArray(const DArray<T, I>& da, const I numParts);
	void setDArray(const DArray<T, I>& da, const I px, const I py, const I pz);
	
	// get entire DArray
	DArray<T, I> getDArray() const;
	// get specific part of DArray
	DArray<T, I>& getDArrayPart(const I partNum);
	DArray<T, I>& getDArrayPart(const I px, I py, I pz);
	
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
	
	I px, py, pz;
	std::vector<DArray<T, I> > dArray;
};

} // namespace rgrid

#endif // DARRAY_CONTAINER_H