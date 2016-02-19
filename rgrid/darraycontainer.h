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
	DArray<T, I>& getDArrayPart(const I i, const I j, const I k) {
		I ind = k * parts[X] * parts[Y] + j * parts[X] + i;
		return dArray[ind];
	}
	
	DArray<T, I>& operator()(const I px, const I py, const I pz) {
		return getDArrayPart(px, py, pz);
	}
	
	I numParts() const { return dArray.size(); }
	I numParts(CartDir dir) const { return parts[dir]; }
	
	// num nodes of entire dArray
	I size() const { return ls[X] * ls[Y] * ls[Z]; }
	// num nodes of entire dArray in direction dir
	I size(CartDir dir) const { return ls[dir]; }
	
	// write line (all nodes on X direction) with coordinates y and z into stream
	void writeLine(std::basic_iostream<char>& stream, I y, I z, rgio::format fmt) const;
	
	// find index in part of darrays
	// dir - direction of indexes and parts
	// contIdx - index of node in container
	inline I locateIndex(const CartDir dir, const I contIdx) {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx % (ml[dir] + 1);
		} else {
			return ml[dir] - 1 - ((ls[dir] - 1 - contIdx) % ml[dir]);
		}
	}
	
	// find what part contains this node
	// dir - direction of parts
	// contIdx - index of node in container
	inline I locatePart(const CartDir dir, const I contIdx) {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx / (ml[dir] + 1);
		} else {
			return parts[dir] - 1 - ((ls[dir] - 1 - contIdx) / ml[dir]);
		}
	}
	
	// return node in position (i, j, k) and number cn
	inline T& getNode(const I i, const I j, const I k, const I cn) {
		I pn[ALL_DIRS] = { locatePart(i), locatePart(j), locatePart(k) };
		I idx[ALL_DIRS] = { locateIndex(i), locateIndex(j), locateIndex(k) };
		return getDArrayPart(pn[X], pn[Y], pn[Z]).val(idx[X], idx[Y], idx[Z], cn);
	}
	
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