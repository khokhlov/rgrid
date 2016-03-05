#ifndef DARRAY_CONTAINER_H
#define DARRAY_CONTAINER_H

#include "rgrid/darray.h"
#include "rgrid/rgcut.h"

#include <vector>

namespace rgrid {

/*
 * DArrayContainer splits DArray to work with small parts
 */
template <typename T, typename I>
class DArrayContainer : public RGCut<I> {
public:
	DArrayContainer() : RGCut<I>() {}
	DArrayContainer(const DArray<T, I>& da, const I numParts) : RGCut<I>() 
	{
		setDArray(da, numParts); 
	}
	DArrayContainer(const DArray<T, I>& da, const I px, const I py, const I pz) : RGCut<I>() 
	{ 
		setDArray(da, px, py, pz); 
	}
		
	void setDArray(const DArray<T, I>& da, const I numParts); 
	void setDArray(const DArray<T, I>& da, const I px, const I py, const I pz);
	
	// get entire DArray
	void getDArray(DArray<T, I>&) const;
	// get specific part of DArray
	DArray<T, I>& getDArrayPart(const I partNum) {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	const DArray<T, I>& getDArrayPart(const I partNum) const {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	DArray<T, I>& getDArrayPart(const I i, const I j, const I k) {
		return dArray[RGCut<I>::linInd(i, j, k)];
	}
	
	DArray<T, I>& operator()(const I px, const I py, const I pz) {
		return getDArrayPart(px, py, pz);
	}
	
	// write line (all nodes on X direction) with coordinates y and z into stream
	void writeLine(std::basic_iostream<char>& stream, I y, I z, rgio::format fmt) const;
	
	// return node in position (i, j, k) and number cn
	inline T& getNode(const I i, const I j, const I k, const I cn) {
		I pn[ALL_DIRS] = { RGCut<I>::locatePart(X, i), RGCut<I>::locatePart(Y, j), RGCut<I>::locatePart(Z, k) };
		I idx[ALL_DIRS] = { RGCut<I>::locateIndex(X, i), RGCut<I>::locateIndex(Y, j), RGCut<I>::locateIndex(Z, k) };
		return getDArrayPart(pn[X], pn[Y], pn[Z]).val(idx[X], idx[Y], idx[Z], cn);
	}
	
	/*
	 * fill fhost nodes of all DArray parts
	 * with values from adjacent DArrays
	 */
	void fillGhost();
#ifdef USE_OPENCL
	// The same as fillGhost(), but without copy to device
	void fillGhostCL();
#endif // USE_OPENCL
	
private:	
	std::vector<DArray<T, I> > dArray;
};

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I numParts) {
	setDArray(da, 1, 1, numParts);
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I px, const I py, const I pz) {
	RGCut<I>::setCutParams(da.localSize(X), da.localSize(Y), da.localSize(Z), px, py, pz);
	dArray.resize(px * py * pz);
	for (I k = 0; k != pz; ++k)
	for (I j = 0; j != py; ++j)
	for (I i = 0; i != px; ++i) {
		I ind = RGCut<I>::linInd(i, j, k);
		I o[ALL_DIRS] = { RGCut<I>::partOrigin(X, i), RGCut<I>::partOrigin(Y, j), RGCut<I>::partOrigin(Z, k) };
		dArray[ind].resize(da.size(X), da.size(Y), da.size(Z),
		              RGCut<I>::partNodes(X, i), RGCut<I>::partNodes(Y, j), RGCut<I>::partNodes(Z, k),
		              da.origin(X) + o[X], da.origin(Y) + o[Y], da.origin(Z) + o[Z],
		              da.ghost(X), da.ghost(Y), da.ghost(Z));
		dArray[ind].alloc(da.getNC());
		for (I cn = 0; cn != da.getNC(); ++cn)
		for (I k2 = 0; k2 != RGCut<I>::partNodes(Z, k); ++k2)
		for (I j2 = 0; j2 != RGCut<I>::partNodes(Y, j); ++j2)
		for (I i2 = 0; i2 != RGCut<I>::partNodes(X, i); ++i2) {
			dArray[ind](i2, j2, k2, cn) = da(o[X] + i2, o[Y] + j2, o[Z] + k2, cn);
		}
	}
}
	
template <typename T, typename I>
void DArrayContainer<T, I>::getDArray(DArray<T, I>& da) const {
	da.resize(
		dArray[0].size(X), dArray[0].size(Y), dArray[0].size(Z),
		RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
		dArray[0].origin(X), dArray[0].origin(Y), dArray[0].origin(Z),
		dArray[0].ghost(X), dArray[0].ghost(Y), dArray[0].ghost(Z));
	I nc = dArray[0].getNC();
	da.alloc(nc);
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		I ind = RGCut<I>::linInd(i, j, k);
		for (I cn = 0; cn != nc; ++cn)
		for (I k2 = 0; k2 != dArray[ind].localSize(Z); ++k2)
		for (I j2 = 0; j2 != dArray[ind].localSize(Y); ++j2)
		for (I i2 = 0; i2 != dArray[ind].localSize(X); ++i2) {
			I orig[ALL_DIRS] = { 
				dArray[ind].origin(X) - da.origin(X),
				dArray[ind].origin(Y) - da.origin(Y),
				dArray[ind].origin(Z) - da.origin(Z)};
			da(orig[X] + i2, orig[Y] + j2, orig[Z] + k2, cn) = dArray[ind](i2, j2, k2, cn);
		}
	}
	
}

template <typename T, typename I>
void DArrayContainer<T, I>::writeLine(std::basic_iostream<char>& stream, I y, I z, rgio::format fmt) const {
	I pn[ALL_DIRS] = { 0, RGCut<I>::locatePart(Y, y), RGCut<I>::locatePart(Z, z) };
	I idx[ALL_DIRS] = { 0, RGCut<I>::locateIndex(Y, y), RGCut<I>::locateIndex(Z, z) };
	for (pn[X] = 0; pn[X] != RGCut<I>::numParts(X); ++pn[X]) {
		dArray[RGCut<I>::linInd(pn[X], pn[Y], pn[Z])].writeLine(stream, idx[Y], idx[Z], fmt);
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::fillGhost() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		DArray<T, I>& da = getDArrayPart(i, j, k);
		if (i != 0)                       da.copyGhost(getDArrayPart(i-1, j, k), X, SIDE_LEFT);
		else                              da.fillGhost(X, SIDE_LEFT);
		if (i != RGCut<I>::numParts(X)-1) da.copyGhost(getDArrayPart(i+1, j, k), X, SIDE_RIGHT);
		else                              da.fillGhost(X, SIDE_RIGHT);
		if (j != 0)                       da.copyGhost(getDArrayPart(i, j-1, k), Y, SIDE_LEFT);
		else                              da.fillGhost(Y, SIDE_LEFT);
		if (j != RGCut<I>::numParts(Y)-1) da.copyGhost(getDArrayPart(i, j+1, k), Y, SIDE_RIGHT);
		else                              da.fillGhost(Y, SIDE_RIGHT);
		if (k != 0)                       da.copyGhost(getDArrayPart(i, j, k-1), Z, SIDE_LEFT);
		else                              da.fillGhost(Z, SIDE_LEFT);
		if (k != RGCut<I>::numParts(Z)-1) da.copyGhost(getDArrayPart(i, j, k+1), Z, SIDE_RIGHT);
		else                              da.fillGhost(Z, SIDE_RIGHT);
	}
}

#ifdef USE_OPENCL
template <typename T, typename I>
void DArrayContainer<T, I>::fillGhostCL() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		DArray<T, I>& da = getDArrayPart(i, j, k);
		if (i != 0)                       da.copyGhostCL(getDArrayPart(i-1, j, k), X, SIDE_LEFT);
		else                              da.fillGhostCL(X, SIDE_LEFT);
		if (i != RGCut<I>::numParts(X)-1) da.copyGhostCL(getDArrayPart(i+1, j, k), X, SIDE_RIGHT);
		else                              da.fillGhostCL(X, SIDE_RIGHT);
		if (j != 0)                       da.copyGhostCL(getDArrayPart(i, j-1, k), Y, SIDE_LEFT);
		else                              da.fillGhostCL(Y, SIDE_LEFT);
		if (j != RGCut<I>::numParts(Y)-1) da.copyGhostCL(getDArrayPart(i, j+1, k), Y, SIDE_RIGHT);
		else                              da.fillGhostCL(Y, SIDE_RIGHT);
		if (k != 0)                       da.copyGhostCL(getDArrayPart(i, j, k-1), Z, SIDE_LEFT);
		else                              da.fillGhostCL(Z, SIDE_LEFT);
		if (k != RGCut<I>::numParts(Z)-1) da.copyGhostCL(getDArrayPart(i, j, k+1), Z, SIDE_RIGHT);
		else                              da.fillGhostCL(Z, SIDE_RIGHT);
	}
}
#endif // USE_OPENCL


} // namespace rgrid

#endif // DARRAY_CONTAINER_H
