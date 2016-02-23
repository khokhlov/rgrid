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
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	const DArray<T, I>& getDArrayPart(const I partNum) const {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	DArray<T, I>& getDArrayPart(const I i, const I j, const I k) {
		return dArray[ind(i, j, k)];
	}
	
	I const ind(const I i, const I j, const I k) const {
		return k * parts[X] * parts[Y] + j * parts[X] + i;
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
	inline I locateIndex(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx % (ml[dir] + 1);
		} else {
			return ml[dir] - 1 - ((ls[dir] - 1 - contIdx) % ml[dir]);
		}
	}
	
	// find what part contains this node
	// dir - direction of parts
	// contIdx - index of node in container
	inline I locatePart(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx / (ml[dir] + 1);
		} else {
			return parts[dir] - 1 - ((ls[dir] - 1 - contIdx) / ml[dir]);
		}
	}
	
	// return node in position (i, j, k) and number cn
	inline T& getNode(const I i, const I j, const I k, const I cn) {
		I pn[ALL_DIRS] = { locatePart(X, i), locatePart(Y, j), locatePart(Z, k) };
		I idx[ALL_DIRS] = { locateIndex(X, i), locateIndex(Y, j), locateIndex(Z, k) };
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

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I numParts) {
	// split first along longest axis from Y or Z
	if (da.localSize[Y] < da.localSize[Z]) setDArray(da, 1, 1, numParts);
	else setDArray(da, 1, numParts, 1);
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I px, const I py, const I pz) {
	parts[X] = px;
	parts[Y] = py;
	parts[Z] = pz;
	for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
		ls[d] = da.localSize(d);
		ml[d] = ls[d] / parts[d];
		iml[d] = ls[d] % parts[d];
	}
	RG_ASSERT(ls[X] >= px || ls[Y] >= py || ls[Z] >= pz, "Not enough nodes to split in this way");
	dArray.resize(px * py * pz);
	
	I nc = da.getNC();
	
	for (I k = 0; k != pz; ++k)
	for (I j = 0; j != py; ++j)
	for (I i = 0; i != px; ++i) {
		I ind = k * px * py + j * px + i;
		I lsp[ALL_DIRS] = {
			i < iml[X] ? ml[X] + 1 : ml[X],
			j < iml[Y] ? ml[Y] + 1 : ml[Y],
			k < iml[Z] ? ml[Z] + 1 : ml[Z]};
		I o[ALL_DIRS] = {
			i <= iml[X] ? (ml[X] + 1) * i: ls[X] - ml[X] * (px - i),
			j <= iml[Y] ? (ml[Y] + 1) * j: ls[Y] - ml[Y] * (py - j),
			k <= iml[Z] ? (ml[Z] + 1) * k: ls[Z] - ml[Z] * (pz - k)};
		dArray[ind].resize(
			da.size(X), da.size(Y), da.size(Z),
			lsp[X], lsp[Y], lsp[Z],
			da.origin(X) + o[X], da.origin(Y) + o[Y], da.origin(Z) + o[Z],
       			da.ghost(X), da.ghost(Y), da.ghost(Z));
		
		dArray[ind].alloc(nc);
		
		for (I cn = 0; cn != nc; ++cn)
		for (I k2 = 0; k2 != lsp[Z]; ++k2)
		for (I j2 = 0; j2 != lsp[Y]; ++j2)
		for (I i2 = 0; i2 != lsp[X]; ++i2) {
			dArray[ind](i2, j2, k2, cn) = da(o[X] + i2, o[Y] + j2, o[Z] + k2, cn);
		}
	}
	
}
	
template <typename T, typename I>
void DArrayContainer<T, I>::getDArray(DArray<T, I>& da) const {
	da.resize(
		dArray[0].size(X), dArray[0].size(Y), dArray[0].size(Z),
		ls[X], ls[Y], ls[Z],
		dArray[0].origin(X), dArray[0].origin(Y), dArray[0].origin(Z),
		dArray[0].ghost(X), dArray[0].ghost(Y), dArray[0].ghost(Z));
	I nc = dArray[0].getNC();
	da.alloc(nc);
	for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		I idx = ind(i, j, k);
		for (I cn = 0; cn != nc; ++cn)
		for (I k2 = 0; k2 != dArray[idx].localSize(Z); ++k2)
		for (I j2 = 0; j2 != dArray[idx].localSize(Y); ++j2)
		for (I i2 = 0; i2 != dArray[idx].localSize(X); ++i2) {
			I orig[ALL_DIRS] = { 
				dArray[idx].origin(X) - da.origin(X),
				dArray[idx].origin(Y) - da.origin(Y),
				dArray[idx].origin(Z) - da.origin(Z)};
			da(orig[X] + i2, orig[Y] + j2, orig[Z] + k2, cn) = dArray[idx](i2, j2, k2, cn);
		}
	}
	
}

template <typename T, typename I>
void DArrayContainer<T, I>::writeLine(std::basic_iostream<char>& stream, I y, I z, rgio::format fmt) const {
	I pn[ALL_DIRS] = { 0, locatePart(Y, y), locatePart(Z, z) };
	I idx[ALL_DIRS] = { 0, locateIndex(Y, y), locateIndex(Z, z) };
	for (pn[X] = 0; pn[X] != parts[X]; ++pn[X]) {
		dArray[ind(pn[X], pn[Y], pn[Z])].writeLine(stream, idx[Y], idx[Z], fmt);
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::synchronize() {
	for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		if (i != 0) 
			(*this)(i, j, k).copyGhost((*this)(i-1, j, k), X, SIDE_LEFT);
		if (i != parts[X]-1) 
			(*this)(i, j, k).copyGhost((*this)(i+1, j, k), X, SIDE_RIGHT);
		if (j != 0) 
			(*this)(i, j, k).copyGhost((*this)(i, j-1, k), Y, SIDE_LEFT);
		if (j != parts[Y]-1) 
			(*this)(i, j, k).copyGhost((*this)(i, j+1, k), Y, SIDE_RIGHT);
		if (k != 0) 
			(*this)(i, j, k).copyGhost((*this)(i, j, k-1), Z, SIDE_LEFT);
		if (k != parts[Z]-1) 
			(*this)(i, j, k).copyGhost((*this)(i, j, k+1), Z, SIDE_RIGHT);
	}
}

#ifdef USE_OPENCL
template <typename T, typename I>
void DArrayContainer<T, I>::clDeviceSynchronize() {
	for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		if (i != 0) 
			(*this)(i, j, k).clCopyGhost((*this)(i-1, j, k), X, SIDE_LEFT);
		if (i != parts[X]-1) 
			(*this)(i, j, k).clCopyGhost((*this)(i+1, j, k), X, SIDE_RIGHT);
		if (j != 0) 
			(*this)(i, j, k).clCopyGhost((*this)(i, j-1, k), Y, SIDE_LEFT);
		if (j != parts[Y]-1) 
			(*this)(i, j, k).clCopyGhost((*this)(i, j+1, k), Y, SIDE_RIGHT);
		if (k != 0) 
			(*this)(i, j, k).clCopyGhost((*this)(i, j, k-1), Z, SIDE_LEFT);
		if (k != parts[Z]-1) 
			(*this)(i, j, k).clCopyGhost((*this)(i, j, k+1), Z, SIDE_RIGHT);
	}
}
#endif // USE_OPENCL


} // namespace rgrid

#endif // DARRAY_CONTAINER_H
