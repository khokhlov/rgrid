#include "rgrid/darraycontainer.h"

namespace rgrid {

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I numParts) {
	// split first along longest axis from Y or Z
	if (da.localSize[Y] < da.localSize[Z]) setDArray(da, 1, 1, numParts);
	else setDArray(da, 1, numParts, 1);
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I>& da, const I px, const I py, const I pz) {
	for (I d = X; d != ALL_DIRS; ++d) {
		ls[d] = da.localSize(d);
		ml[d] = ls[d] / parts[d];
		iml[d] = ls[d] % parts[d];
	}
	RG_ASSERT(ls[X] >= px || ls[Y] >= py || ls[Z] >= pz, "Not enough nodes to split in this way");
	dArray.resize(px * py * pz);
	parts[X] = px;
	parts[Y] = py;
	parts[Z] = pz;
	
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
			dArray[ind](i2, j2, k2, cn) = da(o[X] + i2, o[Y] + j2, o[Z] + k2);
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
	I nc = dArray[0].da.getNC();
	da.alloc(nc);
	for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		I ind = k * parts[X] * parts[Y] + j * parts[X] + i;
		for (I cn = 0; cn != nc; ++cn)
		for (I k2 = 0; k2 != dArray[ind].localSize[Z]; ++k2)
		for (I j2 = 0; j2 != dArray[ind].localSize[Y]; ++j2)
		for (I i2 = 0; i2 != dArray[ind].localSize[X]; ++i2) {
			I orig[ALL_DIRS] = { 
				dArray[ind].origin(X) - da.origin(X),
				dArray[ind].origin(Y) - da.origin(Y),
				dArray[ind].origin(Z) - da.origin(Z)};
			da(orig[X] + i2, orig[Y] + j2, orig[Z] + k2) = dArray[ind](i2, j2, k2, cn);
		}
	}
	
}

template <typename T, typename I>
void DArrayContainer<T, I>::synchronize() {
	for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		DArray<T, I>& da = getDArrayPart(i, j, k);
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
	// TODO
	RG_ASSERT(0, "Not implemented");
}
#endif // USE_OPENCL

} // namespace rgrid

