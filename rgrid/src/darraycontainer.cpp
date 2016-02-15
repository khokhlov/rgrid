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
	I ls[ALL_DIRS] = { da.localSize(X), da.localSize(Y), da.localSize(Z) };
	RG_ASSERT(ls[X] >= px || ls[Y] >= py || ls[Z] >= pz, "Not enough nodes to split in this way");
	dArray.resize(px * py * pz);
	parts[X] = px;
	parts[Y] = py;
	parts[Z] = pz;
	
	I nc = da.getNC();
	
	I ml[ALL_DIRS] = { ls[X] / px, ls[Y] / py, ls[Z] / pz };
	I iml[ALL_DIRS] = { ls[X] % px, ls[Y] % py, ls[Z] % pz };
	
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
DArray<T, I> DArrayContainer<T, I>::getDArray() const {
	// TODO
	RG_ASSERT(0, "Not implemented");
}

template <typename T, typename I>
void DArrayContainer<T, I>::synchronize() {
	// TODO
	RG_ASSERT(0, "Not implemented");
	/*for (I k = 0; k != parts[Z]; ++k)
	for (I j = 0; j != parts[Y]; ++j)
	for (I i = 0; i != parts[X]; ++i) {
		I ind = k * parts[X] * parts[Y] + j * parts[X] + i;
		dArray[ind].copyGhost(dArray[ind]);
	}*/
}

#ifdef USE_OPENCL
template <typename T, typename I>
void DArrayContainer<T, I>::clDeviceSynchronize() {
	// TODO
	RG_ASSERT(0, "Not implemented");
}
#endif // USE_OPENCL

} // namespace rgrid

