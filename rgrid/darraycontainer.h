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
	DArrayContainer(const DArray<T, I> &da, const I numParts) : RGCut<I>() {
		setDArray(da, numParts);
	}
	DArrayContainer(const DArray<T, I> &da, const I px, const I py, const I pz) : RGCut<I>() {
		setDArray(da, px, py, pz);
	}

	void setDArray(const DArray<T, I> &da, const I numParts);
	void setDArray(const DArray<T, I> &da, const I px, const I py, const I pz);

	// get entire DArray
	void getDArray(DArray<T, I> &) const;
	// get specific part of DArray
	DArray<T, I> &getDArrayPart(const I partNum) {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	const DArray<T, I> &getDArrayPart(const I partNum) const {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray[partNum];
	}
	DArray<T, I> &getDArrayPart(const I i, const I j, const I k) {
		return dArray[RGCut<I>::linInd(i, j, k)];
	}

	DArray<T, I> &operator()(const I px, const I py, const I pz) {
		return getDArrayPart(px, py, pz);
	}

	// return node in position (i, j, k) and number cn
	inline T &getNode(const I i, const I j, const I k, const I cn) {
		I pn[ALL_DIRS] = { RGCut<I>::locatePart(X, i), RGCut<I>::locatePart(Y, j), RGCut<I>::locatePart(Z, k) };
		I idx[ALL_DIRS] = { RGCut<I>::locateIndex(X, i), RGCut<I>::locateIndex(Y, j), RGCut<I>::locateIndex(Z, k) };
		return getDArrayPart(pn[X], pn[Y], pn[Z]).val(idx[X], idx[Y], idx[Z], cn);
	}

	/*
	 * copy rect region with size: sz, sy, sz and start: ox, oy, oz to buffer
	 */
	void getSubArray(I ox, I oy, I oz, I sx, I sy, I sz, std::vector<T> &buffer);
	/*
	 * copy rect region with size: sz, sy, sz and start: ox, oy, oz from buffer
	 * buffer must contain sx * sy * sz * cn elements
	 */
	void setSubArray(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer);
	/*
	 * fill boundary ghost nodes with values from nearest nodes
	 */
	void fillGhost();
	/*
	 * fill fhost nodes of all DArray parts with values from adjacent DArrays
	 */
	void sync();

	// write line (all nodes on X direction) with coordinates y and z into stream
	void writeLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) const;

	void saveData(std::iostream &stream, const rgio::format fmt) const;

#ifdef USE_OPENCL
	/*
	 * The same as fillGhost(), but without copy to host
	 */
	void fillGhostCL();
	/*
	 * The same as sync(), but without copy to host
	 */
	void syncCL();
#endif // USE_OPENCL

private:
	std::vector<DArray<T, I> > dArray;
};

template <typename T, typename I>
void DArrayContainer<T, I>::writeLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) const {
	Dim3D<I> pn(0, RGCut<I>::locatePart(Y, y), RGCut<I>::locatePart(Z, z));
	Dim3D<I> idx(0, RGCut<I>::locateIndex(Y, y), RGCut<I>::locateIndex(Z, z));
	for (pn[X] = 0; pn[X] != RGCut<I>::numParts(X); ++pn[X]) {
		dArray.at(RGCut<I>::linInd(pn[X], pn[Y], pn[Z])).writeLine(stream, cn, idx[Y], idx[Z], fmt);
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::saveData(std::iostream &stream, const rgio::format fmt) const {
	Dim3D<I> size(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z));
	I nc = dArray.at(0).getNC();
	rgio::writeHeader(stream, size, nc, fmt);
	for (I cn = 0; cn != nc; ++cn) {
		for (I z = 0; z != RGCut<I>::numNodes(Z); ++z) {
			for (I y = 0; y != RGCut<I>::numNodes(Y); ++y) {
				writeLine(stream, cn, y, z, fmt);
			}
		}
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I> &da, const I numParts) {
	setDArray(da, 1, 1, numParts);
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I> &da, const I px, const I py, const I pz) {
	RGCut<I>::setCutParams(da.localSize(X), da.localSize(Y), da.localSize(Z), px, py, pz);
	dArray.resize(px * py * pz);
	for (I k = 0; k != pz; ++k)
		for (I j = 0; j != py; ++j)
			for (I i = 0; i != px; ++i) {
				I ind = RGCut<I>::linInd(i, j, k);
				I o[ALL_DIRS] = { RGCut<I>::partOrigin(X, i), RGCut<I>::partOrigin(Y, j), RGCut<I>::partOrigin(Z, k) };
				dArray.at(ind).resize(da.globalSize(X), da.globalSize(Y), da.globalSize(Z),
				                   RGCut<I>::partNodes(X, i), RGCut<I>::partNodes(Y, j), RGCut<I>::partNodes(Z, k),
				                   da.origin(X) + o[X], da.origin(Y) + o[Y], da.origin(Z) + o[Z],
				                   da.ghost(X), da.ghost(Y), da.ghost(Z),
				                   da.getNC());
				for (I cn = 0; cn != da.getNC(); ++cn)
					for (I k2 = 0; k2 != RGCut<I>::partNodes(Z, k); ++k2)
						for (I j2 = 0; j2 != RGCut<I>::partNodes(Y, j); ++j2)
							for (I i2 = 0; i2 != RGCut<I>::partNodes(X, i); ++i2) {
								dArray.at(ind)(i2, j2, k2, cn) = da(o[X] + i2, o[Y] + j2, o[Z] + k2, cn);
							}
			}
}

template <typename T, typename I>
void DArrayContainer<T, I>::getDArray(DArray<T, I> &da) const {
	da.resize(
	    dArray.at(0).globalSize(X), dArray.at(0).globalSize(Y), dArray.at(0).globalSize(Z),
	    RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
	    dArray.at(0).origin(X), dArray.at(0).origin(Y), dArray.at(0).origin(Z),
	    dArray.at(0).ghost(X), dArray.at(0).ghost(Y), dArray.at(0).ghost(Z),
	    dArray.at(0).getNC());
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				I ind = RGCut<I>::linInd(i, j, k);
				for (I cn = 0; cn != dArray.at(0).getNC(); ++cn)
					for (I k2 = 0; k2 != dArray.at(ind).localSize(Z); ++k2)
						for (I j2 = 0; j2 != dArray.at(ind).localSize(Y); ++j2)
							for (I i2 = 0; i2 != dArray.at(ind).localSize(X); ++i2) {
								I orig[ALL_DIRS] = {
									dArray.at(ind).origin(X) - da.origin(X),
									dArray.at(ind).origin(Y) - da.origin(Y),
									dArray.at(ind).origin(Z) - da.origin(Z)
								};
								da(orig[X] + i2, orig[Y] + j2, orig[Z] + k2, cn) = dArray.at(ind)(i2, j2, k2, cn);
							}
			}

}

template <typename T, typename I>
void DArrayContainer<T, I>::getSubArray(I ox, I oy, I oz, I sx, I sy, I sz, std::vector<T> &buffer) {
	I nc = dArray.front().getNC();
	buffer.clear();
	buffer.reserve(sx * sy * sz * nc);
	for (I cn = 0; cn != nc; ++cn)
		for (I z = 0; z != sz; ++z)
			for (I y = 0; y != sy; ++y)
				for (I x = 0; x != sx; ++x) {
					I xp = ox + x;
					I yp = oy + y;
					I zp = oz + z;
					I partX = RGCut<I>::locatePart(X, xp);
					I partY = RGCut<I>::locatePart(Y, yp);
					I partZ = RGCut<I>::locatePart(Z, zp);
					I nodeX = RGCut<I>::locateIndex(X, xp);
					I nodeY = RGCut<I>::locateIndex(Y, yp);
					I nodeZ = RGCut<I>::locateIndex(Z, zp);
					buffer.push_back(dArray[RGCut<I>::linInd(partX, partY, partZ)].val(nodeX, nodeY, nodeZ, cn));
				}
}

template <typename T, typename I>
void DArrayContainer<T, I>::setSubArray(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer) {
	I nc = dArray.at(0).getNC();
	RG_ASSERT(buffer.size() == sx * sy * sz * nc, "Wrong buffer size");
	I bufInd = 0;
	for (I cn = 0; cn != nc; ++cn)
		for (I z = 0; z != sz; ++z)
			for (I y = 0; y != sy; ++y)
				for (I x = 0; x != sx; ++x) {
					I xp = ox + x;
					I yp = oy + y;
					I zp = oz + z;
					I partX = RGCut<I>::locatePart(X, xp);
					I partY = RGCut<I>::locatePart(Y, yp);
					I partZ = RGCut<I>::locatePart(Z, zp);
					I nodeX = RGCut<I>::locateIndex(X, xp);
					I nodeY = RGCut<I>::locateIndex(Y, yp);
					I nodeZ = RGCut<I>::locateIndex(Z, zp);
					dArray.at(RGCut<I>::linInd(partX, partY, partZ)).val(nodeX, nodeY, nodeZ, cn) = buffer.at(bufInd++);
				}
}

template <typename T, typename I>
void DArrayContainer<T, I>::sync() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				DArray<T, I> &da = getDArrayPart(i, j, k);
				if (i != 0) {
					da.copyGhost(getDArrayPart(i - 1, j, k), X, SIDE_LEFT);
				}
				if (i != RGCut<I>::numParts(X) - 1) {
					da.copyGhost(getDArrayPart(i + 1, j, k), X, SIDE_RIGHT);
				}
				if (j != 0) {
					da.copyGhost(getDArrayPart(i, j - 1, k), Y, SIDE_LEFT);
				}
				if (j != RGCut<I>::numParts(Y) - 1) {
					da.copyGhost(getDArrayPart(i, j + 1, k), Y, SIDE_RIGHT);
				}
				if (k != 0) {
					da.copyGhost(getDArrayPart(i, j, k - 1), Z, SIDE_LEFT);
				}
				if (k != RGCut<I>::numParts(Z) - 1) {
					da.copyGhost(getDArrayPart(i, j, k + 1), Z, SIDE_RIGHT);
				}
			}
}

template <typename T, typename I>
void DArrayContainer<T, I>::fillGhost() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				DArray<T, I> &da = getDArrayPart(i, j, k);
				if (i == 0) {
					da.fillGhost(X, SIDE_LEFT);
				}
				if (j == RGCut<I>::numParts(X) - 1) {
					da.fillGhost(X, SIDE_RIGHT);
				}
				if (j == 0) {
					da.fillGhost(Y, SIDE_LEFT);
				}
				if (j == RGCut<I>::numParts(Y) - 1) {
					da.fillGhost(Y, SIDE_RIGHT);
				}
				if (k == 0) {
					da.fillGhost(Z, SIDE_LEFT);
				}
				if (k == RGCut<I>::numParts(Z) - 1) {
					da.fillGhost(Z, SIDE_RIGHT);
				}
			}
}

#ifdef USE_OPENCL
template <typename T, typename I>
void DArrayContainer<T, I>::syncCL() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				DArray<T, I> &da = getDArrayPart(i, j, k);
				if (i != 0) {
					da.copyGhostCL(getDArrayPart(i - 1, j, k), X, SIDE_LEFT);
				}
				if (i != RGCut<I>::numParts(X) - 1) {
					da.copyGhostCL(getDArrayPart(i + 1, j, k), X, SIDE_RIGHT);
				}
				if (j != 0) {
					da.copyGhostCL(getDArrayPart(i, j - 1, k), Y, SIDE_LEFT);
				}
				if (j != RGCut<I>::numParts(Y) - 1) {
					da.copyGhostCL(getDArrayPart(i, j + 1, k), Y, SIDE_RIGHT);
				}
				if (k != 0) {
					da.copyGhostCL(getDArrayPart(i, j, k - 1), Z, SIDE_LEFT);
				}
				if (k != RGCut<I>::numParts(Z) - 1) {
					da.copyGhostCL(getDArrayPart(i, j, k + 1), Z, SIDE_RIGHT);
				}
			}
}

template <typename T, typename I>
void DArrayContainer<T, I>::fillGhostCL() {
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				DArray<T, I> &da = getDArrayPart(i, j, k);
				if (i == 0) {
					da.fillGhostCL(X, SIDE_LEFT);
				}
				if (i == RGCut<I>::numParts(X) - 1) {
					da.fillGhostCL(X, SIDE_RIGHT);
				}
				if (j == 0) {
					da.fillGhostCL(Y, SIDE_LEFT);
				}
				if (j == RGCut<I>::numParts(Y) - 1) {
					da.fillGhostCL(Y, SIDE_RIGHT);
				}
				if (k == 0) {
					da.fillGhostCL(Z, SIDE_LEFT);
				}
				if (k == RGCut<I>::numParts(Z) - 1) {
					da.fillGhostCL(Z, SIDE_RIGHT);
				}
			}
}
#endif // USE_OPENCL


} // namespace rgrid

#endif // DARRAY_CONTAINER_H
