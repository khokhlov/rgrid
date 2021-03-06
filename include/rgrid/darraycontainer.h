/**
 * \file
 * \brief One partitioning (local) step for DArray
 */

#ifndef DARRAY_CONTAINER_H
#define DARRAY_CONTAINER_H

#include "rgrid/darray.h"
#include "rgrid/rgcut.h"

#include <vector>

namespace rgrid {

/**
 * \brief DArrayContainer splits DArray to work with small parts
 *
 * One big DArray divided into smaller DArrays, so you can use OpenMP or OpenCL with them
 *
 * \tparam T type of every grid node (i.e. double, float)
 * \tparam I type of grid indexes (i.e. int, long)
 */
template <typename T, typename I>
class DArrayContainer : public RGCut<I> {
public:
	/**
	 * Call rgrid::DArrayContainer<T, I>::setDArray(const DArray<T, I> &da, const I px, const I py, const I pz) if using this constructor for initialization
	 */
	DArrayContainer() : RGCut<I>() {}
	/**
	 * \brief Create DArrayContainer from existing DArray
	 * \sa setDArray
	 * \param[in] da DArray to split
	 * \param[in] numParts number of parts
	 */
	DArrayContainer(const DArray<T, I> &da, const I numParts) : RGCut<I>() {
		setDArray(da, numParts);
	}
	/**
	 * \brief Create DArrayContainer from existing DArray
	 * \param[in] da DArray to split
	 * \param[in] px,py,pz num parts in (X,Y,Z) directions
	 */
	DArrayContainer(const DArray<T, I> &da, const I px, const I py, const I pz) : RGCut<I>() {
		setDArray(da, px, py, pz);
	}
	/**
	 * \brief take DArray and divide into smaller DArrays in one dimension
	 * \param[in] da DArray to split
	 * \param[in] numParts number of parts
	 */
	void setDArray(const DArray<T, I> &da, const I numParts);
	/**
	 * \brief take DArray and divide into smaller DArrays in 3 dims
	 * \param[in] da DArray to split
	 * \param[in] px,py,pz num parts in (X,Y,Z) directions
	 */
	void setDArray(const DArray<T, I> &da, const I px, const I py, const I pz);
	/**
	 * \brief take DArray and divide into smaller DArrays in 3 dims
	 * \param[in] da DArray to split
	 * \param[in] width of blocks
	 */
	void setDArray(const DArray<T, I> &da, const Dim3D<std::vector<I> >& width);
	/**
	 * \brief Allocate memory for DArrays
	 * \param[in] globalSize size of area in which DArrayContainer is placed
	 * \param[in] localSize size of DArrayContainer
	 * \param[in] parts number of DArray's
	 * \param[in] ghost number of ghost nodes
	 * \param[in] origin origin of current DArrayContainer
	 * \param[in] components in each node
	 */
	void setParts(
		const Dim3D<I>& globalSize,
		const Dim3D<I>& localSize,
		const Dim3D<I>& parts,
		const Dim3D<I>& ghost,
		const Dim3D<I>& origin,
		const I components);
	/**
	 * \brief Allocate memory for DArrays
	 * \param[in] globalSize size of area in which DArrayContainer is placed
	 * \param[in] localWidth width of DArray blocks inside DArrayContainer
	 * \param[in] ghost number of ghost nodes
	 * \param[in] origin origin of current DArrayContainer
	 * \param[in] components number of components in each node
	 */
	void setParts(
		const Dim3D<I>& globalSize,
		const Dim3D<std::vector<I> >& localWidth,
		const Dim3D<I>& ghost,
		const Dim3D<I>& origin,
		const I components);
	/**
	 * \brief Get entire DArray from parts in DArrayContainer
	 * \param[out] da entire DArray
	 */
	void getDArray(DArray<T, I> &da) const;
	/**
	 * \brief Get specific part of DArray.
	 *
	 * Use this function if DArray splitted only in one dimension
	 *
	 * \param[in] partNum number of part
	 * \return DArray part
	 */
	DArray<T, I> &getDArrayPart(const I partNum) {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray.at(partNum);
	}
	/**
	 * \sa rgrid::DArray<T, I>::getDArrayPart(const I partNum)
	 */
	const DArray<T, I> &getDArrayPart(const I partNum) const {
		RG_ASSERT(dArray.size() > static_cast<size_t>(partNum), "Out of range");
		return dArray.at(partNum);
	}
	/**
	 * \brief Get specific part of DArray.
	 *
	 * Use this function if DArray splitted in 3 dimensions
	 *
	 * \param[in] i,j,k numbers of parts in dimensions (X,Y,Z)
	 * \return DArray part
	 */
	DArray<T, I> &getDArrayPart(const I i, const I j, const I k) {
		return dArray.at(RGCut<I>::linInd(i, j, k));
	}

	/**
	 * \brief Equivalent to getDArrayPart()
	 */
	DArray<T, I> &operator()(const I px, const I py, const I pz) {
		return getDArrayPart(px, py, pz);
	}

	/**
	 * \brief Check is index present in current container.
	 *
	 * Ghost nodes are not taken in consideration.
	 * \param[in] ind index in container
	 */
	bool isPresent(const Dim3D<I>& ind) const {
		const Dim3D<I>& e = RGCut<I>::numNodes3();
		return
			0 <= ind.x && ind.x < e.x &&
			0 <= ind.y && ind.y < e.y &&
			0 <= ind.z && ind.z < e.z;
	}

	/**
	 * \brief Check is index present in current container.
	 *
	 * Include ghost nodes of container.
	 * \param[in] ind index in container
	 */
	bool isPresentGhost(const Dim3D<I>& ind) const {
		const Dim3D<I>& e = RGCut<I>::numNodes3();
		return
			-getNGhost(X) <= ind.x && ind.x < e.x + getNGhost(X) &&
			-getNGhost(Y) <= ind.y && ind.y < e.y + getNGhost(Y) &&
			-getNGhost(Z) <= ind.z && ind.z < e.z + getNGhost(Z);
	}

	/**
	 * \brief get node in entire DArray when it splitted
	 * \param[in] i,j,k coordinates of nodes in dimensions (X,Y,Z)
	 * \param[in] cn component number
	 * \return node
	 */
	T& getNode(const I i, const I j, const I k, const I cn) {
		return getNode(Dim3D<I>(i, j, k), cn);
	}
	/**
	 * \brief Get Node by index in container.
	 *
	 * Ghost nodes are not taken in consideration.
	 * \param[in] ind index in container
	 * \param[in] cn component number
	 * \return node
	 */
	T &getNode(const Dim3D<I>& ind, const I cn) {
		RG_ASSERT(isPresent(ind), "Trying to get index out of current block");
		Dim3D<I> pn(
			RGCut<I>::locatePart(X, ind.x),
			RGCut<I>::locatePart(Y, ind.y),
			RGCut<I>::locatePart(Z, ind.z));
		Dim3D<I> idx(
			RGCut<I>::locateIndex(X, ind.x),
			RGCut<I>::locateIndex(Y, ind.y),
			RGCut<I>::locateIndex(Z, ind.z));
		return getDArrayPart(pn[X], pn[Y], pn[Z]).val(idx[X], idx[Y], idx[Z], cn);
	}
	/**
	 * \brief get node in entire DArray when it splitted
	 *
	 * Include ghost nodes of container.
	 * \param[in] ind index in container
	 * \param[in] cn component number
	 * \return node
	 */
	T &getNodeGhost(const Dim3D<I>& ind, const I cn) {
		RG_ASSERT(isPresentGhost(ind), "Trying to get index out of current block");
		Dim3D<I> ind2 = ind;
		for (CartDir d = X; d != ALL_DIRS; ++d) {
			if (ind[d] < 0) ind2[d] = 0;
			if (ind[d] >= RGCut<I>::numNodes(d)) ind2[d] = RGCut<I>::numNodes(d) - 1;
		}
		Dim3D<I> pn(
			RGCut<I>::locatePart(X, ind2.x),
			RGCut<I>::locatePart(Y, ind2.y),
			RGCut<I>::locatePart(Z, ind2.z));
		Dim3D<I> idx(
			RGCut<I>::locateIndex(X, ind2.x),
			RGCut<I>::locateIndex(Y, ind2.y),
			RGCut<I>::locateIndex(Z, ind2.z));
		for (CartDir d = X; d != ALL_DIRS; ++d) {
			if (ind[d] < 0) idx[d] -= ind2[d] - ind[d];
			if (ind[d] >= RGCut<I>::numNodes(d)) idx[d] += ind[d] - ind2[d];
		}
		return getDArrayPart(pn[X], pn[Y], pn[Z]).val(idx[X], idx[Y], idx[Z], cn);
	}

	/**
	 * \brief Copy rect region from entire DArrayContainer
	 *
	 * \param[in] ox,oy,oz origins of rect region
	 * \param[in] sx,sy,sz starts of rect region
	 * \param[out] buffer rect region
	 */
	void getSubArray(I ox, I oy, I oz, I sx, I sy, I sz, std::vector<T> &buffer);
	/**
	 * \brief Copy rect region from entire DArrayContainer
	 * \param[in] ox,oy,oz origins of rect region
	 * \param[in] sx,sy,sz starts of rect region
	 * \param[in] buffer rect region
	 * \note buffer must contain sx * sy * sz * cn elements
	 */
	void setSubArray(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer);
	/**
	 * \brief The same as setSubArray(), but indexes can be in ghost nodes of entire DArrayContainer
	 */
	void setSubArrayWithGhost(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer);
	/**
	 * \brief Fill boundary ghost nodes with values from nearest nodes
	 */
	void fillGhost();
	/**
	 * \brief fill fhost nodes of all DArray parts with values from adjacent DArrays
	 */
	void sync();
	/**
	 * \brief write all nodes in X direction into stream
	 *
	 * \param[in] stream output stream
	 * \param[in] cn number of component
	 * \param[in] y,z coordiantes of line
	 * \param[in] fmt writing format
	 */
	void writeLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) const;
	/**
	 * \brief Read all nodes in X direction from stream.
	 * \param[in] stream input stream
	 * \param[in] cn number of component
	 * \param[in] y,z coordiantes of line
	 * \param[in] fmt reading format
	 */
	void readLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt);
	/**
	 * \brief save all DArrayContainer data and header into stream
	 *
	 * \param[in] stream output stream
	 * \param[in] fmt writing format
	 */
	void saveData(std::iostream &stream, const rgio::format fmt) const;
	/**
	 * \brief Append all DArrayContainer data at the end of stream
	 *
	 * \param[in] stream output stream
	 * \param[in] fmt writing format
	 */
	void appendData(std::iostream &stream, const rgio::format fmt) const;
	/**
	 * \brief Load DArrayContainer data from stream.
	 * \param[in] stream input stream
	 * \param[in] fmt reading format
	 */
	void loadData(std::iostream& stream, const rgio::format fmt);
	/**
	 * \brief Get number of ghost nodes
	 */
	Dim3D<I> getNGhost() const {
		if (dArray.empty()) return 0;
		Dim3D<I> ghost;
		for (CartDir d = X; d != ALL_DIRS; ++d) {
			ghost[d] = dArray.at(0).ghost(d);
		}
		return ghost;
	}
	/**
	 * \brief Get number of ghost nodes
	 */
	I getNGhost(CartDir d) const {
		if (dArray.empty()) return 0;
		return dArray.at(0).ghost(d);
	}
	/**
	 * \brief Get number of components
	 */
	I getNC() const {
		if (dArray.empty()) return 0;
		return dArray.at(0).getNC();
	}
#ifdef USE_OPENCL
	/**
	 * \brief The same as fillGhost(), but without copy to host
	 *
	 * When every DArray part on device it will work faster
	 */
	void fillGhostCL();
	/**
	 * \brief The same as sync(), but without copy to host
	 *
	 * When every DArray part on device it will work faster
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
void DArrayContainer<T, I>::readLine(std::iostream &stream, const I cn, const I y, const I z, const rgio::format fmt) {
	Dim3D<I> pn(0, RGCut<I>::locatePart(Y, y), RGCut<I>::locatePart(Z, z));
	Dim3D<I> idx(0, RGCut<I>::locateIndex(Y, y), RGCut<I>::locateIndex(Z, z));
	for (pn[X] = 0; pn[X] != RGCut<I>::numParts(X); ++pn[X]) {
		dArray.at(RGCut<I>::linInd(pn[X], pn[Y], pn[Z])).readLine(stream, cn, idx[Y], idx[Z], fmt);
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::saveData(std::iostream &stream, const rgio::format fmt) const {
	Dim3D<I> size(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z));
	I nc = dArray.at(0).getNC();
	rgio::writeHeader(stream, size, nc, fmt);
	appendData(stream, fmt);
}

template <typename T, typename I>
void DArrayContainer<T, I>::appendData(std::iostream &stream, const rgio::format fmt) const {
	Dim3D<I> size(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z));
	I nc = getNC();
	for (I cn = 0; cn != nc; ++cn) {
		for (I z = 0; z != size.z; ++z) {
			for (I y = 0; y != size.y; ++y) {
				writeLine(stream, cn, y, z, fmt);
			}
		}
	}
}


template <typename T, typename I>
void DArrayContainer<T, I>::loadData(std::iostream &stream, const rgio::format fmt) {
	Dim3D<I> size(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z));
	I nc = getNC();
	for (I cn = 0; cn != nc; ++cn) {
		for (I z = 0; z != size.z; ++z) {
			for (I y = 0; y != size.y; ++y) {
				readLine(stream, cn, y, z, fmt);
			}
		}
	}
}

template <typename T, typename I>
void DArrayContainer<T, I>::setDArray(const DArray<T, I> &da, const I numParts) {
	setDArray(da, 1, 1, numParts);
}

template <typename T, typename I>
void DArrayContainer<T, I>::setParts(
	const Dim3D<I>& globalSize,
	const Dim3D<std::vector<I> >& localWidth,
	const Dim3D<I>& ghost,
	const Dim3D<I>& origin,
	const I components)
{
	RGCut<I>::setCutParams(localWidth);
	dArray.resize(localWidth[X].size() * localWidth[Y].size() * localWidth[Z].size());
	for (typename std::vector<I>::size_type k = 0; k != localWidth[Z].size(); ++k)
		for (typename std::vector<I>::size_type j = 0; j != localWidth[Y].size(); ++j)
			for (typename std::vector<I>::size_type i = 0; i != localWidth[X].size(); ++i) {
				I ind = RGCut<I>::linInd(i, j, k);
				Dim3D<I> o(RGCut<I>::partOrigin(X, i), RGCut<I>::partOrigin(Y, j), RGCut<I>::partOrigin(Z, k));
				dArray.at(ind).resize(
					globalSize.x, globalSize.y, globalSize.z,
					RGCut<I>::partNodes(X, i), RGCut<I>::partNodes(Y, j), RGCut<I>::partNodes(Z, k),
					origin.x + o.x, origin.y + o.y, origin.z + o.z,
					ghost.x, ghost.y, ghost.z,
					components);
			}

}

template <typename T, typename I>
void DArrayContainer<T, I>::setParts(
	const Dim3D<I>& globalSize,
	const Dim3D<I>& localSize,
	const Dim3D<I>& parts,
	const Dim3D<I>& ghost,
	const Dim3D<I>& origin,
	const I components)
{
	Dim3D<std::vector<I> > localWidth;
	for (CartDir d = X; d != ALL_DIRS; ++d) {
		for (int i = 0; i != parts[d]; ++i) {
			I val = localSize[d] / parts[d] + (i < localSize[d] % parts[d] ? 1 : 0);
			if (val > 0)
				localWidth[d].push_back(val);
		}
	}
	setParts(globalSize, localWidth, ghost, origin, components);
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
				dArray.at(ind).resize(
					da.globalSize(X), da.globalSize(Y), da.globalSize(Z),
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
void DArrayContainer<T, I>::setDArray(const DArray<T, I> &da, const Dim3D<std::vector<I> >& width) {
	RGCut<I>::setCutParams(width);
	RG_ASSERT(da.localSize(X) == RGCut<I>::numNodes(X), "Number of nodes mismatch");
	RG_ASSERT(da.localSize(Y) == RGCut<I>::numNodes(Y), "Number of nodes mismatch");
	RG_ASSERT(da.localSize(Z) == RGCut<I>::numNodes(Z), "Number of nodes mismatch");
	dArray.resize(RGCut<I>::numParts());
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
			for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
				I ind = RGCut<I>::linInd(i, j, k);
				I o[ALL_DIRS] = { RGCut<I>::partOrigin(X, i), RGCut<I>::partOrigin(Y, j), RGCut<I>::partOrigin(Z, k) };
				dArray.at(ind).resize(
					da.globalSize(X), da.globalSize(Y), da.globalSize(Z),
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
	RG_ASSERT(!dArray.empty(), "There are no DArrays");
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
					buffer.push_back(dArray.at(RGCut<I>::linInd(partX, partY, partZ)).val(nodeX, nodeY, nodeZ, cn));
				}
}

template <typename T, typename I>
void DArrayContainer<T, I>::setSubArray(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer) {
	I nc = dArray.at(0).getNC();
	RG_ASSERT(buffer.size() == static_cast<typename std::vector<T>::size_type>(sx * sy * sz * nc), "Wrong buffer size");
	I bufInd = 0;
	for (I cn = 0; cn != nc; ++cn)
		for (I z = 0; z != sz; ++z)
			for (I y = 0; y != sy; ++y)
				for (I x = 0; x != sx; ++x) {
					Dim3D<I> pos(ox + x, oy + y, oz + z);
					Dim3D<I> part, node;
					for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
						part[d] = RGCut<I>::locatePart(d, pos[d]);
						node[d] = RGCut<I>::locateIndex(d, pos[d]);
					}
					dArray.at(RGCut<I>::linInd(part.x, part.y, part.z)).val(node.x, node.y, node.z, cn) = buffer.at(bufInd++);
				}
}

template <typename T, typename I>
void DArrayContainer<T, I>::setSubArrayWithGhost(I ox, I oy, I oz, I sx, I sy, I sz, const std::vector<T> &buffer) {
	I nc = dArray.at(0).getNC();
	RG_ASSERT(buffer.size() == static_cast<typename std::vector<T>::size_type>(sx * sy * sz * nc), "Wrong buffer size");
	I bufInd = 0;
	for (I cn = 0; cn != nc; ++cn)
		for (I z = 0; z != sz; ++z)
			for (I y = 0; y != sy; ++y)
				for (I x = 0; x != sx; ++x) {
					Dim3D<I> pos(ox + x, oy + y, oz + z);
					Dim3D<I> part, node;
					for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
						if (pos[d] < 0) {
							part[d] = 0;
							node[d] = pos[d];
						} else if (pos[d] >= RGCut<I>::numNodes(d)) {
							part[d] = RGCut<I>::numParts(d) - 1;
							node[d] = pos[d] - RGCut<I>::numNodes(d) + RGCut<I>::partNodes(d, part[d]);
						} else {
							part[d] = RGCut<I>::locatePart(d, pos[d]);
							node[d] = RGCut<I>::locateIndex(d, pos[d]);
						}
					}
					dArray.at(RGCut<I>::linInd(part.x, part.y, part.z)).val(node.x, node.y, node.z, cn) = buffer.at(bufInd++);
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
				if (i == RGCut<I>::numParts(X) - 1) {
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
