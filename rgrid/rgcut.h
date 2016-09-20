/**
 * \file
 * \brief Partitioning of DArray
 */

#ifndef RG_CUT_H
#define RG_CUT_H

#include "rgrid/types.h"

namespace rgrid {
	
/**
 * \brief Class represents partitioning of rectangular structure (DArray)
 * 
 * It automatically divides 3d rectangular area into specified number of parts
 * with number of nodes that differ not more than one in each direction.
 * 
 * \tparam I type of indexes (i.e. int, long)
 */
template <typename I>
class RGCut {
public:
	/**
	 * Use setCutParams() if RGCut created by this constructor 
	 */
	RGCut() {
		setCutParams(1, 1, 1, 1, 1, 1);
	}
	/**
	 * \param[in] szX,szY,szZ is entire size in each direction
	 * \param[in] ptX,ptY,ptZ is number of parts to devide
	 */
	RGCut(I szX, I szY, I szZ, I ptX, I ptY, I ptZ) { setCutParams(szX, szY, szZ, ptX, ptY, ptZ); }
	/**
	 * \param[in] sz is entire size in nodes in each direction
	 * \param[in] pt is number of parts to devide
	 */
	RGCut(I sz[ALL_DIRS], I pt[ALL_DIRS]) { setCutParams(sz, pt); }
	/**
	 * \param[in] sz is entire size in nodes in each direction
	 * \param[in] pt is number of parts to devide
	 */
	void setCutParams(I const sz[ALL_DIRS], I const pt[ALL_DIRS]) {
		setCutParams(Dim3D<I>(sz), Dim3D<I>(pt));
	}
	/**
	 * \param[in] sz is entire size in nodes in each direction
	 * \param[in] pt is number of parts to devide
	 */
	void setCutParams(const Dim3D<I> sz, const Dim3D<I> pt) {
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			parts[d] = pt[d];
			size[d] = sz[d];
			ml[d] = sz[d] / pt[d];
			iml[d] = sz[d] % pt[d];
		}
		allParts = pt[X] * pt[Y] * pt[Z];
	}
	/**
	 * \brief Set partitioning params
	 * \param[in] szX,szY,szZ is entire size in each direction
	 * \param[in] ptX,ptY,ptZ is number of parts to devide
	 */
	void setCutParams(I szX, I szY, I szZ, I ptX, I ptY, I ptZ) {
		I sz[ALL_DIRS] = {szX, szY, szZ};
		I pt[ALL_DIRS] = {ptX, ptY, ptZ};
		setCutParams(sz, pt);
	}
	/**
	 * \brief number of all parts
	 */
	I numParts() const { return allParts; }
	/**
	 * \brief number of parts in each direction
	 */
	I numParts(CartDir dir) const { return parts[dir]; }
	/**
	 * \brief number of parts in each direction
	 */
	Dim3D<I> numParts3() const { return parts; }
	/** 
	 * \brief Get index of part in linear array
	 */
	I linInd(const I i, const I j, const I k) const {
		return k * parts[X] * parts[Y] + j * parts[X] + i;
	}
	/**
	 * \brief Convert linear index to vector (i,j,k)
	 */
	void vecInd(const I lin, I vec[ALL_DIRS]) const {
		vec[Z] = lin / (parts[X] * parts[Y]);
		const I lin2 = lin % (parts[X] * parts[Y]);
		vec[Y] = lin2 / parts[X];
		vec[X] = lin2 % parts[X];
	}
	/** 
	 * \brief Convert linear index to i, j or k
	 */
	void vecInd(const I lin, CartDir dir) const {
		I vec[ALL_DIRS];
		vecInd(lin, vec);
		return vec[dir]; 
	}
	/**
	 * \brief Number of all nodes
	 */
	I numNodes() const { return size[X] * size[Y] * size[Z]; }
	/** 
	 * \brief Number of nodes in all parts in each direction
	 */
	I numNodes(CartDir dir) const { return size[dir]; }
	/** 
	 * \brief Number of nodes in all parts in each direction
	 */
	Dim3D<I> numNodes3() const { return size; }
	/** 
	 * \brief Locate index of node in its part
	 * \param[in] dir direction to search
	 * \param[in] contIdx index of node
	 */
	I locateIndex(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx % (ml[dir] + 1);
		} else {
			return ml[dir] - 1 - ((size[dir] - 1 - contIdx) % ml[dir]);
		}
	}
	/** 
	 * \brief Locate index of part in specified direction in which specified node contained
	 * \param[in] dir direction of parts
	 * \param[in] contIdx index of node in container
	 */
	I locatePart(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx / (ml[dir] + 1);
		} else {
			return parts[dir] - 1 - ((size[dir] - 1 - contIdx) / ml[dir]);
		}
	}
	/** 
	 * \brief Get origin of part
	 */
	I partOrigin(const CartDir dir, const I partNum) const {
		return partNum <= iml[dir] ? (ml[dir] + 1) * partNum : size[dir] - ml[dir] * (parts[dir] - partNum);
	}
	/**
	 * \brief Get origin of part
	 */
	Dim3D<I> partOrigin(const Dim3D<I>& partNum) const {
		Dim3D<I> r;
		for (CartDir d = X; d < ALL_DIRS; ++d) {
			r[d] = partNum[d] <= iml[d] ? (ml[d] + 1) * partNum[d] : size[d] - ml[d] * (parts[d] - partNum[d]);
		}
		return r;
	}
	/** 
	 * \brief Get part nodes number
	 */
	I partNodes(const CartDir dir, const I partNum) const {
		return partNum < iml[dir] ? ml[dir] + 1 : ml[dir];
	}
	/** 
	 * \brief Get part nodes number
	 * \param partPos coordinate of part
	 */
	Dim3D<I> partNodes(const Dim3D<I>& partPos) const {
		Dim3D<I> r;
		for (CartDir d = X; d < ALL_DIRS; ++d) {
			r[d] = partPos[d] < iml[d] ? ml[d] + 1 : ml[d];
		}
		return r;
	}
	/**
	 * \brief Make the same RGcut as specified
	 */
	RGCut(const RGCut<I>& rhs) {
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			parts[d] = rhs.parts[d];
			size[d] = rhs.size[d];
			ml[d] = rhs.ml[d];
			iml[d] = rhs.iml[d];
		}
		allParts = rhs.allParts;
	}
	/**
	 * \brief Make the same RGcut as specified
	 */
	RGCut<I>& operator=(const RGCut<I>& rhs) {
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			parts[d] = rhs.parts[d];
			size[d] = rhs.size[d];
			ml[d] = rhs.ml[d];
			iml[d] = rhs.iml[d];
		}
		allParts = rhs.allParts;
		return *this;
	}
private:
	/* number of parts in each direction */
	Dim3D<I> parts;
	/* total number of parts */
	I allParts;
	/* number of nodes in each direction in ALL parts */
	Dim3D<I> size;
	/* minimal length - length in nodes of last parts in each direction */
	Dim3D<I> ml;
	/* index of part before that number of nodes equals ml[dir] + 1 */
	Dim3D<I> iml;
};

} /* namespace rgrid */

#endif /* RG_CUT_H */
