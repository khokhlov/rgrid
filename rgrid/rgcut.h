/**
 * \file
 * \brief Partitioning of DArray
 */

#ifndef RG_CUT_H
#define RG_CUT_H

#include "rgrid/types.h"

#include <iostream>
#include <vector>
#include <map>

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
	 * \param[in] sz is entire size in nodes in each direction
	 * \param[in] pt is number of parts to devide
	 */
	void setCutParams(const Dim3D<I> sz, const Dim3D<I> pt) {
		Dim3D<std::vector<I> > width;
		for (CartDir d = X; d != ALL_DIRS; ++d) {
			for (I i = 0; i != pt[d]; ++i) {
				I w = sz[d] / pt[d] + (i < sz[d] % pt[d] ? 1 : 0);
				if (w > 0) width[d].push_back(w);
			}
		}
		setCutParams(width);
	}
	/**
	 * \param[in] width is width of each part
	 */
	void setCutParams(const Dim3D<std::vector<I> >& width) {
		this->width = width;
		for (CartDir d = X; d != ALL_DIRS; ++d) {
			locator[d].clear();
			RG_ASSERT(width[d].size() > 0, "Empty dimension in setCutParams");
			parts[d] = width[d].size();
			start[d].resize(width[d].size());
			start[d].at(0) = 0;
			locator[d][width[d].at(0)] = 0;
			for (typename std::vector<I>::size_type i = 1; i < width[d].size(); ++i) {
				start[d].at(i) = start[d].at(i-1) + width[d].at(i-1);
				locator[d][start[d].at(i) + width[d].at(i)] = i;
			}
			size[d] = start[d].at(width[d].size()-1) + width[d].at(width[d].size()-1);
		}
		allParts = parts[X] * parts[Y] * parts[Z];
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
	I vecInd(const I lin, CartDir dir) const {
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
	 * \brief Get local index from global
	 * \param[in] dir direction to search
	 * \param[in] contIdx global index of node
	 */
	I locateIndex(const CartDir dir, const I contIdx) const {
		I partNum = locatePart(dir, contIdx);
		return contIdx - partOrigin(dir, partNum);
	}
	/**
	 * \brief Locate index of part in specified direction in which specified node contained
	 * \param[in] dir direction of parts
	 * \param[in] contIdx global index of node
	 */
	I locatePart(const CartDir dir, const I contIdx) const {
		return locator[dir].upper_bound(contIdx)->second;
	}
	/**
	 * \brief Get origin of part
	 */
	I partOrigin(const CartDir dir, const I partNum) const {
		return start[dir].at(partNum);
	}
	/**
	 * \brief Get origin of part
	 */
	Dim3D<I> partOrigin(const Dim3D<I>& partNum) const {
		return Dim3D<I>(start[X].at(partNum[X]), start[Y].at(partNum[Y]), start[Z].at(partNum[Z]));
	}
	/**
	 * \brief Get part nodes number
	 */
	I partNodes(const CartDir dir, const I partNum) const {
		return width[dir].at(partNum);
	}
	/**
	 * \brief Get part nodes number
	 * \param partPos coordinate of part
	 */
	Dim3D<I> partNodes(const Dim3D<I>& partPos) const {
		return Dim3D<I>(width[X].at(partPos[X]), width[Y].at(partPos[Y]), width[Z].at(partPos[Z]));
	}
	/**
	 * \brief Get new copy of width array
	 */
	Dim3D<std::vector<I> > getWidth() const {
		return width;
	}
private:
	/* number of parts in each direction */
	Dim3D<I> parts;
	/* total number of parts */
	I allParts;
	/* number of nodes in each direction in ALL parts */
	Dim3D<I> size;
	/* origins of blocks */
	Dim3D<std::vector<I> > start;
	/* width of blocks */
	Dim3D<std::vector<I> > width;
	/* key - first index of next block, val - current block num */
	Dim3D<std::map<I, I> > locator;
};

} /* namespace rgrid */

#endif /* RG_CUT_H */
