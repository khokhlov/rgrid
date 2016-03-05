#ifndef RG_CUT_H
#define RG_CUT_H

#include "rgrid/types.h"

namespace rgrid {
	
/*
 * divide 3d rectangular area into specified number of parts
 * with number of nodes that differ not more than one in each direction
 */
template <typename I>
class RGCut {
public:
	/* Use setCutParams() if RGCut created by this constructor */
	RGCut() {
		setCutParams(1, 1, 1, 1, 1, 1);
	}
	/*
	 * szX, szY, szZ is entire size in each direction
	 * ptX, ptY, ptZ is number of parts to devide
	 */
	RGCut(I szX, I szY, I szZ, I ptX, I ptY, I ptZ) { setCutParams(szX, szY, szZ, ptX, ptY, ptZ); }
	/*
	 * sz is entire size in nodes in each direction
	 * pt is number of parts to devide
	 */
	RGCut(I sz[ALL_DIRS], I pt[ALL_DIRS]) { setCutParams(sz, pt); }
	/*
	 * sz is entire size in nodes in each direction
	 * pt is number of parts to devide
	 */
	void setCutParams(I sz[ALL_DIRS], I pt[ALL_DIRS]) {
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			parts[d] = pt[d];
			size[d] = sz[d];
			ml[d] = sz[d] / pt[d];
			iml[d] = sz[d] % pt[d];
		}
		allParts = pt[X] * pt[Y] * pt[Z];
	}
	/*
	 * szX, szY, szZ is entire size in each direction
	 * ptX, ptY, ptZ is number of parts to devide
	 */
	void setCutParams(I szX, I szY, I szZ, I ptX, I ptY, I ptZ) {
		I sz[ALL_DIRS] = {szX, szY, szZ};
		I pt[ALL_DIRS] = {ptX, ptY, ptZ};
		setCutParams(sz, pt);
	}
	/* number of all parts */
	I numParts() const { return allParts; }
	/* number of parts in each direction */
	I numParts(CartDir dir) const { return parts[dir]; }
	/* index of part in linear array */
	I linInd(const I i, const I j, const I k) const {
		return k * parts[X] * parts[Y] + j * parts[X] + i;
	}
	/* number of all nodes */
	I numNodes() const { return size[X] * size[Y] * size[Z]; }
	/* number of nodes in all parts in each direction */
	I numNodes(CartDir dir) const { return size[dir]; }
	/* 
	 * locate index of node in its part
	 * dir - direction to search
	 * contIdx - index of node
	 */
	I locateIndex(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx % (ml[dir] + 1);
		} else {
			return ml[dir] - 1 - ((size[dir] - 1 - contIdx) % ml[dir]);
		}
	}
	/* 
	 * locate index of part in specified direction in which specified node contained
	 * dir - direction of parts
	 * contIdx - index of node in container
	 */
	I locatePart(const CartDir dir, const I contIdx) const {
		if (contIdx < iml[dir] * (ml[dir] + 1)) { 
			return contIdx / (ml[dir] + 1);
		} else {
			return parts[dir] - 1 - ((size[dir] - 1 - contIdx) / ml[dir]);
		}
	}
	/* get origin of part */
	I partOrigin(const CartDir dir, const I partNum) const {
		return partNum <= iml[dir] ? (ml[dir] + 1) * partNum : size[dir] - ml[dir] * (parts[dir] - partNum);
	}
	/* get part nodes number */
	I partNodes(const CartDir dir, const I partNum) const {
		return partNum < iml[dir] ? ml[dir] + 1 : ml[dir];
	}
private:
	/* number of parts in each direction */
	I parts[ALL_DIRS];
	/* total number of parts */
	I allParts;
	/* number of nodes in each direction in ALL parts */
	I size[ALL_DIRS];
	/* minimal length - length in nodes of last parts in each direction */
	I ml[ALL_DIRS];
	/* index of part before that number of nodes equals ml[dir] + 1 */
	I iml[ALL_DIRS];
};

} /* namespace rgrid */

#endif /* RG_CUT_H */
