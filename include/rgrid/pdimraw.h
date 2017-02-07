/**
 * \file
 * \brief Sizes of DArray
 */


#ifndef RG_PDIM_RAW_H
#define RG_PDIM_RAW_H

#include "rgrid/types.h"
#include "rgrid/range.h"

namespace rgrid {

/**
 * \brief This class responsible for sizes of DArray's
 *
 * This class is helpful because it can be used inside OpenCL kernels (When T is POD type, PDimRaw also POD type)
 * \tparam T type of indexes (i.e. int, long)
 * \sa DArray
 */
template <typename T>
struct PDimRaw {
	/**
	 * \brief convert 3D indexes to index in linear array
	 * \param[in] x,y,z coordinates
	 * \param[in] cn component number
	 * \return index in linear array
	 */
	T ind(const T x, const T y, const T z, const T cn) const;
	/**
	 * \sa rgrid::PDimRaw< T >::ind(const T x, const T y, const T z, const T cn)
	 */
	T ind(const T x, const T y, const T cn) const { return ind(x, y, static_cast<T>(0), cn); }
	/**
	 * \sa rgrid::PDimRaw< T >::ind(const T x, const T y, const T z, const T cn) const
	 */
	T ind(const T x, const T cn) const { return ind(x, static_cast<T>(0), cn); }

	/**
	 * \brief convert 3D indexes to index in linear array
	 * \param[in] d1,d2,d3 different dimensions in any order
	 * \param[in] i1,i2,i3 coordinates in the same order as d1,d2,d3
	 * \param[in] cn component number
	 * \return index in linear array
	 */
	T ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3, const T cn) const;
	/**
	 * \sa rgrid::PDimRaw< T >::ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3, const T cn) const
	 */
	T ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T cn) const { return ind(i1, d1, i2, d2, static_cast<T>(0), Z, cn); }
	/**
	 * \sa rgrid::PDimRaw< T >::ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3, const T cn) const
	 */
	T ind(const T i1, const CartDir d1, const T cn) const { return ind(i1, d1, static_cast<T>(0), Y, cn); }

	/**
	 * \brief Get number of ghost nodes on specific side in specific direction where other directions has single node
	 * \param[in] d direction
	 * \return Number of ghost nodes
	 */
	T ghost(const CartDir d = X) const {
		return m_ghost_size[d];
	}
	/**
	 * \brief Get number of all nodes in bigger rectangular structure in specific direction
	 * \param[in] d direction
	 * \return number of nodes
	 */
	T globalSize(const CartDir d) const {
		return m_global_size[d];
	}
	/**
	 * \brief Get number of all nodes in bigger rectangular structure
	 * \return Number of nodes
	 */
	T globalSize() const {
		return m_global_size_all;
	}
	/**
	 * \brief Get number of all nodes in current rectangular structure in specific direction
	 * \param[in] d direction
	 * \return Number of nodes
	 */
	T localSize(const CartDir d) const {
		return m_local_size[d];
	}
	/**
	 * \brief Get number of all nodes in current rectangular structure
	 * \return Number of nodes
	 */
	T localSize() const {
		return m_local_size_all;
	}

	/**
	 * \brief Get origin of current rectangular structure inside bigger rectangular structure in specific direction
	 * \param[in] d direction
	 * \return origin
	 */
	T origin(const CartDir d) const {
		return m_origin[d];
	}
	
	Dim3D<T> origin() const {
		return m_origin;
	}

	/**
	 * \brief Get number of all nodes in current rectangular structure in specific direction includding ghost nodes
	 * \param[in] d direction
	 * \return Number of nodes
	 */
	T localGhostSize(const CartDir d) const {
		return m_local_ghost_size[d];
	}
	/**
	 * \brief Get number of all nodes in local rect struct including ghost nodes
	 * \return Number of nodes
	 */
	T localGhostSize() const {
		return m_local_ghost_size_all;
	}
	/**
	 * \brief Distance in linear array between two nodes with difference only in one index
	 * \note suitable for MPI ralated stuff
	 * \param[in] d direction where indexes differs
	 * \return Distance in nodes
	 */
	T localStride(const CartDir d) const {
		return m_local_stride[d];
	}

	/**
	 * \brief Convert index in linear array to 3D indexes
	 * \param[in] index linear index
	 * \param[in] dir direction of desired index
	 * \return One of indexes in 3D array
	 */
	T localToIJK(const T index, const CartDir dir);
	/**
	 * \brief The same as ind(), but assuming that there are no ghost nodes
	 */
	T indNoGhost(const T x, const T y, const T z, const T cn);

	/**
	 * \brief Convert index in local rect struct to index in bigger rect struct
	 * \param[in] index in local rect struct
	 * \param[in] dir direction
	 * \return index in bigger rect struct
	 */
	T localToGlobal(const T index, const CartDir dir);
	/**
	 * \brief Convert index in global rect struct to index in linear array of local rect struct
	 * \param x,y,z 3D index
	 * \param cn component number
	 * \return Index in local linear array
	 */
	T localFromGlobal(const T x, const T y, const T z, const T cn);
	/**
	 * \sa rgrid::PDimRaw< T >::localFromGlobal(const T x, const T y, const T z, const T cn)
	 */
	T localFromGlobal(const T x, const T y, const T cn) {
		return localFromGlobal(x, y, static_cast<T>(0), cn);
	}
	/**
	 * \sa rgrid::PDimRaw< T >::PDimRaw< T >::localFromGlobal(const T x, const T y, const T z, const T cn)
	 */
	T localFromGlobal(const T x, const T cn) {
		return localFromGlobal(x, static_cast<T>(0), cn);
	}

	/**
	 * \brief Equivalent to localFromGlobal
	 */
	T indGlobal(const T x, const T y, const T z, const T cn) {
		return localFromGlobal(x, y, z, cn);
	}
	/**
	 * \brief Equivalent to localFromGlobal
	 */
	T indGlobal(const T x, const T y, const T cn) {
		return localFromGlobal(x, y, cn);
	}
	/**
	 * \brief Equivalent to localFromGlobal
	 */
	T indGlobal(const T x, const T cn) {
		return localFromGlobal(x, cn);
	}

	/**
	 * \brief Check is global index inside local rect struct
	 * \param[in] index in global rect struct
	 * \param[in] dir direction of index
	 */
	bool check(const T index, const CartDir dir) const;
	bool check(const T i, const T j, const T k, const T cn) const {
		return check(i, X) && check(j, Y) && check(k, Z) && (0 <= cn && cn < m_nc);
	}
	/**
	 * \brief The same as check(), but index can be inside ghost nodes of local rect struct
	 */
	bool checkg(const T index, const CartDir dir) const;
	/**
	 * \brief The same as check(), but index can be inside ghost nodes of local rect struct
	 */
	bool checkg(const T i, const T j, const T k, const T cn) const {
		return checkg(i, X) && checkg(j, Y) && checkg(k, Z) && (0 <= cn && cn < m_nc);
	}
	/**
	 * \brief Check is index in local rect struct in one of corners of global rect struct
	 * \param[in] d1,d2,d3 different directions in any order
	 * \param[in] i1,i2,i3 coordinates in local rect struct in the same order as d1,d2,d3
	 */
	bool isCorner(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3);
	/**
	 * \brief Check is index in local rect struct in one of edges of global rect struct
	 * \param[in] d1,d2,d3 different directions in any order
	 * \param[in] i1,i2,i3 coordinates in local rect struct in the same order as d1,d2,d3
	 */
	bool isEdge(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3);
	/**
	 * \brief Check is local index inside ghost nodes of local rect struct
	 * \param[in] i,j,k indexes in local rect struct
	 */
	bool isGhost(const T i, const T j, const T k);
	/**
	 * \sa rgrid::PDimRaw< I >::isGhost(const T i, const T j, const T k)
	 */
	bool isGhost(const T i, const T j) {
		return isGhost(i, j, static_cast<T>(0));
	}
	/**
	 * \sa rgrid::PDimRaw< I >::isGhost(const T i, const T j, const T k)
	 */
	bool isGhost(const T i) {
		return isGhost(i, static_cast<T>(0));
	}
	/**
	 * \brief Check is local rect struct placed on the edge of global rect structure
	 * \param[in] d direction
	 * \param[in] s side
	 */
	bool isOnFace(const CartDir d, const CartSide s);

	/**
	 * \brief Get number of dimensions of rect struct
	 */
	CartDim dim() const;

	/**
	 * \brief Get number of components inside each node
	 */
	T getNC() const {
		return m_nc;
	};

	/**
	 * \brief Check is this rect structure has the same sizes as the other
	 */
	bool isEqual(const PDimRaw<T>& rhs) const {
		return isEqualNoGhost(rhs) && (m_ghost_size == rhs.m_ghost_size);
	}

	/**
	 * \brief Check is this rect structure has the same sizes as the other,
	 * but not check size of ghost nodes
	 */
	bool isEqualNoGhost(const PDimRaw<T>& rhs) const {
		const PDimRaw<T>& lhs = *this;
		const bool is1 = lhs.m_global_size == rhs.m_global_size;
		const bool is2 = lhs.m_local_size == rhs.m_local_size;
		const bool is3 = lhs.m_origin == rhs.m_origin;
		const bool is4 = lhs.m_nc == rhs.m_nc;
		return is1 && is2 && is3 && is4;
	}
	
	Range<T> getRange() const {
		return Range<T>(m_origin, m_local_size);
	}
	
	Range<T> getGhostRange() const {
		return Range<T>(m_origin-m_ghost_size, m_origin+m_local_ghost_size);
	}
	
	/// size of bigger rect structure (global size)
	Dim3D<T> m_global_size;
	/// multiplied global sizes in all directions
	T m_global_size_all;
	/// Distance in linear array between two nodes with difference only in one index for global rect struct
	Dim3D<T> m_global_stride;

	/// size of local rect struct
	Dim3D<T> m_local_size;
	/// number of all local nodes
	T m_local_size_all;
	/// Distance in linear array between two nodes with difference only in one index for local rect struct
	Dim3D<T> m_local_stride;

	/// number of ghost sizes in each direction
	Dim3D<T> m_ghost_size;
	/** m_local_ghost_size = m_local_size + 2 * m_ghost_size */
	Dim3D<T> m_local_ghost_size;
	/// multiplied ghost sizes in all directions
	T m_local_ghost_size_all;
	/// origins of local rect struct in global rect struct
	Dim3D<T> m_origin;

	/// Number of components
	T m_nc;
}; // PDimRaw

template <typename T>
CartDim PDimRaw<T>::dim() const {
	CartDim d = DIM_3D;
	if (m_global_size[Z] == 1) {
		d = DIM_2D;
	}
	if (m_global_size[Y] == 1) {
		d = DIM_1D;
	}
	return d;
}

template <typename T>
bool PDimRaw<T>::isOnFace(const CartDir d, const CartSide s) {
	bool ret = false;
	if (s == SIDE_LEFT && origin(d) == 0) {
		ret = true;
	}
	if (s == SIDE_RIGHT && origin(d) + localSize(d) == globalSize(d)) {
		ret = true;
	}
	return ret;
}

template <typename T>
bool PDimRaw<T>::isGhost(const T i, const T j, const T k) {
	return i < 0 || i >= localSize(X) || j < 0 || j >= localSize(Y) || k < 0 || k >= localSize(Z);
}

/// Check is two PDimRaw equal
template <typename T>
bool operator==(const PDimRaw<T> &lhs, const PDimRaw<T> &rhs) {
	return lhs.isEqual(rhs);
}

/// Check is two PDimRaw not equal
template <typename T>
bool operator!=(const PDimRaw<T> &lhs, const PDimRaw<T> &rhs) {
	return !(lhs == rhs);
}

template <typename T>
bool PDimRaw<T>::isCorner(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3) {
	bool a = (localToGlobal(i1, d1) == static_cast<T>(0)) || (localToGlobal(i1, d1) == m_global_size[d1] - 1);
	bool b = (localToGlobal(i2, d2) == static_cast<T>(0)) || (localToGlobal(i2, d2) == m_global_size[d2] - 1);
	bool c = (localToGlobal(i3, d3) == static_cast<T>(0)) || (localToGlobal(i3, d3) == m_global_size[d3] - 1);
	return a && b && c;
}

template <typename T>
bool PDimRaw<T>::isEdge(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3) {
	bool a = (localToGlobal(i1, d1) == static_cast<T>(0)) || (localToGlobal(i1, d1) == m_global_size[d1] - 1);
	bool b = (localToGlobal(i2, d2) == static_cast<T>(0)) || (localToGlobal(i2, d2) == m_global_size[d2] - 1);
	bool c = (localToGlobal(i3, d3) == static_cast<T>(0)) || (localToGlobal(i3, d3) == m_global_size[d3] - 1);
	return (a && b) || (b && c) || (c && a);
}

template <typename T>
bool PDimRaw<T>::check(const T index, const CartDir dir) const {
	return index >= m_origin[dir] && index < m_origin[dir] + m_local_size[dir];
}

template <typename T>
bool PDimRaw<T>::checkg(const T index, const CartDir dir) const {
	return index >= -m_ghost_size[dir] + m_origin[dir] && index < m_origin[dir] + m_local_size[dir] + m_ghost_size[dir];
}


template <typename T>
T PDimRaw<T>::localToGlobal(const T index, const CartDir dir) {
	return index + m_origin[dir];
}

template <typename T>
T PDimRaw<T>::indNoGhost(const T x, const T y, const T z, const T cn) {
	return x + y * m_local_size[X] + z * m_local_size[X] * m_local_size[Y] + cn * m_local_size_all;
}

template <typename T>
T PDimRaw<T>::localFromGlobal(const T x, const T y, const T z, const T cn) {
	return ind(x - m_origin[X], y - m_origin[Y], z - m_origin[Z], cn);
}

template <typename T>
T PDimRaw<T>::ind(const T x, const T y, const T z, const T cn) const {
	return (x + m_ghost_size[X]) +
	       (y + m_ghost_size[Y]) * m_local_stride[Y] +
	       (z + m_ghost_size[Z]) * m_local_stride[Z] +
	       cn * m_local_ghost_size_all;
}

template <typename T>
T PDimRaw<T>::ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3, const T cn) const {
	return (i1 + m_ghost_size[d1]) * m_local_stride[d1] + (i2 + m_ghost_size[d2]) * m_local_stride[d2] + (i3 + m_ghost_size[d3]) * m_local_stride[d3] + cn * m_local_ghost_size_all;
}

template <typename T>
T PDimRaw<T>::localToIJK(const T index, const CartDir dir) {
	T k = index / m_local_stride[Z] - ghost(Z);
	T j = (index - (k + ghost(Z)) * m_local_stride[Z]) / m_local_stride[Y] - ghost(Y);
	T i = (index - (k + ghost(Z)) * m_local_stride[Z] - (j + ghost(Y)) * m_local_stride[Y]) / m_local_stride[X] - ghost(X);
	if (dir == Z) {
		return k;
	} else if (dir == Y) {
		return j;
	} else if (dir == X) {
		return i;
	}
	return -1;
}

} // namesace rgrid

#endif // RG_PDIM_RAW_H

