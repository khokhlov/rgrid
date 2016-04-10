#ifndef RG_PDIM_RAW_H
#define RG_PDIM_RAW_H

namespace rgrid {

#include "rgrid/types.h"

template <typename T>
struct PDimRaw {
	T ind(const T x, const T y, const T z, const T cn) const;
	T ind(const T x, const T y, const T cn) const { return ind(x, y, static_cast<T>(0), cn); }
	T ind(const T x, const T cn) const { return ind(x, static_cast<T>(0), cn); }

	T ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3, const T cn) const;
	T ind(const T i1, const CartDir d1, const T i2, const CartDir d2, const T cn) const { return ind(i1, d1, i2, d2, static_cast<T>(0), Z, cn); }
	T ind(const T i1, const CartDir d1, const T cn) const { return ind(i1, d1, static_cast<T>(0), Y, cn); }

	T ghost(const CartDir d) const {
		return m_ghost_size[d];
	}
	T ghost() const {
		return ghost(X);
	}

	T globalSize(const CartDir d) const {
		return m_global_size[d];
	}
	T globalSize() const {
		return m_global_size_all;
	}

	T localSize(const CartDir d) const {
		return m_local_size[d];
	}
	T localSize() const {
		return m_local_size_all;
	}

	T origin(const CartDir d) const {
		return m_origin[d];
	}
	T origin() const {
		return origin(X);
	}

	T localGhostSize(const CartDir d) const {
		return m_local_ghost_size[d];
	}
	T localGhostSize() const {
		return m_local_ghost_size_all;
	}
	T localStride(const CartDir d) const {
		return m_local_stride[d];
	}

	T localToIJK(const T index, const CartDir dir);
	T indNoGhost(const T x, const T y, const T z, const T cn);

	T localToGlobal(const T index, const CartDir dir);
	T localFromGlobal(const T x, const T y, const T z, const T cn);
	T localFromGlobal(const T x, const T y, const T cn) {
		return localFromGlobal(x, y, static_cast<T>(0), cn);
	}
	T localFromGlobal(const T x, const T cn) {
		return localFromGlobal(x, static_cast<T>(0), cn);
	}

	T indGlobal(const T x, const T y, const T z, const T cn) {
		return localFromGlobal(x, y, z, cn);
	}
	T indGlobal(const T x, const T y, const T cn) {
		return localFromGlobal(x, y, cn);
	}
	T indGlobal(const T x, const T cn) {
		return localFromGlobal(x, cn);
	}

	bool check(const T index, const CartDir dir);
	bool check(const T i) {
		return check(i, X);
	}
	bool check(const T i, const T j) {
		return check(i, X) && check(j, Y);
	}
	bool check(const T i, const T j, const T k) {
		return check(i, X) && check(j, Y) && check(k, Z);
	}
	bool checkg(const T index, const CartDir dir);
	bool checkg(const T i) {
		return checkg(i, X);
	}
	bool checkg(const T i, const T j) {
		return checkg(i, X) && checkg(j, Y);
	}
	bool checkg(const T i, const T j, const T k) {
		return checkg(i, X) && checkg(j, Y) && checkg(k, Z);
	}

	bool isCorner(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3);
	bool isEdge(const T i1, const CartDir d1, const T i2, const CartDir d2, const T i3, const CartDir d3);
	bool isGhost(const T i, const T j, const T k);
	bool isGhost(const T i, const T j) {
		return isGhost(i, j, static_cast<T>(0));
	}
	bool isGhost(const T i) {
		return isGhost(i, static_cast<T>(0));
	}
	bool isOnFace(const CartDir d, const CartSide s);

	CartDim dim() const;

	/* get number of components */
	T getNC() const {
		return m_nc;
	};
	
	bool isEqual(const PDimRaw<T>& rhs) const {
		return isEqualNoGhost(rhs) && (m_ghost_size == rhs.m_ghost_size);
	}
	
	// is equal without ghost check
	bool isEqualNoGhost(const PDimRaw<T>& rhs) const {
		const PDimRaw<T>& lhs = *this;
		const bool is1 = lhs.m_global_size == rhs.m_global_size;
		const bool is2 = lhs.m_local_size == rhs.m_local_size;
		const bool is3 = lhs.m_origin == rhs.m_origin;
		const bool is4 = lhs.m_nc == rhs.m_nc;
		return is1 && is2 && is3 && is4;
	}

	Dim3D<T> m_global_size;
	T m_global_size_all;
	Dim3D<T> m_global_stride;

	Dim3D<T> m_local_size;
	T m_local_size_all;
	Dim3D<T> m_local_stride;

	Dim3D<T> m_ghost_size;
	Dim3D<T> m_local_ghost_size; /* m_local_ghost_size = m_local_size + 2 * m_ghost_size */
	T m_local_ghost_size_all;
	Dim3D<T> m_origin;

	T m_nc; /* number of components */
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

template <typename T>
bool operator==(const PDimRaw<T> &lhs, const PDimRaw<T> &rhs) {
	return lhs.isEqual(rhs);
}

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
bool PDimRaw<T>::check(const T index, const CartDir dir) {
	return index >= m_origin[dir] && index < m_origin[dir] + m_local_size[dir];
}

template <typename T>
bool PDimRaw<T>::checkg(const T index, const CartDir dir) {
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
