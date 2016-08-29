/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_TYPES_H
#define RGRID_TYPES_H

namespace rgrid {
enum CartDir {
	X = 0,
	Y,
	Z,
	ALL_DIRS,
	DIR_UNDEFINED
};

enum CartSide {
	SIDE_LEFT = 0,
	SIDE_RIGHT,
	SIDE_ALL,
	SIDE_UNDEFINED
};

enum CartDim {
	DIM_1D = 0,
	DIM_2D,
	DIM_3D,
	DIM_ALL
};

template <typename T>
struct Dim3D {
	T x, y, z;
	Dim3D(T i /*= T()*/, T j /*= T()*/, T k /*= T()*/)
		: x(i), y(j), z(k) {
	}
	Dim3D() {
		
	}
	T &operator[](const CartDir d) {
		return (d == X) ? x : ((d == Y) ? y : z);
	}
	T operator[](const CartDir d) const {
		return (d == X) ? x : ((d == Y) ? y : z);
	}
};

template <typename T>
bool operator==(const Dim3D<T>& lhs, const Dim3D<T>& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

template <typename T>
bool operator!=(const Dim3D<T>& lhs, const Dim3D<T>& rhs) {
	return !(lhs == rhs);
}

namespace rgio {

enum format { BINARY = 0, 
              TEXT = 1 };

} // namespace rgio

} // rgrid

#endif // RGRID_TYPES_H

