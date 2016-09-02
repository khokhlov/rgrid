/**
 * \file
 * \brief Types common for all classes
 */

/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_TYPES_H
#define RGRID_TYPES_H

namespace rgrid {
/**
 * CartDir represents direction in 3D space
 */
enum CartDir {
	X = 0,
	Y,
	Z,
	ALL_DIRS,
	DIR_UNDEFINED
};

/**
 * CartSide represents side of something: left or right
 */
enum CartSide {
	SIDE_LEFT = 0,
	SIDE_RIGHT,
	SIDE_ALL,
	SIDE_UNDEFINED
};

/**
 * CartDim represents number of dimensions
 */
enum CartDim {
	DIM_1D = 0,
	DIM_2D,
	DIM_3D,
	DIM_ALL
};

/**
 * \brief 3D coordinates
 * \tparam T type of coordinates
 */
template <typename T>
struct Dim3D {
	/**@{*/ 
	/** coordinates */
	T x, y, z;
	/**@}*/
	/**
	 * \brief Init by some values
	 */
	Dim3D(T i /*= T()*/, T j /*= T()*/, T k /*= T()*/)
		: x(i), y(j), z(k) {
	}
	Dim3D() {
		
	}
	/**
	 * \brief Get coord by index X, Y or Z
	 */
	T &operator[](const CartDir d) {
		return (d == X) ? x : ((d == Y) ? y : z);
	}
	
	/**
	 * \brief Get coord by index X, Y or Z
	 */
	T operator[](const CartDir d) const {
		return (d == X) ? x : ((d == Y) ? y : z);
	}
};

/// Check is two Dim3D equal
template <typename T>
bool operator==(const Dim3D<T>& lhs, const Dim3D<T>& rhs) {
	return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

/// Check is two Dim3D not equal
template <typename T>
bool operator!=(const Dim3D<T>& lhs, const Dim3D<T>& rhs) {
	return !(lhs == rhs);
}

namespace rgio {

/**
 * \brief Formats in which data can be stored in file
 */
enum format { BINARY = 0, 
              TEXT = 1 };

} // namespace rgio

} // rgrid

#endif // RGRID_TYPES_H

