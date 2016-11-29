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

inline const CartDir& operator++(CartDir& d) {
	d = static_cast<CartDir>(d + 1);
	return d;
}

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
	Dim3D(const T a[ALL_DIRS])
		: x(a[X]), y(a[Y]), z(a[Z])
	{
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

	template <typename U>
	operator Dim3D<U>() const {
		Dim3D<U> t;
		t.x = static_cast<U>(x);
		t.y = static_cast<U>(y);
		t.z = static_cast<U>(z);
		return t;
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
 *
 * Header of DARRAY data format described here
 * \code{.unparsed}
 * # DARRAY DATA FORMAT
 * SIZE: X Y Z
 * COMPONENTS: C
 * FORMAT: format
 * DATA START
 * \endcode
 * 'SIZE' is 3 numbers - 'X Y Z'. 'COMPONENTS' is the number of values inside node. 'FORMAT' can be 'binary' or 'text'.
 * After 'DATA START' placed all the data in binary format or in text format(values divided by spaces).
 * Number of values is X*Y*Z*C. Values placed in the next order: line of X, plane of Y, volume of Z, and after that in the same order volumes for next components.
 *
 * Example:
 * \code{.unparsed}
 * # DARRAY DATA FORMAT
 * SIZE: 5 10 1
 * COMPONENTS: 3
 * FORMAT: binary
 * DATA START
 * \endcode
 */
enum format {
	BINARY, ///< own DARRAY binary format
	TEXT, ///< own DARRAY text format
	CUSTOM_HEADER ///< binary format with custom header
};

} // namespace rgio

} // rgrid

#endif // RGRID_TYPES_H

