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
} // rgrid

#endif // RGRID_TYPES_H

