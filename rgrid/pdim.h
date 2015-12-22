/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_PDIM_H
#define RGRID_PDIM_H

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

	void ortDirs(const CartDir &axis, CartDir &dir1, CartDir &dir2)
	{
		if (axis == X) {
			dir1 = Y;
			dir2 = Z;
		} else if (axis == Y) {
			dir1 = Z;
			dir2 = X;
		} else if (axis == Z) {
			dir1 = X;
			dir2 = Y;
		} else {
			dir1 = dir2 = DIR_UNDEFINED;
		}
	}


template <typename T>
class PDim {
	public:
		PDim() { resize(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); }
		PDim(const PDim &p) { resize(p); }
		void resize(const PDim &p);
		void resize(const T &x, const T &y, const T&z,
		            const T &px, const T &py, const T &pz,
		            const T &ox, const T &oy, const T &oz,
		            const T &gx, const T &gy, const T &gz);
		void resize(const T &x, const T &y, const T &z);
		void resize(const T &x, const T &y, const T &z, const T &gx, const T &gy, const T &gz);
		void resize(const T &x, const T &y) { resize(x, y, static_cast<T>(1)); }
		void resize(const T &x) { resize(x, static_cast<T>(1)); }
		void resize_d(const CartDir &d1, const T &s1, const CartDir &d2, const T &s2, const CartDir &d3, const T &s3);

		T ind(const T &x, const T &y, const T &z);
		T ind(const T &x, const T &y) { return ind(x, y, static_cast<T>(0)); }
		T ind(const T &x) { return ind(x, static_cast<T>(0)); }

		T ind(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3);
		T ind(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2) { return ind(i1, d1, i2, d2, static_cast<T>(0), Z); }
		T ind(const T &i1, const CartDir &d1) { return ind(i1, d1, static_cast<T>(0), Y); }

		T ghost(const CartDir &d) const { return ghost_size[d]; }
		T ghost() const { return ghost(X); }

		T size(const CartDir &d) const { return all_ext[d]; }
		T size() const { return allocalSize; }

		T localSize(const CartDir &d) const { return local_ext[d]; }
		T localSize() const { return local_size; }

		T origin(const CartDir &d) const { return o[d]; }
		T origin() const { return origin(X); }

		T localSizeGhost(const CartDir &d) const { return local_ghost_ext[d]; }
		T localSizeGhost() const { return local_ghost_size; }
		T stride(const CartDir &d) const { return local_stride[d]; }

		T localToIJK(const T &index, const CartDir &dir);
		T indNoGhost(const T &x, const T &y, const T &z);

		T localToGlobal(const T &index, const CartDir &dir);
		T localFromGlobal(const T &x, const T &y, const T &z);
		T localFromGlobal(const T &x, const T &y) { return localFromGlobal(x, y, static_cast<T>(0)); }
		T localFromGlobal(const T &x) { return localFromGlobal(x, static_cast<T>(0)); }

		T indGlobal(const T &x, const T &y, const T &z) { return localFromGlobal(x, y, z); }
		T indGlobal(const T &x, const T &y) { return localFromGlobal(x, y); }
		T indGlobal(const T &x) { return localFromGlobal(x); }

		bool check(const T &index, const CartDir &dir);
		bool check(const T &i) { return check(i, X); }
		bool check(const T &i, const T &j) { return check(i, X) && check(j, Y); }
		bool check(const T &i, const T &j, const T &k) { return check(i, X) && check(j, Y) && check(k, Z); }
		bool checkg(const T &index, const CartDir &dir);
		bool checkg(const T &i) { return checkg(i, X); }
		bool checkg(const T &i, const T &j) { return checkg(i, X) && checkg(j, Y); }
		bool checkg(const T &i, const T &j, const T &k) { return checkg(i, X) && checkg(j, Y) && checkg(k, Z); }

		bool isCorner(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3);
		bool isEdge(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3);
		bool isGhost(const T &i, const T &j, const T &k);
		bool isGhost(const T &i, const T &j) { return isGhost(i, j, static_cast<T>(0)); }
		bool isGhost(const T &i) { return isGhost(i, static_cast<T>(0)); }
		bool isOnFace(const CartDir &d, const CartSide &s);

		bool operator == (PDim<T> a);
		CartDim dim();

		const static int GHOST_SIZE = 2;
		int all_ext[ALL_DIRS];
		int allocalSize;
		int all_stride[ALL_DIRS];

		int local_ext[ALL_DIRS];
		int local_size;
		int local_stride[ALL_DIRS];
	
		int ghost_size[ALL_DIRS];
		int local_ghost_ext[ALL_DIRS];
		int local_ghost_size;

		int o[ALL_DIRS];
	}; // PDim

template <typename T>
CartDim PDim<T>::dim()
{
	CartDim d = DIM_3D;
	if (size(Z) == 1) d = DIM_2D;
	if (size(Y) == 1) d = DIM_1D;
	return d;
}

template <typename T>
bool PDim<T>::isOnFace(const CartDir &d, const CartSide &s)
{
	bool ret = false;
	if (s == SIDE_LEFT && origin(d) == 0) ret = true;
	if (s == SIDE_RIGHT && origin(d) + localSize(d) == size(d)) ret = true;
	return ret;
}

template <typename T>
bool PDim<T>::isGhost(const T &i, const T &j, const T &k)
{
	return i < 0 || i >= localSize(X) || j < 0 || j >= localSize(Y) || k < 0 || k >= localSize(Z);
}

template <typename T>
bool PDim<T>::operator == (PDim<T> a)
{
	bool is1 = (all_ext[X] == a.all_ext[X]) && (all_ext[Y] == a.all_ext[Y]) && (all_ext[Z] == a.all_ext[Z]);
	bool is2 = (local_ext[X] == a.local_ext[X]) && (local_ext[Y] == a.local_ext[Y]) && (local_ext[Z] == a.local_ext[Z]);
	bool is3 = (ghost_size[X] == a.ghost_size[X]) && (ghost_size[Y] == a.ghost_size[Y]) && (ghost_size[Z] == a.ghost_size[Z]);
	bool is4 = (o[X] == a.o[X]) && (o[Y] == a.o[Y]) && (o[Z] == a.o[Z]);
	return is1 && is2 && is3 && is4;
}

template <typename T>
bool PDim<T>::isCorner(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3)
{
	bool a = (localToGlobal(i1, d1) == static_cast<T>(0)) || (localToGlobal(i1, d1) == all_ext[d1] - 1);
	bool b = (localToGlobal(i2, d2) == static_cast<T>(0)) || (localToGlobal(i2, d2) == all_ext[d2] - 1);
	bool c = (localToGlobal(i3, d3) == static_cast<T>(0)) || (localToGlobal(i3, d3) == all_ext[d3] - 1);
	return a && b && c;
}

template <typename T>
bool PDim<T>::isEdge(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3)
{
	bool a = (localToGlobal(i1, d1) == static_cast<T>(0)) || (localToGlobal(i1, d1) == all_ext[d1] - 1);
	bool b = (localToGlobal(i2, d2) == static_cast<T>(0)) || (localToGlobal(i2, d2) == all_ext[d2] - 1);
	bool c = (localToGlobal(i3, d3) == static_cast<T>(0)) || (localToGlobal(i3, d3) == all_ext[d3] - 1);
	return (a && b) || (b && c) || (c && a);
}

template <typename T>
bool PDim<T>::check(const T& index, const CartDir& dir)
{
	return index >= o[dir] && index < o[dir] + local_ext[dir];
}

template <typename T>
bool PDim<T>::checkg(const T& index, const CartDir& dir)
{
	return index >= -ghost_size[dir] + o[dir] && index < o[dir] + local_ext[dir] + ghost_size[dir];
}


template <typename T>
T PDim<T>::localToGlobal(const T& index, const CartDir& dir)
{
	return index + o[dir];
}

template <typename T>
T PDim<T>::indNoGhost(const T& x, const T& y, const T& z)
{
	return x + y * local_ext[X] + z * local_ext[X] * local_ext[Y];
}

template <typename T>
T PDim<T>::localFromGlobal(const T& x, const T& y, const T& z)
{
	return ind(x - o[X], y - o[Y], z - o[Z]);
}


template <typename T>
T PDim<T>::ind(const T& x, const T& y, const T& z)
{
	return x + ghost_size[X] + (y + ghost_size[Y]) * local_stride[Y] + (z + ghost_size[Z]) * local_stride[Z];
}

template <typename T>
T PDim<T>::ind(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3)
{
	return (i1 + ghost_size[d1]) * local_stride[d1] + (i2 + ghost_size[d2]) * local_stride[d2] + (i3 + ghost_size[d3]) * local_stride[d3];
}

template <typename T>
void PDim<T>::resize_d(const CartDir &d1, const T &s1, const CartDir &d2, const T &s2, const CartDir &d3, const T &s3)
{
	all_ext[d1] = s1;
	all_ext[d2] = s2;
	all_ext[d3] = s3;
	resize(all_ext[X], all_ext[Y], all_ext[Z]);
}

template <typename T>
void PDim<T>::resize(const T &x, const T &y, const T &z, const T &gx, const T &gy, const T &gz)
{
	resize(x, y, z,
	       x, y, z,
	       0, 0, 0,
	       gx, gy, gz);
}

template <typename T>
void PDim<T>::resize(const T &x, const T &y, const T &z)
{
	resize(x, y, z,
	       x > 1 ? GHOST_SIZE : 0, y > 1 ? GHOST_SIZE : 0, z > 1 ? GHOST_SIZE : 0);
}

template <typename T>
void PDim<T>::resize(const PDim<T> &p)
{
	resize(p.all_ext[X], p.all_ext[Y], p.all_ext[Z],
	       p.local_ext[X], p.local_ext[Y], p.local_ext[Z],
	       p.o[X], p.o[Y], p.o[Z],
	       p.ghost_size[X], p.ghost_size[Y], p.ghost_size[Z]
	);
}


template <typename T>
void PDim<T>::resize(const T& x, const T& y, const T& z, const T& px, const T& py, const T& pz, const T& ox, const T& oy, const T& oz, const T& gx, const T& gy, const T& gz)
{
	// Set base variables.
	all_ext[X] = x;
	all_ext[Y] = y;
	all_ext[Z] = z;

	local_ext[X] = px;
	local_ext[Y] = py;
	local_ext[Z] = pz;

	o[X] = ox;
	o[Y] = oy;
	o[Z] = oz;

	ghost_size[X] = gx;
	ghost_size[Y] = gy;
	ghost_size[Z] = gz;

	// Calculate auxilary variables.
	allocalSize = all_ext[X] * all_ext[Y] * all_ext[Z];

	all_stride[X] = 1;
	all_stride[Y] = all_ext[X];
	all_stride[Z] = all_ext[X] * all_ext[Y];

	local_size = local_ext[X] * local_ext[Y] * local_ext[Z];

	local_ghost_ext[X] = local_ext[X] + ghost_size[X] * 2;
	local_ghost_ext[Y] = local_ext[Y] + ghost_size[Y] * 2;
	local_ghost_ext[Z] = local_ext[Z] + ghost_size[Z] * 2;

	local_ghost_size = local_ghost_ext[X] * local_ghost_ext[Y] * local_ghost_ext[Z];

	local_stride[X] = 1;
	local_stride[Y] = local_ghost_ext[X];
	local_stride[Z] = local_ghost_ext[X] * local_ghost_ext[Y];
}

template <typename T>
T PDim<T>::localToIJK(const T& index, const CartDir& dir)
{
	T k = index / local_stride[Z] - ghost(Z);
	T j = (index - (k + ghost(Z)) * local_stride[Z]) / local_stride[Y] - ghost(Y);
	T i = (index - (k + ghost(Z)) * local_stride[Z] - (j + ghost(Y)) * local_stride[Y]) / local_stride[X] - ghost(X);
	if (dir == Z) return k;
	else if (dir == Y) return j;
	else if (dir == X) return i;
	return -1;
}

}; // rgrid

#endif // RGRID_PDIM_H
