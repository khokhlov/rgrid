/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_PDIM_H
#define RGRID_PDIM_H

#include "rgrid/types.h"
#include "rgrid/utils.h"
#include "rgrid/pdimraw.h"

namespace rgrid {

template <typename T>
struct PDim : public PDimRaw<T> {
	virtual void resize(const T x, const T y, const T z,
	                    const T px, const T py, const T pz,
	                    const T ox, const T oy, const T oz,
	                    const T gx, const T gy, const T gz,
	                    const T cn);
	virtual void resize(const PDim<T> &p);
	virtual void resize(const T x, const T y, const T z, const T cn);
	virtual void resize(const T x, const T y, const T z,
	                    const T gx, const T gy, const T gz,
	                    const T cn);
	virtual void resize_d(const CartDir d1, const T s1, const CartDir d2, const T s2, const CartDir d3, const T s3, const T cn);

	const static T s_ghost_size = 2;

	PDim(): PDimRaw<T>() {
	}
	
	PDim(const T x, const T y = 1, const T z = 1, const T cn = 1) : PDimRaw<T>() {
		resize(x, y, z, x, y, z, 0, 0, 0, 0, 0, 0, cn);
	}
	PDim(const T x, const T y, const T z,
	     const T px, const T py, const T pz,
	     const T ox, const T oy, const T oz,
	     const T gx, const T gy, const T gz,
	     const T cn) : PDimRaw<T>() {
		resize(x, y, z, px, py, pz, ox, oy, oz, gx, gy, gz, cn);
	}
	PDim(const PDim<T> &rhs) : PDimRaw<T>()  {
		resize(rhs);
	}
	PDim<T> &operator=(const PDim<T> &rhs) {
		resize(rhs);
		return *this;
	}
	
	virtual ~PDim() {}
}; // PDim

template <typename T>
void PDim<T>::resize_d(const CartDir d1, const T s1, const CartDir d2, const T s2, const CartDir d3, const T s3, const T cn) {
	PDimRaw<T>::m_global_size[d1] = s1;
	PDimRaw<T>::m_global_size[d2] = s2;
	PDimRaw<T>::m_global_size[d3] = s3;
	resize(PDimRaw<T>::m_global_size[X], PDimRaw<T>::m_global_size[Y], PDimRaw<T>::m_global_size[Z], cn);
}

template <typename T>
void PDim<T>::resize(const T x, const T y, const T z, const T gx, const T gy, const T gz, const T cn) {
	resize(x, y, z,
	       x, y, z,
	       0, 0, 0,
	       gx, gy, gz,
	       cn);
}

template <typename T>
void PDim<T>::resize(const T x, const T y, const T z, const T cn) {
	resize(x, y, z,
	       x > 1 ? s_ghost_size : 0,
	       y > 1 ? s_ghost_size : 0,
	       z > 1 ? s_ghost_size : 0,
	       cn);
}

template <typename T>
void PDim<T>::resize(const PDim<T> &p) {
	resize(p.m_global_size[X], p.m_global_size[Y], p.m_global_size[Z],
	       p.m_local_size[X], p.m_local_size[Y], p.m_local_size[Z],
	       p.m_origin[X], p.m_origin[Y], p.m_origin[Z],
	       p.m_ghost_size[X], p.m_ghost_size[Y], p.m_ghost_size[Z],
	       p.m_nc);
}


template <typename T>
void PDim<T>::resize(const T x, const T y, const T z, const T px, const T py, const T pz, const T ox, const T oy, const T oz, const T gx, const T gy, const T gz, const T cn) {
	// Set base variables.
	PDimRaw<T>::m_global_size[X] = x;
	PDimRaw<T>::m_global_size[Y] = y;
	PDimRaw<T>::m_global_size[Z] = z;

	PDimRaw<T>::m_local_size[X] = px;
	PDimRaw<T>::m_local_size[Y] = py;
	PDimRaw<T>::m_local_size[Z] = pz;

	PDimRaw<T>::m_origin[X] = ox;
	PDimRaw<T>::m_origin[Y] = oy;
	PDimRaw<T>::m_origin[Z] = oz;

	PDimRaw<T>::m_ghost_size[X] = gx;
	PDimRaw<T>::m_ghost_size[Y] = gy;
	PDimRaw<T>::m_ghost_size[Z] = gz;

	PDimRaw<T>::m_nc = cn;

	// Calculate auxilary variables.
	PDimRaw<T>::m_global_size_all = PDimRaw<T>::m_global_size[X] * PDimRaw<T>::m_global_size[Y] * PDimRaw<T>::m_global_size[Z];

	PDimRaw<T>::m_global_stride[X] = 1;
	PDimRaw<T>::m_global_stride[Y] = PDimRaw<T>::m_global_size[X];
	PDimRaw<T>::m_global_stride[Z] = PDimRaw<T>::m_global_size[X] * PDimRaw<T>::m_global_size[Y];

	PDimRaw<T>::m_local_size_all = PDimRaw<T>::m_local_size[X] * PDimRaw<T>::m_local_size[Y] * PDimRaw<T>::m_local_size[Z];

	PDimRaw<T>::m_local_ghost_size[X] = PDimRaw<T>::m_local_size[X] + PDimRaw<T>::m_ghost_size[X] * 2;
	PDimRaw<T>::m_local_ghost_size[Y] = PDimRaw<T>::m_local_size[Y] + PDimRaw<T>::m_ghost_size[Y] * 2;
	PDimRaw<T>::m_local_ghost_size[Z] = PDimRaw<T>::m_local_size[Z] + PDimRaw<T>::m_ghost_size[Z] * 2;

	PDimRaw<T>::m_local_ghost_size_all = PDimRaw<T>::m_local_ghost_size[X] * PDimRaw<T>::m_local_ghost_size[Y] * PDimRaw<T>::m_local_ghost_size[Z];

	PDimRaw<T>::m_local_stride[X] = 1;
	PDimRaw<T>::m_local_stride[Y] = PDimRaw<T>::m_local_ghost_size[X];
	PDimRaw<T>::m_local_stride[Z] = PDimRaw<T>::m_local_ghost_size[X] * PDimRaw<T>::m_local_ghost_size[Y];
}

}; // rgrid

#endif // RGRID_PDIM_H



