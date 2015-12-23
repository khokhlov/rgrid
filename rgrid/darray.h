/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_DARRAY_H
#define RGRID_DARRAY_H

#include <vector>
#include <iostream>

#include "rgrid/pdim.h"
#include "rgrid/utils.h"

namespace rgrid {

template <typename T, typename I>
class DArray : public PDim<I> {
	public:
		enum StorageType {
			STORAGE_SOA = 0,
			STORAGE_AOS
		};
		DArray() : PDim<I>() {
			st = STORAGE_SOA;
			nc = 1;
		}
		void alloc() { alloc(static_cast<I>(1)); }
		void alloc(const I &nc);

		void fillGhost(const CartDir &d, const CartSide &s);
		void fillGhost(const CartDir &d) { fillGhost(d, SIDE_LEFT); fillGhost(d, SIDE_RIGHT); }
		void fillGhost() { fillGhost(X); fillGhost(Y); fillGhost(Z); }

		void fill(const T &v);

		T &operator[](const I index) { return data[index]; }
		const T &operator[](const I index) const { data[index]; }

		T &val(const I i, const I j, const I k) { return data[this->ind(i, j, k)]; }
		T &val(const I i, const I j) { return val(i, j, static_cast<T>(0)); }
		T &val(const I i) { return val(i, static_cast<T>(0)); }
		T &val(const T &i1, const CartDir &d1, const T &i2, const CartDir &d2, const T &i3, const CartDir &d3) { return data[this->ind(i1, d1, i2, d2, i3, d3)]; }

		T &operator()(const I i, const I j, const I k) { return val(i, j, k); }
		T &operator()(const I i, const I j) { return (*this)(i, j, static_cast<I>(0)); }
		T &operator()(const I i) { return (*this)(i, static_cast<I>(0)); }

		const T &operator()(const I i, const I j, const I k) const { return data[this->ind(i, j, k)]; }
		const T &operator()(const I i, const I j) const { return (*this)(i, j, static_cast<I>(0)); }
		const T &operator()(const I i) const { return (*this)(i, static_cast<I>(0)); }
		
		/*void convert_to_storage(const StorageType &st);*/

		/* Data. */
		std::vector<T> data;

		/* Storage type. */
		StorageType st;

		/* Number of components. */
		I nc;
}; // DArray

template <typename T, typename I>
void DArray<T, I>::alloc(const I &nc)
{
	this->nc = nc;
	data.resize(this->localSizeGhost() * nc);
}

template <typename T, typename I>
void DArray<T, I>::fill(const T &v)
{
	for (I i = 0; i < this->localSizeGhost(); i++) { data[i] = v; }
}

template <typename T, typename I>
void DArray<T, I>::fillGhost(const CartDir &d, const CartSide &s)
{
	CartDir ort1, ort2;
	ortDirs(d, ort1, ort2);
	I k = 0;
	I sign;
	if (this->isOnFace(d, s)) {
		if (s == SIDE_LEFT) {
			k = 0;
			sign = 1;
		} else {
			k = this->localSize(d) - 1;
			sign = -1;
		}
	} else {
		return;
	}
	for (I i = 0; i < this->localSize(ort1); i++) {
		for (I j = 0; j < this->localSize(ort2); j++) {
			for (I gs = 1; gs <= this->ghost(d); gs++) {
				this->val(i, ort1, j, ort2, k - sign * gs, d) = this->val(i, ort1, j, ort2, k, d);
			}
		}
	}
}
}; // rgrid

#endif // RGRID_DARRAY_H
