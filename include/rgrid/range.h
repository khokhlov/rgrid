#ifndef RGRID_RANGE_H
#define RGRID_RANGE_H

#include <algorithm>

#include "rgrid/types.h"

namespace rgrid {

template <typename I>
struct Range {
	static const int L = SIDE_LEFT;
	static const int R = SIDE_RIGHT;
	
	Range() : m_r({Dim3D<I>(0,0,0), Dim3D<I>(0,0,0)}) {}
	Range(Dim3D<I> l, Dim3D<I> r) : m_r({l, r}) {}
	Range(
		I li, I lj, I lk,
		I ri, I rj, I rk) : Range(Dim3D<I>(li, lj, lk), Dim3D<I>(ri, rj, rk))
	{
	}
	
	void set(I v, CartDir d, CartSide s) { 
		m_r[s][d] = v;
	}
	
	void set(const Dim3D<I>& v, CartSide s) { 
		m_r[s] = v;
	}
	
	I get(CartDir d, CartSide s) const { 
		return m_r[s][d]; 
	}
	
	Dim3D<I> get(CartSide s) const { return m_r[s]; }
	
	bool isValid() { 
		return 
			!( m_r[L].x > m_r[R].x
			|| m_r[L].y > m_r[R].y
			|| m_r[L].z > m_r[R].z);
	}
	
	bool isEmpty() { 
		return 
			m_r[L][X] == m_r[R][X] ||
			m_r[L][Y] == m_r[R][Y] ||
			m_r[L][Z] == m_r[R][Z];
	}
	
	bool isInside(I i, I j, I k) const {
		return 
			m_r[L][X] <= i && i < m_r[R][X] &&
			m_r[L][Y] <= j && j < m_r[R][Y] &&
			m_r[L][Z] <= k && k < m_r[R][Z];
	}
	
	bool isInside(const Dim3D<I>& v) const {
		return isInside(v.x, v.y, v.z);
	}
	
	void clear() {
		for (CartSide s = SIDE_LEFT; s != SIDE_ALL; ++s)
			for (CartDir d = X; d != ALL_DIRS; ++d)
				m_r[s][d] = 0;
	}
	
	Range<I> shift(const Dim3D<I>& v) const {
		return Range<I>(m_r[L] + v, m_r[R] + v);
	}
	
private:
	
	Dim3D<I> m_r[SIDE_ALL];
};

template <typename I>
inline Range<I> intersect(const Range<I>& r1, const Range<I>& r2) {
	Range<I> r;
	r.set(std::max(r1.get(X, SIDE_LEFT), r2.get(X, SIDE_LEFT)), X, SIDE_LEFT);
	r.set(std::min(r1.get(X, SIDE_RIGHT), r2.get(X, SIDE_RIGHT)), X, SIDE_RIGHT);
	r.set(std::max(r1.get(Y, SIDE_LEFT), r2.get(Y, SIDE_LEFT)), Y, SIDE_LEFT);
	r.set(std::min(r1.get(Y, SIDE_RIGHT), r2.get(Y, SIDE_RIGHT)), Y, SIDE_RIGHT);
	r.set(std::max(r1.get(Z, SIDE_LEFT), r2.get(Z, SIDE_LEFT)), Z, SIDE_LEFT);
	r.set(std::min(r1.get(Z, SIDE_RIGHT), r2.get(Z, SIDE_RIGHT)), Z, SIDE_RIGHT);
	if (!r.isValid()) r.clear();
	return r;
}

template <typename I>
inline bool isIntersect(const Range<I>& r1, const Range<I>& r2) {
	return !intersect(r1, r2).isEmpty();
}

} // namespace rgrid

#endif
