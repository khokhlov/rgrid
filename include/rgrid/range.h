#ifndef RGRID_RANGE_H
#define RGRID_RANGE_H

#include "rgrid/types.h"

namespace rgrid {

template <typename I>
struct Range {
	static const int L = SIDE_LEFT;
	static const int R = SIDE_RIGHT;
	
	Range() : m_r({Dim3D<I>(0,0,0), Dim3D<I>(0,0,0)}) {}
	Range(Dim3D<I> l, Dim3D<I> r) : m_r({l, r}) {
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
	
	Dim3D<I> get(CartSide s) const { m_r[s]; }
	
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
	r.set(max(r1.getLeft(X), r2.getLeft(X)), X, Range<I>::L);
	r.set(max(r1.getRight(X), r2.getRight(X)), X, Range<I>::R);
	r.set(max(r1.getLeft(Y), r2.getLeft(Y)), Y, Range<I>::L);
	r.set(min(r1.getRight(Y), r2.getRight(Y)), Y, Range<I>::R);
	r.set(min(r1.getLeft(Z), r2.getLeft(Z)), Z, Range<I>::L);
	r.set(min(r1.getRight(Z), r2.getRight(Z)), Z, Range<I>::R);
	if (!r.isValid()) r.clear();
	return r;
}

template <typename I>
inline bool isIntersect(const Range<I>& r1, const Range<I>& r2) {
	return !intersect(r1, r2).isEmpty();
}

} // namespace rgrid

#endif
