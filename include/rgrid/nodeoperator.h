#ifndef NODE_OPERATOR_H
#define NODE_OPERATOR_H

#include <iostream>

#include "rgrid/darrayscatter.h"
#include "rgrid/fdc.h"
#include "rgrid/debug.h"

namespace rgrid {

namespace operators {

/*
 * Leaf operators
 */
	
template <typename T, typename I>
struct ConstOp {
	typedef T DataType;
	typedef I IndexType;
	ConstOp(const T& val) : m_val(val) {}
	T eval(I, I, I) const {
		return m_val;
	}
	void attachInput(const std::string&, DArray<T, I>&) {
	}
private:
	const T m_val;
};

template <typename T, typename I>
struct VarOp {
	typedef T DataType;
	typedef I IndexType;
	VarOp(const std::string& var_name) : m_var_name(var_name), m_da(NULL) {}
	void attachInput(const std::string& var_name, DArray<T, I>& da) {
		if (m_var_name == var_name) {
			m_da = &da;
		}
	}
	T eval(I i, I j, I k) const {
		DEBUG_ASSERT(m_da != NULL, "VarOp is not initialized");
		return m_da->val(i, j, k, 0);
	}
private:
	
	const std::string m_var_name;
	DArray<T, I>* m_da;
};

template <
	CartDir D, typename E>
struct FD1Op {
public:
	typedef typename E::DataType DataType;
	typedef typename E::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
	
public:
	FD1Op(DataType space_step, IndexType half_order, const E& e)
	: m_space_step(space_step), m_half_order(half_order), m_e(e)
	{
		STATIC_ASSERT(D == X || D == Y || D == Z, Wrong_CartDir_param);
		DEBUG_ASSERT(half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
	}
	void attachInput(const std::string& var_name, DArray<T, I>& da) {
		m_e.attachInput(var_name, da);
	}
	T eval(I i, I j, I k) const {
		T result = 0;
		for (int n = 0; n < m_half_order; ++n) {
			if (D == X) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i+n+1, j, k) - m_e.eval(i-n, j, k));
			} else if (D == Y) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i, j+n+1, k) - m_e.eval(i, j-n, k));
			} else if (D == Z) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i, j, k+n+1) - m_e.eval(i, j, k-n));
			}
		}
		result /= m_space_step;
		return result;
	}
	
private:
	
	const T m_space_step;
	const I m_half_order;
	E m_e;
};

template <
	CartDir D, typename E>
struct BD1Op {
public:
	typedef typename E::DataType DataType;
	typedef typename E::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
public:
	BD1Op(T space_step, I half_order, const E& e)
	: m_space_step(space_step), m_half_order(half_order), m_e(e) 
	{
		STATIC_ASSERT(D == X || D == Y || D == Z, Wrong_CartDir_param);
		DEBUG_ASSERT(half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
	}
	void attachInput(const std::string& var_name, DArray<T, I>& da) {
		m_e.attachInput(var_name, da);
	}
	T eval(I i, I j, I k) const {
		T result = 0;
		for (int n = 0; n < m_half_order; ++n) {
			if (D == X) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i+n, j, k) - m_e.eval(i-n-1, j, k));
			} else if (D == Y) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i, j+n, k) - m_e.eval(i, j-n-1, k));
			} else if (D == Z) {
				result += fdc::SC1[m_half_order-1][n] *
					(m_e.eval(i, j, k+n) - m_e.eval(i, j, k-n-1));
			}
		}
		result /= m_space_step;
		return result;
	}
private:
	
	const T m_space_step;
	const I m_half_order;
	E m_e;
};

/*
 * Derived operators
 */

template <typename E1, typename E2, template <typename> class BinOp>
struct BinaryOp {
public:
	typedef typename E1::DataType DataType;
	typedef typename E1::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
public:
	BinaryOp(const E1& e1, const E2& e2) 
	: m_e1(e1), m_e2(e2)
	{
		STATIC_ASSERT((IsEqualTypes<typename E1::DataType, typename E2::DataType>::result), Different_data_types);
		STATIC_ASSERT((IsEqualTypes<typename E1::IndexType, typename E2::IndexType>::result), Different_index_types);
	}
	void attachInput(const std::string& var_name, DArray<T, I>& da) {
		m_e1.attachInput(var_name, da);
		m_e2.attachInput(var_name, da);
	}
	T eval(I i, I j, I k) const {
		T result = BinOp<T>::f(m_e1.eval(i, j, k), m_e2.eval(i, j, k));
		return result;
	}
private:
	
	E1 m_e1;
	E2 m_e2;
};

/*
 * Subtypes for BinaryOp
 */
template <typename T>
struct AddOp {
	static T f(const T& e1, const T& e2) { return e1 + e2; }
};

template <typename T>
struct SubOp {
	static T f(const T& e1, const T& e2) { return e1 - e2; }
};

template <typename T>
struct MulOp {
	static T f(const T& e1, const T& e2) { return e1 * e2; }
};

template <typename T>
struct DivOp {
	static T f(const T& e1, const T& e2) { return e1 / e2; }
};



template <typename E, template <typename> class UnOp>
struct UnaryOp {
public:
	typedef typename E::DataType DataType;
	typedef typename E::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
public:
	UnaryOp(const E& e) 
	: m_e(e) {}
	void attachInput(const std::string& var_name, DArray<T, I>& da) {
		m_e.attachInput(var_name, da);
	}
	T eval(I i, I j, I k) const {
		return UnOp<T>::f(m_e.eval(i, j, k));
	}
private:
	
	E m_e;
};

/*
 * Subtypes for UnaryOp
 */
template <typename T>
struct NegOp {
	NegOp() {}
	T operator()(const T& e) { return -e; }
};

/*
 * types creators
 */

template <typename E1, typename E2>
inline BinaryOp<E1, E2, MulOp> operator*(const E1& e1, const E2& e2) {
	return BinaryOp<E1, E2, MulOp>(e1, e2);
}

template <typename E1, typename E2>
inline BinaryOp<E1, E2, DivOp> operator/(const E1& e1, const E2& e2) {
	return BinaryOp<E1, E2, DivOp>(e1, e2);
}

template <typename E1, typename E2>
inline BinaryOp<E1, E2, AddOp> operator+(const E1& e1, const E2& e2) {
	return BinaryOp<E1, E2, AddOp>(e1, e2);
}

template <typename E1, typename E2>
inline BinaryOp<E1, E2, SubOp> operator-(const E1& e1, const E2& e2) {
	return BinaryOp<E1, E2, SubOp>(e1, e2);
}

template <CartDir D, typename E>
inline FD1Op<D, E> fd1(
	typename E::DataType space_step,
	typename E::IndexType half_order,
	const E& expr)
{
	return FD1Op<D, E>(space_step, half_order, expr);
}

template <CartDir D, typename E>
inline BD1Op<D, E> bd1(
	typename E::DataType space_step, 
	typename E::IndexType half_order, 
	const E& expr)
{
	return BD1Op<D, E>(space_step, half_order, expr);
}

} // namespace operators

} // namespace rgrid

#endif

