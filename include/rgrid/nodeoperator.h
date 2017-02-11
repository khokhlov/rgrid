#ifndef NODE_OPERATOR_H
#define NODE_OPERATOR_H

#include <iostream>

#include "rgrid/darrayscatter.h"
#include "rgrid/fdc.h"
#include "rgrid/debug.h"

namespace rgrid {

namespace operators {

struct B {
	double val;
};
	
template <int N>
struct A {
	B val;
	A<N-1> prev;
	template <int Num>
	B& getStruct() {
		return Num == N ? val : prev.getStruct<Num>();
	}
	enum { number = N };
};

template <>
struct A<1> {
	B val;
	template <int Num>
	B& getStruct() {
		STATIC_ASSERT(1 <= Num, Wrong_index);
		return val;
	}
};

void foo() {
	A<5> a5;
	(void) a5.getStruct<1>();
}

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
private:
	ConstOp(const ConstOp&);
	ConstOp& operator=(const ConstOp&);
	const T m_val;
};

template <typename T, typename I>
struct VarOp {
	typedef T DataType;
	typedef I IndexType;
	VarOp() : m_da(NULL) {}
	void attachInput(DArray<T, I>& da) {
		m_da = &da;
	}
	T eval(I i, I j, I k) const {
		DEBUG_ASSERT(m_da != NULL, "VarOp is not initialized");
		return m_da->val(i, j, k, 0);
	}
private:
	VarOp(const VarOp&);
	VarOp& operator=(const VarOp&);
	
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
	
	FD1Op(const FD1Op&);
	FD1Op& operator=(const FD1Op&);
	
	const T m_space_step;
	const I m_half_order;
	const E& m_e;
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
	BD1Op(const BD1Op&);
	BD1Op& operator=(const BD1Op&);
	
	const T m_space_step;
	const I m_half_order;
	const E& m_e;
};

/*
 * Derived operators
 */

template <typename E1, typename E2, typename BinOp>
struct BinaryOp {
public:
	typedef typename E1::DataType DataType;
	typedef typename E1::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
public:
	BinaryOp(const E1& e1, const E2& e2, const BinOp& op = BinOp()) 
	: m_e1(e1), m_e2(e2), m_op(op) 
	{
		STATIC_ASSERT((IsEqualTypes<typename E1::DataType, typename E2::DataType>::result), Different_data_types);
		STATIC_ASSERT((IsEqualTypes<typename E1::IndexType, typename E2::IndexType>::result), Different_index_types);
	}
	T eval(I i, I j, I k) const {
		T result = m_op(m_e1.eval(i, j, k), m_e2.eval(i, j, k));
		return result;
	}
private:
	BinaryOp(const BinaryOp&);
	BinaryOp& operator=(const BinaryOp&);
	
	const E1& m_e1;
	const E2& m_e2;
	const BinOp& m_op;
};

/*
 * Subtypes for BinaryOp
 */
template <typename T>
struct AddOp {
	AddOp() {}
	T operator()(const T& e1, const T& e2) const { return e1 + e2; }
};

template <typename T>
struct SubOp {
	SubOp() {}
	T operator()(const T& e1, const T& e2) const { return e1 - e2; }
};

template <typename T>
struct MulOp {
	MulOp() {}
	T operator()(const T& e1, const T& e2) const { return e1 * e2; }
};



template <typename E, typename UnOp>
struct UnaryOp {
public:
	typedef typename E::DataType DataType;
	typedef typename E::IndexType IndexType;
private:
	typedef DataType T;
	typedef IndexType I;
public:
	UnaryOp(const E& e, const UnOp& op = UnOp()) 
	: m_e(e), m_op(op) {}
	T eval(I i, I j, I k) const {
		return m_op(m_e.eval(i, j, k));
	}
private:
	UnaryOp(const UnaryOp&);
	UnaryOp& operator=(const UnaryOp&);
	
	const E& m_e;
	const UnOp& m_op;
};

/*
 * Subtypes for UnaryOp
 */
template <typename T>
struct NegOp {
	NegOp() {}
	T operator()(const T& e) { return -e; }
};

} // namespace operators

} // namespace rgrid

#endif

