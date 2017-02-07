#ifndef NODE_OPERATOR_H
#define NODE_OPERATOR_H

#include "rgrid/darrayscatter.h"
#include "rgrid/fdc.h"
#include "rgrid/debug.h"
#include "rgrid/taskparams.h"

namespace rgrid {

namespace operators {

/*
 * Grid operators
 */

/**
 * \brief First derivative staggered forward operator
 */
template <typename T, typename I, CartDir D>
struct FD1 {
	typedef DArray<T, I> DA;
	static inline T apply(
		const TaskParams<T, I>& tp,
		I i, I j, I k,
		const DA& in)
	{
		STATIC_ASSERT(D == X || D == Y || D == Z, BD_Wrong_CartDir_param);
		DEBUG_ASSERT(tp.half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		T result = 0;
		for (int n = 0; n < tp.half_order; ++n) {
			if (D == X) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i+n+1, j, k, 0) - in(i-n, j, k, 0));
			} else if (D == Y) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i, j+n+1, k, 0) - in(i, j-n, k, 0));
			} else if (D == Z) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i, j, k+n+1, 0) - in(i, j, k-n, 0));
			}
		}
		result /= tp.space_step[D];
		return result;
	}
};

/**
 * \brief First derivative staggered backward operator
 */
template <typename T, typename I, CartDir D>
struct BD1 {
	typedef DArray<T, I> DA;
	static inline T apply(
		const TaskParams<T, I>& tp,
		I i, I j, I k,
		const DA& in)
	{
		STATIC_ASSERT(D == X || D == Y || D == Z, BD_Wrong_CartDir_param);
		DEBUG_ASSERT(tp.half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		T result = 0;
		for (int n = 0; n < tp.half_order; ++n) {
			if (D == X) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i+n, j, k, 0) - in(i-n-1, j, k, 0));
			} else if (D == Y) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i, j+n, k, 0) - in(i, j-n-1, k, 0));
			} else if (D == Z) {
				result += fdc::SC1[tp.half_order-1][n] *
					(in(i, j, k+n, 0) - in(i, j, k-n-1, 0));
			}
		}
		result /= tp.space_step[D];
		return result;
	}
};

/**
 * \brief Second derivative operator
 */
template <typename T, typename I, CartDir D>
struct D2 {
	typedef DArray<T, I> DA;
	static inline T apply(
		const TaskParams<T, I>& tp,
		I i, I j, I k,
		const DA& in)
	{
		STATIC_ASSERT(D == X || D == Y && D == Z, D2_Wrong_CartDir_param);
		DEBUG_ASSERT(tp.half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		T result = fdc::C2[tp.half_order-1][0] * in(i, j, k, 0);
		for (int n = 0; n < tp.half_order; ++n) {
			if (D == X) {
				result += fdc::C2[tp.half_order-1][n+1] *
					(in(i+n+1, j, k, 0) + in(i-n-1, j, k, 0));
			} else if (D == Y) {
				result += fdc::C2[tp.half_order-1][n+1] *
					(in(i, j+n+1, k, 0) + in(i, j-n-1, k, 0));
			} else if (D == Z) {
				result += fdc::C2[tp.half_order-1][n+1] *
					(in(i, j, k+n+1, 0) + in(i, j, k-n-1, 0));
			}
		}
		result /= tp.space_step[D] * tp.space_step[D];
		return result;
	}
};

/**
 * \brief First derivative operator
 */
template <typename T, typename I, CartDir D>
struct D1 {
	typedef DArray<T, I> DA;
	static inline T apply(
		const TaskParams<T, I>& tp,
		I i, I j, I k,
		const DA& in)
	{
		STATIC_ASSERT(D == X || D == Y && D == Z, D1_Wrong_CartDir_param);
		DEBUG_ASSERT(tp.half_order < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		T result = 0;
		for (int n = 0; n < tp.half_order; ++n) {
			if (D == X) {
				result += fdc::C1[tp.half_order-1][n] *
					(in(i+n+1, j, k, 0) + in(i-n-1, j, k, 0));
			} else if (D == Y) {
				result += fdc::C1[tp.half_order-1][n] *
					(in(i, j+n+1, k, 0) + in(i, j-n-1, k, 0));
			} else if (D == Z) {
				result += fdc::C1[tp.half_order-1][n] *
					(in(i, j, k+n+1, 0) + in(i, j, k-n-1, 0));
			}
		}
		result /= tp.space_step[D];
		return result;		
	}
};

/*
 * Next types as params to grid operators
 */

template <typename T, typename I> struct VAL {
	struct Type {
		typedef DArray<T, I> DA;
		static inline T apply(
			const TaskParams<T, I>& tp,
			I i, I j, I k,
			const DA& in)
		{
			(void) tp;
			(void) i; (void) j; (void) k;
			(void) in;
			return m_val;
		}
	};
	void setVal(T val) { m_val = val; }
private:
	static T m_val;
};

template <typename T, typename I> struct ID {
	struct Type {
		typedef DArray<T, I> DA;
		static inline T apply(
			const TaskParams<T, I>& tp,
			I i, I j, I k,
			const DA& in)
		{
			(void) tp;
			return in(i, j, k, 0);
		}
	};
};

template <typename T, typename I> struct D1X { typedef D1<T, I, X> Type; };
template <typename T, typename I> struct D1Y { typedef D1<T, I, Y> Type; };
template <typename T, typename I> struct D1Z { typedef D1<T, I, Z> Type; };

template <typename T, typename I> struct D2X { typedef D2<T, I, X> Type; };
template <typename T, typename I> struct D2Y { typedef D2<T, I, Y> Type; };
template <typename T, typename I> struct D2Z { typedef D2<T, I, Z> Type; };

template <typename T, typename I> struct FD1X { typedef FD1<T, I, X> Type; };
template <typename T, typename I> struct FD1Y { typedef FD1<T, I, Y> Type; };
template <typename T, typename I> struct FD1Z { typedef FD1<T, I, Z> Type; };

template <typename T, typename I> struct BD1X { typedef BD1<T, I, X> Type; };
template <typename T, typename I> struct BD1Y { typedef BD1<T, I, Y> Type; };
template <typename T, typename I> struct BD1Z { typedef BD1<T, I, Z> Type; };

/*
 * Func
 */

template <typename T>
struct FuncId {
	static inline void f(T& out, const T& in) { out = in; }
};

template <typename T>
struct FuncAdd {
	static inline void f(T& out, const T& in) { out = out + in; }
};

template <typename T>
struct FuncSub {
	static inline void f(T& out, const T& in) { out = out - in; }
};

template <typename T>
struct FuncDiv {
	static inline void f(T& out, const T& in) { out = out / in; }
};

template <typename T>
struct FuncSubR {
	static inline void f(T& out, const T& in) { out = in - out; }
};

template <typename T>
struct FuncDivR {
	static inline void f(T& out, const T& in) { out = in / out; }
};

template <typename T>
struct FuncMul {
	static inline void f(T& out, const T& in) { out = out * in; }
};

} // namespace operators

} // namespace rgrid

#endif

