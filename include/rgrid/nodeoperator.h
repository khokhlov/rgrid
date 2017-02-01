#ifndef NODE_OPERATOR_H
#define NODE_OPERATOR_H

#include "rgrid/darrayscatter.h"
#include "rgrid/fdc.h"

namespace rgrid {

namespace operators {

/**
 * \brief First derivative staggered forward operator
 */
template <typename T, typename I, int HalfOrder>
struct DxF {
	typedef DArray<T, I> DA;
	static inline void apply(
		I i, I j, I k,
		const DA& in, DA& out)
	{
		DEBUG_ASSERT(HalfOrder < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		out(i, j, k, 0) = 0;
		for (int n = 0; n < HalfOrder; ++n) {
			out(i, j, k, 0) += fdc::SC1[HalfOrder-1][n] *
				(in(i+n+1, j, k, 0) - in(i-n, j, k, 0));
		}
	}
};

/**
 * \brief First derivative staggered backward operator
 */
template <typename T, typename I, int HalfOrder>
struct DxB {
	typedef DArray<T, I> DA;
	static inline void apply(
		I i, I j, I k,
		const DA& in, DA& out)
	{
		DEBUG_ASSERT(HalfOrder < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		out(i, j, k, 0) = 0;
		for (int n = 0; n < HalfOrder; ++n) {
			out(i, j, k, 0) += fdc::SC1[HalfOrder-1][n] *
				(in(i+n, j, k, 0) - in(i-n-1, j, k, 0));
		}
	}
};

/**
 * \brief Second derivative operator
 */
template <typename T, typename I, int HalfOrder>
struct D2x {
	typedef DArray<T, I> DA;
	static inline void apply(
		I i, I j, I k,
		const DA& in, DA& out)
	{
		DEBUG_ASSERT(HalfOrder < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		out(i, j, k, 0) = fdc::C2[HalfOrder-1][0] * in(i, j, k, 0);
		for (int n = 0; n < HalfOrder; ++n) {
			out(i, j, k, 0) += fdc::C2[HalfOrder-1][n+1] *
				(in(i+n+1, j, k, 0) + in(i-n-1, j, k, 0));
		}
	}
};

/**
 * \brief First derivative operator
 */
template <typename T, typename I, int HalfOrder>
struct Dx {
	typedef DArray<T, I> DA;
	static inline void apply(
		I i, I j, I k,
		const DA& in, DA& out)
	{
		DEBUG_ASSERT(HalfOrder < fdc::MAX_HALF_ORDER,
			"There are no such precomputed coefficients");
		out(i, j, k, 0) = 0;
		for (int n = 0; n < HalfOrder; ++n) {
			out(i, j, k, 0) += fdc::C1[HalfOrder-1][n] *
				(in(i+n+1, j, k, 0) + in(i-n-1, j, k, 0));
		}
	}
};

} // namespace operators

} // namespace rgrid

#endif
