#ifndef GRID_OPERATOR_H
#define GRID_OPERATOR_H

#include "rgrid/darrayscatter.h"

namespace rgrid {

namespace operators {

template <typename T, typename I, typename NodeOperator>
struct GridOperator {
	typedef DArrayScatter<T, I> DAS;
	typedef DArray<T, I> DA;
	static void apply(
		const Range<I>& r, const DAS& in, DAS& out)
	{
		for (I gk = 0; gk != in.numParts(Z); ++gk)
		for (I gj = 0; gj != in.numParts(Y); ++gj)
		for (I gi = 0; gi != in.numParts(X); ++gi) {
			const DA& u = in.getDArrayPart(gi, gj, gk);
			DA& nu = out.getDArrayPart(gi, gj, gk);
			
			Range<I> glob = intersect(r, u.getRange());
			Range<I> loc = glob.shift(-u.origin());
			
			for (I k = loc.get(SIDE_LEFT, Z); k != loc.get(SIDE_RIGHT, Z); ++k)
			for (I j = loc.get(SIDE_LEFT, Y); j != loc.get(SIDE_RIGHT, Y); ++j)
			for (I i = loc.get(SIDE_LEFT, X); i != loc.get(SIDE_RIGHT, X); ++i) {
				NodeOperator::apply(i, j, k, in, out);
			}
		}
	}
};

} // namespace operators

} // namespace rgrid

#endif
