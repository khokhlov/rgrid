#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"
#include "rgrid/rgio.h"
#include "rgrid/rgmpi.h"
#include "rgrid/darrayscatter.h"

#ifdef USE_MPI

TEST_CASE(
		"DArrayScatter",
		"scatter and gather"
		)
{
	rgmpi::forceInit();
	
	if (rgmpi::worldSize() == 3 * 2 * 2) {
	
		rgrid::DArrayScatter<int, int> das;
		int const size[3] = { 100, 200, 150 };
		int const gp[3] = { 3, 2, 2 };
		int const lp[3] = { 3, 4, 5 };
		int const ghost[3] = { 2, 2, 2 };
		das.setSizes(size, gp, lp, ghost, 3);
			
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2;
			d1.resize(100, 200, 150);
			d1.alloc(3);
			d1.fill(2);
			d1(3, 3, 3, 0) = 7;
			d1(3, 3, 3, 1) = 8;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.gatherAndGet(0, d2);
			REQUIRE(d1 == d2);
		} else {
			rgrid::DArray<int, int> d1, d2;
			das.setAndScatter(0, d1);
			das.gatherAndGet(0, d2);
		}
	}
	
	rgmpi::forceFinalize();
}

#endif /* USE_MPI */
