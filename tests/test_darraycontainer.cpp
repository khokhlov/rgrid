#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"
#include "rgrid/rgio.h"

TEST_CASE(
		"DArrayContainer",
		"splitAndCombine"
	 )
{
	rgrid::DArray<int, int> d1, d2;
	d1.resize(100, 200, 150);
	d1.alloc(3);
	d1.fill(2);
	d1(3, 3, 3, 0) = 7;
	d1(3, 3, 3, 1) = 8;
	d1.fillGhost();
	rgrid::DArrayContainer<int, int> dac(d1, 13, 2, 7);
	dac.getDArray(d2);
	REQUIRE(d1 == d2);
	dac.getNode(35, 50, 51, 2) = 777;
	dac.getDArray(d2);
	d1(35, 50, 51, 2) = 777;
	REQUIRE(d1 == d2);
}
