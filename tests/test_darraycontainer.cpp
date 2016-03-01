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

TEST_CASE(
		"DArrayContainer",
		"fillGhost"
	 )
{
	rgrid::DArray<int, int> d;
	d.resize(30, 20, 10);
	d.alloc(2);
	d.fill(3);
	d.val(0, 0, 0, 0) = 5;
	d.fillGhost(rgrid::Y, rgrid::SIDE_LEFT);
	REQUIRE(d.val(0, -1, 0, 0) == 5);
	rgrid::DArrayContainer<int, int> dac(d, 5, 5, 3);
	dac.getNode(5, 2, 3, 0) = 17;
	dac.fillGhost();
	REQUIRE(dac.getDArrayPart(0, 0, 0).val(-2, 0, 0, 0) == 5);
	REQUIRE(dac.getDArrayPart(0, 0, 0).val(5, 2, 3, 0) == 17);
	REQUIRE(dac.getDArrayPart(1, 0, 0).val(-1, 2, 3, 0) == 17);
	REQUIRE(dac.getDArrayPart(0, 1, 0).val(5, -2, 3, 0) == 17);
	REQUIRE(dac.getDArrayPart(0, 0, 1).val(5, 2, -1, 0) == 17);
}
