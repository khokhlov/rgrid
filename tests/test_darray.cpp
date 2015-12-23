#define CATCH_CONFIG_MAIN


#include "catch.hpp"
#include "rgrid/darray.h"

#include <iostream>

using namespace std;

typedef rgrid::DArray<int, int> DArray;

TEST_CASE(
		"DArray",
		"alloc"
	 )
{
	DArray d;
	REQUIRE(d.size() == 0);
	REQUIRE(d.localSize() == 0);
	REQUIRE(d.data.size() == 0);
	d.resize(1, 2, 3);
	REQUIRE(d.data.size() != d.localSizeGhost());
	d.alloc();
	REQUIRE(d.data.size() == d.localSizeGhost());
}

TEST_CASE(
		"DArray",
		"fill"
	 )
{
	DArray d;
	d.resize(2, 2, 3);
	d.alloc();
	d.fill(1);
	bool ok = true;
	for (int i = 0; i < d.localSizeGhost(); i++) { ok &= (d[i] == 1); }
	REQUIRE(ok);
	d.fill(2);
	ok = true;
	for (int i = 0; i < d.localSizeGhost(); i++) { ok &= (d[i] == 2); }
	REQUIRE(ok);
}

TEST_CASE(
		"DArray",
		"fillGhost1"
	 )
{
	DArray d;
	d.resize(10);
	d.alloc();
	d.fill(1);
	d(0) = 2;
	d(9) = 3;
	d.fillGhost(rgrid::X);
	REQUIRE(d(-1) == 2);
	REQUIRE(d(-2) == 2);
	REQUIRE(d(10) == 3);
	REQUIRE(d(11) == 3);
}

TEST_CASE(
		"DArray",
		"fillGhost2"
	 )
{
	DArray d;
	d.resize(10, 1, 1, 3, 0, 0);
	d.alloc();
	d.fill(1);
	d(0) = 2;
	d(9) = 3;
	d.fillGhost(rgrid::X);
	REQUIRE(d(-1) == 2);
	REQUIRE(d(-2) == 2);
	REQUIRE(d(-3) == 2);
	REQUIRE(d(10) == 3);
	REQUIRE(d(11) == 3);
	REQUIRE(d(12) == 3);
}
