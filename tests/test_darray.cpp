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
	d(0, 0) = 2;
	d(9, 0) = 3;
	d.fillGhost(rgrid::X);
	REQUIRE(d(-1, 0) == 2);
	REQUIRE(d(-2, 0) == 2);
	REQUIRE(d(10, 0) == 3);
	REQUIRE(d(11, 0) == 3);
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
	d(0, 0) = 2;
	d(9, 0) = 3;
	d.fillGhost(rgrid::X);
	REQUIRE(d(-1, 0) == 2);
	REQUIRE(d(-2, 0) == 2);
	REQUIRE(d(-3, 0) == 2);
	REQUIRE(d(10, 0) == 3);
	REQUIRE(d(11, 0) == 3);
	REQUIRE(d(12, 0) == 3);
}

TEST_CASE(
		"DArray",
		"components"
	 )
{
	DArray d;
	d.resize(2, 3, 4);
	d.alloc(5);
	d.fill(18);
	for (int i = 0; i != 5; ++i) {
		d(1, 0, 2, i) = 77+i;
	}
	d.fillGhost();
	for (int i = 0; i != 5; ++i) {
		REQUIRE(d(1, -1, 2, i) == d(1, 0, 2, i));
		REQUIRE(d(1, -2, 2, i) == d(1, 0, 2, i));
	}
}