#define CATCH_CONFIG_MAIN


#include "catch.hpp"
#include "rgrid/pdim.h"

#include <iostream>

using namespace std;

typedef rgrid::PDim<int> PDim;

TEST_CASE("pdim creation")
{
	PDim d;
	d.resize(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	REQUIRE(d.size() == 0);
	REQUIRE(d.localSize() == 0);
	REQUIRE(d.origin(rgrid::X) == 0);
	REQUIRE(d.origin(rgrid::Y) == 0);
	REQUIRE(d.origin(rgrid::Z) == 0);
}

TEST_CASE("prim resize")
{
	PDim d;
	d.resize(1, 2, 3);
	REQUIRE(d.size() == 6);
	REQUIRE(d.localSize() == 6);

	REQUIRE(d.origin(rgrid::X) == 0);
	REQUIRE(d.origin(rgrid::Y) == 0);
	REQUIRE(d.origin(rgrid::Z) == 0);

	REQUIRE(d.localSize(rgrid::X) == 1);
	REQUIRE(d.localSize(rgrid::Y) == 2);
	REQUIRE(d.localSize(rgrid::Z) == 3);

	REQUIRE(d.size(rgrid::X) == 1);
	REQUIRE(d.size(rgrid::Y) == 2);
	REQUIRE(d.size(rgrid::Z) == 3);

	PDim d1(d);
	REQUIRE(d1 == d);
	d.resize(3, 2, 1, 1, 2, 3);
	REQUIRE_FALSE(d1 == d);

	REQUIRE(d.localSize(rgrid::X) == 3);
	REQUIRE(d.localSize(rgrid::Y) == 2);
	REQUIRE(d.localSize(rgrid::Z) == 1);

	REQUIRE(d.size(rgrid::X) == 3);
	REQUIRE(d.size(rgrid::Y) == 2);
	REQUIRE(d.size(rgrid::Z) == 1);

	REQUIRE(d.ghost(rgrid::X) == 1);
	REQUIRE(d.ghost(rgrid::Y) == 2);
	REQUIRE(d.ghost(rgrid::Z) == 3);
}
