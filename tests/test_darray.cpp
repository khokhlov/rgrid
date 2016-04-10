#define CATCH_CONFIG_MAIN


#include "catch.hpp"
#include "rgrid/darray.h"

#include <iostream>

using namespace std;

typedef rgrid::DArray<int, int> DArray;

TEST_CASE("DArray resize & alloc")
{
	DArray d;
	REQUIRE(d.globalSize() == 0);
	REQUIRE(d.localSize() == 0);
	REQUIRE(d.dataSize() == 0);
	d.resize(1, 2, 3, 1);
	REQUIRE(d.dataSize() == d.localGhostSize());
}

TEST_CASE("DArray fill")
{
	DArray d;
	d.resize(2, 2, 3, 1);
	d.fill(1);
	bool ok = true;
	for (int i = 0; i < d.localGhostSize(); i++) { ok &= (d[i] == 1); }
	REQUIRE(ok);
	d.fill(2);
	ok = true;
	for (int i = 0; i < d.localGhostSize(); i++) { ok &= (d[i] == 2); }
	REQUIRE(ok);
}

TEST_CASE("DArray fillGhost", "fillGhost1")
{
	SECTION("resize 1") {
		DArray d;
		d.resize(10, 1, 1, 1);
		d.fill(1);
		d(0, 0) = 2;
		d(9, 0) = 3;
		d.fillGhost(rgrid::X);
		REQUIRE(d(-1, 0) == 2);
		REQUIRE(d(-2, 0) == 2);
		REQUIRE(d(10, 0) == 3);
		REQUIRE(d(11, 0) == 3);
	}
	SECTION("resize 2") {
		DArray d;
		d.resize(10, 1, 1, 3, 0, 0, 1);
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
}

TEST_CASE("DArray components")
{
	DArray d;
	d.resize(2, 3, 4, 5);
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

TEST_CASE("DArray IO in string")
{
	rgrid::DArray<int, int> d1, d2, d3;
	d1.resize(7, 6, 5, 5);
	d1.fill(17);
	d1(4, 5, 1, 0) = 1;
	d1(4, 5, 1, 2) = 2;
	d1.fillGhost();
	std::stringstream ss1, ss2;
	d1.saveData(ss1, rgrid::rgio::BINARY);
	d2.loadData(ss1);
	REQUIRE(d1 == d2);
	d2.saveData(ss2, rgrid::rgio::TEXT);
	d3.loadData(ss2);
	REQUIRE(d1 == d3);
}

TEST_CASE("DArray IO in file")
{
	rgrid::DArray<int, int> d1, d2, d3;
	d1.resize(5, 10, 1, 3);
	d1.fill(17);
	d1(4, 7, 0, 0) = 4;
	d1(4, 5, 0, 2) = 5;
	d1.fillGhost();
	std::fstream fs("testfile1.txt", std::ios_base::out | std::ios_base::trunc);
	d1.saveData(fs, rgrid::rgio::BINARY);
	fs.close();
	fs.open("testfile1.txt", std::ios_base::in);
	d2.loadData(fs);
	fs.close();
	REQUIRE(d1 == d2);
	fs.open("testfile2.txt", std::ios_base::out | std::ios_base::trunc);
	d2.saveData(fs, rgrid::rgio::TEXT);
	fs.close();
	fs.open("testfile2.txt", std::ios_base::in);
	d3.loadData(fs);
	fs.close();
	REQUIRE(d1 == d3);
}
