#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "rgrid/darray.h"
#include "rgrid/rgio.h"

TEST_CASE(
		"rgio",
		"DArrayInString"
	 )
{
	rgrid::DArray<int, int> d1, d2, d3;
	d1.resize(7, 6, 5);
	d1.alloc(5);
	d1.fill(17);
	d1(4, 5, 1, 0) = 1;
	d1(4, 5, 1, 2) = 2;
	d1.fillGhost();
	std::stringstream ss1, ss2;
	rgrid::rgio::saveData(ss1, d1, rgrid::rgio::BINARY);
	rgrid::rgio::loadData(ss1, d2);
	REQUIRE(d1 == d2);
	rgrid::rgio::saveData(ss2, d2, rgrid::rgio::TEXT);
	rgrid::rgio::loadData(ss2, d3);
	REQUIRE(d1 == d3);
}

TEST_CASE(
		"rgio",
		"DArrayInFile"
	 )
{
	rgrid::DArray<int, int> d1, d2, d3;
	d1.resize(5, 10, 1);
	d1.alloc(3);
	d1.fill(17);
	d1(4, 7, 0, 0) = 4;
	d1(4, 5, 0, 2) = 5;
	d1.fillGhost();
	std::fstream fs("testfile1.txt", std::ios_base::out | std::ios_base::trunc);
	rgrid::rgio::saveData(fs, d1, rgrid::rgio::BINARY);
	fs.close();
	fs.open("testfile1.txt", std::ios_base::in);
	rgrid::rgio::loadData(fs, d2);
	fs.close();
	REQUIRE(d1 == d2);
	fs.open("testfile2.txt", std::ios_base::out | std::ios_base::trunc);
	rgrid::rgio::saveData(fs, d2, rgrid::rgio::TEXT);
	fs.close();
	fs.open("testfile2.txt", std::ios_base::in);
	rgrid::rgio::loadData(fs, d3);
	fs.close();
	REQUIRE(d1 == d3);
}
