// prowide main
#define CATCH_CONFIG_RUNNER
// provide cout, cerr
#define CATCH_CONFIG_NOSTDOUT

#include "catch.hpp"

#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"
#include "rgrid/rgio.h"
#include "rgrid/rgmpi.h"
#include "rgrid/darrayscatter.h"

// disable output from processes with rank != 0
namespace Catch {
	std::ostringstream emptyoss;
	std::ostream& cout() {
		if (rgmpi::worldRank() == 0) {
			return std::cout;
		} else {
			return emptyoss;
		}
	}
	std::ostream& cerr() {
		if (rgmpi::worldRank() == 0) {
			return std::cerr;
		} else {
			return emptyoss;
		}
	}
}

int main(int argc, char** argv)
{
	rgmpi::init(&argc, &argv);
	
	int result = Catch::Session().run(argc, argv);
	
	rgmpi::forceFinalize();
	
	return result;
}

#ifdef USE_MPI

TEST_CASE("DArrayScatter scatter and gather")
{	
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
}

TEST_CASE("DArrayScatter input/output") {	
	if (rgmpi::worldSize() == 2 * 1 * 1) {
	
		rgrid::DArrayScatter<int, int> das;
		int const size[3] = { 4, 1, 1 };
		int const gp[3] = { 2, 1, 1 };
		int const lp[3] = { 1, 1, 1 };
		int const ghost[3] = { 0, 0, 0 };
		das.setSizes(size, gp, lp, ghost, 1);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2;
			d1.resize(4, 1, 1, 4, 1, 1, 0, 0, 0, 0, 0, 0);
			d1.alloc(1);
			d1.fill(2);
			d1(0, 0, 0, 0) = 7;
			d1(3, 0, 0, 0) = 8;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d2);
			REQUIRE(d1 == d2);
			std::fstream fs;
			fs.open("test2.txt", std::ios_base::out);
			rgrid::rgio::saveData(fs, d2, rgrid::rgio::BINARY);
			fs.close();
		} else {
			rgrid::DArray<int, int> d1, d2;
			das.setAndScatter(0, d1);
			das.saveDataBegin("test.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d2);
		}
	}
}

#endif /* USE_MPI */
