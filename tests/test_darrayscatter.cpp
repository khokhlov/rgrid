// provide main
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
		int rank = 0;
#ifdef USE_MPI		
		rank = rgmpi::worldRank();
#endif
		if (rank == 0) {
			return std::cout;
		} else {
			return emptyoss;
		}
	}
	std::ostream& cerr() {
		int rank = 0;
#ifdef USE_MPI		
		rank = rgmpi::worldRank();
#endif
		if (rank == 0) {
			return std::cerr;
		} else {
			return emptyoss;
		}
	}
}

int main(int argc, char** argv)
{
#ifdef USE_MPI
	rgmpi::init(&argc, &argv);
#endif
	
	int result = Catch::Session().run(argc, argv);
	
#ifdef USE_MPI
	rgmpi::forceFinalize();
#endif
	
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
			d1.resize(100, 200, 150, 3);
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

TEST_CASE("DArrayScatter IO 1 - MPI save") {	
	if (rgmpi::worldSize() == 1 * 1 * 1) {
	
		rgrid::DArrayScatter<int, int> das;
		int const size[3] = { 3, 1, 5 };
		int const gp[3] = { 1, 1, 1 };
		int const lp[3] = { 2, 1, 2 };
		int const ghost[3] = { 0, 0, 0 };
		das.setSizes(size, gp, lp, ghost, 2);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2;
			d1.resize(3, 1, 5, 3, 1, 5, 0, 0, 0, 0, 0, 0, 2);
			d1.fill(2);
			d1(0, 0, 0, 0) = 7;
			d1(2, 0, 1, 0) = 3;
			d1(2, 0, 1, 1) = 5;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io_simple_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d1);
			std::fstream fs;
			fs.open("test_das_io_simple_2.txt", std::ios_base::out);
			d1.saveData(fs, rgrid::rgio::BINARY);
			fs.close();
			fs.open("test_das_io_simple_1.txt", std::ios_base::in);
			d2.loadData(fs);
			fs.close();
			REQUIRE(d1 == d2);
		} else {
			rgrid::DArray<int, int> d1, d2;
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io_simple_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d2);
		}
	}
}

TEST_CASE("DArrayScatter IO 2 - MPI save") {	
	if (rgmpi::worldSize() == 2 * 1 * 1) {
	
		rgrid::DArrayScatter<int, int> das;
		int const size[3] = { 4, 1, 1 };
		int const gp[3] = { 2, 1, 1 };
		int const lp[3] = { 2, 1, 1 };
		int const ghost[3] = { 0, 1, 0 };
		das.setSizes(size, gp, lp, ghost, 1);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2;
			d1.resize(4, 1, 1, 4, 1, 1, 0, 0, 0, 0, 1, 0, 1);
			d1.fill(2);
			d1(0, 0, 0, 0) = 7;
			d1(3, 0, 0, 0) = 3;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d1);
			std::fstream fs;
			fs.open("test_das_io_2.txt", std::ios_base::out);
			d1.saveData(fs, rgrid::rgio::BINARY);
			fs.close();
			fs.open("test_das_io_1.txt", std::ios_base::in);
			d2.loadData(fs);
			fs.close();
			REQUIRE(d1 == d2);
		} else {
			rgrid::DArray<int, int> d1, d2;
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d2);
		}
	}
}

TEST_CASE("DArrayScatter IO 3 - MPI save") {	
	if (rgmpi::worldSize() == 3 * 4 * 5) {
	
		rgrid::DArrayScatter<int, int> das;
		int const size[3] = { 17, 50, 6 };
		int const gp[3] = { 3, 4, 5 };
		int const lp[3] = { 2, 1, 1 };
		int const ghost[3] = { 3, 0, 1 };
		das.setSizes(size, gp, lp, ghost, 2);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2;
			d1.resize(17, 50, 6, 17, 50, 6, 0, 0, 0, 3, 0, 1, 2);
			d1.fill(4);
			d1(12, 23, 0, 0) = 7;
			d1(3, 3, 5, 1) = 3;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io3_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d1);
			std::fstream fs;
			fs.open("test_das_io3_2.txt", std::ios_base::out);
			d1.saveData(fs, rgrid::rgio::BINARY);
			fs.close();
			fs.open("test_das_io3_1.txt", std::ios_base::in);
			d2.loadData(fs);
			fs.close();
			REQUIRE(d1 == d2);
		} else {
			rgrid::DArray<int, int> d1, d2;
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io3_1.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			das.gatherAndGet(0, d2);
		}
	}
}

TEST_CASE("DArrayScatter IO 4 - MPI save/load") {	
	if (rgmpi::worldSize() == 1 * 1 * 1) {
	
		rgrid::DArrayScatter<int, int> das, das2;
		int const size[3] = { 17, 50, 6 };
		int const gp[3] = { 1, 1, 1 };
		int const lp[3] = { 1, 1, 1 };
		int const ghost[3] = { 0, 0, 0 };
		das.setSizes(size, gp, lp, ghost, 1);
		das2.setSizes(size, gp, lp, ghost, 1);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2, d3;
			d1.resize(17, 50, 6, 17, 50, 6, 0, 0, 0, 0, 0, 0, 1);
			d1.fill(4);
			d1(12, 23, 0, 0) = 7;
			d1(3, 3, 5, 0) = 3;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io4.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			std::fstream fs;
			fs.open("test_das_io4.txt", std::ios_base::out);
			d1.saveData(fs, rgrid::rgio::BINARY);
			fs.close();
			
			das.gatherAndGet(0, d3);
						
			das2.loadDataBegin("test_das_io4.txt");
			das2.loadDataEnd();
			das2.gatherAndGet(0, d2);
			REQUIRE(d3 == d2);
		} else {
			rgrid::DArray<int, int> d1;
			das.setAndScatter(0, d1);
			
			das.saveDataBegin("test_das_io4.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			das.gatherAndGet(0, d1);
			
			das2.loadDataBegin("test_das_io4.txt");
			das2.loadDataEnd();
			
 			das2.gatherAndGet(0, d1);
		}
	}
}

TEST_CASE("DArrayScatter IO 5 - MPI save/load") {	
	if (rgmpi::worldSize() == 1 * 2 * 1) {
	
		rgrid::DArrayScatter<int, int> das, das2;
		int const size[3] = { 17, 50, 6 };
		int const gp[3] = { 1, 2, 1 };
		int const lp[3] = { 1, 1, 1 };
		int const ghost[3] = { 0, 0, 0 };
		das.setSizes(size, gp, lp, ghost, 1);
		das2.setSizes(size, gp, lp, ghost, 1);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2, d3;
			d1.resize(17, 50, 6, 17, 50, 6, 0, 0, 0, 0, 0, 0, 1);
			d1.fill(4);
			d1(12, 23, 0, 0) = 7;
			d1(3, 3, 5, 0) = 3;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io5.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			das.gatherAndGet(0, d3);
			
			das2.loadDataBegin("test_das_io5.txt");
			das2.loadDataEnd();
			das2.gatherAndGet(0, d2);
			REQUIRE(d3 == d2);
		} else {
			rgrid::DArray<int, int> d1;
			das.setAndScatter(0, d1);
			
			das.saveDataBegin("test_das_io5.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			das.gatherAndGet(0, d1);
			
			das2.loadDataBegin("test_das_io5.txt");
			das2.loadDataEnd();
 			das2.gatherAndGet(0, d1);			
		}
	}
}

TEST_CASE("DArrayScatter IO - MPI save/load") {	
	if (rgmpi::worldSize() == 10 * 2 * 3) {
	
		rgrid::DArrayScatter<int, int> das, das2;
		int const size[3] = { 17, 50, 6 };
		int const gp[3] = { 10, 2, 3 };
		int const lp[3] = { 1, 1, 1 };
		int const ghost[3] = { 3, 0, 1 };
		das.setSizes(size, gp, lp, ghost, 2);
		das2.setParts(gp, lp, ghost);
		
		if (das.getInternalRank() == 0) {
			rgrid::DArray<int, int> d1, d2, d3;
			d1.resize(17, 50, 6, 17, 50, 6, 0, 0, 0, 3, 0, 1, 2);
			d1.fill(4);
			d1(12, 23, 0, 0) = 7;
			d1(3, 3, 5, 1) = 3;
			d1.fillGhost();
			das.setAndScatter(0, d1);
			das.saveDataBegin("test_das_io_mpisl.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			das.gatherAndGet(0, d3);
			
			das2.loadDataBegin("test_das_io_mpisl.txt");
			das2.loadDataEnd();
			das2.gatherAndGet(0, d2);
			REQUIRE(d3 == d2);
		} else {
			rgrid::DArray<int, int> d1;
			das.setAndScatter(0, d1);
			
			das.saveDataBegin("test_das_io_mpisl.txt", rgrid::rgio::BINARY);
			das.saveDataEnd();
			
			das.gatherAndGet(0, d1);
			
			das2.loadDataBegin("test_das_io_mpisl.txt");
			das2.loadDataEnd();
 			das2.gatherAndGet(0, d1);			
		}
	}
}

TEST_CASE("DArrayScatter external sync") {
	if (rgmpi::worldSize() == 2 * 2 * 2) {
	
		rgrid::DArrayScatter<int, int> das;
		int const gp[3] = { 2, 2, 2 };
		int const lp[3] = { 1, 1, 1 };
		int const ghost[3] = { 1, 1, 1 };
		das.setParts(gp, lp, ghost);
		
		rgrid::DArray<int, int> d1;
		
		if (das.getInternalRank() == 0) {
			d1.resize(2, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1);
			d1.fill(4);
			d1(0, 0, 0, 0) = 7;
			d1(1, 0, 0, 0) = 3;
		}	
		
		das.setAndScatter(0, d1);
			
		rgrid::DArrayContainer<int, int>& dac = das.getLocalContainer();
			
		das.externalSyncStart();
		das.externalSyncEnd();
		
		if (das.getInternalRank() == 0) {
			REQUIRE(dac.getDArrayPart(0).val(1, 0, 0, 0) == 3);
		}
		if (das.getInternalRank() == 1) {
			REQUIRE(dac.getDArrayPart(0).val(-1, 0, 0, 0) == 7);
		}
		if (das.getInternalRank() == 3) {
			REQUIRE(dac.getDArrayPart(0).val(0, 0, 1, 0) == 4);
		}
		
		das.gatherAndGet(0, d1);
	}
}

TEST_CASE("DArrayScatter external sync 2") {
	if (rgmpi::worldSize() == 3 * 4 * 5) {
	
		rgrid::DArrayScatter<int, int> das;
		int const gp[3] = { 3, 4, 5 };
		int const lp[3] = { 2, 1, 1 };
		int const ghost[3] = { 3, 0, 1 };
		das.setParts(gp, lp, ghost);
		
		rgrid::DArray<int, int> d1;
		
		if (das.getInternalRank() == 0) {
			d1.resize(17, 50, 6, 17, 50, 6, 0, 0, 0, 3, 0, 1, 2);
			d1.fill(4);
			
			d1(6, 0, 0, 1) = 7;
			d1(8, 0, 0, 0) = 5;
			
			d1(6, 13, 1, 0) = 10;
		}
		
		das.setAndScatter(0, d1);
			
		rgrid::DArrayContainer<int, int>& dac = das.getLocalContainer();

 		das.externalSyncStart();
 		das.externalSyncEnd();

		if (das.getInternalRank() == 0) {
			REQUIRE(dac.getDArrayPart(1).val(3, 0, 0, 1) == d1(6, 0, 0, 1));
			REQUIRE(dac.getDArrayPart(1).val(4, 0, 0, 1) == d1(7, 0, 0, 1));
			REQUIRE(dac.getDArrayPart(1).val(5, 0, 0, 1) == d1(8, 0, 0, 1));
			REQUIRE(dac.getDArrayPart(1).val(5, 0, 0, 0) == d1(8, 0, 0, 0));
		}
		if (das.getInternalPos(rgrid::X) == 0 && das.getInternalPos(rgrid::Y) == 1 && das.getInternalPos(rgrid::Z) == 0) {
			REQUIRE(dac.getDArrayPart(1).val(3, 0, 1, 0) == 10);
		}
		if (das.getInternalPos(rgrid::X) == 1 && das.getInternalPos(rgrid::Y) == 1 && das.getInternalPos(rgrid::Z) == 1) {
			REQUIRE(dac.getDArrayPart(0).val(0, 0, -1, 0) == 10);
		}
	}
}

#endif /* USE_MPI */
