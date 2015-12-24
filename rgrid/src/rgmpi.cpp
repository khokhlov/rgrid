#include "rgrid/rgmpi.h"

#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>

#ifdef USE_MPI

using namespace std;

namespace rgmpi
{
std::string getError(const int rc)
{
	char error_string[BUFSIZ];
	int length_of_error_string, error_class;
	MPI_Error_class(rc, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	ostringstream ss;
	ss << " rgmpi::rc: " << rc;
	ss << ", rank: " << rgmpi::worldRank();
	ss << ", error class: \"" << error_string;
	MPI_Error_string(rc, error_string, &length_of_error_string);
	ss << "\", error msg: \"" << error_string << "\"";
	return ss.str();
}

void init(int *argc, char ***argv)
{
	MPI_CHECK(MPI_Init(argc, argv));
}

void init()
{
	init(NULL, NULL);
}

bool forceInit()
{
	int flag;
	MPI_CHECK(MPI_Initialized(&flag));
	if (flag == 0) {
		init();
	}
	return flag == 1;
}

int worldSize()
{
	return commSize(MPI_COMM_WORLD);
}

void barrier()
{
	MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
}

int worldRank()
{
	return commRank(MPI_COMM_WORLD);
}

int commSize(MPI_Comm comm)
{
	int size;
	RG_ASSERT(comm != MPI_COMM_NULL, "rgmpi::commRank: Invalid communicator.");
	MPI_CHECK(MPI_Comm_size(comm, &size));
	return size;
}

int commRank(MPI_Comm comm)
{
	int rank;
	RG_ASSERT(comm != MPI_COMM_NULL, "rgmpi::commRank: Invalid communicator.");
	MPI_CHECK(MPI_Comm_rank(comm, &rank));
	return rank;
}

int groupSize(MPI_Group g)
{
	int size;
	MPI_CHECK(MPI_Group_size(g, &size));
	return size;
}

int groupRank(MPI_Group g)
{
	int rank;
	MPI_CHECK(MPI_Group_rank(g, &rank));
	return rank;
}

void forceFinalize()
{
	int flag;
	MPI_CHECK(MPI_Finalized(&flag));
	if (flag == 0) {
		MPI_CHECK(MPI_Finalize());
	}
}
};
#endif // USE_MPI
