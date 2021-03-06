#include "rgrid/rgmpi.h"

#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

namespace rgmpi
{

#ifdef USE_MPI
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
#endif

void init(int *argc, char ***argv)
{
#ifdef USE_MPI
	MPI_CHECK(MPI_Init(argc, argv));
#else
	(void)argc;
	(void)argv;
#endif
}

void init()
{
	init(NULL, NULL);
}

bool forceInit()
{
#ifdef USE_MPI
	int flag;
	MPI_CHECK(MPI_Initialized(&flag));
	if (flag == 0) {
		init();
	}
	return flag == 1;
#else
	return 0;
#endif
}

int worldSize()
{
#ifdef USE_MPI
	return commSize(MPI_COMM_WORLD);
#else
	return 1;
#endif
}

void barrier()
{
#ifdef USE_MPI
	MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
#endif
}

int worldRank()
{
#ifdef USE_MPI
	return commRank(MPI_COMM_WORLD);
#else
	return 0;
#endif
}

#ifdef USE_MPI
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

void commFree(MPI_Comm& comm)
{
	if (comm == MPI_COMM_NULL) return;
	MPI_CHECK(MPI_Comm_free(&comm));
	comm = MPI_COMM_NULL;
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
#endif

void forceFinalize()
{
#ifdef USE_MPI
	int flag;
	MPI_CHECK(MPI_Finalized(&flag));
	if (flag == 0) {
		MPI_CHECK(MPI_Finalize());
	}
#endif
}

#ifdef USE_MPI
void cartCreate(MPI_Comm& cartComm, int const parts[3])
{
	int periods[3] = { false, false, false };
	int iParts[3] = { parts[2], parts[1], parts[0] };
	MPI_CHECK(MPI_Cart_create(MPI_COMM_WORLD, 3, iParts, periods, false, &cartComm));
}

void cartCoords(MPI_Comm cartComm, int const rank, int coords[3])
{
	int iCoords[3];
	MPI_CHECK(MPI_Cart_coords(cartComm, rank, 3, iCoords));
	for (int i = 0; i != 3; ++i) {
		coords[2 - i] = iCoords[i];
	}
}

rgrid::Dim3D<int> cartCoords(MPI_Comm const cartComm, int const rank) {
	int iCoords[3];
	rgrid::Dim3D<int> coords;
	MPI_CHECK(MPI_Cart_coords(cartComm, rank, 3, iCoords));
	for (int i = 0; i != 3; ++i) {
		coords[static_cast<rgrid::CartDir>(2 - i)] = iCoords[static_cast<rgrid::CartDir>(i)];
	}
	return coords;
}

int cartRank(MPI_Comm const cartComm, int const coords[3])
{
	int rank;
	int iCoords[3] = { coords[2], coords[1], coords[0] };
	MPI_CHECK(MPI_Cart_rank(cartComm, iCoords, &rank));
	return rank;
}

int cartRank(MPI_Comm const cartComm, const rgrid::Dim3D<int>& coords) {
	int rank;
	int iCoords[3] = { coords[rgrid::Z], coords[rgrid::Y], coords[rgrid::X] };
	MPI_CHECK(MPI_Cart_rank(cartComm, iCoords, &rank));
	return rank;
}

template <> MPI_Datatype getMPItype<int>() { return MPI_INT; }
template <> MPI_Datatype getMPItype<float>() { return MPI_FLOAT; }
template <> MPI_Datatype getMPItype<double>() { return MPI_DOUBLE; }

void freeSubarrayType(MPI_Datatype& dt) {
	MPI_CHECK(MPI_Type_free(&dt));
}
#endif

} /* namespace rgmpi */
