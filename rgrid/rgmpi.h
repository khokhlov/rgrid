/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_MPI_H
#define RGRID_MPI_H

#include <string>

#include "rgrid/utils.h"

#ifdef USE_MPI
#include <mpi.h>

#define MPI_CHECK(x) \
	{ int rc = x; \
	if (rc != MPI_SUCCESS) RG_ASSERT(rc == MPI_SUCCESS, STR(#x) + ": " + rgmpi::getError(rc)); }


namespace rgmpi
{
	
std::string getError(const int rc);
void init(int *argc, char ***argv);
void init();
bool forceInit();
int worldRank();
int worldSize();
int commSize(MPI_Comm comm);
int commRank(MPI_Comm comm);
int groupSize(MPI_Group g);
int groupRank(MPI_Group s);
void forceFinalize();
void barrier();
void cartCreate(MPI_Comm& cartComm, int const parts[3]);
void cartCoords(MPI_Comm const cartComm, int const rank, int coords[3]);
int cartRank(MPI_Comm const cartComm, int const coords[3]);

template <typename T>
MPI_Datatype getMPItype();

template <typename T>
MPI_Datatype createSubarrayType(int gsize[4], int lsize[4], int start[4]) {
	int igsize[4] = { gsize[3], gsize[2], gsize[1], gsize[0] };
	int ilsize[4] = { lsize[3], lsize[2], lsize[1], lsize[0] };
	int istart[4] = { start[3], start[2], start[1], start[0] };
	MPI_Datatype dt;
// 	igsize[0] *= sizeof(T);
// 	ilsize[0] *= sizeof(T);
// 	istart[0] *= sizeof(T);
	MPI_CHECK(MPI_Type_create_subarray(4, igsize, ilsize, istart, MPI_ORDER_C, getMPItype<T>()/*MPI_BYTE*/, &dt));
	MPI_CHECK(MPI_Type_commit(&dt));
	return dt;
}

void freeSubarrayType(MPI_Datatype& dt);

};
#endif // USE_MPI

#endif // RGRID_MPI_H
