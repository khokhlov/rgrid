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
};
#endif // USE_MPI

#endif // RGRID_MPI_H
