/**
 * \file
 * \brief MPI related functions
 */

/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_MPI_H
#define RGRID_MPI_H

#include <string>

#include "rgrid/utils.h"

#ifdef USE_MPI
#include <mpi.h>

/// Check for errors in execution of MPI functions
#define MPI_CHECK(x) \
	{ int rc = x; \
	if (rc != MPI_SUCCESS) RG_ASSERT(rc == MPI_SUCCESS, STR(#x) + ": " + rgmpi::getError(rc)); }

#endif // USE_MPI

namespace rgmpi
{

#ifdef USE_MPI
/**
 * \brief Get MPI error by code
 */
std::string getError(const int rc);
#endif
/**
 * \brief Init MPI
 * \param[out] argc get number of command line args
 * \param[out] argv get command line args
 */
void init(int *argc, char ***argv);
/**
 * \brief Init MPI
 */
void init();
/**
 * \brief Init MPI
 * \return Is MPI already initialized
 */
bool forceInit();
/**
 * \brief Get rank of process in MPI_COMM_WORLD
 */
int worldRank();
/**
 * \brief Get size of MPI_COMM_WORLD
 */
int worldSize();
#ifdef USE_MPI
/**
 * \brief Get size of communicator
 */
int commSize(MPI_Comm comm);
/**
 * \brief Get rank of current process in communicator
 */
int commRank(MPI_Comm comm);
/**
 * \brief Free communicator
 */
void commFree(MPI_Comm& comm);
/**
 * \brief Get group size
 */
int groupSize(MPI_Group g);
/**
 * \brief Get group rank
 */
int groupRank(MPI_Group s);
#endif
/**
 * \brief Finilize MPI
 */
void forceFinalize();
/**
 * \brief Synchronize all processes in MPI_COMM_WORLD
 */
void barrier();
/**
 * \brief Create cart comm from MPI_COMM_WORLD communicator
 * \param[out] cartComm
 * \param[in] parts number of parts in each direction
 */
#ifdef USE_MPI
void cartCreate(MPI_Comm& cartComm, int const parts[3]);
/**
 * \brief Create cart comm from MPI_COMM_WORLD communicator
 * \param[out] cartComm
 * \param[in] parts number of parts in each direction
 */
template <typename I>
inline void cartCreate(MPI_Comm& cartComm, rgrid::Dim3D<I> parts) {
	int p[3];
	p[0] = parts.x;
	p[1] = parts.y;
	p[2] = parts.z;
	cartCreate(cartComm, p);
}
/**
 * \brief Get coords of process with specific rank in cart comm
 * \param[in] cartComm
 * \param[in] rank
 * \param[out] coords
 */
void cartCoords(MPI_Comm const cartComm, int const rank, int coords[3]);
/**
 * \brief Get rank of process with specific coords in cart comm
 * \param[in] cartComm
 * \param[in] coords
 * \return rank
 */
int cartRank(MPI_Comm const cartComm, int const coords[3]);

/**
 * \brief get MPI type corresponding to specific type
 * \tparam T specific type
 */
template <typename T>
MPI_Datatype getMPItype();

/**
 * \brief Create subarray type (array inside bigger array)
 * \param[in] gsize global size (size of bigger array)
 * \param[in] lsize local size (size of internal array)
 * \param[in] start origins of local array in global array
 * \tparam T types for each element of array
 */
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

/**
 * \brief Release subarray type
 */
void freeSubarrayType(MPI_Datatype& dt);

#endif // USE_MPI

} // namespace rgmpi


#endif // RGRID_MPI_H
