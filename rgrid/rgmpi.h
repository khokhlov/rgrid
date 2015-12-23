/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_MPI_H
#define RGRID_MPI_H

#include <string>

#ifdef USE_MPI
#include <mpi.h>

namespace rgmpi
{
std::string getError(const int rc);
};
#endif // USE_MPI

#endif // RGRID_MPI_H
