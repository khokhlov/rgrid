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
	//ss << ", rank: " << mpiutils::commRank(MPI_COMM_WORLD);
	ss << ", error class: \"" << error_string;
	MPI_Error_string(rc, error_string, &length_of_error_string);
	ss << "\", error msg: \"" << error_string << "\"";
	return ss.str();
}
};
#endif // USE_MPI
