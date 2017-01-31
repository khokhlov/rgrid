#include "rgrid/rgio.h"

namespace rgrid {
	
namespace rgio {
	
#ifdef USE_MPI

long getLineMPI(MPI_File fh, std::string& str, char delim) {
	str.clear();
	char c;
	long count = 0;
	while (1) {
		MPI_CHECK(MPI_File_read(fh, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE));
		++count;
		if (c == delim) break;
		str += c;
	}
	return count;
}

#endif // USE_MPI

} // namespace rgio

} // namespace rgrid