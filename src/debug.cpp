#include "rgrid/debug.h"

#ifdef __GLIBC__
	#include <execinfo.h>
#endif

#include <iostream>
#include <cstdlib>

using namespace std;

namespace dassert {

void stackTrace() {
#ifdef __GLIBC__
	void * array[250];
	int nSize = backtrace(array, 250);
	char ** symbols = backtrace_symbols(array, nSize);
	for (int i = 0; i < nSize; i++)
		cerr << symbols[i] << endl;
	free(symbols);
#endif
}

void assert(const char *file, const int line, const char *msg, const char *exp)
{
	cerr << "\033[91mASSERT! " << file << ":" << line << ": " << exp << "\033[0m";
	cerr << " : " << msg << endl;
	stackTrace();
	exit(-1);
}

} // namespace dassert
