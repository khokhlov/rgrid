#include "rgrid/utils.h"

#include <iostream>
#include <stdexcept>

using namespace std;

namespace rgrid
{
void rgassert(const char *file, const int line, const char *msg, const char *exp)
{
	cerr << "\033[91mERROR! " << file << ":" << line << ": " << exp << "\033[0m" << endl;
	throw logic_error(msg);
}
};
