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

void ortDirs(const CartDir &axis, CartDir &dir1, CartDir &dir2)
{
	if (axis == X) {
		dir1 = Y;
		dir2 = Z;
	} else if (axis == Y) {
		dir1 = Z;
		dir2 = X;
	} else if (axis == Z) {
		dir1 = X;
		dir2 = Y;
	} else {
		dir1 = dir2 = DIR_UNDEFINED;
	}
}
};
