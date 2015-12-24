/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_UTILS_H
#define RGRID_UTILS_H

#include <string>
#include <sstream>

#include "rgrid/types.h"

#define STR(x) rgrid::toString(x)

#define RG_ASSERT(exp, message) \
	if (!(exp)) rgrid::rgassert(__FILE__, __LINE__, STR(message).c_str(), #exp);


namespace rgrid
{
void rgassert(const char *file, const int line, const char *msg, const char *exp);
void ortDirs(const CartDir &axis, CartDir &dir1, CartDir &dir2);

template <typename T> std::string toString(const T& x)
{
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

};

#endif // RGRID_UTILS_H
