/**
 * \file
 * \brief Helpful utils
 */

/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_UTILS_H
#define RGRID_UTILS_H

#include <string>
#include <sstream>

#include "rgrid/types.h"

/// Try to convert any object to string
#define STR(x) rgrid::toString(x)

/**
 * \brief Check is expression is true, overwise show message and throw exception
 * \param[in] exp
 * \param[in] message
 */
#define RG_ASSERT(exp, message) \
	if (!(exp)) rgrid::rgassert(__FILE__, __LINE__, STR(message).c_str(), #exp);

namespace rgrid
{
/**
 * \brief Show error message and throw exception
 */
void rgassert(const char *file, const int line, const char *msg, const char *exp);
/**
 * \brief Get two different directions ortogonal to first
 * \param[in] axis first direction
 * \param[out] dir1 direction ortogonal to first
 * \param[out] dir2 direction ortogonal to first
 */
void ortDirs(const CartDir &axis, CartDir &dir1, CartDir &dir2);

/**
 * Try to convert any object to string
 */
template <typename T> std::string toString(const T& x)
{
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

}

#endif // RGRID_UTILS_H
