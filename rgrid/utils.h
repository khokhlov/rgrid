/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, (C) 2015
 */

#ifndef RGRID_UTILS_H
#define RGRID_UTILS_H

#define RG_ASSERT(exp, message) \
	if (!(exp)) rgrid::rgssert(__FILE__, __LINE__, STR(message).c_str(), #exp);

namespace rgrig
{
void rgassert(const char *file, const int line, const char *msg, const char *exp);
};

#endif // RGRID_UTILS_H
