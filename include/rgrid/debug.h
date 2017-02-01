#ifndef DEBUG_H
#define DEBUG_H

#include <string>
#include <sstream>

#ifdef DEBUG
	#define DEBUG_ASSERT(exp, msg) \
		if (!(exp)) \
			dassert::assert(__FILE__, __LINE__, dassert::toString(msg).c_str(), #exp);
#else
	#define DEBUG_ASSERT(exp, msg)
#endif

namespace dassert
{
/**
 * \brief Show error message and stop operation
 */
void assert(const char *file, const int line, const char *msg, const char *exp);

/**
 * Try to convert any object to string
 */
template <typename T> std::string toString(const T& x)
{
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

} // namespace dassert

#endif
