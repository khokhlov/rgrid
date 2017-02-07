#ifndef DEBUG_H
#define DEBUG_H

#include <string>
#include <sstream>

// static assert

namespace dassert {

template <bool b> struct CompileTimeError;
template <> struct CompileTimeError<false> {};

} // namespace dassert

#define SA_SUB_CAT(X,Y) X##Y
#define SA_CAT(X,Y) SA_SUB_CAT(X,Y)

#define STATIC_ASSERT(expr, msg) { \
	dassert::CompileTimeError<!(expr)> SA_CAT(ERROR__##msg##__LINE_,__LINE__); \
	(void) SA_CAT(ERROR__##msg##__LINE_,__LINE__); }

// debug assert

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
