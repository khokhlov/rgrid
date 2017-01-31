#include "rgrid/vtksaver.h"

namespace rgrid {

template <> std::string getVTKType<unsigned char>() { return "unsigned_char"; }
template <> std::string getVTKType<char>() { return "char"; }
template <> std::string getVTKType<unsigned short>() { return "unsigned_short"; }
template <> std::string getVTKType<short>() { return "short"; }
template <> std::string getVTKType<unsigned int>() { return "unsigned_int"; }
template <> std::string getVTKType<int>() { return "int"; }
template <> std::string getVTKType<unsigned long>() { return "unsigned_long"; }
template <> std::string getVTKType<long>() { return "long"; }
template <> std::string getVTKType<float>() { return "float"; } 
template <> std::string getVTKType<double>() { return "double"; }

} // namespace rgrid