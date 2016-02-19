#include "rgrid/rgio.h"

using namespace rgrid;

template <typename T, typename I> 
void saveData(std::basic_iostream<char>& stream, const DArray<T, I>& dArray, rgio::format fmt) {
	for (I z = 0; z != dArray.localSize(Z); ++z) {
		for (I y = 0; y != dArray.localSize(Y); ++y) {
			dArray.writeLine(stream, y, z, fmt);
		}
	}
}

template <typename T, typename I> 
void saveData(std::basic_iostream<char>& stream, const DArrayContainer<T, I>& dArrayContainer, rgio::format fmt) {
	for (I z = 0; z != dArrayContainer.size(Z); ++z) {
		for (I y = 0; y != dArrayContainer.size(Y); ++y) {
			dArrayContainer.writeLine(stream, y, z, fmt);
		}
	}
}