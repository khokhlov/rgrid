#ifndef RG_IO_H
#define RG_IO_H

#include <iostream>

#include "rgrid/types.h"
#include "rgrid/darray.h"
#include "rgrid/darraycontainer.h"

namespace rgrid {

namespace rgio {

template <typename T, typename I> 
void saveData(std::basic_iostream<char>& stream, const DArray<T, I>& dArray, format fmt);
template <typename T, typename I>
void saveData(std::basic_iostream<char>& stream, const DArrayContainer<T, I>& dArrayContainer, format fmt);

template <typename T, typename I>
void loadData(std::basic_iostream<char>& stream, DArray<T, I>& dArray);

template <typename T, typename I>
static void writeHeader(std::basic_iostream<char>& stream, const DArray<T, I>& dArray);
template <typename T, typename I>
static void writeHeader(std::basic_iostream<char>& stream, const DArrayContainer<T, I>& dArrayContainer);

template <typename T, typename I>
static void loadHeader(std::basic_iostream<char>& stream, DArray<T, I>& dArray);

} // namespace rgio

} // namespace rgrid

#endif // RG_IO_H
