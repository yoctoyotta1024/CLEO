// Author: Clara Bayley
// File: "observer_singlevariable.cpp"
/* Classes and helper functions in a namespace
useful for using storage clases with buffers to
write values of 1D data into chunks of arrays
in a zarr store */

#include "observer_singlevariable.hpp" 

template class SingleVarStorage<double>;
template class SingleVarStorage<unsigned int>;
template class SingleVarStorage<size_t>;

template <typename T>
SingleVarStorage<T>::SingleVarStorage(FSStore &store, const unsigned int maxcsize, const std::string name,
                              const std::string dtype, const std::string units,
                              const double scale_factor)
      : store(store), name(name), units(units), scale_factor(scale_factor),
        buffer(maxcsize), chunksize(maxcsize),
        chunkcount(0), bufferfill(0), ndata(0), dtype(dtype)
        {
          init_buffer(maxcsize);
        }

template <>
void SingleVarStorage<double>::init_buffer(const unsigned int maxcsize)
{
  buffer = std::vector<double>(maxcsize, std::nan(""));
}

template <>
void SingleVarStorage<unsigned int>::init_buffer(const unsigned int maxcsize)
{
  buffer = std::vector<unsigned int>(maxcsize, -1);
}

template <>
void SingleVarStorage<size_t>::init_buffer(const unsigned int maxcsize)
{
  buffer = std::vector<size_t>(maxcsize, -1);
}