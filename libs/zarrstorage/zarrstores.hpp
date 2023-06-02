// Author: Clara Bayley and Tobias Kölling
// File: zarr_stores.hpp
/* objects that can be used as stores obyeying the
zarr storage specification version 2 (e.g. see FSStore)
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef ZARRSTORES_HPP
#define ZARRSTORES_HPP

#include <string>
#include <string_view>
#include <span>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "sdmgridboxes/gridbox.hpp"
#include "superdrop_solver/superdrop.hpp"

template <typename Store>
struct StoreAccessor
/* functions for converting types (e.g. vectors of
unsigned integers or doubles) into vectors of single bytes to
write to store under a given key. Store can be anything that
satisfies the zarr storage specifcaiton version 2 */
{
  Store &store;
  std::string_view key;

  StoreAccessor &operator=(std::span<const uint8_t> buffer)
  /* write range of memory representing uint8_ts to store */
  {
    store.write(key, buffer);
    return *this;
  }

  StoreAccessor &operator=(std::string_view buffer)
  /* reinterpret range of memory representing string as
  a range of memory representing uint8_ts, then write to store */
  {
    return operator=(std::span<const uint8_t>(
        reinterpret_cast<const uint8_t *>(buffer.data()),
        buffer.size()));
  }

  template <typename T>
  StoreAccessor &operator=(std::span<const T> buffer)
  /* re-interpret range of memory representing vector of type T as
  a range of memory representing uint8_ts, then write to store */
  {
    return operator=(std::span<const uint8_t>(
        reinterpret_cast<const uint8_t *>(buffer.data()),
        buffer.size() * sizeof(T)));
  }
};

class FSStore
/* A file system (with root in 'basedir' directory) obeying Zarr
version 2 requirements for a Store. Store contins a series
of key, values where values may be data arrays or groups in the store.
data for a given key is written to the store via the functions
in StoreAccessor */
{
private:
  const std::filesystem::path basedir;

public:
  FSStore(std::filesystem::path basedir) : basedir(basedir)
  {
    // initialize a zarr group (i.e. dataset)
    const std::string zarr_format("2"); // storage spec. version 2
    const std::string zgroupjson("{\"zarr_format\": " +
                                   zarr_format + "}");
    (*this)[".zgroup"] = zgroupjson;

    // global metadata (optional)
    (*this)[".zattrs"] = "{\"creator\": \"Clara Bayley\", "
                         "\"title\": \"store for output of coupled SDM\"}";
  }

  StoreAccessor<FSStore> operator[](std::string_view key)
  {
    return {*this, key};
  }

  bool write(std::string_view key, std::span<const uint8_t> buffer);
  /* write function called by StoreAccessor once data has been
  converted into a vector of unsigned integer types */
};

namespace storagehelper
/* namespace for generic helper functions used to
write a double to a buffer, a buffer to a chunk of an
array in a store, and an array's metadata to a store */
{

  template <typename T>
  inline unsigned int val2buffer(const T val,
                                 std::vector<T> &buffer,
                                 unsigned int j)
  /* copy a type T (e.g. a double) called
  'val', to buffer at index j */
  {
    buffer.at(j) = val;

    return ++j;
  }

  template <typename T>
  inline unsigned int vec2buffer(const std::vector<T> &vec,
                                 std::vector<T> &buffer,
                                 unsigned int j)
  /* copy vector of type T (e.g. a double) called
  'vec', to buffer at index j. Function is equivalent to 
  std::copy(vec.begin(), vec.end(), buffer.begin()+j);
  but faster for copying a large vector (not iterative) */
  {
    size_t nvalues(vec.size());
 
    buffer.erase(buffer.end() - nvalues, buffer.end());
    buffer.insert(buffer.begin()+j, vec.begin(), vec.end());

    // return j + nvalues;
    return 0; 
  }

  template <typename T>
  std::pair<unsigned int, unsigned int>
  writebuffer2chunk(FSStore &store,
                    std::vector<T> &buffer,
                    const std::string name,
                    const std::string chunknum,
                    unsigned int chunkcount)
  /* write buffer vector into attr's store at chunk no. 'kk', then
  replace contents of buffer with max numeric limit of type.
  Return incremented value of chunkcount */
  {
    store[name + "/" + chunknum].operator=<T>(buffer);
    buffer.assign(buffer.size(), std::numeric_limits<T>::max());

    return std::pair(++chunkcount, 0);
  }

  template <typename T>
  std::pair<unsigned int, unsigned int>
  writebuffer2chunk(FSStore &store,
                    std::vector<T> &buffer,
                    const std::string name,
                    unsigned int chunkcount)
  /* write buffer vector into attr's store at 'chunkcount' and then
  return incremented chunkcount */
  {
    const std::string chunknum = std::to_string(chunkcount);

    return storagehelper::writebuffer2chunk(store, buffer, name,
                                            chunknum, chunkcount);
  }

  inline void writezarrjsons(FSStore &store,
                             const std::string name,
                             const std::string &metadata,
                             const std::string &arrayattrs)
  /* write .zarray and .zattr json files into store for the
  metadata of an array of a variable called 'name' */
  {
    // strictly required metadata to decode chunks (MUST)
    store[name + "/.zarray"] = metadata;

    // define dimension names of this array, to make xarray and netCDF happy
    // (not a MUST, ie. not strictly required, by zarr)
    // e.g. "{\"_ARRAY_DIMENSIONS\": [\"x\"]}";
    store[name + "/.zattrs"] = arrayattrs;
  }

  inline std::string metadata(const char zarr_format,
                              const char order,
                              const std::string &shape,
                              const std::string &chunks,
                              const std::string &dtype,
                              const std::string &compressor,
                              const std::string &fill_value,
                              const std::string &filters)
  /* make string of metadata for an array in a zarr store */
  {
    const std::string metadata("{"
                               "\"shape\": " +
                               shape + ", "
                                       "\"chunks\": " +
                               chunks + ", "
                                        "\"dtype\": \"" +
                               dtype + "\", "
                                       "\"order\": \"" +
                               order + "\", "
                                       "\"compressor\": " +
                               compressor + ", "
                                            "\"fill_value\": " +
                               fill_value + ", "
                                            "\"filters\": " +
                               filters + ", "
                                         "\"zarr_format\": " +
                               zarr_format + "}");
    return metadata;
  }

  inline std::string metadata(const char zarr_format,
                              const char order,
                              const unsigned int ndata,
                              const size_t chunksize,
                              const std::string &dtype,
                              const std::string &compressor,
                              const std::string &fill_value,
                              const std::string &filters)
  /* make string of metadata for an array in a zarr store */
  {
    const auto shape("[" + std::to_string(ndata) + "]");
    const auto chunks("[" + std::to_string(chunksize) + "]");

    return storagehelper::metadata(zarr_format, order, shape, chunks,
                                   dtype, compressor, fill_value, filters);
  }

  inline std::string arrayattrs(const std::string &dims,
                                const std::string &units = " ",
                                const double scale_factor = 1)
  /* make string of zattrs attribute information for an array in a zarr store */
  {
    std::ostringstream sfstr;
    sfstr << std::scientific << scale_factor;

    const std::string arrayattrs = "{\"_ARRAY_DIMENSIONS\": " + dims +
                                   ", \"units\": " + "\"" + units +
                                   "\", \"scale_factor\": " +
                                   sfstr.str() + "}";
    return arrayattrs;
  }
};

#endif // ZARRSTORES_HPP