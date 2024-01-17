/*
 * ----- CLEO -----
 * File: storehelpers.hpp
 * Project: zarr
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * some helper functions for writng to buffers and writing
 * buffers to arrays in stores that obey the
 * zarr storage specification version 2 (e.g. see FSStore)
 * https://zarr.readthedocs.io/en/stable/spec/v2.html
 */

#ifndef STOREHELPERS_HPP
#define STOREHELPERS_HPP

#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <limits>
#include <sstream>

#include "./fsstore.hpp"

namespace storehelpers
/* namespace for generic helper functions used to
write a double to a buffer, a buffer to a chunk of an
array in a store, and an array's metadata to a store */
{
  inline unsigned int good2Dchunk(const unsigned int maxchunk,
                                  const size_t ndim1)
  /* given max chunksize, returns the (largest)
  suitable chunksize such that chunks are always an
  integer multiple of ndim1 (which should be the
  length of the 2nd dimension of 2D data) */
  {
    return std::floor(maxchunk / ndim1) * ndim1;
  }

  template <typename T>
  inline std::pair<unsigned int, unsigned int>
  val2buffer(const T val, std::vector<T> &buffer,
             const unsigned int ndata, const unsigned int j)
  /* copy a type T (e.g. a double) called
  'val', to buffer at index j */
  {
    buffer.at(j) = val;

    return std::pair(ndata + 1, j + 1); // updated {ndata, buffersfill}
  }

  template <typename T>
  inline std::pair<unsigned int, unsigned int>
  vec2buffer(const std::vector<T> &vec, std::vector<T> &buffer,
             const unsigned int ndata, const unsigned int j)
  /* copy vector of type T (e.g. a double) called
  'vec', to buffer at index j. Function is equivalent to
  std::copy(vec.begin(), vec.end(), buffer.begin()+j);
  but faster for copying a large vector (not iterative) */
  {
    size_t nvalues(vec.size());

    buffer.erase(buffer.end() - nvalues, buffer.end());
    buffer.insert(buffer.begin() + j, vec.begin(), vec.end());

    return std::pair(ndata + nvalues, j + nvalues);
  }

  template <typename T>
  inline std::pair<unsigned int, unsigned int>
  writebuffer2chunk(FSStore &store,
                    std::vector<T> &buffer,
                    const std::string &name,
                    const std::string &chunknum,
                    const unsigned int chunkcount)
  /* write buffer vector into attr's store at chunk no. 'kk', then
  replace contents of buffer with max numeric limit of type.
  Return incremented value of chunkcount */
  {
    store[name + "/" + chunknum].operator= <T>(buffer);
    buffer.assign(buffer.size(), std::numeric_limits<T>::max());

    return std::pair(chunkcount + 1, 0); // updated {chunkcount, bufferfill}
  }

  template <typename T>
  inline std::pair<unsigned int, unsigned int>
  writebuffer2chunk(FSStore &store,
                    std::vector<T> &buffer,
                    const std::string &name,
                    const unsigned int chunkcount)
  /* write buffer vector into attr's store at 'chunkcount' and then
  return incremented chunkcount */
  {
    const std::string chunknum = std::to_string(chunkcount);

    return storehelpers::writebuffer2chunk(store, buffer, name,
                                           chunknum, chunkcount);
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

    return storehelpers::metadata(zarr_format, order, shape, chunks,
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

  inline void writezattrsjson(FSStore &store,
                              const std::string &name,
                              const std::string &dims,
                              const std::string &units = " ",
                              const double scale_factor = 1)
  /* write .zattr json file into store for an array
  of a variable called 'name' */
  {
    // define dimension names of this array, to make xarray and netCDF
    // happy (not a MUST, ie. not strictly required, by zarr)
    // e.g. "{\"_ARRAY_DIMENSIONS\": [\"x\"]}";
    store[name + "/.zattrs"] = arrayattrs(dims, units, scale_factor);
  }

  inline void writejsons(FSStore &store,
                         const std::string &name,
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

  inline void writejsons(FSStore &store,
                         const std::string &name,
                         const std::string &metadata,
                         const std::string &dims,
                         const std::string &units,
                         const double scale_factor)
  /* make arrayattrs then write it and array's
  metadata to .json files */
  {
    const std::string arrayattrs(
        storehelpers::arrayattrs(dims, units, scale_factor));

    storehelpers::writejsons(store, name, metadata, arrayattrs);
  }
};

#endif // STOREHELPERS_HPP
