/*
 * ----- CLEO -----
 * File: singlevarstorage.hpp
 * Project: zarr
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Incomplete templated class for writing a single variable
 * (could be 1-D or 2-D variable or a coordinate etc.)
 * via a buffer into chunks of arrays in a fsstore obeying
 * the Zarr Storage specification version 2
 * */

#ifndef SINGLEVARSTORAGE_HPP
#define SINGLEVARSTORAGE_HPP

#include <limits>
#include <vector>
#include <string>
#include <stdexcept>

#include "./fsstore.hpp"
#include "./storehelpers.hpp"

template <typename T>
class SingleVarStorage
{
protected:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2
  std::vector<T> buffer;     // buffer to store values in until writing to array chunk
  const std::string name;    // name to call variable being stored
  const std::string units;   // units of coordinate being stored (for arrayattrs json)
  const double scale_factor; // scale_factor of data (for array .zattrs json)

  const size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed
 
  void zarrayjsons(const std::string shape,
                   const std::string chunks,
                   const std::string dims)
  /* write array's metadata to .json files */
  {
    const std::string metadata = storehelpers::
        metadata(zarr_format, order, shape, chunks, dtype,
                 compressor, fill_value, filters);

    const std::string arrayattrs = storehelpers::
        arrayattrs(dims, units, scale_factor);

    storehelpers::writezarrjsons(store, name, metadata, arrayattrs);
  }

private:
  virtual void writechunk() = 0;
  virtual void writejsons() = 0;

  const char zarr_format = '2';          // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype;               // datatype stored in arrays

  void copy2buffer(const T val)
  /* copy value 'val' to buffer */
  {
    std::tie(ndata, bufferfill) =
        storehelpers::val2buffer<T>(val, buffer, ndata, bufferfill);
  }

  void copy2buffer(const std::vector<T> &vec)
  /* copy values of type T in vector 'vec' to buffer */
  {
    std::tie(ndata, bufferfill) =
        storehelpers::vec2buffer<T>(vec, buffer, ndata, bufferfill);
  }

public:
  SingleVarStorage(FSStore &store, const unsigned int maxchunk,
                   const std::string name, const std::string dtype,
                   const std::string units, const double scale_factor)
      : store(store), buffer(maxchunk, std::numeric_limits<T>::max()),
        name(name), units(units), scale_factor(scale_factor),
        chunksize(maxchunk), chunkcount(0), bufferfill(0),
        ndata(0), dtype(dtype) {}

  virtual ~SingleVarStorage(){};

  unsigned int get_ndata() const { return ndata; }

  void is_name(const std::string &goodname) const
  {
    if (name != goodname)
    {
      const std::string errmsg("name of storage is " + name +
                               ", but should be " + goodname);
      throw std::invalid_argument(errmsg);
    }
  }

  void value_to_storage(const T val)
  /* write val in the zarr store. First copy it to a buffer,
  then write buffer to a chunk in the store when the number
  of values in the buffer reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      writechunk();
    }

    copy2buffer(val);
  }

  void value_to_storage(const std::vector<T> &vec)
  /* write 'vec' vector of type T in the zarr store.
  First copy vector to a buffer, then write buffer to a
  chunk in the store when the number of values in
  the buffer reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      writechunk();
    }

    copy2buffer(vec);
  }
};

#endif // SINGLEVARSTORAGE_HPP