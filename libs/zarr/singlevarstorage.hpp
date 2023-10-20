/*
 * ----- CLEO -----
 * File: singlevarstorage.hpp
 * Project: zarr
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
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
private:
  virtual void writechunk() = 0;
  virtual void writejsons() = 0;

  size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)

protected:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2
  const std::string name;    // name to call variable being stored
  const std::string units;   // units of coordinate being stored (for arrayattrs json)
  const double scale_factor; // scale_factor of data (for array .zattrs json)
  std::vector<T> buffer;     // buffer to store values in until writing to array chunk

  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed

  const char zarr_format = '2';          // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype;               // datatype stored in arrays

  unsigned int get_chunksize() const { return chunksize; }

  void set_buffer_chunksize(const unsigned int i_chunksize)
  {
    if (chunksize != storehelpers::NOTSETCHUNKSIZE)
    {
      const std::string err("chunksize already set; it cannot be changed");
      throw std::invalid_argument(err);
    }
    chunksize = i_chunksize;
    buffer.assign(chunksize, std::numeric_limits<T>::max());
  }

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

  void copy2buffer(const T val)
  /* copy value 'val' to buffer */
  {
    bufferfill = storehelpers::val2buffer<T>(val, buffer, bufferfill);
    ++ndata;
  }

  void copy2buffer(const std::vector<T> &vec)
  /* copy values of type T in vector 'vec' to buffer */
  {
    bufferfill = storehelpers::vec2buffer<T>(vec, buffer, bufferfill);
    ndata += vec.size();
  }

public:
  SingleVarStorage(FSStore &store, const unsigned int chunksize,
                   const std::string name, const std::string dtype,
                   const std::string units, const double scale_factor)
      : chunksize(chunksize), store(store),
        name(name), units(units), scale_factor(scale_factor),
        buffer(chunksize, std::numeric_limits<T>::max()),
        chunkcount(0), bufferfill(0), ndata(0), dtype(dtype) {}

  virtual ~SingleVarStorage(){};

  int get_ndata() const { return ndata; }

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