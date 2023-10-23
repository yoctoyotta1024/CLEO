/*
 * ----- CLEO -----
 * File: massmomentsstorage.hpp
 * Project: zarr
 * Created Date: Sunday 22nd October 2023
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
 * Storage similar to twoDstorage for writing
 * variables to a fsstore via buffers and occording
 * to the Zarr storage specification version 2.0,
 * but extended to more than one variable and with
 * metadata written specifically for the 0th, 1st
 * and 2nd moments of the (real) droplet mass
 * distribution
 */

#ifndef MASSMOMENTSSTORAGE_HPP
#define MASSMOMENTSSTORAGE_HPP

#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include <cassert>
#include <stdexcept>

#include "./fsstore.hpp"
#include "./storehelpers.hpp"
#include "../cleoconstants.hpp"

template <typename T>
struct MassMomentsBuffers
{
private:
  const std::string endname; // name to add to end of massmom[X] being stored
  std::vector<T> mom0;       // buffer for 0th mass moment data until writing to array chunk
  std::vector<T> mom1;       // buffer for 1st mass moment data until writing to array chunk
  std::vector<T> mom2;       // buffer for 2nd mass moment data until writing to array chunk

  std::string get_name(const std::string mom) const
  {
    return "massmom" + mom + endname;
  }

  void writezarrjsons(FSStore &store,
                      const std::string &name,
                      const std::string &metadata,
                      const std::string &dims,
                      const std::string &units,
                      const double scale_factor) const
  /* write array's metadata to .json files */
  {
    const std::string arrayattrs = storehelpers::
        arrayattrs(dims, units, scale_factor);

    storehelpers::writezarrjsons(store, name, metadata, arrayattrs);
  }

public:
  MassMomentsBuffers(const std::string endname,
                     const unsigned int chunksize)
      : endname(endname),
        mom0(chunksize, std::numeric_limits<T>::max()),
        mom1(chunksize, std::numeric_limits<T>::max()),
        mom2(chunksize, std::numeric_limits<T>::max()) {}

  std::pair<unsigned int, unsigned int>
  copy2buffer(const std::array<T, 3> moms,
              const unsigned int ndata,
              const unsigned int buffersfill)
  /* copy value to mass moments to their respective buffers */
  {
    storehelpers::val2buffer<T>(moms.at(0), mom0, ndata, buffersfill);
    storehelpers::val2buffer<T>(moms.at(1), mom1, ndata, buffersfill);
    storehelpers::val2buffer<T>(moms.at(2), mom2, ndata, buffersfill);

    return std::pair(ndata + 1, buffersfill + 1); // updated {ndata, buffersfill}
  }

  std::pair<unsigned int, unsigned int>
  writechunks(FSStore &store, const unsigned int chunkcount)
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    const std::string chunknum = std::to_string(chunkcount) + ".0";

    storehelpers::writebuffer2chunk(store, mom0, get_name("0"),
                                    chunknum, chunkcount);

    storehelpers::writebuffer2chunk(store, mom1, get_name("1"),
                                    chunknum, chunkcount);

    storehelpers::writebuffer2chunk(store, mom2, get_name("2"),
                                    chunknum, chunkcount);

    return std::pair(chunkcount + 1, 0); // updated {chunkcount, bufferfill}
  }

  void writejsons(FSStore &store,
                  const std::string &metadata) const
  /* write array's metadata to .json files */
  {
    const std::string dims = "[\"time\", \"gbxindex\"]";

    const std::string units0 = " ";
    constexpr double scale_factor0 = 1.0;
    writezarrjsons(store, get_name("0"), metadata, dims,
                   units0, scale_factor0);

    const std::string units1 = "g";
    constexpr double scale_factor1 = dlc::MASS0grams; // grams
    writezarrjsons(store, get_name("1"), metadata, dims,
                   units1, scale_factor1);

    const std::string units2 = "g^2";
    constexpr double scale_factor2 = dlc::MASS0grams * dlc::MASS0grams; // grams squared
    writezarrjsons(store, get_name("2"), metadata, dims,
                   units2, scale_factor2);
  }
};

template <typename T>
struct MassMomentsStorage
/* 2D storage with dimensions [time, gbxindex] for 0th, 1st and 
2nd momments of the (real) droplet mass distirbution.
nobs is number of observation events (no. time outputs)
and ngbxs is the number of elements in 1st dimension
of 2-D data i.e. no. gridboxes observed for each time */
{
private:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2

  const size_t chunksize;           // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount;  // number of chunks of array so far written to store
  unsigned int buffersfill; // number of datapoints so far copied into buffer
  unsigned int ndata;       // number of data points that have been observed

  const char zarr_format = '2';          // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype;               // datatype stored in arrays

  MassMomentsBuffers<T> buffers;
  const size_t ngbxs;               // number elements in 1st dimensin (e.g. number of gridboxes that are observed)

  void writejsons() const
  /* write strictly required metadata to decode chunks (MUST).
  Assert also check 2D data dimensions is as expected */
  {
    assert((ndata == nobs * ngbxs) &&
           "1D data length matches 2D array size");
    assert((chunksize % ngbxs == 0.0) &&
           "chunks are integer multiple of 1st dimension of 2-D data");

    const auto n1str = std::to_string(ngbxs);
    const auto nobstr = std::to_string(nobs);
    const auto nchstr = std::to_string(chunksize / ngbxs);
    
    const auto shape("[" + nobstr + ", " + n1str + "]");
    const auto chunks("[" + nchstr + ", " + n1str + "]");
    
    const std::string metadata = storehelpers::
        metadata(zarr_format, order, shape, chunks, dtype,
                 compressor, fill_value, filters);

    buffers.writejsons(store, metadata);
  }

  void writechunks()
  /* write data from buffers into chunks in store,
  then reset buffersfill and write associated metadata */
  {
    std::tie(chunkcount, buffersfill) =
        buffers.writechunks(store, chunkcount);
      
    writejsons();
  }

  void copy2buffers(const std::array<T, 3> moms)
  /* copy data from thermostate to buffers */
  {
    std::tie(ndata, buffersfill) =
        buffers.copy2buffer(moms, ndata, buffersfill);
  }

public:
  unsigned int nobs; // number of output times that have been observed

  MassMomentsStorage(FSStore &store, const unsigned int maxchunk,
                     const std::string dtype, const size_t ngbxs,
                     const std::string endname)
      : store(store),
        chunksize(storehelpers::good2Dchunk(maxchunk, ngbxs)),
        chunkcount(0), buffersfill(0), ndata(0), dtype(dtype), 
        buffers(endname, chunksize), ngbxs(ngbxs), nobs(0) {}

  ~MassMomentsStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (buffersfill != 0)
    {
      writechunks();
    }
  }

  void massmoments_to_storage(const T mom0, const T mom1, const T mom2)
  /* write val in the zarr store. First copy it to a buffer,
  then write buffer to a chunk in the store when the number
  of values in the buffer reaches the chunksize */
  {
    if (buffersfill == chunksize)
    {
      writechunks();
    }

    copy2buffers({mom0, mom1, mom2});
  }
};

#endif // MASSMOMENTSSTORAGE_HPP  