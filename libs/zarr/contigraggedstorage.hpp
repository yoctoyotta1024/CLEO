/*
 * ----- CLEO -----
 * File: contigraggedstorage.hpp
 * Project: zarr
 * Created Date: Monday 23rd October 2023
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
 * File for Contiguous Ragged Array Storage
 * used to store superdroplet attributes
 * (see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
 * in a FFStore obeying zarr storage specification verion 2:
 * https://zarr.readthedocs.io/en/stable/spec/v2.html */


#ifndef CONTIGRAGGEDSTORAGE
#define CONTIGRAGGEDSTORAGE 

#include <concepts>
#include <vector>
#include <string>
#include <utility>
#include <tuple>

#include "./fsstore.hpp"
#include "./storehelpers.hpp"
#include "./superdropsbuffers.hpp"
#include "../cleoconstants.hpp"
#include "superdrops/superdrop.hpp"

template <SuperdropsBuffers Buffers>
class ContigRaggedStorage
/* Class for outputting Superdrop's data into zarr storage in
arrays of contigous ragged representation with 'chunkcount' number
of chunks that have a fixed chunksize. Works by filling buffers in
sdbuffers with superdrop data and then writing these buffers
into chunks in their corresponding array stores when number of
datapoints copied to the buffers reaches chunksize. */
{
private:
  FSStore &store;         // file system store satisfying zarr store specificaiton v2
  const size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)

  std::vector<size_t> rgdcount;                      // count variable for contiguous ragged representation of arrays
  unsigned int rgdcount_chunkcount;                  // number of chunks of rgdcount array so far written to store
  unsigned int rgdcount_bufferfill;                  // number of rgdcount values so far copied into its buffer
  unsigned int rgdcount_ndata;                       // number of rgdcount values observed so far
  const std::string rgdcount_name = "rgdtotnsupers"; // name of raggedcount zarray in store
  const std::string rgdcount_dtype = "<u8";          // datatype of raggedcount variable

  Buffers buffers;         // buffers and their handler functions for wrting SD data to store
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int buffersfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed (= size of array written to store)

  const char zarr_format = '2';          // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression

  template <typename T>
  void copy2sdbuffers(const T &value)
  /* copy data from superdrop to buffer(s) and
  increment required counting variables */
  {
    std::tie(ndata, buffersfill) =
        buffers.copy2buffer(value, ndata, buffersfill);
  }

  void copy2rgdcount(const size_t raggedn)
  /* write raggedn into rgdcount buffer */
  {
    std::tie(rgdcount_ndata, rgdcount_bufferfill) =
        storehelpers::val2buffer<size_t>(raggedn, rgdcount,
                                         rgdcount_ndata,
                                         rgdcount_bufferfill);
  }

public:
  ContiguousRaggedSDStorage(FSStore &store,
                            const Buffers buffers,
                            const size_t maxchunk)
      : store(store), chunksize(maxchunk), rgdcount(maxchunk),
        rgdcount_chunkcount(0), rgdcount_bufferfill(0),
        rgdcount_ndata(0), buffers(buffers),
        chunkcount(0), buffersfill(0), ndata(0)
  {
    buffers.set_buffersize(chunksize);
  }

  ~ContiguousRaggedSDStorage()
  {
    if (buffersfill != 0)
    {
      sdbuffers_writechunk();
    }

    if (rgdcount_bufferfill != 0)
    {
      rgdcount_writechunk();
    }
  }

  template <typename T>
  void data_to_raggedstorage(const T &value)
  /* write 'value' in contiguous ragged representation of an array
  in the zarr store. First copy data to buffer(s), then write buffer(s)
  to chunks in the store when the number of datapoints they contain
  reaches the chunksize */
  {
    if (buffersfill == chunksize)
    {
      buffers_writechunk();
    }

    copy2buffers(value);
  }

  void raggedarray_count(const size_t raggedn)
  /* add element 'raggedn' to rgdcount. 'raggedn' should be
  number of datapoints written to sdbuffer(s) during one event.
  rgdcount is then count variable for contiguous ragged
  representation of arrays written to store via sdbuffer(s). */
  {
    if (rgdcount_bufferfill == chunksize)
    {
      rgdcount_writechunk();
    }

    copy2rgdcount(raggedn);
  }
}

#endif // CONTIGRAGGEDSTORAGE