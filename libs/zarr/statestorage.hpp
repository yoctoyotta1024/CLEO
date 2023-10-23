/*
 * ----- CLEO -----
 * File: statestorage.hpp
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
 * metadata written specifically for variables
 * in the state of a gridbox
 */

#ifndef STATESTORAGE_HPP 
#define STATESTORAGE_HPP 

#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include <cassert>
#include <stdexcept>

#include "./fsstore.hpp"
#include "./storehelpers.hpp"

template <typename T>
struct StateStorage
/* 2D storage with dimensions [time, gbxindex] for
variables in the state of each gridbox over time.
nobs is number of observation events (no. time outputs)
and ngbxs is the number of elements in 1st dimension
of 2-D data i.e. no. gridboxes observed for each time */
{
private:
  const size_t chunksize;           // fixed size of array chunks (=max no. datapoints in buffer before writing)
  const size_t ngbxs;               // number elements in 1st dimensin (e.g. number of gridboxes that are observed)

  void writechunk()
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    const std::string chunknum = std::to_string(this->chunkcount) + ".0";

    storehelpers::writebuffer2chunk(this->store, this->buffers.mom0,
                                    get_name("0"), chunknum,
                                    this->chunkcount);

    storehelpers::writebuffer2chunk(this->store, this->buffers.mom1,
                                    get_name("1"), chunknum,
                                    this->chunkcount);

    std::tie(this->chunkcount, this->buffersfill) = storehelpers::
        writebuffer2chunk(this->store, this->buffers.mom2,
                          get_name("2"), chunknum,
                          this->chunkcount);

    writejsons();
  }
};

#endif //  STATESTORAGE_HPP  