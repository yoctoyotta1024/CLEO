/*
 * ----- CLEO -----
 * File: onedstorage.hpp
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
 * templated class for writing a single
 * 1-D variable via inherited features from
 * SingleVarStorage
 */

#ifndef ONEDSTORAGE_HPP
#define ONEDSTORAGE_HPP

#include <tuple>
#include <string>

#include "./storehelpers.hpp"
#include "./singlevarstorage.hpp"

template <typename T>
struct OneDStorage : SingleVarStorage<T>
/* storage of a 1-D variable with
dimensions 'time' in .zattrs metadata,
ie variable is one value output at
each time */
{
private:
  void writechunk()
  /* write data in buffer to a chunk in store */
  {
    std::tie(this->chunkcount, this->bufferfill) =
        storehelpers::
            writebuffer2chunk(this->store, this->buffer,
                              this->name, this->chunkcount);

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST) */
  {
    const auto shape("[" + std::to_string(this->ndata) + "]");
    const auto chunks("[" + std::to_string(this->chunksize) + "]");
    const std::string dims = "[\"time\"]";

    this->zarrayjsons(shape, chunks, dims);
  }

public:
  OneDStorage(FSStore &store, const unsigned int chunksize,
              const std::string name, const std::string dtype,
              const std::string units, const double scale_factor)
      : SingleVarStorage<T>(store, chunksize, name, dtype,
                            units, scale_factor) {}

  ~OneDStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      writechunk();
    }
  }
};

#endif // ONEDSTORAGE_HPP