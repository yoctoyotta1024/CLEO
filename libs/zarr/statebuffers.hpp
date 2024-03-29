/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: statebuffers.hpp
 * Project: zarr
 * Created Date: Sunday 22nd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Storage similar to TwoDstorage for writing
 * variables to a fsstore via buffers and occording
 * to the Zarr storage specification version 2.0,
 * but extended to more than one variable and with
 * metadata written specifically for variables
 * in the state of a gridbox
 */

#ifndef LIBS_ZARR_STATEBUFFERS_HPP_
#define LIBS_ZARR_STATEBUFFERS_HPP_

#include <array>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../cleoconstants.hpp"
#include "superdrops/state.hpp"
#include "zarr/fsstore.hpp"
#include "zarr/storehelpers.hpp"

namespace dlc = dimless_constants;

template <typename T>
struct StateBuffers {
 private:
  std::vector<T> press;
  std::vector<T> temp;
  std::vector<T> qvap;
  std::vector<T> qcond;

 public:
  StateBuffers(const std::string endname, const unsigned int chunksize)
      : press(chunksize, std::numeric_limits<T>::max()),
        temp(chunksize, std::numeric_limits<T>::max()),
        qvap(chunksize, std::numeric_limits<T>::max()),
        qcond(chunksize, std::numeric_limits<T>::max()) {}

  /* copy value to mass moments to their respective buffers */
  std::pair<unsigned int, unsigned int> copy2buffer(const State &state, const unsigned int ndata,
                                                    const unsigned int buffersfill) {
    storehelpers::val2buffer<T>(state.press, press, ndata, buffersfill);
    storehelpers::val2buffer<T>(state.temp, temp, ndata, buffersfill);
    storehelpers::val2buffer<T>(state.qvap, qvap, ndata, buffersfill);
    storehelpers::val2buffer<T>(state.qcond, qcond, ndata, buffersfill);

    return std::pair(ndata + 1, buffersfill + 1);  // updated {ndata, buffersfill}
  }

  /* write data in buffer to a chunk in store alongside metadata jsons */
  std::pair<unsigned int, unsigned int> writechunks(FSStore &store, const unsigned int chunkcount) {
    const std::string chunknum = std::to_string(chunkcount) + ".0";

    storehelpers::writebuffer2chunk(store, press, "press", chunknum, chunkcount);

    storehelpers::writebuffer2chunk(store, temp, "temp", chunknum, chunkcount);

    storehelpers::writebuffer2chunk(store, qvap, "qvap", chunknum, chunkcount);

    storehelpers::writebuffer2chunk(store, qcond, "qcond", chunknum, chunkcount);

    return std::pair(chunkcount + 1, 0);  // updated {chunkcount, bufferfill}
  }

  /* write array's metadata to .json files */
  void writejsons(FSStore &store, const std::string &metadata) const {
    const std::string dims = "[\"time\", \"gbxindex\"]";

    storehelpers::writejsons(store, "press", metadata, dims, "hPa", dlc::P0 / 100);

    storehelpers::writejsons(store, "temp", metadata, dims, "K", dlc::TEMP0);

    storehelpers::writejsons(store, "qvap", metadata, dims, " ", 1.0);

    storehelpers::writejsons(store, "qcond", metadata, dims, " ", 1.0);
  }
};

#endif  // LIBS_ZARR_STATEBUFFERS_HPP_
