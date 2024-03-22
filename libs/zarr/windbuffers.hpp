/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: windbuffers.hpp
 * Project: zarr
 * Created Date: Friday 22nd March 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 22nd March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Storage similar to TwoDstorage for writing variables to a fsstore via buffers
 * and occording to the Zarr storage specification version 2.0, but extended to
 * more than one variable for the wvel, vvel and uvel winds in each gridbox with
 * metadata written specifically for these variables
 */

#ifndef LIBS_ZARR_WINDBUFFERS_HPP_
#define LIBS_ZARR_WINDBUFFERS_HPP_

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
struct WindBuffers {
 private:
  std::vector<T> wvel;
  std::vector<T> uvel;
  std::vector<T> vvel;

 public:
  WindBuffers(const std::string endname, const unsigned int chunksize)
      : wvel(chunksize, std::numeric_limits<T>::max()),
        vvel(chunksize, std::numeric_limits<T>::max()),
        uvel(chunksize, std::numeric_limits<T>::max()) {}

  /* copy value to mass moments to their respective buffers */
  std::pair<unsigned int, unsigned int> copy2buffer(const State &state, const unsigned int ndata,
                                                    const unsigned int buffersfill) {
    storehelpers::val2buffer<T>(state.wvelcentre(), wvel, ndata, buffersfill);
    storehelpers::val2buffer<T>(state.uvelcentre(), uvel, ndata, buffersfill);
    storehelpers::val2buffer<T>(state.vvelcentre(), vvel, ndata, buffersfill);

    return std::pair(ndata + 1, buffersfill + 1);  // updated {ndata, buffersfill}
  }

  /* write data in buffer to a chunk in store alongside metadata jsons */
  std::pair<unsigned int, unsigned int> writechunks(FSStore &store, const unsigned int chunkcount) {
    const std::string chunknum = std::to_string(chunkcount) + ".0";

    storehelpers::writebuffer2chunk(store, wvel, "wvel", chunknum, chunkcount);
    storehelpers::writebuffer2chunk(store, uvel, "uvel", chunknum, chunkcount);
    storehelpers::writebuffer2chunk(store, vvel, "vvel", chunknum, chunkcount);

    return std::pair(chunkcount + 1, 0);  // updated {chunkcount, bufferfill}
  }

  /* write array's metadata to .json files */
  void writejsons(FSStore &store, const std::string &metadata) const {
    const std::string dims = "[\"time\", \"gbxindex\"]";

    storehelpers::writejsons(store, "wvel", metadata, dims, "m/s", dlc::W0);
    storehelpers::writejsons(store, "wvel", metadata, dims, "m/s", dlc::W0);
    storehelpers::writejsons(store, "wvel", metadata, dims, "m/s", dlc::W0);
  }
};

#endif  // LIBS_ZARR_WINDBUFFERS_HPP_
