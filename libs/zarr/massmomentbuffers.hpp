/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: massmomentbuffers.hpp
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
 * Buffers to use in TwoDMultiVarstorage for writing
 * 0th, 1st and 2nd mass moments to a fsstore via
 * buffers and according to the Zarr storage specification
 * version 2.0, with metadata written specifically for
 * these moments of the (real) droplet mass distribution
 */

#ifndef MASSMOMENTBUFFERS_HPP
#define MASSMOMENTBUFFERS_HPP

#include <string>
#include <vector>
#include <limits>
#include <array>
#include <tuple>
#include <utility>

#include "../cleoconstants.hpp"
#include "./fsstore.hpp"
#include "./storehelpers.hpp"

namespace dlc = dimless_constants;

template <typename T>
struct MassMomentBuffers
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

public:
  MassMomentBuffers(const std::string endname,
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
    storehelpers::writejsons(store, get_name("0"), metadata,
                             dims, units0, scale_factor0);

    const std::string units1 = "g";
    constexpr double scale_factor1 = dlc::MASS0grams; // grams
    storehelpers::writejsons(store, get_name("1"), metadata,
                             dims, units1, scale_factor1);

    const std::string units2 = "g^2";
    constexpr double scale_factor2 = dlc::MASS0grams * dlc::MASS0grams; // grams squared
    storehelpers::writejsons(store, get_name("2"), metadata,
                             dims, units2, scale_factor2);
  }
};

#endif // MASSMOMENTBUFFERS_HPP
