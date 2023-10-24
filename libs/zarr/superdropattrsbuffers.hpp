/*
 * ----- CLEO -----
 * File: superdropattrsbuffers.hpp
 * Project: zarr
 * Created Date: Monday 23rd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 24th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * structs obeying the superdropsbuffers concept in order to
 * write out superdrop attributes into a ragged array in a 
 * fsstore via a buffer
 */


#ifndef SUPERDROPATTRSBUFFERS_HPP 
#define SUPERDROPATTRSBUFFERS_HPP 

#include <concepts>
#include <string>
#include <vector>
#include <limits>
#include <utility>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "../cleoconstants.hpp"
#include "./fsstore.hpp"
#include "./storehelpers.hpp"
#include "./superdropsbuffers.hpp"
#include "superdrops/superdrop.hpp"

namespace dlc = dimless_constants;

template <typename T>
struct SuperdropAttrBuffer
/* generic struct satisfies the SuperdropletsBuffer concept.
Can be inherited and then used to put a single attribute called
'attr' into a buffer (given an implementation of the
'copy2buffer' function) and write the accompanying metadata
obeying the zarr storage specificatino version 2.0 via
'writechunk' and 'writemetadata' functions */
{
  const std::string attr;  // name of attribute in fsstore
  const std::string dtype; // datatype stored in arrays
  std::vector<T> buffer;   // buffer to fill before writing to store

  SuperdropAttrBuffer(const std::string attr,
                      const std::string dtype)
      : attr(attr), dtype(dtype), buffer(0) {}

  virtual ~SuperdropAttrBuffer(){};

  virtual std::pair<unsigned int, unsigned int> 
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j) = 0;
  /* virtual void function placeholding function for
  copying superdrop's data into a buffer vector at j'th index */

  std::pair<unsigned int, unsigned int>
  writechunk(FSStore &store, unsigned int chunkcount)
  /* write buffer vector into attr's store at chunkcount
  and then replace contents of buffer with numeric limit */
  {
    return storehelpers::writebuffer2chunk(store, buffer,
                                           attr, chunkcount);
  }

  void writejsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    const std::string metadata = storehelpers::
        metadata(md.zarr_format, md.order, md.shape,
                 md.chunks, dtype, md.compressor,
                 md.fill_value, md.filters);

    const std::string arrayattrs = "{\"_ARRAY_DIMENSIONS\": " + md.dims + "}";

    storehelpers::writezarrjsons(store, attr, metadata, arrayattrs);
  }

  void set_buffer(const size_t maxchunk)
  {
    buffer = std::vector<T>(maxchunk, std::numeric_limits<T>::max());
  }
};

struct SdIdBuffer : SuperdropAttrBuffer<size_t>
{
  SdIdBuffer() : SuperdropAttrBuffer("sdId", "<u8"){};

  std::pair<unsigned int, unsigned int> 
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j)
  {
    return storehelpers::val2buffer<size_t>(superdrop.sdId.value,
                                             buffer, ndata, j);
  }
};

struct XiBuffer : SuperdropAttrBuffer<unsigned long long>
{
  XiBuffer()
      : SuperdropAttrBuffer("xi", "<u8"){};

  std::pair<unsigned int, unsigned int> 
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j)
  {
    return storehelpers::
        val2buffer<unsigned long long>(superdrop.get_xi(),
                                       buffer, ndata, j);
  }
};

struct RadiusBuffer : SuperdropAttrBuffer<double>
{
  RadiusBuffer()
      : SuperdropAttrBuffer("radius", "<f8"){};

  std::pair<unsigned int, unsigned int> 
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j)
  {
    return storehelpers::val2buffer<double>(superdrop.get_radius(),
                                            buffer, ndata, j);
  }

  void writejsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    /* write array metadata (and array .zattrs) json */
    SuperdropAttrBuffer::writejsons(store, md);

    /* rewrite array .zattrs json */
    constexpr double sf = dlc::R0 * 1e6; // scale factor to convert dimless radius to microns
    const std::string arrayattrs(storehelpers::
                                     arrayattrs(md.dims, "micro m", sf));
    store[attr + "/.zattrs"] = arrayattrs;
  }
};

#endif // SUPDROPRSATTRSBUFFERS_HPP 
