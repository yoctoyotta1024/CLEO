/*
 * ----- CLEO -----
 * File: superdropsbuffers.hpp
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
 * File for concept of buffers used by
 * Contiguous Ragged Array Storage
 * to store superdroplet attributes
 * (see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
 * in a FFStore obeying zarr storage specification verion 2:
 * https://zarr.readthedocs.io/en/stable/spec/v2.html */


#ifndef SUPERDROPSBUFFERS_HPP 
#define SUPERDROPSBUFFERS_HPP 

#include <concepts>
#include <string>
#include <utility>

#include "./fsstore.hpp
#include "superdrops/superdrop.hpp"

struct SomeMetadata
{
  const char zarr_format;
  const char order;
  const std::string shape;
  const std::string chunks;
  const std::string compressor;
  const std::string fill_value;
  const std::string filters;
  const std::string dims;

  SomeMetadata(const char zarr_format, const char order,
               const unsigned int ndata, const size_t chunksize,
               const std::string compressor,
               const std::string fill_value,
               const std::string filters,
               const std::string dims)
      : zarr_format(zarr_format),
        order(order),
        shape("[" + std::to_string(ndata) + "]"),
        chunks("[" + std::to_string(chunksize) + "]"),
        compressor(compressor),
        fill_value(fill_value),
        filters(filters),
        dims(dims) {}
};

template <typename B>
concept SuperdropsBuffers = requires(B b,
                                    Superdrop &superdrop,
                                    FSStore &store,
                                    std::string &str,
                                    unsigned int u,
                                    SomeMetadata &md)
/* concept for all types that have correct signatures
for functions which provide way of copying some superdroplet's
data into a buffer, writing buffer to a chunk of array in the
store, and writing array metadata and attribute .json files
using the contigraggedstorage structure */
{
  {
    b.copy2buffer(superdrop, u, u)
  } -> std::same_as<std::pair<unsigned int, unsigned int>>;

  {
    b.writechunk(store, u)
  } -> std::same_as<std::pair<unsigned int, unsigned int>>;

  {
    b.writejsons(store, md)
  } -> std::same_as<void>;

  {
    b.set_buffersize(u)
  } -> std::same_as<void>;
};

template <SuperdropsBuffers B1, SuperdropsBuffers B2>
struct CombinedSuperdropsBuffers
/* combination of two types (B1, B2) that satisfiy
SuperdropsBuffers is B1 followed by B2 */
{
  B1 a;
  B2 b;

  CombinedSuperdropsBuffers(B1 a, B2 b)
      : b1(a), b2(b) {}

  std::pair<unsigned int, unsigned int>
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j)
  {
    a.copy2buffer(superdrop, j);
    b.copy2buffer(superdrop, j);

    return std::pair(ndata + 1, j + 1); // TODO update this to std:pair with bufferfill and ndata
  }

  std::pair<unsigned int, unsigned int>
  writechunk(FSStore &store, const unsigned int chunkcount)
  {
    a.writechunk(store, chunkcount);
    b.writechunk(store, chunkcount);

    return std::pair(chunkcount + 1, 0); // updated {chunkcount, bufferfill}
  }

  void writejsons(FSStore &store, const SomeMetadata &md)

  {
    a.writejsons(store, md);
    b.writejsons(store, md);
  }

  void set_buffersize(const size_t maxchunk)
  {
    a.set_buffersize(maxchunk);
    b.set_buffersize(maxchunk);
  }
};

auto operator>>(SuperdropsBuffers auto a, SuperdropsBuffers auto b)
/* define ">>" operator that combines two
SuperdropsBuffers types */
{
  return CombinedSuperdropsBuffers(a, b);
}

struct NullSuperdropsBuffer
/* Null does nothing at all (is defined for
completeness of a Monoid Structure) */
{
  std::pair<unsigned int, unsigned int>
  copy2buffer(const Superdrop &superdrop,
              const unsigned int ndata, const unsigned int j) const
  {
    return std::pair(ndata, j);
  }

  std::pair<unsigned int, unsigned int>
  writechunk(FSStore &store, const unsigned int chunkcount) const
  {
    return std::pair(chunkcount, 0);
  }

  void writejsons(FSStore &store, const SomeMetadata &md) const {}
  
  void set_buffersize(const size_t maxchunk) const {}
};

#endif //SUPERDROPSBUFFERS_HPP 