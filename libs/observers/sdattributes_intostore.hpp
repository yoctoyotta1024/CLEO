// Author: Clara Bayley
// File: sdattributes_intostore.hpp
/* structs to create a SDsAttributeObserver that writes 
data for (any chosen combination of) a superdroplet's 
attribute(s) into ragged contiguous array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef SDATTRIBUTES_INTOSTORE_HPP 
#define SDATTRIBUTES_INTOSTORE_HPP 

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <concepts>
#include <limits>

#include "../claras_SDconstants.hpp"
#include "./zarrstores.hpp"
#include "./contigraggedsdstorage.hpp"

namespace dlc = dimless_constants;

template <typename T>
struct AttributeIntoStoreViaBuffer
/* generic structure satisfying SuperdropIntoStoreViaBuffer concept.
Useful for putting a superdrop's attribute called attr into a buffer
(given an implementation of the copy2buffer function) and then writing
buffer into an array in a Zarr store using writechunk and writemetadata */
{
  const std::string attr;                // name of attribute in fsstore
  const std::string dtype;               // datatype stored in arrays
  std::vector<T> buffer;                 // buffer to fill before writing to store

  AttributeIntoStoreViaBuffer(const std::string attr, const std::string dtype)
      : attr(attr), dtype(dtype), buffer(0) {}

  virtual ~AttributeIntoStoreViaBuffer(){};
  
  virtual void copy2buffer(const Superdrop &superdrop, const int j) = 0;
  /* virtual void function placeholding function for
  copying superdrop's data into a buffer vector at j'th index */

  void writechunk(FSStore &store, const int chunkcount)
  /* write buffer vector into attr's store at chunkcount
  and then replace contents of buffer with std::nans */
  {
    const std::string chunknum = std::to_string(chunkcount);
    storagehelper::writebuffer2chunk(store, buffer, attr, chunknum);
  }

  void zarrayjsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    const std::string metadata = storagehelper::metadata(md.zarr_format, md.order,
                                                         md.shape, md.chunks, dtype,
                                                         md.compressor, md.fill_value,
                                                         md.filters);
    
    const std::string arrayattrs = "{\"_ARRAY_DIMENSIONS\": "+md.dims+"}";
    
    storagehelper::write_zarrarrayjsons(store, attr, metadata, arrayattrs);
  }

  void set_buffersize(const size_t csize)
  {
    if (buffer.size() != csize)
    {
      buffer = std::vector<T>(csize, std::numeric_limits<T>::max());
    }
  } 
};

struct IdIntoStore : AttributeIntoStoreViaBuffer<size_t>
{
  IdIntoStore()
      : AttributeIntoStoreViaBuffer("sdindex", "<u8"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<size_t>(superdrop.id.value, buffer, j);
  }
};

struct EpsIntoStore : AttributeIntoStoreViaBuffer<size_t>
{
  EpsIntoStore()
      : AttributeIntoStoreViaBuffer("eps", "<u8"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<size_t>(superdrop.eps, buffer, j);
  }
};

struct RadiusIntoStore : AttributeIntoStoreViaBuffer<double>
{
  RadiusIntoStore()
      : AttributeIntoStoreViaBuffer("radius", "<f8"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.radius, buffer, j);
  }

  void zarrayjsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    /* write array metadata (and array attrs) json */
    AttributeIntoStoreViaBuffer::zarrayjsons(store, md);

    /* rewrite array attrs json */
    const double scale_factor = dlc::R0 * 1e6; // convert radius to microns
    const std::string arrayattrs = storagehelper::arrayattrs(md.dims, "micro m",
                                                             scale_factor);
    store[attr + "/.zattrs"] = arrayattrs; 
  }
};

struct M_solIntoStore : AttributeIntoStoreViaBuffer<double>
{
  M_solIntoStore()
      : AttributeIntoStoreViaBuffer("m_sol", "<f8"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.m_sol, buffer, j);
  }

  void zarrayjsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    /* write array metadata (and array attrs) json */
    AttributeIntoStoreViaBuffer::zarrayjsons(store, md);

    /* rewrite array attrs json */
    const double scale_factor = pow(dlc::R0, 3.0) * dlc::RHO0 * 1000; // convert mass to grams
    const std::string arrayattrs = storagehelper::arrayattrs(md.dims, "g",
                                                             scale_factor);
    store[attr + "/.zattrs"] = arrayattrs; 
  }
};

struct SdCoordIntoStore : AttributeIntoStoreViaBuffer<double>
{
  SdCoordIntoStore(const std::string attr)
      : AttributeIntoStoreViaBuffer(attr, "<f8"){};

  virtual void copy2buffer(const Superdrop &superdrop, const int j) = 0;

  void zarrayjsons(FSStore &store, const SomeMetadata &md) const
  /* write metadata for attr's array into store */
  {
    /* write array metadata (and array attrs) json */
    AttributeIntoStoreViaBuffer::zarrayjsons(store, md);

    /* rewrite array attrs json */
    const std::string arrayattrs = storagehelper::arrayattrs(md.dims, "m",
                                                             dlc::COORD0);
    store[attr + "/.zattrs"] = arrayattrs; 
  }
};

struct Coord3IntoStore : SdCoordIntoStore
{
  Coord3IntoStore() : SdCoordIntoStore("coord3"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.coord3, buffer, j);
  }
};

struct Coord1IntoStore : SdCoordIntoStore
{
  Coord1IntoStore() : SdCoordIntoStore("coord1"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.coord1, buffer, j);
  }
};

struct Coord2IntoStore : SdCoordIntoStore
{
  Coord2IntoStore() : SdCoordIntoStore("coord2"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.coord2, buffer, j);
  }
};

#endif // SDATTRIBUTES_INTOSTORE_HPP 