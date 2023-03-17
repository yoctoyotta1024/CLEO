// Author: Clara Bayley
// File: observer_superdropletattributes.hpp
/* structs to create a SDsAttributeObserver that writes 
data for (any chosen combination of) a superdroplet's 
attribute(s) into ragged contiguous array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef OBSERVER_SUPERDROPLETATTRIBUTES_HPP
#define OBSERVER_SUPERDROPLETATTRIBUTES_HPP

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <concepts>

#include "zarrstores.hpp"

struct SomeMetadata
{
  unsigned int zarr_format;
  char order;
  std::string shape;
  std::string chunks;
  std::string compressor;
  std::string fill_value;
  std::string filters;
  std::string dims;
};

template <typename Aah>
concept SuperdropIntoStoreViaBuffer = requires(Aah aah, const Superdrop &superdrop,
                                               FSStore &store, const std::string &str,
                                               const int j, const unsigned int u,
                                               const SomeMetadata &md)
/* concept SuperdropIntoStoreViaBuffer is all types that have correct
signatures for these 3 void functions. The motivation is that these functions
provide way of copying some superdroplet's data into a buffer, writing buffer to
a chunk of array in the store, and writing array metadata and attribute .json files */
{
  {
    aah.copy2buffer(superdrop, j)
    } -> std::same_as<void>;

  {
    aah.writechunk(store, j)
    } -> std::same_as<void>;

  {
    aah.zarrayjsons(store, md)
    } -> std::same_as<void>;
  
  {
    aah.set_buffersize(u)
    } -> std::same_as<void>;
};

template <SuperdropIntoStoreViaBuffer A1, SuperdropIntoStoreViaBuffer A2>
struct CombinedSuperdropIntoStoreViaBuffer
/* combination of two types (A1, A2) that satisfiy 
SuperdropIntoStoreViaBuffer is A1 followed by A2 */
{
  A1 aah1;
  A2 aah2;

  CombinedSuperdropIntoStoreViaBuffer(A1 aah1, A2 aah2)
      : aah1(aah1), aah2(aah2) {}

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    aah1.copy2buffer(superdrop, j);
    aah2.copy2buffer(superdrop, j);
  }

  void writechunk(FSStore &store, const int chunkcount)
  {
    aah1.writechunk(store, chunkcount);
    aah2.writechunk(store, chunkcount);
  }

  void zarrayjsons(FSStore &store,
                     const SomeMetadata &md)
                     
  {
    aah1.zarrayjsons(store, md);
    aah2.zarrayjsons(store, md);
  }

  void set_buffersize(const size_t csize)
  {
    aah1.set_buffersize(csize);
    aah2.set_buffersize(csize);
  }
};

auto operator>>(SuperdropIntoStoreViaBuffer auto aah1,
                SuperdropIntoStoreViaBuffer auto aah2)
/* define ">>" operator that combines two 
SuperdropIntoStoreViaBuffer types */
{
  return CombinedSuperdropIntoStoreViaBuffer{aah1, aah2};
}

struct NullSuperdropIntoStoreViaBuffer
/* Null does nothing at all (is defined for 
completeness of a Monoid Structure) */
{
  void copy2buffer(const Superdrop &superdrop, const int j) const {}
  void writechunk(FSStore &store, const int chunkcount) const {}
  void zarrayjsons(FSStore &store, const SomeMetadata &md) const {}
  void set_buffersize(const size_t csize) const {}
};

template <SuperdropIntoStoreViaBuffer SDIntoStore>
class ContiguousRaggedSuperdropStorage
/* Class for outputting Superdrop's data into zarr storage in
arrays of contigous ragged representation with 'chunkcount' number
of chunks that have a fixed chunksize. Works by filling buffers in
sdbuffers with superdrop data and then writing these buffers
into chunks in their corresponding array stores when number of
datapoints copied to the buffers reaches chunksize. */
{
private:
  FSStore &store;        // file system store satisfying zarr store specificaiton v2
  SDIntoStore sdbuffers; // buffers and their handler functions for wrting SD data to store
  std::vector<size_t> raggedcount; // count variable for contiguous ragged representation of arrays

  const size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount;      // number of chunks of array so far written to store
  unsigned int bufferfill;      // number of datapoints so far copied into buffer
  unsigned int ndata;           // number of data points that have been observed (= size of array written to store)

  unsigned int raggedcount_chunkcount;      // number of chunks of raggedcount array so far written to store
  unsigned int raggedcount_bufferfill;      // number of raggedcount values so far copied into its buffer
  unsigned int raggedcount_ndata;           // number of raggedcount values observed so far

  const unsigned int zarr_format = 2;    // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
 
  void raggedcount_zarrayjsons()
  /* store count variable array 'raggedcount', 
  in 1 chunk in store under 'count_ragged'*/
  {
    const std::string count_name("raggedcount");
    const std::string count_dtype("<u8");
    const auto count_shape("[" + std::to_string(raggedcount_ndata) + "]");
    const auto count_chunks("[" + std::to_string(chunksize) + "]");

    const std::string count_metadata = storagehelper::metadata(zarr_format, order,
                                                         count_shape, count_chunks,
                                                         count_dtype, compressor,
                                                         fill_value, filters);
    
    const std::string count_arrayattrs = "{\"_ARRAY_DIMENSIONS\": [\"time\"],"
                                         "\"sample_dimension\": \"superdroplets\"}";
    
    storagehelper::write_zarrarrayjsons(store, count_name, count_metadata, count_arrayattrs); 
  }

public:
  ContiguousRaggedSuperdropStorage(FSStore &store, const SDIntoStore sdbuffers_i,
                                  const size_t csize)
      : store(store), sdbuffers(sdbuffers_i), raggedcount(csize),
        chunksize(csize), chunkcount(0), bufferfill(0), ndata(0),
        raggedcount_chunkcount(0), raggedcount_bufferfill(0), raggedcount_ndata(0)
  {
    // initialise buffer(s) to size 'chunksize' (filled with nans)
    sdbuffers.set_buffersize(chunksize);                                                   
  }

  ~ContiguousRaggedSuperdropStorage()
  {
    if (bufferfill != 0)
    {
      // write data in buffer to a chunk in store
      sdbuffers.writechunk(store, chunkcount);
      ++chunkcount;
    }

    if (raggedcount_bufferfill != 0)
    {
      // write data in buffer to a chunk in store
      const std::string chunknum = std::to_string(raggedcount_chunkcount);
      storagehelper::writebuffer2chunk(store, raggedcount, "raggedcount", chunknum);
      ++raggedcount_chunkcount;
    } 

    // write strictly required metadata to decode chunks (MUST)
    const auto shape("[" + std::to_string(ndata) + "]");
    const auto chunks("[" + std::to_string(chunksize) + "]");
    const std::string dims = "[\"sdindex\"]";
    const SomeMetadata md = {zarr_format, order, shape, chunks,
                            compressor, fill_value, filters, dims};
    sdbuffers.zarrayjsons(store, md);
    
    // write count variable array for contiguous ragged representation of arrays
    raggedcount_zarrayjsons();
  }

  void data_to_contigraggedarray(const Superdrop &superdrop)
  /* write superdrop's data in contoguous ragged representation of an array
  in the zarr store. First copy data to buffer(s), then write buffer(s) to chunks
  in the store when the number of datapoints they contain reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      // write data in buffer to a chunk in store
      sdbuffers.writechunk(store, chunkcount);
      ++chunkcount;
      bufferfill = 0;
    }

    // copy data from superdrop to buffer(s)
    sdbuffers.copy2buffer(superdrop, bufferfill);
    ++bufferfill;

    ++ndata;
  }

  void contigraggedarray_count(const size_t n)
  /* add element to raggedcount that is number of datapoints
  written to buffer(s) during one event. This is count variable 
  for contiguous ragged representation */
  {
    if (raggedcount_bufferfill == chunksize)
    {
      // write data in buffer to a chunk in store
      const std::string chunknum = std::to_string(raggedcount_chunkcount);
      storagehelper::writebuffer2chunk(store, raggedcount, "raggedcount", chunknum);
      ++raggedcount_chunkcount;
      raggedcount_bufferfill = 0;
    }

  // copy double to buffer
  storagehelper::val2buffer<size_t>(n, raggedcount, raggedcount_bufferfill);
  ++raggedcount_bufferfill;
  ++raggedcount_ndata;
  }
};

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
      buffer = std::vector<T>(csize, std::nan(""));
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

  void set_buffersize(const size_t csize)
  {
    if (buffer.size() != csize)
    {
      buffer = std::vector<size_t>(csize, -1);
    }
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

  void set_buffersize(const size_t csize)
  {
    if (buffer.size() != csize)
    {
      buffer = std::vector<size_t>(csize, -1);
    }
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

struct Coord3IntoStore : AttributeIntoStoreViaBuffer<double>
{
  Coord3IntoStore()
      : AttributeIntoStoreViaBuffer("coord3", "<f8"){};

  void copy2buffer(const Superdrop &superdrop, const int j)
  {
    storagehelper::val2buffer<double>(superdrop.coord3, buffer, j);
  }

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

#endif // OBSERVER_SUPERDROPLETATTRIBUTES_HPP