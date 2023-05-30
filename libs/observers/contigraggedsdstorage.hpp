// Author: Clara Bayley
// File: contigraggedsdstorage.hpp
/* File for ContiguousRaggedSDStorage
used to store superdroplet attributes
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef CONTIGRAGGEDSDSTORAGE_HPP
#define CONTIGRAGGEDSDSTORAGE_HPP

#include <vector>
#include <string>

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

  SomeMetadata(const unsigned int zarr_format, const char order,
               const unsigned int ndata, const size_t chunksize,
               const std::string compressor, const std::string fill_value,
               const std::string filters, const std::string dims)
      : zarr_format(zarr_format),
        order(order),
        shape("[" + std::to_string(ndata) + "]"),
        chunks("[" + std::to_string(chunksize) + "]"),
        compressor(compressor),
        fill_value(fill_value),
        filters(filters),
        dims(dims)
  {
  }
};

template <typename Aah>
concept SuperdropIntoStoreViaBuffer = requires(Aah aah, const Superdrop &superdrop,
                                               FSStore &store, const std::string &str,
                                               unsigned int j, const unsigned int u,
                                               const SomeMetadata &md)
/* concept SuperdropIntoStoreViaBuffer is all types that have correct
signatures for these 3 void functions. The motivation is that these functions
provide way of copying some superdroplet's data into a buffer, writing buffer to
a chunk of array in the store, and writing array metadata and attribute .json files */
{
  {
    aah.copy2buffer(superdrop, j)
  } -> std::same_as<unsigned int>;

  {
    aah.writechunk(store, j)
  } -> std::same_as<unsigned int>;

  {
    aah.writejsons(store, md)
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

  unsigned int copy2buffer(const Superdrop &superdrop, unsigned int j)
  {
    aah1.copy2buffer(superdrop, j);
    aah2.copy2buffer(superdrop, j);

    return ++j;
  }

  unsigned int writechunk(FSStore &store, unsigned int chunkcount)
  {
    aah1.writechunk(store, chunkcount);
    aah2.writechunk(store, chunkcount);

    return ++chunkcount;
  }

  void writejsons(FSStore &store,
                   const SomeMetadata &md)

  {
    aah1.writejsons(store, md);
    aah2.writejsons(store, md);
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
  unsigned int copy2buffer(const Superdrop &superdrop,
                           const unsigned int j) const
  {
    return j;
  }

  unsigned int writechunk(FSStore &store, const unsigned int chunkcount) const
  {
    return chunkcount;
  }

  void writejsons(FSStore &store, const SomeMetadata &md) const {}
  
  void set_buffersize(const size_t csize) const {}
};

template <SuperdropIntoStoreViaBuffer SDIntoStore>
class ContiguousRaggedSDStorage
/* Class for outputting Superdrop's data into zarr storage in
arrays of contigous ragged representation with 'chunkcount' number
of chunks that have a fixed chunksize. Works by filling buffers in
sdbuffers with superdrop data and then writing these buffers
into chunks in their corresponding array stores when number of
datapoints copied to the buffers reaches chunksize. */
{
private:
  FSStore &store;               // file system store satisfying zarr store specificaiton v2
  SDIntoStore sdbuffers;        // buffers and their handler functions for wrting SD data to store
  std::vector<size_t> rgdcount; // count variable for contiguous ragged representation of arrays

  const size_t chunksize;  // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed (= size of array written to store)

  unsigned int rgdcount_chunkcount; // number of chunks of rgdcount array so far written to store
  unsigned int rgdcount_bufferfill; // number of rgdcount values so far copied into its buffer
  unsigned int rgdcount_ndata;      // number of rgdcount values observed so far

  constexpr unsigned int zarr_format = 2; // storage spec. version 2
  constexpr char order = 'C';                 // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  constexpr std::string compressor = "null";  // compression of data when writing to store
  constexpr std::string fill_value = "null";  // fill value for empty datapoints in array
  constexpr std::string filters = "null";     // codec configurations for compression

  constexpr std::string rgdcount_name = "raggedcount"; // name of rgdcount zarray in store
  constexpr std::string rgdcount_dtype = "<u8";        // datatype of rgdcount variable

  void sdbuffers_writejsons()
  {
    // write strictly required metadata to decode chunks (MUST)
    const std::string dims = "[\"sdindex\"]";
    const SomeMetadata md(zarr_format, order, ndata, chunksize,
                          compressor, fill_value, filters, dims);
    sdbuffers.writejsons(store, md);
  }

  void rgdcount_writejsons()
  /* write zarray jsons for array of rgdcount variable in store */
  {
    constexpr std::string
        count_arrayattrs = "{\"_ARRAY_DIMENSIONS\": [\"time\"],"
                           "\"sample_dimension\": \"superdroplets\"}";

    const std::string
        count_metadata = storagehelper::
            metadata(zarr_format, order, rgdcount_ndata, chunksize,
                     rgdcount_dtype, compressor, fill_value, filters);

    storagehelper::writezarrjsons(store, rgdcount_name,
                                        count_metadata,
                                        count_arrayattrs);
  }

  void sdbuffers_writechunk()
  /* write data in sdbuffers to chunks of zarrays in store
  and (re)write associated metadata for zarrays */
  {
    chunkcount = sdbuffers.writechunk(store, chunkcount);
    bufferfill = 0; // reset bufferfill

    sdbuffers_writejsons();
  }

  void rgdcount_writechunk()
  /* write rgdcount data in buffers to a chunk of its
  zarray in store and (re)write its associated metadata */
  {
    rgdcount_chunkcount = storagehelper::writebuffer2chunk(store, rgdcount,
                                                           rgdcount_name,
                                                           rgdcount_chunkcount);
    rgdcount_bufferfill = 0; // reset rgdcount bufferfill

    rgdcount_writejsons();
  }

  template <typename T>
  void copy2sdbuffers(const T &value)
  /* copy data from superdrop to buffer(s) and
  increment required counting variables */
  {
    bufferfill = sdbuffers.copy2buffer(value, bufferfill);
    ++ndata;
  }

  void copy2rgdcount(const size_t raggedn)
  /* write raggedn into rgdcount buffer */
  {
    rgdcount_bufferfill = storagehelper::
        val2buffer<size_t>(raggedn, rgdcount, rgdcount_bufferfill);
    ++rgdcount_ndata;
  }

public:
  ContiguousRaggedSDStorage(FSStore &store,
                            const SDIntoStore sdbuffers_i,
                            const size_t csize)
      : store(store), sdbuffers(sdbuffers_i), rgdcount(csize),
        chunksize(csize), chunkcount(0), bufferfill(0), ndata(0),
        rgdcount_chunkcount(0), rgdcount_bufferfill(0),
        rgdcount_ndata(0)
  {
    // initialise buffer(s) to size 'chunksize' (filled with numeric limit)
    sdbuffers.set_buffersize(chunksize);
  }

  ~ContiguousRaggedSDStorage()
  {
    if (bufferfill != 0)
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
    if (bufferfill == chunksize)
    {
      sdbuffers_writechunk();
    }

    copy2sdbuffers(value);
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
};

#endif // CONTIGRAGGEDSDSTORAGE_HPP