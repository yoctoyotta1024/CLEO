// Author: Clara Bayley
// File: "singlevarstorage.hpp"
/* Classes and helper functions in a namespace
useful for using storage clases with buffers to
write values of 1D data into chunks of arrays
in a zarr store */

#ifndef SINGLEVARSTORAGE_HPP 
#define SINGLEVARSTORAGE_HPP 

#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>

#include "./zarrstores.hpp"

template <typename T>
class SingleVarStorage
{
private:
  virtual unsigned int writechunk() = 0;

protected:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2
  const std::string name;    // name to call variable being stored
  const std::string units;   // units of coordinate being stored (for arrayattrs json)
  const double scale_factor; // scale_factor of data (for array .zattrs json)
  std::vector<T> buffer;     // buffer to store values in until writing to array chunk

  const size_t chunksize;  // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed

  const unsigned int zarr_format = 2;    // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype;               // datatype stored in arrays

  void zarrayjsons(const std::string shape,
                   const std::string chunks,
                   const std::string dims)
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    const std::string metadata = storagehelper::metadata(zarr_format, order,
                                                         shape, chunks, dtype,
                                                         compressor, fill_value,
                                                         filters);

    const std::string arrayattrs = storagehelper::arrayattrs(dims, units,
                                                             scale_factor);
    
    storagehelper::write_zarrarrayjsons(store, name,
                                        metadata, arrayattrs);
  }

public:
  SingleVarStorage(FSStore &store, const unsigned int maxchunk,
                   const std::string name, const std::string dtype,
                   const std::string units, const double scale_factor)
      : store(store), name(name), units(units), scale_factor(scale_factor),
        buffer(maxchunk, std::numeric_limits<T>::max()), chunksize(maxchunk),
        chunkcount(0), bufferfill(0), ndata(0), dtype(dtype) {}

  virtual ~SingleVarStorage(){};
  
  std::string get_name() const { return name; }

  int get_ndata() const { return ndata; }

  template <typename T1>
  void value_to_storage(const T1 val)
  /* write val in the zarr store. First copy it to a buffer,
  then write buffer to a chunk in the store when the number
  of values in the buffer reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      chunkcount = writechunk();
      bufferfill = 0;
    }

    // copy double to buffer
    storagehelper::val2buffer<T>(val, buffer, bufferfill);
    ++(bufferfill);
    ++(ndata);
  }
};

template <typename T>
struct CoordinateStorage : SingleVarStorage<T>
/* storage of a 1D variable with 'dims' in .zattrs metadata
equal to name of variable (ie. variable is an xarray coord)*/
{
private:

  unsigned int writechunk()
  /* write data in buffer to a chunk in store */
  {
    const std::string chunknum = std::to_string(this->chunkcount);
    storagehelper::writebuffer2chunk(this->store, this->buffer,
                                     this->name, chunknum);

    return ++(this->chunkcount);
  }

public:
  CoordinateStorage(FSStore &store, const unsigned int maxchunk, const std::string name,
               const std::string dtype, const std::string units,
               const double scale_factor)
      : SingleVarStorage<T>(store, maxchunk, name,
                            dtype, units, scale_factor) {}

  ~CoordinateStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      this->chunkcount = writechunk(); 
    }

    // write strictly required metadata to decode chunks (MUST)
    const auto shape("[" + std::to_string(this->ndata) + "]");
    const auto chunks("[" + std::to_string(this->chunksize) + "]");
    const std::string dims = "[\"" + this->name + "\"]";

    this->zarrayjsons(shape, chunks, dims); 
  }
};

template <typename T>
struct TwoDStorage : SingleVarStorage<T>
{
private:
  const unsigned int ngridboxes; // number of output times that have been observed 

  unsigned int writechunk()
  /* write data in buffer to a chunk in store */
  {
    const std::string chunknum = std::to_string(this->chunkcount)+".0";
    storagehelper::writebuffer2chunk(this->store, this->buffer,
                                     this->name, chunknum);

    return ++(this->chunkcount);
  }

public:
  unsigned int nobs; // number of output times that have been observed

  TwoDStorage(FSStore &store, const unsigned int maxchunk, const std::string name,
              const std::string dtype, const std::string units,
              const double scale_factor, const unsigned int ngrid)
      : SingleVarStorage<T>(store, floor(maxchunk / ngrid) * ngrid, name,
                            dtype, units, scale_factor),
        ngridboxes(ngrid), nobs(0) {}

  ~TwoDStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      this->chunkcount = writechunk(); 
    }

    // write strictly required metadata to decode chunks (MUST)
    assert((this->ndata == nobs*ngridboxes) && "1D data length matches 2D array size");
    assert((this->chunksize % ngridboxes == 0.0) && "chunks are integer multple of number of gridboxes");
    const auto ngstr = std::to_string(ngridboxes);
    const auto shape("[" + std::to_string(nobs) + ", " + ngstr + "]");
    const auto chunks("[" + std::to_string(this->chunksize/ngridboxes) + ", " + ngstr + "]");
    const std::string dims = "[\"time\", \"gbxindex\"]";
    
    this->zarrayjsons(shape, chunks, dims); 
  }
};

#endif // SINGLEVARSTORAGE_HPP 