// Author: Clara Bayley
// File: thermostatestorage.hpp
/* structs/classes to create a ThermoStateObserver that writes
data from thermostate into orthogonal multidimensional array(s)
(see: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#_contiguous_ragged_array_representation)
in a FFStore obeying zarr storage specification verion 2:
https://zarr.readthedocs.io/en/stable/spec/v2.html */

#ifndef THERMOSTATESTORAGE_HPP
#define THERMOSTATESTORAGE_HPP

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

#include "claras_SDconstants.hpp"
#include "./zarrstores.hpp"
#include "superdrop_solver/thermostate.hpp"

namespace dlc = dimless_constants

struct ThermoIntoStore
{
  std::vector<double> pressbuffer;
  std::vector<double> tempbuffer;
  std::vector<double> qvapbuffer;
  std::vector<double> qcondbuffer;

  ThermoIntoStore(const unsigned int buffersize)
      : pressbuffer(buffersize, std::nan("")), tempbuffer(buffersize, std::nan("")),
        qvapbuffer(buffersize, std::nan("")), qcondbuffer(buffersize, std::nan("")) {}

  void copy2buffers(const ThermoState &state, const int j);
  /* copy press, temp, qvap and qcond data in the state to buffers at index j */

  void writechunks(FSStore &store, const int chunkcount);
  /* write buffer vector into attr's store at chunkcount
  and then replace contents of buffer with std::nans */

  void zarrayjsons(FSStore &store, const std::string &metadata) const;
  /* write same .zarray metadata to a json file for each thermostate array
  in store alongside distinct .zattrs json files */

};

class ThermoStateStorage
{
private:
  FSStore &store;           // file system store satisfying zarr store specificaiton v2
  ThermoIntoStore buffers; // buffers and their handler functions for wrting SD data to store

  const size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount;      // number of chunks of array so far written to store
  unsigned int bufferfill;      // number of datapoints so far copied into buffer
  unsigned int ndata; // number of data points that have been observed should = nobs * ngridboxes

  const unsigned int zarr_format = 2;    // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype = "<f8";       // datatype stored in arrays

  const unsigned int ngridboxes; // number of output times that have been observed 

public:
  unsigned int nobs; // number of output times that have been observed 
  
  ThermoStateStorage(FSStore &store, const unsigned int maxcsize,
                     const unsigned int ngrid)
      : store(store), buffers(floor(maxcsize / ngrid)*ngrid), 
        chunksize(floor(maxcsize / ngrid)*ngrid), chunkcount(0),
        bufferfill(0), ndata(0), ngridboxes(ngrid), nobs(0) {}

  ~ThermoStateStorage()
  /* upon destruction write any data leftover in buffers 
  to chunks and write arrays' metadata to .json files */
  {
    if (bufferfill != 0)
    {
      // write data in buffer to a chunk in store
      buffers.writechunks(store, chunkcount);
      ++chunkcount;
    }

    // write strictly required metadata to decode chunks (MUST)
    assert((ndata == nobs*ngridboxes) && "1D data length matches 2D array size");
    assert((chunksize % ngridboxes == 0.0) && "chunks are integer multple of number of gridboxes");
    const auto ngstr = std::to_string(ngridboxes);
    const auto shape("[" + std::to_string(nobs) + ", " + ngstr + "]");
    const auto chunks("[" + std::to_string(chunksize/ngridboxes) + ", " + ngstr + "]");
    const std::string metadata = storagehelper::metadata(zarr_format, order,
                                                         shape, chunks, dtype,
                                                         compressor, fill_value,
                                                         filters);
    buffers.zarrayjsons(store, metadata);
  }

  void thermodata_to_storage(const ThermoState &state);
  /* write thermo variables from a thermostate in arrays in the zarr store. 
  First copy data to buffers, then write buffers to chunks in the store 
  when the number of datapoints they contain reaches the chunksize */

};

#endif // THERMOSTATESTORAGE_HPP 