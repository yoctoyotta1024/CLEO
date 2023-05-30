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

#include "../claras_SDconstants.hpp"
#include "./zarrstores.hpp"
#include "superdrop_solver/thermostate.hpp"

namespace dlc = dimless_constants;

struct ThermoIntoStore
{
  std::vector<double> pressbuffer;
  std::vector<double> tempbuffer;
  std::vector<double> qvapbuffer;
  std::vector<double> qcondbuffer;

  ThermoIntoStore(const unsigned int buffersize)
      : pressbuffer(buffersize, std::numeric_limits<double>::max()),
        tempbuffer(buffersize, std::numeric_limits<double>::max()), 
        qvapbuffer(buffersize, std::numeric_limits<double>::max()), 
        qcondbuffer(buffersize, std::numeric_limits<double>::max()) {}

  unsigned int copy2buffers(const ThermoState &state, unsigned int j);
  /* copy press, temp, qvap and qcond data in the state to buffers at index j */

  unsigned int writechunks(FSStore &store, unsigned int chunkcount);
  /* write buffer vector into attr's store at chunkcount
  and then replace contents of buffer with numeric limit */

  void zarrayjsons(FSStore &store, const std::string &metadata) const;
  /* write same .zarray metadata to a json file for each thermostate array
  in store alongside distinct .zattrs json files */
};

class ThermoStateStorage
{
private:
  FSStore &store;          // file system store satisfying zarr store specificaiton v2
  ThermoIntoStore buffers; // buffers and their handler functions for wrting SD data to store

  const size_t chunksize;  // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed should = nobs * ngridboxes

  const unsigned int zarr_format = 2;    // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype = "<f8";       // datatype stored in arrays

  const unsigned int ngridboxes; // number of output times that have been observed 

  void copy2buffers(const ThermoState &state)
  /* copy data from thermostate to buffers */
  {
    bufferfill = buffers.copy2buffers(state, bufferfill);
    ++ndata;
  }

  void writechunks()
  /* write data from thermo buffers into chunks in store,
  then reset bufferfill and write associated metadata */
  {
    chunkcount = buffers.writechunks(store, chunkcount);
    bufferfill = 0;

    zarrayjsons();
  }

  void zarrayjsons()
  /* write strictly required metadata to decode chunks (MUST) */
  {
    assert((ndata == nobs*ngridboxes) && "1D data length matches 2D array size");
    assert((chunksize % ngridboxes == 0.0) && "chunks are integer multple of number of gridboxes");
    
    const auto ngstr = std::to_string(ngridboxes);
    const auto nobstr = std::to_string(nobs);
    const auto nchstr = std::to_string(chunksize / ngridboxes);

    const auto shape("[" + nobstr + ", " + ngstr + "]");
    const auto chunks("[" + nchstr+ ", " + ngstr + "]");
    const auto metadata(storagehelper::
                            metadata(zarr_format, order, shape,
                                     chunks, dtype, compressor,
                                     fill_value, filters));
    buffers.zarrayjsons(store, metadata);
  }

public:
  unsigned int nobs; // number of output times that have been observed 
  
  ThermoStateStorage(FSStore &store, const unsigned int maxchunk,
                     const unsigned int ngrid)
      : store(store), buffers(floor(maxchunk / ngrid)*ngrid), 
        chunksize(floor(maxchunk / ngrid)*ngrid), chunkcount(0),
        bufferfill(0), ndata(0), ngridboxes(ngrid), nobs(0) {}

  ~ThermoStateStorage()
  /* upon destruction write any data leftover in buffers 
  to chunks and write arrays' metadata to .json files */
  {
    if (bufferfill != 0)
    {
      writechunks();
    }
  }

  void thermodata_to_storage(const ThermoState &state)
  /* write thermo variables from a thermostate in arrays in the zarr store.
  First copy data to buffers, then write buffers to chunks in the store
  when the number of datapoints they contain reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      writechunks()
    }

    copy2buffers(state);
  }
};

#endif // THERMOSTATESTORAGE_HPP 