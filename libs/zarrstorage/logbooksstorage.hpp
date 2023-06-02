// Author: Clara Bayley
// File: "logbooksstorage.hpp"
/* Classes and helper functions
writing a records from logbooks
via a buffer into chunks of arrays
in a zarr store */

#ifndef LOGBOOKSSTORAGE_HPP
#define LOGBOOKSSTORAGE_HPP

#include <string>

#include "./zarrstores.hpp"
#include "./singlevarstorage.hpp"

template <typename T>
struct LogbooksStorage : TwoDStorage<T>
{
private:
  const unsigned int maxchunk;

public:
  LogbooksStorage(FSStore &store, const unsigned int maxchunk,
                  const std::string name, const std::string dtype,
                  const std::string units, const double scale_factor,
                  const std::string i_dim1name)
      : TwoDStorage<T>(store, 0, name, dtype, units,
                       scale_factor, i_dim1name, 0),
        maxchunk(maxchunk) {}

  void set_chunks(const unsigned int ndim1)
  /* given 'size' (number of entries in logbook)
  change ndims1 and chunksize of zarr storage */
  {
    this->set_ndim1(ndim1);
    this->set_chunksize(storagehelper::good2Dchunk(maxchunk, ndim1));
  }
};

#endif // LOGBOOKSSTORAGE_HPP