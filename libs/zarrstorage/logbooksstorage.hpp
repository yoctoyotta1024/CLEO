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

  void set_chunks(const unsigned int ndim1)
  /* given 'size' (number of entries in logbook)
  change ndims1 and chunksize of zarr storage */
  {
    this->set_ndim1(ndim1);
    this->set_buffer_chunksize(storagehelper::good2Dchunk(maxchunk, ndim1));
  }

public:
  LogbooksStorage(FSStore &store, const unsigned int maxchunk,
                  const std::string name, const std::string dtype,
                  const std::string units, const double scale_factor,
                  const std::string i_dim1name)
      : TwoDStorage<T>(store, NOTSETCHUNKSIZE(), name, dtype, units,
                       scale_factor, i_dim1name, NOTSETVALUE()),
        maxchunk(maxchunk) {}

  void prepare(const size_t chunksize)
  {
    set_chunks(chunksize);
  }
};

LogbooksStorage<double> make_logbookszarr(FSStore &store,
                                          const unsigned int maxchunk)
{
  const std::string name("surfpp");
  const std::string dtype("<f8");
  const std::string units("mm");
  constexpr double r0cubed(dlc::R0 * dlc::R0 * dlc::R0);
  constexpr double c0sqrd(dlc::COORD0 * dlc::COORD0);
  constexpr double scale_factor(r0cubed / c0sqrd * 1000); // conversion of precip detection to mm
  const std::string dim1name("logbooktags");

  return LogbooksStorage<double>(store, maxchunk, name, dtype,
                                 units, scale_factor, dim1name);
}

#endif // LOGBOOKSSTORAGE_HPP