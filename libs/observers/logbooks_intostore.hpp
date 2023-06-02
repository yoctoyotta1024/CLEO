// Author: Clara Bayley
// File: "logbooks_intostore.hpp"
/* structures obeying the ObserveLbks
concept for various ways of observing
logbooks which ends up writing data
into a (zarr) store on disk */

#ifndef LOGBOOKS_INTOSTORE_HPP
#define LOGBOOKS_INTOSTORE_HPP

#include <vector>

#include "sdmgridboxes/logbooks.hpp"
#include "zarrstorage/singlevarstorage.hpp"

template <typename T>
struct LogbookStorage : SingleVarStorage<T>
{
private:
  const unsigned int ndim0;

  void writechunk()
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    // const std::string chunknum = std::to_string(this->chunkcount) + ".0";
    // this->chunkcount = storagehelper::
    //     writebuffer2chunk(this->store, this->buffer,
    //                       this->name, chunknum,
    //                       this->chunkcount);

    // writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST).
  Assert also check 2D data dimensions is as expected */
  {
    // assert((this->ndata == nobs * ndim0) &&
    //        "1D data length matches 2D array size");
    // assert((this->chunksize % ndim0 == 0.0) &&
    //        "chunks are integer multple of 0th dimension");

    // const auto ngstr = std::to_string(ndim0);
    // const auto nobstr = std::to_string(nobs);
    // const auto nchstr = std::to_string(this->chunksize / ndim0);

    // const auto shape("[" + nobstr + ", " + ngstr + "]");
    // const auto chunks("[" + nchstr + ", " + ngstr + "]");
    // // const std::string dims = "[\"time\", \"gbxindex\"]";
    // // this->zarrayjsons(shape, chunks, dims);
  }

public:
  unsigned int nobs; // number of output times that have been observed

  LogbookStorage(FSStore &store, const unsigned int maxchunk,
                  const std::string name, const std::string dtype,
                  const std::string units, const double scale_factor,
                  const unsigned int ndim0)
      : SingleVarStorage<T>(store, floor(maxchunk / ngrid) * ngrid,
                            name, dtype, units, scale_factor),
        ndim0(ndim0), nobs(0) {}

  ~LogbookStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    // if (this->bufferfill != 0)
    // {
    //   writechunk();
    // }
  }
};

struct ObservePrecip
/* satisfies ObserveLbks concept and
writes precipation data to zarr storage */
{
private:
  LogbookStorage &zarr;

public:  
  void observe_accumprecip(const std::shared_ptr<Logbook<double>> logbook) const
  {
    // std::vector<double> record = logbook.get_and_reset_record();
    // zarr.value_to_storage(record);
  }

  void operator()(const DetectorLogbooks &logbooks) const
  {
    observe_accumprecip(logbooks.accumprecip); 
  }
};

// double surface_precipitation(const GridBox &gbx, const double coord3lim)
// /* calculates mm of precipitation in a gridbox
// from mass of all superdrops which have
// radius >= rlim and coord3 <= zlim  */
// {
//   constexpr double rlim(40e-6 / dlc::R0);   // dimless minimum radius of precip

//   double precip(0.0);
//   const double area = gbx.area;
//   {
//     if r >= rlim && coord3 <= zlim
//     {
//       precip += alknca / area * COORD0 etc. // dimless
//     }
//   }
  
//   return precip;
// }

#endif // LOGBOOKS_INTOSTORE_HPP