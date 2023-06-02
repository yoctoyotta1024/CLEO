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

struct ObservePrecip
/* satisfies ObserveLbks concept and
writes precipation data to zarr storage */
{
private:
  TwoDStorage<double> &zarr;

public:  
  ObservePrecip(TwoDStorage<double> &zarr) : zarr(zarr)
  {
    zarr.is_name("accumprecip");
    zarr.is_dim1(0, "logbooktags");
  }

  void observe_accumprecip(const std::shared_ptr<Logbook<double>> logbook) const
  {
    std::vector<double> record = logbook -> get_and_reset_record(0.0);
    zarr.value_to_storage(record);
  }

  void operator()(const DetectorLogbooks &logbooks) const
  {
    observe_accumprecip(logbooks.accumprecip); 
    zarr.is_dim1(logbooks.accumprecip -> get_size(), "logbooktags");
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