// Author: Clara Bayley
// File: "logbooks_intostore.hpp"
/* structures obeying the ObserveLbks
concept for various ways of observing
logbooks which ends up writing data
into a (zarr) store on disk */

#ifndef LOGBOOKS_INTOSTORE_HPP
#define LOGBOOKS_INTOSTORE_HPP

struct ObserveAccumPrecip
/* satisfies Observer concept and
prints out details about gridboxes'
thermodynamic states and superdroplets */
{
  void printprecip(
      const std::shared_ptr<Logbook<double>> logbook) const
  {
    double totaccumpp(0.0);
    for (size_t idx = 0; idx < logbook  -> get_size(); ++idx)
    {
      totaccumpp += logbook ->get_from_record(idx);
    }
    
    constexpr int printprec(4); // precision to print data with
    std::cout << std::scientific
            << std::setprecision(printprec)
            << "tot accum. precip = "
            << totaccumpp << '\n';
  }

  void operator()(const DetectorLogbooks &logbooks) const
  {
    printprecip(logbooks.accumprecip); 
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