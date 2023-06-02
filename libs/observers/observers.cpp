// Author: Clara Bayley
// File: "observers.cpp"
/* Observer Concept and related structures 
for various ways of observing gridboxes 
and logbooks of the superdroplet model. 
An example of an observer is printing some data
from a gridbox's thermostate to the terminal */

#include "observers.hpp"

void print_thermostate_with_precision(const ThermoState &state)
/* prints to terminal a datavalue followed
by "lineend" string with precision "prec" */
{
  constexpr int printprec(4); // precision to print data with
 
  std::cout << std::scientific
            << std::setprecision(printprec)
            << "[P,T,qv,qc]=["
            << state.press << ", "
            << state.temp << ", "
            << state.qvap << ", "
            << state.qcond << "]\n";
}

void PrintObserver::observe_gridboxes(const size_t ngbxs,
                                  const Kokkos::View<GridBox *> h_gridboxes) const
/* print t, kinematic data (p, temp, qv, qc) and total
number of sueprdrops to terminal */
{
  constexpr int printprec(4); // precision to print data with

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    const auto &gbx = h_gridboxes(ii);
    std::cout << "GBx " << gbx.gbxindex << ": "
              << "t=" << std::fixed
              << std::setprecision(printprec)
              << gbx.state.time * dlc::TIME0 << "s, "
              << "nsupers=" << gbx.span4SDsinGBx.size() << ", ";

    print_thermostate_with_precision(gbx.state);
  }
}

void PrintObserver::observe_logbooks(const DetectorLogbooks &logbooks) const
{
  size_t size(logbooks.accumprecip->get_size());
  std::cout << "no. logbook entries: "
               "[accumprecip: " << size << "]\n";
}