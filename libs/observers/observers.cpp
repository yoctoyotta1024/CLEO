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

double PrintObserver::
    sum_surfpp(const std::shared_ptr<Logbook<double>> &logbook) const
{
  double totpp(0.0);
  for (size_t idx = 0; idx < logbook->get_size(); ++idx)
  {
    totpp += logbook->get_entry(idx);
  }

  return totpp * dlc::MASS0grams; // [grams]
}

void PrintObserver::observe_logbooks(const DetectorLogbooks &logbooks) const
{
  constexpr int printprec(4); // precision to print data with

  size_t size(logbooks.surfpp->get_size());
  double totpp(sum_surfpp(logbooks.surfpp));

  std::cout << "logbooks: [surfpp (" << size << "), "
            << std::scientific << std::setprecision(printprec)
            << totpp << "g]\n";
}