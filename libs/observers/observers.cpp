// Author: Clara Bayley
// File: "observers.cpp"
/* functionality for observer structures. Each
observer is a way of observing a gridbox of the
superdroplet model, gridbox contained the
thermodynamic state and the vector of
superdroplets' state. Observation is
for example printing some thermodynamic data
to terminal or writing them to a csv file */

#include "observers.hpp"

void print_thermostate_with_precision(const ThermoState &state,
                                      const int prec)
/* prints to terminal a datavalue followed
by "lineend" string with precision "prec" */
{

  std::cout << std::scientific
            << std::setprecision(prec)
            << "[P,T,qv,qc]=["
            << state.press << ", "
            << state.temp << ", "
            << state.qvap << ", "
            << state.qcond << "]\n";
}

void PrintObserver::observe_state(const Kokkos::View<GridBox*> h_gridboxes) const
/* print t, kinematic data (p, temp, qv, qc) and total
number of sueprdrops to terminal */
{
  const size_t Ngrid = h_gridboxes.size();
  for (size_t ii(0); ii < Ngrid; ++ii)
  {
    const auto &gbx = h_gridboxes(ii);
    std::cout << "GBx " << gbx.gbxindex << ": "
              << "t=" << std::fixed
              << std::setprecision(printprec)
              << gbx.state.time * dlc::TIME0 << "s, "
              << "nsupers=" << gbx.span4SDsinGBx.size() << ", ";
    print_thermostate_with_precision(gbx.state, printprec);
  }
}