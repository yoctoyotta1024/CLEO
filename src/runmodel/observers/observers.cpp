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

void print_with_precision(const double datavalue,
                          const std::string lineend,
                          const int prec)
/* prints to terminal a datavalue followed
by "lineend" string with precision "prec" */
{
  std::cout << std::scientific
            << std::setprecision(prec)
            << datavalue << lineend;
}

void PrintObserver::observe_state(const std::vector<GridBox> &gridboxes) const
/* print t, kinematic data (p, temp, qv, qc) and total
number of sueprdrops to terminal */
{

  for (auto &gbx : gridboxes)
  {
    std::cout << "t=" << std::fixed
              << std::setprecision(printprec)
              << gbx.state.time*dlc::TIME0 << "s, y=[";
    print_with_precision(gbx.state.press, ", ", printprec);
    print_with_precision(gbx.state.temp, ", ", printprec);
    print_with_precision(gbx.state.qvap, ", ", printprec);
    print_with_precision(gbx.state.qcond, "], ", printprec);
    std::cout << "nsupers = " << gbx.span4SDsinGBx.size() << '\n';
  }
}