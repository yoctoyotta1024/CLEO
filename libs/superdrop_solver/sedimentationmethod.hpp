// Author: Clara Bayley
// File: sedimentationmethod.hpp
/* Header file for method for
sedimentation of superdroplets */

#ifndef SEDIMENTATIONMETHOD_HPP
#define SEDIMENTATIONMETHOD_HPP

#include <span>
#include <random>
#include <functional>
#include <concepts>

#include "./terminalvelocity.hpp"
#include "./superdrop.hpp"
#include "./thermostate.hpp"
#include "./sdmprocess.hpp"

template <VelocityFormula TerminalVelocity>
class SedimentationMethod
/* class for implementing superdroplet sedimentation in SDM */
{
private:
  const double delt;
  const TerminalVelocity terminal_velocity; // returns terminal velocity given a superdroplet

  void sediment_drop(Superdrop &drop) const
  /* enacts sedimentation by changing coord3
  (z coord) of superdroplet */
  {
    drop.coord3 -= terminal_velocity(drop) * delt;
  }

public:
  SedimentationMethod(const double delt, const TerminalVelocity v)
      : delt(delt),
        terminal_velocity(v) {}


  void sediment_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx) const
  /* sediment all superdroplets stored in some span of contigous memory.
  Here the span points to some subsection of a vector containing
  superdroplet in gridbox instances 'SDinGBx' */
  {
    for (auto &SDinGBx : span4SDsinGBx)
    {
      sediment_drop(SDinGBx.superdrop);
    }
  }

  inline void operator()(const int currenttimestep,
                         std::span<SuperdropWithGbxindex> span4SDsinGBx,
                         ThermoState &state,
                         std::mt19937 &gen) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call sediment_superdroplets. (*hint* run_step
  usually found within a type that satisfies the SdmProcess concept) */
  {
    sediment_superdroplets(span4SDsinGBx);
  }
};

template <VelocityFormula TerminalVelocity>
SdmProcess auto SedimentationProcess(const int interval,
                                     const std::function<double(int)> int2time,
                                     const TerminalVelocity v)
/* constructs SdmProcess for sedimentation with constant timestep 'interval'
given a function to convert the interval to a (dimensionless) time
and a terminal velocity formula */
{
  const double dimlesststep = int2time(interval);
  return ConstTstepProcess{interval, SedimentationMethod(dimlesststep, v)};
}

#endif // SEDIMENTATIONMETHOD_HPP