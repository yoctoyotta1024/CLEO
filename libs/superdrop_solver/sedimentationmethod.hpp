// Author: Clara Bayley
// File: sedimentationmethod.hpp
/* Header file for method for
sedimentation of superdroplets */

#ifndef SEDIMENTATIONMETHOD_HPP
#define SEDIMENTATIONMETHOD_HPP

#include <span>
#include <random>

#include "./terminalvelocity.hpp"
#include "./superdrop.hpp"
#include "./thermostate.hpp"

template <VelocityFormula TerminalVelocity>
class SedimentationMethod
/* class for implementing superdroplet sedimentation in SDM */
{
private:
  const double delt;
  TerminalVelocity terminal_velocity;
  /* returns terminal velocity of a given superdroplet */

  void sediment_drop(Superdrop &drop) const
  /* enacts sedimentation by changing coord3
  (z coord) of superdroplet */
  {
    const double vertical_velocity = -1.0 * terminal_velocity(drop); 
    drop.coord3 += vertical_velocity * delt;
  }

public:
  SedimentationMethod(const double delt, TerminalVelocity v)
      : delt(delt),
        terminal_velocity(v) {}

  void sediment_superdroplets(std::span<SuperdropWithGridbox> span4SDsinGBx) const
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
                         std::span<SuperdropWithGridbox> span4SDsinGBx,
                         ThermoState &state,
                         std::mt19937 &gen) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call sediment_superdroplets. (*hint* run_step
  usually found within a type that satisfies the SdmProcess concept) */
  {
    sediment_superdroplets(span4SDsinGBx);
  }
};

#endif // SEDIMENTATIONMETHOD_HPP