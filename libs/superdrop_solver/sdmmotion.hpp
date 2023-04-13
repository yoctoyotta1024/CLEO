// Author: Clara Bayley
// File: "sdmmotion.hpp"
/* Header file for functions related to
moving superdroplets (updating their
coordinates) */

#ifndef SDMMOTION_HPP
#define SDMMOTION_HPP

#include <concepts>

#include "./superdrop.hpp"
#include "./sedimentationmethod.hpp"
#include "./sdmprocess.hpp"
#include "./terminalvelocity.hpp"
#include "./thermostate.hpp"

template <typename P>
concept SdmMotion = requires(P p, const ThermoState &state,
                             Superdrop &superdrop)
/* concept SdmMotion is all types that meet requirements
(constraints) of void function called "move_superdroplet"
which takes a ThermoState and Superdrop as arguments */
{
  {
    p.move_superdroplet(state, superdrop)
  };
};

struct NullMovement
{
  NullMovement(){};

  void move_superdroplet(const ThermoState &state,
                         Superdrop &superdrop) const {}
};

template <VelocityFormula TerminalVelocity>
struct MoveWithSedimentation
{
  const double delt; //dimensionless delta time durign which motion occurs
  TerminalVelocity terminal_velocity; // returns terminal velocity given a superdroplet
  
  MoveWithSedimentation(const double delt, TerminalVelocity v)
      : delt(delt), terminal_velocity(v){};

  void move_superdroplet(const ThermoState &state,
                         Superdrop &superdrop) const
  {
    // const double vel3 = state.wvel; // w component of wind velocity (z=3)
    // const double vel1 = state.uvel; // u component of wind velocity (x=1)
    // const double vel2 = state.vvel; // v component of wind velocity (y=2)
  }
};

#endif // SDMMOTION_HPP