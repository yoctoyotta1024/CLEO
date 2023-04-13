// Author: Clara Bayley
// File: "superdropmotion.hpp"
/* Header file for functions related to
moving superdroplets (updating their
coordinates according to equations of motion) */

#ifndef SDMMOTION_HPP
#define SDMMOTION_HPP

#include <concepts>

#include "./superdrop.hpp"
#include "./sedimentationmethod.hpp"
#include "./sdmprocess.hpp"
#include "./terminalvelocity.hpp"
#include "./thermostate.hpp"

template <typename P>
concept SuperdropMotion = requires(P p, const ThermoState &state,
                                   Superdrop &superdrop)
/* concept SuperdropMotion is all types that meet requirements
(constraints) of void function called "move_superdroplet"
which takes a ThermoState and Superdrop as arguments */
{
  {
    p(state, superdrop)
  };
};

struct NullMotion
{
  NullMotion(){};

  void operator()(const ThermoState &state,
                Superdrop &superdrop) const {}
};

template <VelocityFormula TerminalVelocity>
class MoveWithSedimentation
{
private:
  const double delt;                  // dimensionless delta time durign which motion occurs
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

public:
  void operator()(const ThermoState &state,
                Superdrop &superdrop) const
  {
    move_superdroplet(state, superdrop);
  }
};

#endif // SDMMOTION_HPP