// Author: Clara Bayley
// File: "sdmotion.hpp"
/* Header file for functions related to
moving superdroplets (updating their
coordinates according to equations of motion) */

#ifndef SDMOTION_HPP
#define SDMOTION_HPP

#include <concepts>
#include <functional>
#include <limits>

#include "./superdrop.hpp"
#include "./sedimentationmethod.hpp"
#include "./sdmprocess.hpp"
#include "./terminalvelocity.hpp"
#include "./thermostate.hpp"

template <typename M>
concept SdMotion = requires(M m, const int currenttimestep,
                            const ThermoState &state,
                            Superdrop &superdrop)
/* concept SdMotion is all types that meet requirements
(constraints) of void function called "move_superdroplet"
which takes a ThermoState and Superdrop as arguments */
{
  {
    m.next_move(currenttimestep)
    } -> std::convertible_to<int>;
  {
    m.on_move(currenttimestep)
    } -> std::convertible_to<bool>;
  {
    m.change_superdroplet_coords(state, superdrop)
  };
};

struct NullMotion
{
  int next_move(const int currenttimestep) const
  {
    return std::numeric_limits<int>::max();
  }

  bool on_move(const int currenttimestep) const
  {
    return false;
  }

  void change_superdroplet_coords(const ThermoState &state,
                                    Superdrop &superdrop) const {}
};

template <VelocityFormula TerminalVelocity>
class MoveWithSedimentation
{
private:
  const int interval;                 // integer timestep for movement
  const double delt;                  // equivalent of interval as dimensionless time
  
  TerminalVelocity terminal_velocity; // returns terminal velocity given a superdroplet

public:
  MoveWithSedimentation(const int interval,
                        const std::function<double(int)> int2time,
                        const TerminalVelocity v)
      : interval(interval),
        delt(int2time(interval)),
        terminal_velocity(v) {}

  int next_move(const int t) const
  {
    return ((t / interval) + 1) * interval;
  }

  bool on_move(const int t) const
  {
    return t % interval == 0;
  }

  void change_superdroplet_coords(const ThermoState &state,
                                  Superdrop &superdrop) const
  {
    // const double vel3 = state.wvel; // w component of wind velocity (z=3)
    // const double vel1 = state.uvel; // u component of wind velocity (x=1)
    // const double vel2 = state.vvel; // v component of wind velocity (y=2)
  }
};

#endif // SDMOTION_HPP