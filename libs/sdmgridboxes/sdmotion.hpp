// Author: Clara Bayley
// File: "sdmotion.hpp"
/* Header file for functions related to
updatings superdroplets positions 
(updating their
coordinates according to equations of motion) */

#ifndef SDMOTION_HPP
#define SDMOTION_HPP

#include <concepts>
#include <functional>
#include <limits>

#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/terminalvelocity.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "./gridbox.hpp"

bool cfl_criterion

template <typename M>
concept SdMotion = requires(M m, const int currenttimestep,
                            const GridBox &gbx,
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
    m.change_superdroplet_coords(gbx, superdrop)
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

  void change_superdroplet_coords(const GridBox &gbx,
                                  Superdrop &superdrop) const {}
};

template <VelocityFormula TerminalVelocity>
class MoveWithSedimentation
{
private:
  const int interval;                 // integer timestep for movement
  const double delt;                  // equivalent of interval as dimensionless time
  
  TerminalVelocity terminalv; // returns terminal velocity given a superdroplet

  double deltacoord(const double vel) const
  /* returns change in a coord given a velocity component 'vel' */
  {
    return vel * delt;
  }

public:
  MoveWithSedimentation(const int interval,
                        const std::function<double(int)> int2time,
                        const TerminalVelocity terminalv)
      : interval(interval),
        delt(int2time(interval)),
        terminalv(terminalv) {}

  int next_move(const int t) const
  {
    return ((t / interval) + 1) * interval;
  }

  bool on_move(const int t) const
  {
    return t % interval == 0;
  }

  void change_superdroplet_coords(const ThermoState &state,
                                  Superdrop &drop) const
  {
    const double vel3 = state.wvel - terminalv(drop); // w wind + terminal velocity
    drop.coord3 += deltacoord(vel3);
    
    const double vel1 = state.uvel; // u component of wind velocity
    drop.coord1 += deltacoord(vel1);

    const double vel2 = state.vvel; // v component of wind velocity (y=2)
    drop.coord2 += deltacoord(vel2);
  }
};

#endif // SDMOTION_HPP