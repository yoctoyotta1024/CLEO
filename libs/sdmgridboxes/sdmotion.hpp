// Author: Clara Bayley
// File: "sdmotion.hpp"
/* Header file for functions related to
updatings superdroplets positions 
(updating their
coordinates according to equations of motion) */

#ifndef SDMOTION_HPP
#define SDMOTION_HPP

#include <concepts>
#include <limits>
#include <stdexcept>
#include <utility>
#include <functional>

#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/terminalvelocity.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"

bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                  const unsigned int gbxindex,
                  const double delta3, const double delta1,
                  const double delta2);
/* returns false if any of z, x or y (3,1,2) directions
  do not meet their cfl criterion. For each direction,
  Criterion is C = delta[X] / gridstep =< 1 where the
  gridstep is calculated from the gridbox boundaries map */

inline bool cfl_criterion(const double gridstep,
                          const double sdstep)
/* sdstep = change in superdroplet coordinate position.
returns false if cfl criterion, C = sdstep / gridstep, > 1 */
{
  return (sdstep <= gridstep);
}

struct WindsAtCoord
/* struct containing method to interpolate w, u, and v
wind velocities defined on faces on gridbox to a
superdroplet's (z,x,y) coordinates at
(coord3, coord1, coord2) */
{
private:
  const Maps4GridBoxes &gbxmaps;
  const ThermoState &state;
  const unsigned int gbxindex;
  const double coord3;
  const double coord1;
  const double coord2;

  double WindsAtCoord::interpolate_wind(const std::pair<double, double> bounds,
                                        const std::pair<double, double> vel,
                                        const double coord);
  /* Given [X = z,x or y] wind velocity component, vel, that is
  defined on the faces of a gridbox at {lower, upper} [X] bounds,
  return wind at [X] coord. Method is 'simple' linear interpolation
  from Grabowski et al. (2018) */

public:
  double interp_wvel() const;
  /* returns w wind velocity at z=coord3 for gridbox gbxindex */

  double interp_uvel() const;
  /* returns u wind velocity at x=coord1 for gridbox gbxindex */

  double interp_vvel() const;
  /* returns v wind velocity at y=coord2 for gridbox gbxindex */
};

template <typename M>
concept SdMotion = requires(M m, const int currenttimestep,
                            const GridBox &gbx,
                            const Maps4GridBoxes &gbxmaps,
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
    m.change_superdroplet_coords(gbxmaps, gbx, superdrop)
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

  void change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                                  const GridBox &gbx,
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

  void change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                                  const GridBox &gbx,
                                  Superdrop &drop) const
  /* very crude method to forward timestep the velocity
  using the velocity from the gridbox thermostate, ie.
  without interpolation to the SD position and using
  single step forward euler method to integrate dx/dt */
  {
    const WindsAtCoord winds{gbxmaps, gbx.state, gbx.gbxindex,
                             drop.coord3, drop.coord1, drop.coord2};



    const double delta3 = deltacoord(winds.interp_wvel() - terminalv(drop)); // w wind + terminal velocity
    const double delta1 = deltacoord(winds.interp_uvel());                   // u component of wind velocity
    const double delta2 = deltacoord(winds.interp_vvel());                   // v component of wind velocity (y=2)

    cfl_criteria(gbxmaps, gbx.gbxindex, delta3, delta1, delta2);

    drop.coord3 += delta3;
    drop.coord1 += delta1;
    drop.coord2 += delta2;
  }
};

#endif // SDMOTION_HPP