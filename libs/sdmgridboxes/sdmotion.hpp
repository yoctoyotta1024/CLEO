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

  double interpolate_wind(const std::pair<double, double> bounds,
                          const std::pair<double, double> vel,
                          const double coord) const;
  /* Given [X = z,x or y] wind velocity component, vel, that is
  defined on the faces of a gridbox at {lower, upper} [X] bounds,
  return wind at [X] coord. Method is 'simple' linear interpolation
  from Grabowski et al. (2018) */

public:
  const Maps4GridBoxes &gbxmaps;
  const ThermoState &state;
  unsigned int gbxindex;
  double coord3;
  double coord1;
  double coord2;

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
  const int interval; // integer timestep for movement
  const double delt;  // equivalent of interval as dimensionless time

  TerminalVelocity terminalv; // returns terminal velocity given a superdroplet

  struct Deltas
  {
    double delta3;
    double delta1;
    double delta2;
  };

  Deltas predictor_corrector(const Maps4GridBoxes &gbxmaps,
                             const GridBox &gbx,
                             const Superdrop &drop) const
  {
    const double terminal = terminalv(drop);

    WindsAtCoord winds{gbxmaps, gbx.state, gbx.gbxindex,
                       drop.coord3, drop.coord1, drop.coord2};

    /* corrector velocities based on predicted coords */
    const double vel3 = winds.interp_wvel() - terminal;
    const double vel1 = winds.interp_uvel();
    const double vel2 = winds.interp_vvel();

    /* predictor coords given velocity at previous coords */
    winds.coord3 += vel3 * delt; // move by w wind + terminal velocity
    winds.coord1 += vel1 * delt; // move by u wind
    winds.coord2 += vel2 * delt; // move by v wind

    /* corrector velocities based on predicted coords */
    const double corrvel3 = winds.interp_wvel() - terminal;
    const double corrvel1 = winds.interp_uvel();
    const double corrvel2 = winds.interp_vvel();

    /* predicted-corrected change to superdrop coords */
    const double delta3((vel3 + corrvel3) * (delt / 2));
    const double delta1((vel1 + corrvel1) * (delt / 2));
    const double delta2((vel2 + corrvel2) * (delt / 2));

    return Deltas{delta3, delta1, delta2};
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
  /* Uses predictor-corrector method to forward timestep 
  a superdroplet's coordinates using the interpolated 
  wind velocity from a gridbox's thermostate */
  {
    /* USe predictor-corrector method to get change in SD coords */
    Deltas d{predictor_corrector(gbxmaps, gbx, drop)};

    /* CFL check predicted change to SD coords */
    cfl_criteria(gbxmaps, gbx.gbxindex, d.delta3, d.delta1, d.delta2);

    /* update SD coords */
    drop.coord3 += d.delta3;
    drop.coord1 += d.delta1;
    drop.coord2 += d.delta2;
  }
};

inline double WindsAtCoord::interpolate_wind(const std::pair<double, double> bounds,
                                      const std::pair<double, double> vel,
                                      const double coord) const
/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018) */
{
  const double alpha((coord - bounds.first) / (bounds.second - bounds.first));

  const double interpolated_vel(alpha*vel.second + (1-alpha)*vel.first);
  
  return interpolated_vel;
} 

inline double WindsAtCoord::interp_wvel() const
/* returns w wind velocity at z=coord3 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_z(gbxindex),
                          state.wvel, coord3);
}

inline double WindsAtCoord::interp_uvel() const
/* returns u wind velocity at x=coord1 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_x(gbxindex),
                          state.uvel, coord1); 
}

inline double WindsAtCoord::interp_vvel() const
/* returns v wind velocity at y=coord2 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_y(gbxindex),
                          state.vvel, coord2); 
}

#endif // SDMOTION_HPP