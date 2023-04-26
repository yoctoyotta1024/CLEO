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
#include <functional>
#include <stdexcept>
#include <cmath>
#include <numbers>

#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/terminalvelocity.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"

bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                    const unsigned int gbxindex,
                    const double delt,const double wvel,
                    const double uvel, const double vvel);

inline bool cfl_criterion(const double gridstep,
                          const double speed,
                          const double delt)
/* returns false if cfl criterion, C = speed*delt / gridstep, > 1 */
{
  return (speed*delt <= gridstep);
}

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
class NoInterpMoveWithSedimentation
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
  NoInterpMoveWithSedimentation(const int interval,
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
                                  Superdrop &drop) const;
  /* very crude method to forward timestep the velocity
  using the velocity from the gridbox thermostate, ie.
  without interpolation to the SD position and using 
  forward euler method to integrate dx/dt instead of
  a better method e.g. a predictor-corrector scheme */
};

class Prescribed2DFlow
{
  const double ztilda;
  const double xtilda;
  const double wamp;
  const double rhotilda;

  Prescribed2DFlow(const double zlength,
                   const double xlength,
                   const double wmax,
                   const double rhotilda)
      /* Fixed 2D flow with constant density from
      Arabas et al. 2015 with lengthscales
      xlength = 2*pi*xtilda and zlength = pi*ztilda */
      : ztilda(zlength / std::numbers::pi),         // 1/wavenumber given dimensionless wavelength
        xtilda(xlength / (2.0 * std::numbers::pi)), // 1/wavenumber given dimensionless wavelength
        wamp(2.0 * wmax),                           // amplitude of velocioty variations
        rhotilda(rhotilda) {} // normaliseddry air density

  double prescribed_wvel(const double zcoord, const double xcoord) const
  {
    return wamp / rhotilda * std::sin(zcoord/ztilda) * std::sin(xcoord/xtilda);
  } 

  double prescribed_uvel(const double zcoord, const double xcoord) const
  {
    const double u_amplitude = wamp / rhotilda * xtilda / ztilda; 
    return u_amplitude * std::cos(zcoord/ztilda) * std::cos(xcoord/xtilda);
  }
}

template <VelocityFormula TerminalVelocity>
class MoveWith2DFixedFlow
{
private:
  const int interval;                 // integer timestep for movement
  const double delt;                  // equivalent of interval as dimensionless time

  const Prescribed2DFlow flow2d;            // method to get wvel and uvel from 2D flow field

public:
  MoveWith2DFixedFlow(const int interval,
                      const std::function<double(int)> int2time,
                      const Prescribed2DFlow flow2d)
      : interval(interval),
        delt(int2time(interval)),
        flow2d(flow2d) {}

  MoveWith2DFixedFlow(const int interval,
                      const std::function<double(int)> int2time,
                      const double zlength,
                      const double xlength,
                      const double wmax,
                      const double rhotilda)
      : interval(interval),
        delt(int2time(interval)),
        flow2d(Prescribed2DFlow(zlength, xlength, wmax, rhotilda)) {}

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
                                  Superdrop &drop) const;
  /* Use predictor-corrector scheme from Grabowksi et al. 2018
  (similar to Arabas et al. 2015) to update a SD position.
  The velocity required for this scheme is determined
  from the PrescribedFlow2D instance */
};

#endif // SDMOTION_HPP