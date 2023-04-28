// Author: Clara Bayley
// File: "prescribedmotion.hpp"
/* Header file for struct obeying 
sdmotion concept to update 
superdroplet positions 
according to a prescribed 2D flow */

#ifndef PRESCRIBEDMOTION_HPP
#define PRESCRIBEDMOTION_HPP

#include <functional>
#include <cmath>
#include <numbers>
#include <utility>

#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/terminalvelocity.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmotion.hpp"

class Prescribed2DFlow
/* Fixed 2D flow with constant density from
Arabas et al. 2015 with lengthscales
xlength = 2*pi*xtilda and zlength = pi*ztilda */
{
private:
  const double ztilda;
  const double xtilda;
  const double wamp;
  const std::function<double(ThermoState)> rhotilda; // function for normalised rho(z)

public:
  Prescribed2DFlow(const double zlength,
                   const double xlength,
                   const double wmax,
                   const std::function<double(ThermoState)> rhotilda);

  double prescribed_wvel(const ThermoState &state, const double zcoord,
                         const double xcoord) const;

  double prescribed_uvel(const ThermoState &state, const double zcoord,
                         const double xcoord) const;
};

class MoveWith2DPrescribedFlow
{
private:
  const int interval;                 // integer timestep for movement
  const double delt;                  // equivalent of interval as dimensionless time

  const Prescribed2DFlow flow2d; // method to get wvel and uvel from 2D flow field

  std::pair<double, double> predictor_corrector(const ThermoState &state,
                                                const double coord3,
                                                const double coord1) const;
  /* returns change in (z,x) coordinates = (delta3, delta1)
  obtained using predictor-corrector method and velocities
  calculated from a Prescribed2DFlow */
  
  std::pair<double, double> leapfrog(const ThermoState &state,
                              const double coord3,
                              const double coord1) const;
  /* returns change in (z,x) coordinates = (delta3, delta1)
  obtained using a simple leapfrog method and velocities
  calculated from a Prescribed2DFlow */

public:
  MoveWith2DPrescribedFlow(const int interval,
                      const std::function<double(int)> int2time,
                      const Prescribed2DFlow flow2d)
      : interval(interval),
        delt(int2time(interval)),
        flow2d(flow2d) {}

  MoveWith2DPrescribedFlow(const int interval,
                      const std::function<double(int)> int2time,
                      const double zlength,
                      const double xlength,
                      const double wmax,
                      const std::function<double(ThermoState)> rhotilda)
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


#endif // PRESCRIBEDMOTION_HPP