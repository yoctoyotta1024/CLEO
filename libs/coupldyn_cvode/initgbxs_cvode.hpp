/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: initgbxs_cvode.hpp
 * Project: coupldyn_cvode
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 14th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for griboxes' initial conditions
 * for CLEO SDM (e.g. thermodynamics)
 * when coupled to cvode thermodynamics solver.
 * Struct can be used by InitConds as
 * GbxInitConds type
 */

#ifndef INITGBXS_CVODE_HPP
#define INITGBXS_CVODE_HPP

#include <vector>
#include <utility>

#include "../cleoconstants.hpp"
#include "./differentialfuncs.hpp"
#include "initialise/config.hpp"

namespace dlc = dimless_constants;

struct InitGbxsCvode
/* struct containing functions which return data
for the initial conditions needed to create
gridboxes e.g. via the create_gbxs function where
all gridboxes are initially the same */
{
private:
  size_t ngbxs;
  double press_i;   // initial pressure [Pa]
  double temp_i;    // initial temperature [Pa]
  double relh_init; // initial relative humidity (%)
  double qcond_i;   // initial liquid water mixing ratio

public:
  InitGbxsCvode(const Config &config)
      : ngbxs(config.ngbxs),
        press_i(config.P_INIT / dlc::P0),
        temp_i(config.TEMP_INIT / dlc::TEMP0),
        relh_init(config.relh_init),
        qcond_i(0.0) {}

  size_t get_ngbxs() const { return ngbxs; }

  std::vector<double> press() const
  /* pressure for all gbxs is same initial
  (dimless) 'press_i' given by PRESS_INIT */
  {
    return std::vector<double>(ngbxs, press_i);
  }

  std::vector<double> temp() const
  /* temperature for all gbxs is same initial
  (dimless) 'temp' given by TEMP_INIT */
  {
    return std::vector<double>(ngbxs, temp_i);
  }

  std::vector<double> qvap() const
  /* vapour mass mixing ratio, 'qvap' for all
  gbxs is same as given by temp_i, press_i and
  relh_init */
  {
    const auto psat = double{cvode_saturationpressure(temp_i)};
    const auto vapp = double{psat * relh_init / 100.0}; // initial vapour pressure
    const auto qvap_i = double{cvode_massmixingratio(vapp, press_i)};

    return std::vector<double>(ngbxs, qvap_i);
  }

  std::vector<double> qcond() const
  /* liquid mass mixing ratio, 'qcond'
  for all gbxs is qcond_i (= 0.0 by default) */
  {
    return std::vector<double>(ngbxs, qcond_i);
  }

  std::vector<std::pair<double, double>> wvel() const
  /* wind velocity for al gbxs = 0.0 */
  {
    std::pair<double, double> wvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, wvel_i);
  }

  std::vector<std::pair<double, double>> uvel() const
  /* wind velocity for al gbxs = 0.0 */
  {
    std::pair<double, double> uvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, uvel_i);
  }

  std::vector<std::pair<double, double>> vvel() const
  /* wind velocity for al gbxs = 0.0 */
  {
    std::pair<double, double> vvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, vvel_i);
  }
};

#endif // INITGBXS_CVODE_HPP
