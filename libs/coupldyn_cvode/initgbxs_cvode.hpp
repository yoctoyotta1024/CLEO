/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: initgbxs_cvode.hpp
 * Project: coupldyn_cvode
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for griboxes' initial conditions for CLEO SDM (e.g. thermodynamics)
 * when coupled to cvode thermodynamics solver. Struct can be used by InitConds as
 * GbxInitConds type
 */

#ifndef LIBS_COUPLDYN_CVODE_INITGBXS_CVODE_HPP_
#define LIBS_COUPLDYN_CVODE_INITGBXS_CVODE_HPP_

#include <utility>
#include <vector>

#include "../cleoconstants.hpp"
#include "configuration/optional_config_params.hpp"
#include "coupldyn_cvode/differentialfuncs.hpp"

namespace dlc = dimless_constants;

/* struct containing functions which return data
for the initial conditions needed to create
gridboxes e.g. via the create_gbxs function where
all gridboxes are initially the same */
struct InitGbxsCvode {
 private:
  size_t ngbxs;
  double press_i;    // initial pressure [Pa]
  double temp_i;     // initial temperature [Pa]
  double relh_init;  // initial relative humidity (%)
  double qcond_i;    // initial liquid water mixing ratio

 public:
  explicit InitGbxsCvode(const OptionalConfigParams::CvodeDynamicsParams &config)
      : ngbxs(config.ngbxs),
        press_i(config.P_init / dlc::P0),
        temp_i(config.TEMP_init / dlc::TEMP0),
        relh_init(config.relh_init),
        qcond_i(0.0) {}

  size_t get_ngbxs() const { return ngbxs; }

  /* pressure for all gbxs is same initial
  (dimless) 'press_i' given by PRESS_init */
  std::vector<double> press() const { return std::vector<double>(ngbxs, press_i); }

  /* temperature for all gbxs is same initial
  (dimless) 'temp' given by TEMP_init */
  std::vector<double> temp() const { return std::vector<double>(ngbxs, temp_i); }

  /* vapour mass mixing ratio, 'qvap' for all
  gbxs is same as given by temp_i, press_i and
  relh_init */
  std::vector<double> qvap() const {
    const auto psat = double{cvode_saturationpressure(temp_i)};
    const auto vapp = double{psat * relh_init / 100.0};  // initial vapour pressure
    const auto qvap_i = double{cvode_massmixingratio(vapp, press_i)};

    return std::vector<double>(ngbxs, qvap_i);
  }

  /* liquid mass mixing ratio, 'qcond'
  for all gbxs is qcond_i (= 0.0 by default) */
  std::vector<double> qcond() const { return std::vector<double>(ngbxs, qcond_i); }

  /* wind velocity for al gbxs = 0.0 */
  std::vector<std::pair<double, double>> wvel() const {
    std::pair<double, double> wvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, wvel_i);
  }

  /* wind velocity for al gbxs = 0.0 */
  std::vector<std::pair<double, double>> uvel() const {
    std::pair<double, double> uvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, uvel_i);
  }

  /* wind velocity for al gbxs = 0.0 */
  std::vector<std::pair<double, double>> vvel() const {
    std::pair<double, double> vvel_i(0.0, 0.0);
    return std::vector<std::pair<double, double>>(ngbxs, vvel_i);
  }
};

#endif  // LIBS_COUPLDYN_CVODE_INITGBXS_CVODE_HPP_
