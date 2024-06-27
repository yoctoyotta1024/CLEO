/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: optional_config_params.cpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions involved in reading optional configuration parameters from a config file.
 */

#include "initialise/optional_config_params.hpp"

/* read configuration file given by config_filename to set members of required configuration */
OptionalConfigParams::OptionalConfigParams(const std::filesystem::path config_filename) {
  const YAML::Node config = YAML::LoadFile(std::string{config_filename});

  if (config["microphysics"]) {
    set_microphysics(config);
  }

  if (config["initsupers"]) {
    set_initsupers(config);
  }

  if (config["coupled_dynamics"]) {
    set_coupled_dynamics(config);
  }

  if (config["boundary_conditions"]) {
    set_boundary_conditions(config);
  }
}

void OptionalConfigParams::set_microphysics(const YAML::Node &config) {
  const YAML::Node node = config["microphysics"];

  if (node["condensation"]) {
    condensation.set_params(config);
    condensation.print_params();
  }

  if (node["breakup"]) {
    breakup.set_params(config);
    breakup.print_params();
  }
}

void OptionalConfigParams::set_initsupers(const YAML::Node &config) {
  const auto type = config["initsupers"]["type"].as<std::string>();

  if (type == "frombinary") {
    initsupersfrombinary.set_params(config);
    initsupersfrombinary.print_params();
  } else {
    throw std::invalid_argument("unknown initsupers 'type': " + type);
  }
}

void OptionalConfigParams::set_coupled_dynamics(const YAML::Node &config) {
  const auto type = config["coupled_dynamics"]["type"].as<std::string>();

  if (type == "fromfile") {
    fromfiledynamics.set_params(config);
    fromfiledynamics.print_params();
  } else if (type == "cvode") {
    cvodedynamics.set_params(config);
    cvodedynamics.print_params();
  } else if (type == "yac") {
    yac_dynamics.set_params(config);
    yac_dynamics.print_params();
  } else {
    throw std::invalid_argument("unknown coupled_dynamics 'type': " + type);
  }
}

void OptionalConfigParams::set_boundary_conditions(const YAML::Node &config) {
  const auto type = config["boundary_conditions"]["type"].as<std::string>();

  if (type == "addsupersatdomaintop") {
    addsupersatdomaintop.set_params(config);
    addsupersatdomaintop.print_params();
  } else {
    throw std::invalid_argument("unknown boundary_conditions 'type': " + type);
  }
}

void OptionalConfigParams::CondensationParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["microphysics"]["condensation"];

  do_alter_thermo = node["do_alter_thermo"].as<bool>();
  maxniters = node["maxniters"].as<size_t>();
  MINSUBTSTEP = node["MINSUBTSTEP"].as<double>();
  rtol = node["rtol"].as<double>();
  atol = node["atol"].as<double>();
}

void OptionalConfigParams::CondensationParams::print_params() const {
  std::cout << "\n-------- Condensation Configuration Parameters --------------"
            << "\ndo_alter_thermo: " << do_alter_thermo << "\nmaxniters: " << maxniters
            << "\nMINSUBSTEP: " << MINSUBTSTEP << "\nrtol: " << rtol << "\natol: " << atol
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::BreakupParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["microphysics"]["breakup"]["constnfrags"];
  constnfrags.nfrags = node["nfrags"].as<double>();
}

void OptionalConfigParams::BreakupParams::print_params() const {
  std::cout << "\n-------- Breakup Configuration Parameters --------------"
            << "\nConstNFrags nfrags: " << constnfrags.nfrags
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::InitSupersFromBinaryParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["initsupers"];

  assert((node["type"].as<std::string>() == "frombinary"));

  maxnsupers = config["domain"]["maxnsupers"].as<size_t>();
  initsupers_filename = std::filesystem::path(node["initsupers_filename"].as<std::string>());
  nspacedims = config["domain"]["nspacedims"].as<unsigned int>();
  if (node["initnsupers"]) {
    initnsupers = node["initnsupers"].as<size_t>();
  } else {
    initnsupers = maxnsupers;
  }
}

void OptionalConfigParams::InitSupersFromBinaryParams::print_params() const {
  std::cout << "\n-------- InitSupersFromBinary Configuration Parameters --------------"
            << "\nmaxnsupers: " << maxnsupers << "\nnspacedims: " << nspacedims
            << "\ninitsupers_filename: " << initsupers_filename << "\ninitnsupers: " << initnsupers
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::FromFileDynamicsParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["coupled_dynamics"];

  assert((node["type"].as<std::string>() == "fromfile"));

  /* convert string to std::filesystem::path type */
  auto fspath_from_yaml = [&node](const std::string &key) {
    return std::filesystem::path(node[key].as<std::string>());
  };

  nspacedims = config["domain"]["nspacedims"].as<unsigned int>();
  press = fspath_from_yaml("press");
  temp = fspath_from_yaml("temp");
  qvap = fspath_from_yaml("qvap");
  qcond = fspath_from_yaml("qcond");
  switch (nspacedims) {
    case 3:  // 3-D model
      vvel = fspath_from_yaml("vvel");
      [[fallthrough]];
    case 2:  // 3-D or 2-D model
      uvel = fspath_from_yaml("uvel");
      [[fallthrough]];
    case 1:  // 3-D, 2-D or 1-D model
      wvel = fspath_from_yaml("wvel");
  }
}

void OptionalConfigParams::FromFileDynamicsParams::print_params() const {
  std::cout << "\n-------- FromFileDynamics Configuration Parameters --------------"
            << "\nnspacedims: " << nspacedims << "\npress: " << press << "\ntemp: " << temp
            << "\nqvap: " << qvap << "\nqcond: " << qcond << "\nwvel: " << wvel
            << "\nuvel: " << uvel << "\nvvel: " << vvel
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::CvodeDynamicsParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["coupled_dynamics"];

  assert((node["type"].as<std::string>() == "cvode"));

  ngbxs = config["domain"]["ngbxs"].as<unsigned int>();
  P_init = node["P_init"].as<double>();
  TEMP_init = node["TEMP_init"].as<double>();
  relh_init = node["relh_init"].as<double>();
  W_avg = node["W_avg"].as<double>();
  TAU_half = node["TAU_half"].as<double>();
  rtol = node["rtol"].as<double>();
  atol = node["atol"].as<double>();
}

void OptionalConfigParams::CvodeDynamicsParams::print_params() const {
  std::cout << "\n-------- CvodeDynamics Configuration Parameters --------------"
            << "\nngbxs: " << ngbxs << "\nP_init: " << P_init << "\nTEMP_init: " << TEMP_init
            << "\nrelh_init: " << relh_init << "\nW_avg: " << W_avg << "\nTAU_half: " << TAU_half
            << "\nrtol: " << rtol << "\natol: " << atol
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::YacDynamicsParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["coupled_dynamics"];

  assert((node["type"].as<std::string>() == "yac"));

  if (node["lower_longitude"])
    lower_longitude = node["lower_longitude"].as<double>();

  if (node["upper_longitude"])
    upper_longitude = node["upper_longitude"].as<double>();

  if (node["lower_latitude"])
    lower_latitude = node["lower_latitude"].as<double>();

  if (node["upper_latitude"])
    upper_latitude = node["upper_latitude"].as<double>();
}

void OptionalConfigParams::YacDynamicsParams::print_params() const {
  std::cout << "\n-------- YacDynamics Configuration Parameters --------------"
            << "\nlower_longitude: " << lower_longitude << "\nupper_longitude: " << upper_longitude
            << "\nlower_latitude: " << lower_latitude << "\nupper_latitude: " << upper_latitude
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::AddSupersAtDomainTopParams::set_params(const YAML::Node &config) {
  const YAML::Node node = config["boundary_conditions"];

  if (config["initsupers"] && config["initsupers"]["initnsupers"]) {
    initnsupers = config["initsupers"]["initnsupers"].as<size_t>();
  } else {
    initnsupers = config["domain"]["maxnsupers"].as<size_t>();
  }
  newnsupers = node["newnsupers"].as<size_t>();
  COORD3LIM = node["COORD3LIM"].as<double>();
  DRYRADIUS = node["DRYRADIUS"].as<double>();
  MINRADIUS = node["MINRADIUS"].as<double>();
  MAXRADIUS = node["MAXRADIUS"].as<double>();
  NUMCONC_a = node["NUMCONC_a"].as<double>();
  GEOMEAN_a = node["GEOMEAN_a"].as<double>();
  geosigma_a = node["geosigma_a"].as<double>();
  NUMCONC_b = node["NUMCONC_b"].as<double>();
  GEOMEAN_b = node["GEOMEAN_b"].as<double>();
  geosigma_b = node["geosigma_b"].as<double>();
}

void OptionalConfigParams::AddSupersAtDomainTopParams::print_params() const {
  std::cout << "\n-------- AddSupersAtDomainTop Configuration Parameters --------------"
            << "\ninitnsupers: " << initnsupers << "\nnewnsupers: " << newnsupers
            << "\nCOORD3LIM: " << COORD3LIM << "\nDRYRADIUS: " << DRYRADIUS
            << "\nMINRADIUS: " << MINRADIUS << "\nMAXRADIUS: " << MAXRADIUS
            << "\nNUMCONC_a: " << NUMCONC_a << "\nGEOMEAN_a: " << GEOMEAN_a
            << "\ngeosigma_a: " << geosigma_a << "\nNUMCONC_b: " << NUMCONC_b
            << "\nGEOMEAN_b: " << GEOMEAN_b << "\ngeosigma_b: " << geosigma_b
            << "\n---------------------------------------------------------\n";
}
