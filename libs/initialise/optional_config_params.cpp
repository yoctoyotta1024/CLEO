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
 * Last Modified: Thursday 18th April 2024
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
  const YAML::Node yaml = config["microphysics"];

  if (yaml["condensation"]) {
    condensation.set_params(config);
    condensation.print_params();
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
  const YAML::Node yaml = config["microphysics"]["condensation"];

  do_alter_thermo = yaml["do_alter_thermo"].as<bool>();
  niters = yaml["niters"].as<unsigned int>();
  SUBTSTEP = yaml["SUBTSTEP"].as<double>();
  rtol = yaml["rtol"].as<double>();
  atol = yaml["atol"].as<double>();
}

void OptionalConfigParams::CondensationParams::print_params() const {
  std::cout << "\n-------- Condensation Configuration Parameters --------------"
            << "\ndo_alter_thermo: " << do_alter_thermo << "\nniters: " << niters
            << "\nSUBSTEP: " << SUBTSTEP << "\nrtol: " << rtol << "\natol: " << atol
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::InitSupersFromBinaryParams::set_params(const YAML::Node &config) {
  const YAML::Node yaml = config["initsupers"];

  assert((yaml["type"].as<std::string>() == "frombinary"));

  initsupers_filename = std::filesystem::path(yaml["initsupers_filename"].as<std::string>());
  initnsupers = yaml["initnsupers"].as<size_t>();
  nspacedims = config["domain"]["nspacedims"].as<unsigned int>();
}

void OptionalConfigParams::InitSupersFromBinaryParams::print_params() const {
  std::cout << "\n-------- InitSupersFromBinary Configuration Parameters --------------"
            << "\nnspacedims: " << nspacedims << "\ninitsupers_filename: " << initsupers_filename
            << "\ninitnsupers: " << initnsupers
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::FromFileDynamicsParams::set_params(const YAML::Node &config) {
  const YAML::Node yaml = config["coupled_dynamics"];

  assert((yaml["type"].as<std::string>() == "fromfile"));

  /* convert string to std::filesystem::path type */
  auto fspath_from_yaml = [&yaml](const std::string &key) {
    return std::filesystem::path(yaml[key].as<std::string>());
  };

  nspacedims = config["domain"]["nspacedims"].as<unsigned int>();
  press = fspath_from_yaml("press");
  temp = fspath_from_yaml("temp");
  qvap = fspath_from_yaml("qvap");
  qcond = fspath_from_yaml("qcond");
  switch (nspacedims) {
    case 3:  // 3-D model
      vvel = fspath_from_yaml("vvel");
    case 2:  // 3-D or 2-D model
      uvel = fspath_from_yaml("uvel");
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
  const YAML::Node yaml = config["coupled_dynamics"];

  assert((yaml["type"].as<std::string>() == "cvode"));

  ngbxs = config["domain"]["ngbxs"].as<unsigned int>();
  P_init = yaml["P_init"].as<double>();
  TEMP_init = yaml["TEMP_init"].as<double>();
  relh_init = yaml["relh_init"].as<double>();
  W_avg = yaml["W_avg"].as<double>();
  TAU_half = yaml["TAU_half"].as<double>();
  rtol = yaml["rtol"].as<double>();
  atol = yaml["atol"].as<double>();
}

void OptionalConfigParams::CvodeDynamicsParams::print_params() const {
  std::cout << "\n-------- CvodeDynamics Configuration Parameters --------------"
            << "\nngbxs: " << ngbxs << "\nP_init: " << P_init << "\nTEMP_init: " << TEMP_init
            << "\nrelh_init: " << relh_init << "\nW_avg: " << W_avg << "\nTAU_half: " << TAU_half
            << "\nrtol: " << rtol << "\natol: " << atol
            << "\n---------------------------------------------------------\n";
}

void OptionalConfigParams::AddSupersAtDomainTopParams::set_params(const YAML::Node &config) {
  const YAML::Node yaml = config["boundary_conditions"];

  COORD3LIM = yaml["COORD3LIM"].as<double>();
  newnsupers = yaml["newnsupers"].as<size_t>();
}

void OptionalConfigParams::AddSupersAtDomainTopParams::print_params() const {
  std::cout << "\n-------- AddSupersAtDomainTop Configuration Parameters --------------"
            << "\nCOORD3LIM: " << COORD3LIM << "\nnewnsupers: " << newnsupers
            << "\n---------------------------------------------------------\n";
}
