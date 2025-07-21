/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: required_config_params.cpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions involved in reading required configuration parameters from a config file.
 */

#include "configuration/required_config_params.hpp"

/* read configuration file given by config_filename to set members of required configuration */
RequiredConfigParams::RequiredConfigParams(const std::filesystem::path config_filename) {
  const YAML::Node config = YAML::LoadFile(std::string{config_filename});

  /* convert string to std::filesystem::path type */
  auto fspath_from_yaml = [](YAML::Node& node, const std::string& key) {
    return std::filesystem::path(node[key].as<std::string>());
  };

  YAML::Node node = config["inputfiles"];
  inputfiles.constants_filename = fspath_from_yaml(node, "constants_filename");
  inputfiles.grid_filename = fspath_from_yaml(node, "grid_filename");

  node = config["outputdata"];
  outputdata.setup_filename = fspath_from_yaml(node, "setup_filename");
  outputdata.zarrbasedir = fspath_from_yaml(node, "zarrbasedir");
  outputdata.maxchunk = node["maxchunk"].as<size_t>();

  node = config["domain"];
  domain.nspacedims = node["nspacedims"].as<unsigned int>();
  domain.ngbxs = node["ngbxs"].as<size_t>();
  domain.maxnsupers = node["maxnsupers"].as<size_t>();

  node = config["timesteps"];
  timesteps.CONDTSTEP = node["CONDTSTEP"].as<double>();
  timesteps.COLLTSTEP = node["COLLTSTEP"].as<double>();
  timesteps.MOTIONTSTEP = node["MOTIONTSTEP"].as<double>();
  timesteps.COUPLTSTEP = node["COUPLTSTEP"].as<double>();
  timesteps.OBSTSTEP = node["OBSTSTEP"].as<double>();
  timesteps.T_END = node["T_END"].as<double>();

  print_params();
}

void RequiredConfigParams::print_params() const {
  std::cout << "\n-------- Required Configuration Parameters --------------"
            << "\nconstants_filename : " << inputfiles.constants_filename
            << "\ngrid_filename : " << inputfiles.grid_filename
            << "\nsetup_filename : " << outputdata.setup_filename
            << "\nzarrbasedir : " << outputdata.zarrbasedir
            << "\nmaxchunk : " << outputdata.maxchunk << "\nnspacedims : " << domain.nspacedims
            << "\nngbxs : " << domain.ngbxs << "\nmaxnsupers : " << domain.maxnsupers
            << "\nCONDTSTEP : " << timesteps.CONDTSTEP << "\nCOLLTSTEP : " << timesteps.COLLTSTEP
            << "\nMOTIONTSTEP : " << timesteps.MOTIONTSTEP
            << "\nCOUPLTSTEP : " << timesteps.COUPLTSTEP << "\nOBSTSTEP : " << timesteps.OBSTSTEP
            << "\nT_END : " << timesteps.T_END
            << "\n---------------------------------------------------------\n";
}
