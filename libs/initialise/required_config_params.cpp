/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: required_config_params.cpp
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
 * Functions involved in reading required configuration parameters from a config file.
 */

#include "initialise/required_config_params.hpp"

/* read configuration file given by config_filename to set members of required configuration */
RequiredConfigParams::RequiredConfigParams(const std::filesystem::path config_filename) {
  const YAML::Node config = YAML::LoadFile(std::string{config_filename});

  /* convert string to std::filesystem::path type */
  auto fspath_from_yaml = [](YAML::Node& yaml, const std::string& key) {
    return std::filesystem::path(yaml[key].as<std::string>());
  };

  YAML::Node yaml = config["inputfiles"];
  inputfiles.constants_filename = fspath_from_yaml(yaml, "constants_filename");
  inputfiles.grid_filename = fspath_from_yaml(yaml, "grid_filename");

  yaml = config["outputdata"];
  outputdata.setup_filename = fspath_from_yaml(yaml, "setup_filename");
  outputdata.stats_filename = fspath_from_yaml(yaml, "stats_filename");
  outputdata.zarrbasedir = fspath_from_yaml(yaml, "zarrbasedir");
  outputdata.maxchunk = yaml["maxchunk"].as<size_t>();

  yaml = config["domain"];
  domain.nspacedims = yaml["nspacedims"].as<unsigned int>();
  domain.ngbxs = yaml["ngbxs"].as<size_t>();
  domain.maxnsupers = yaml["maxnsupers"].as<size_t>();

  yaml = config["timesteps"];
  timesteps.CONDTSTEP = yaml["CONDTSTEP"].as<double>();
  timesteps.COLLTSTEP = yaml["COLLTSTEP"].as<double>();
  timesteps.MOTIONTSTEP = yaml["MOTIONTSTEP"].as<double>();
  timesteps.COUPLTSTEP = yaml["COUPLTSTEP"].as<double>();
  timesteps.OBSTSTEP = yaml["OBSTSTEP"].as<double>();
  timesteps.T_END = yaml["T_END"].as<double>();

  print_params();
}

void RequiredConfigParams::print_params() const {
  std::cout << "\n-------- Required Configuration Parameters --------------"
            << "\nconstants_filename : " << inputfiles.constants_filename
            << "\ngrid_filename : " << inputfiles.grid_filename
            << "\nsetup_filename : " << outputdata.setup_filename
            << "\nstats_filename : " << outputdata.stats_filename
            << "\nzarrbasedir : " << outputdata.zarrbasedir
            << "\nmaxchunk : " << outputdata.maxchunk << "\nnspacedims : " << domain.nspacedims
            << "\nngbxs : " << domain.ngbxs << "\nmaxnsupers : " << domain.maxnsupers
            << "\nCONDTSTEP : " << timesteps.CONDTSTEP << "\nCOLLTSTEP : " << timesteps.COLLTSTEP
            << "\nMOTIONTSTEP : " << timesteps.MOTIONTSTEP
            << "\nCOUPLTSTEP : " << timesteps.COUPLTSTEP << "\nOBSTSTEP : " << timesteps.OBSTSTEP
            << "\nT_END : " << timesteps.T_END
            << "\n---------------------------------------------------------\n";
}
