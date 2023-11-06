/*
 * ----- CLEO -----
 * File: fromfile_cartesian_dynamics.cpp
 * Project: coupldyn_fromfile
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 6th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for dynamics solver in CLEO
 * where coupling is one-way and dynamics
 * are read from file
 */

#include "./fromfile_cartesian_dynamics.hpp"

std::vector<double>
thermodynamicvar_from_binary(std::string_view filename)
{
  /* open file and read in the metatdata
  for all the variables in that file */
  std::ifstream file(open_binary(filename));
  std::vector<VarMetadata> meta(metadata_from_binary(file));

  /* read in the data for the 1st variable in the file */
  std::vector<double>
      thermovar(vector_from_binary<double>(file, meta.at(0)));

  return thermovar;
}

CartesianDynamics::
    CartesianDynamics(const Config &config,
                      const std::array<size_t, 3> i_ndims)
    : ndims(i_ndims),
      pos(0),
      pos_zface(0),
      pos_xface(0),
      pos_yface(0)
      press(thermodynamicvar_from_binary(config.press_filename)),
      temp(thermodynamicvar_from_binary(config.temp_filename)),
      qvap(thermodynamicvar_from_binary(config.qvap_filename)),
      qcond(thermodynamicvar_from_binary(config.qcond_filename)),
      // wvel(), uvel(), vvel(),
{
  std::cout << "\nFinished reading thermodynamics from binaries for:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n";

  set_windvelocities(config);

  // check_thermodyanmics_vectorsizes(config.SDnspace, ndims, nsteps);
}

void CartesianDynamics::set_windvelocities(const Config &config)
/* depending on SDnspace, read in data
for wind velocity components (or not)*/
{
  switch (config.nspacedims)
  {
  case 0:
    std::cout << "0-D model has no wind data\n";
    break;
  case 1: //TODO 

  default:
    throw std::invalid_argument("nspacedims for wind data is invalid");
  }

  else if (SDnspace <= 3) // means 1 <= SDnspace < 4
  {
    const std::string windstr(set_windvelocities_frombinaries(config));
    std::cout << windstr;
  }

}

void CartesianDynamics::increment_position()
/* updates positions to gbx0 in vector (for
acessing value at next timestep). Assumes domain
is decomposed into cartesian C grid with dimensions
(ie. number of gridboxes in each dimension) ndims */
{
  pos += ndims[0] * ndims[1] * ndims[2];
  pos_zface += (ndims[0] + 1) * ndims[1] * ndims[2];
  pos_xface += ndims[0] * (ndims[1] + 1) * ndims[2];
  pos_yface += ndims[0] * ndims[1] * (ndims[2] + 1);
}