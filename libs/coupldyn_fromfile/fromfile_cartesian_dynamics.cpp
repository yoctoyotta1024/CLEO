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

std::vector<double> wvel_from_binary(std::string_view filename);
/* returns vector of wvel retrieved from binary
file called 'filename' where wvel is defined on
the z-faces (coord3) of gridboxes */

std::vector<double> uvel_from_binary(std::string_view filename);
/* returns vector of yvel retrieved from binary
file called 'filename' where uvel is defined on
the x-faces (coord1) of gridboxes */

std::vector<double> vvel_from_binary(std::string_view filename);
/* returns vector of vvel retrieved from binary
file called 'filename' where vvel is defined on
the y-faces (coord2) of gridboxes */

std::vector<double>
thermodynamicvar_from_binary(std::string_view filename)
/* open file called 'filename' and return vector
of doubles for first variable in that file */
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

std::function<std::pair<double, double>(const unsigned int)> nullwinds()
/* nullwinds retuns an empty function 'func' that returns
{0.0, 0.0}. Useful for setting get_[X]vel[Y]faces functions
in case of non-existent wind component e.g. get_uvelyface
when setup is 2-D model (x and z only) */
{
  const auto func = [](const unsigned int ii)
  { return std::pair<double, double>{0.0, 0.0}; };

  return func;
}

CartesianDynamics::
    CartesianDynamics(const Config &config,
                      const std::array<size_t, 3> i_ndims)
    : wvel_zfaces(0), uvel_xfaces(0), vvel_yfaces(0),
      ndims(i_ndims),
      pos(0),
      pos_zface(0),
      pos_xface(0),
      pos_yface(0),
      press(thermodynamicvar_from_binary(config.press_filename)),
      temp(thermodynamicvar_from_binary(config.temp_filename)),
      qvap(thermodynamicvar_from_binary(config.qvap_filename)),
      qcond(thermodynamicvar_from_binary(config.qcond_filename)),
      get_wvel(nullwinds()), get_uvel(nullwinds()), get_vvel(nullwinds())
{
  std::cout << "\nFinished reading thermodynamics from binaries for:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n";

  set_winds(config);

  // check_thermodyanmics_vectorsizes(config.nspacedims, ndims, nsteps);
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

void CartesianDynamics::set_winds(const Config &config)
/* depending on nspacedims, read in data
for 1-D, 2-D or 3-D wind velocity components */
{
  const unsigned int nspacedims(config.nspacedims);

  switch (nspacedims)
  {
  case 0:
    std::cout << "0-D model has no wind data\n";
    break;

  case 1:
  case 2:
  case 3: // 1-D, 2-D or 3-D model
    const std::string windstr(
        set_winds_from_binaries(nspacedims,
                                config.wvel_filename,
                                config.uvel_filename,
                                config.vvel_filename));
    std::cout << windstr;

  default:
    throw std::invalid_argument("nspacedims for wind data is invalid");
  }
}

std::string CartesianDynamics::
    set_winds_from_binaries(const unsigned int nspacedims,
                            std::string_view wvel_filename,
                            std::string_view uvel_filename,
                            std::string_view vvel_filename)
/* Read in data from binary files for wind
velocity components in 1D, 2D or 3D model
and check they have correct size */
{
  std::string info(std::to_string(nspacedims) +
                   "-D model, wind velocity");

  std::string vel;
  switch (nspacedims)
  {
  case 3: // 3-D model
    vvel_yfaces = get_vvel_from_binary();
    vvel_yfaces = thermodynamicvar_from_binary(vvel_filename);
    vel = ", u";
  case 2: // 3-D or 2-D model
    uvel_xfaces = get_uvel_from_binary();
    uvel_xfaces = thermodynamicvar_from_binary(uvel_filename); 
    vel = ", v" + vel;
  case 1: // 3-D, 2-D or 1-D model
    get_wvel = get_wvel_from_binary();
    wvel_zfaces = thermodynamicvar_from_binary(wvel_filename);
    vel = "w" + vel;
  }

  info += " = [" + vel + "]\n";
  return info;
}