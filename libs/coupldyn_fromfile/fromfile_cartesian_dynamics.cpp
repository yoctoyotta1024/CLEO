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

void FromFileDynamics::prepare_to_timestep() const
{
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

CartesianDynamics::CartesianDynamics(const Config &config)
    : ndims(ndims),
      atpos(0),
      atpos_zface(0),
      atpos_xface(0),
      atpos_yface(0),
      press(thermodynamicvar_from_binary(config.press_filename)),
      temp(thermodynamicvar_from_binary(config.temp_filename)),
      qvap(thermodynamicvar_from_binary(config.qvap_filename)),
      qcond(thermodynamicvar_from_binary(config.qcond_filename)),
      wvel(), uvel(), vvel(),
      get_wvelzfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; }),
      get_uvelxfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; }),
      get_vvelyfaces([](const unsigned int ii)
                     { return std::pair<double, double>{0.0, 0.0}; })
{
  std::string windstr = set_windvelocities(config);
  std::cout << "\nFinished reading thermodynamics from binaries for:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n  "
            << windstr << '\n';

  check_thermodyanmics_vectorsizes(config.SDnspace, ndims, nsteps);
}