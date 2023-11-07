/*
 * ----- CLEO -----
 * File: fromfile_cartesian_dynamics.hpp
 * Project: coupldyn_fromfile
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct obeying coupleddyanmics concept for
 * dynamics solver in CLEO where coupling is
 * one-way and dynamics are read from file
 */

#ifndef FROMFILE_CARTESIAN_DYNAMICS_HPP 
#define FROMFILE_CARTESIAN_DYNAMICS_HPP 

#include <array>
#include <vector>
#include <memory>
#include <string_view>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include <utility>
#include <stdexcept>

#include "initialise/config.hpp"
#include "initialise/readbinary.hpp"

struct CartesianDynamics
/* contains 1-D vector for each (thermo)dynamic
variable which is ordered by gridbox at every timestep
e.g. press = [p_gbx0(t0), p_gbx1(t0), ,... , p_gbxN(t0),
p_gbx0(t1), p_gbx1(t1), ..., p_gbxN(t1), ..., p_gbxN(t_end)]
"pos[_X]" gives position of variable in a vector to read
current timestep from for the first gridbox (gbx0)  */
{
private:
  using get_winds_func = std::function<std::pair<double,double>(const unsigned int)>;

  std::vector<double> wvel_zfaces; // w velocity defined on coord3 faces of gridboxes
  std::vector<double> uvel_xfaces; // u velocity defined on coord1 faces of gridboxes
  std::vector<double> vvel_yfaces; // v velocity defined on coord2 faces of gridboxes

  void set_winds(const Config &config);
  /* depending on nspacedims, read in data
  for 1-D, 2-D or 3-D wind velocity components */

  std::string set_winds_from_binaries(const unsigned int nspacedims,
                                      std::string_view wvel_filename,
                                      std::string_view uvel_filename,
                                      std::string_view vvel_filename);
  /* Read in data from binary files for wind
  velocity components in 1D, 2D or 3D model
  and check they have correct size */

  get_winds_func nullwinds() const;
  /* nullwinds retuns an empty function 'func' that returns
  {0.0, 0.0}. Useful for setting get_[X]vel[Y]faces functions
  in case of non-existent wind component e.g. get_uvelyface
  when setup is 2-D model (x and z only) */

  get_winds_func get_wvel_from_binary() const;
  /* returns vector of wvel retrieved from binary
  file called 'filename' where wvel is defined on
  the z-faces (coord3) of gridboxes */

  get_winds_func get_uvel_from_binary() const;
  /* returns vector of yvel retrieved from binary
  file called 'filename' where uvel is defined on
  the x-faces (coord1) of gridboxes */

  get_winds_func get_vvel_from_binary() const;
  /* returns vector of vvel retrieved from binary
  file called 'filename' where vvel is defined on
  the y-faces (coord2) of gridboxes */

  void check_thermodyanmics_vectorsizes(const unsigned int nspacedims,
                                        const std::array<size_t, 3> &ndims,
                                        const unsigned int nsteps) const;

public:
  /* position in vector for 0th gridbox at current timestep  */
  const std::array<size_t, 3> ndims; // number of (centres of) gridboxes in [coord3, coord1, coord2] directions
  size_t pos;                        // for variable defined at gridbox centres
  size_t pos_zface;                  // for variable defined at gridbox coord3 faces
  size_t pos_xface;                  // for variable defined at gridbox coord1 faces
  size_t pos_yface;                  // for variable defined at gridbox coord2 faces

  /* (thermo)dynamic variables read from file */
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;
  
  get_winds_func get_wvel; // funcs to get velocity defined in construction of class 
  get_winds_func get_uvel; // warning: these functions are not const member funcs by default
  get_winds_func get_vvel;

  CartesianDynamics(const Config &config,
                    const std::array<size_t, 3> i_ndims,
                    const unsigned int nsteps);

  void increment_position();
  /* updates positions to gbx0 in vector (for
  acessing value at next timestep). Assumes domain
  is decomposed into cartesian C grid with dimensions
  (ie. number of gridboxes in each dimension) ndims */
};

struct FromFileDynamics
/* type satisfying CoupledDyanmics solver concept
specifically for thermodynamics and wind velocities
that are read from binary files */
{
private:
  const unsigned int interval;
  std::unique_ptr<CartesianDynamics> dynvars; // pointer to (thermo)dyanmic variables

  void run_dynamics(const unsigned int t_mdl) const
  /* increment position of thermodata for 0th gridbox
  to positon at next timestep (ie. ngridbox_faces
  further along vector) */
  {
    dynvars->increment_position();
  }

public:
  FromFileDynamics(const Config &config,
                   const unsigned int couplstep,
                   const std::array<size_t, 3> ndims,
                   const unsigned int nsteps)
      : interval(couplstep),
        dynvars(std::make_unique<CartesianDynamics>(config,
                                                    ndims,
                                                    nsteps)) {}

  auto get_couplstep() const
  {
    return interval;
  }

  void prepare_to_timestep() const {}

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void run_step(const unsigned int t_mdl,
                const unsigned int t_next) const
  {
    if (on_step(t_mdl))
    {
      run_dynamics(t_mdl);
    }
  }

};

#endif // FROMFILE_CARTESIAN_DYNAMICS_HPP  