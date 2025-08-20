/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: yac_cartesian_dynamics.hpp
 * Project: coupldyn_yac
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupleddynamics concept for dynamics solver in CLEO where coupling is
 * one-way and dynamics are read from file
 */

#ifndef LIBS_COUPLDYN_YAC_YAC_CARTESIAN_DYNAMICS_HPP_
#define LIBS_COUPLDYN_YAC_YAC_CARTESIAN_DYNAMICS_HPP_

#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "configuration/communicator.hpp"
#include "configuration/config.hpp"

/* contains 1-D vector for each (thermo)dynamic
variable which is ordered by gridbox at every timestep
e.g. press = [p_gbx0(t0), p_gbx1(t0), ,... , p_gbxN(t0),
p_gbx0(t1), p_gbx1(t1), ..., p_gbxN(t1), ..., p_gbxN(t_end)]
"pos[_X]" gives position of variable in a vector to read
current timestep from for the first gridbox (gbx0)  */
struct CartesianDynamics {
 private:
  using get_winds_func = std::function<std::pair<double, double>(const unsigned int)>;

  // number of (centres of) gridboxes in [coord3, coord1, coord2] directions
  const std::array<size_t, 3> ndims;
  const Config &config;

  /* --- (thermo)dynamic variables received from YAC --- */

  // Containers for cell-centered fields
  std::vector<double> press, temp, qvap, qcond;

  // Target for lon and lat edge data respectively
  // (these are copied from united_edge_data after receiving from YAC)
  std::vector<double> vvel, uvel;

  // Container for cell-centered vertical wind velocities
  std::vector<double> wvel;

  // YAC field ids
  int pressure_yac_id;
  int temp_yac_id;
  int qvap_yac_id;
  int qcond_yac_id;
  int eastward_wind_yac_id;
  int northward_wind_yac_id;
  int vertical_wind_yac_id;

  // Containers to receive data from YAC
  double **yac_raw_cell_data;
  double **yac_raw_edge_data;
  double **yac_raw_vertical_wind_data;

  /* --- Private functions --- */

  /* depending on nspacedims, read in data
  for 1-D, 2-D or 3-D wind velocity components */
  void set_winds(const Config &config);

  /* Read in data from YAC coupling for wind
  velocity components in 1D, 2D or 3D model */
  std::string set_winds_from_yac(const unsigned int nspacedims);

  /* nullwinds retuns an empty function 'func' that returns
  {0.0, 0.0}. Useful for setting get_[X]vel[Y]faces functions
  in case of non-existent wind component e.g. get_uvelyface
  when setup is 2-D model (x and z only) */
  get_winds_func nullwinds() const;

  /* returns vector of wvel, uvel and vvel retrieved from the YAC coupling
   * where wvel is defined on the z-faces (coord3), uvel is defined on the
   * x-faces (coord1) and vvel is defined on the y-faces (coord2) of gridboxes */
  get_winds_func get_wvel_from_yac() const;
  get_winds_func get_uvel_from_yac() const;
  get_winds_func get_vvel_from_yac() const;

  /* Receives an horizontal slice from yac,
   meaning a 2D set of grid boxes along u and w directions */
  void receive_hor_slice_from_yac(int cell_offset, int u_edges_offset, int w_edges_offset);

 public:
  CartesianDynamics(const Config &config, const std::array<size_t, 3> i_ndims,
                    const unsigned int nsteps);
  ~CartesianDynamics();

  get_winds_func get_wvel;  // funcs to get velocity defined in construction of class
  get_winds_func get_uvel;  // warning: these functions are not const member funcs by default
  get_winds_func get_vvel;

  double get_press(const size_t ii) const { return press.at(ii); }

  double get_temp(const size_t ii) const { return temp.at(ii); }

  double get_qvap(const size_t ii) const { return qvap.at(ii); }

  double get_qcond(const size_t ii) const { return qcond.at(ii); }

  /* Public call to receive data from YAC
   * If the problem is 2D turns into a wrapper for receive_hor_slice_from_yac */
  void receive_fields_from_yac();
  void receive_yac_field(unsigned int yac_field_id, double **yac_raw_data,
                         std::vector<double> &target_array, const size_t ndims_north,
                         const size_t ndims_east, const size_t vertical_levels,
                         double conversion_factor) const;
};

/* type satisfying CoupledDyanmics solver concept
specifically for thermodynamics and wind velocities
that are received from YAC */
struct YacDynamics {
 private:
  const unsigned int interval;
  const unsigned int end_time;
  std::shared_ptr<CartesianDynamics> dynvars;  // pointer to (thermo)dynamic variables

  /* Calls the get operations to receive data from YAC for each of the fields of interest */
  void run_dynamics(const unsigned int t_mdl) const { dynvars->receive_fields_from_yac(); }

 public:
  YacDynamics(const Config &config, const unsigned int couplstep, const std::array<size_t, 3> ndims,
              const unsigned int nsteps)
      : interval(couplstep),
        end_time(config.get_timesteps().T_END),
        dynvars(std::make_shared<CartesianDynamics>(config, ndims, nsteps)) {}

  auto get_couplstep() const { return interval; }

  void prepare_to_timestep() const {}

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void run_step(const unsigned int t_mdl, const unsigned int t_next) const {
    // Temporary simple solution to prevent a 4th coupling with yac from happening
    if (on_step(t_mdl) && t_mdl != end_time * 100) {
      run_dynamics(t_mdl);
    }
  }

  double get_press(const size_t ii) const { return dynvars->get_press(ii); }

  double get_temp(const size_t ii) const { return dynvars->get_temp(ii); }

  double get_qvap(const size_t ii) const { return dynvars->get_qvap(ii); }

  double get_qcond(const size_t ii) const { return dynvars->get_qcond(ii); }

  std::pair<double, double> get_wvel(const size_t ii) const { return dynvars->get_wvel(ii); }

  std::pair<double, double> get_uvel(const size_t ii) const { return dynvars->get_uvel(ii); }

  std::pair<double, double> get_vvel(const size_t ii) const { return dynvars->get_vvel(ii); }
};

#endif  // LIBS_COUPLDYN_YAC_YAC_CARTESIAN_DYNAMICS_HPP_
