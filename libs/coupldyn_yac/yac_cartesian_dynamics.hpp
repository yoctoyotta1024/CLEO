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

#include <Kokkos_Core.hpp>
#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cartesiandomain/cartesian_decomposition.hpp"
#include "configuration/communicator.hpp"
#include "configuration/config.hpp"
#include "superdrops/state.hpp"

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
  const Config& config;
  int yac_coupling_flag;

  /* --- (thermo)dynamic variables sent/received via YAC --- */

  // Containers for sending cell-centered fields
  Kokkos::View<double*, Kokkos::HostSpace> delta_temp_send;
  Kokkos::View<double*, Kokkos::HostSpace> delta_qvap_send;
  Kokkos::View<double*, Kokkos::HostSpace> delta_qcloud_send;
  Kokkos::View<double*, Kokkos::HostSpace> delta_qrain_send;

  // Containers for receiving cell-centered fields
  std::vector<double> press_recv, temp_recv, qvap_recv, qcloud_recv, qrain_recv;

  // Containers for receiving edge datat on lon and lat edges respectively
  // (these are copied from united_edge_data after receiving from YAC)
  std::vector<double> vvel_recv, uvel_recv;

  // Container for receiving cell-centered vertical wind velocities
  std::vector<double> wvel_recv;

  // YAC field ids
  int pressure_yac_id_recv;
  int temp_yac_id_recv;
  int temp_yac_id_send;
  int qvap_yac_id_recv;
  int qvap_yac_id_send;
  int qcloud_yac_id_recv;
  int qcloud_yac_id_send;
  int qrain_yac_id_send;
  int qrain_yac_id_recv;
  int eastward_wind_yac_id_recv;
  int northward_wind_yac_id_recv;
  int vertical_wind_yac_id_recv;

  // Containers to receive data via YAC
  double** yac_raw_cell_data;
  double** yac_raw_edge_data;
  double** yac_raw_vertical_wind_data;

  // Containers to send data via YAC
  double*** yac_raw_cell_data_send;

  std::array<size_t, 3> partition_origin;
  std::array<size_t, 3> partition_size;
  std::vector<std::vector<double>> gridbox_bounds;
  std::array<std::array<double, 3>, 2> domain_bounds;

  /* --- Private functions --- */

  /* depending on nspacedims, read in data
  for 1-D, 2-D or 3-D wind velocity components */
  void set_winds(const Config& config);

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

  /* functions for handling YAC field communication to/from Cleo */
  void receive_yac_cell_field(unsigned int yac_field_id, double** yac_raw_cell_data,
                              std::vector<double>& target_array, const size_t vertical_levels,
                              double conversion_factor) const;
  void receive_yac_edge_field(unsigned int yac_field_id, double** yac_raw_edge_data,
                              std::vector<double>& target_array, double conversion_factor,
                              bool eastward_edge) const;
  void send_yac_cell_field(int field_id, double* field_data, double conversion_factor);

 public:
  CartesianDynamics(const Config& config, const std::array<size_t, 3> i_ndims,
                    const unsigned int nsteps, const CartesianDecomposition& decomp);
  ~CartesianDynamics();

  get_winds_func get_wvel;  // funcs to get velocity defined in construction of class
  get_winds_func get_uvel;  // warning: these functions are not const member funcs by default
  get_winds_func get_vvel;

  int get_yac_coupling_flag() const { return yac_coupling_flag; }

  double get_press(const size_t ii) const { return press_recv.at(ii); }

  double get_temp(const size_t ii) const { return temp_recv.at(ii); }

  double get_qvap(const size_t ii) const { return qvap_recv.at(ii); }

  void set_temp_delta(const size_t ii, const double new_temp) const {
    delta_temp_send(ii) = new_temp - temp_recv.at(ii);
  }

  void set_qvap_delta(const size_t ii, const double new_qvap) const {
    delta_qvap_send(ii) = new_qvap - qvap_recv.at(ii);
  }

  void set_qcloud_delta(const size_t ii, const double new_qcloud) const {
    delta_qcloud_send(ii) = 0.0;  // new_qcloud - qcloud_recv.at(ii);
  }

  void set_qrain_delta(const size_t ii, const double new_qrain) const {
    delta_qrain_send(ii) = 0.0;  // new_qrain - qrain_recv.at(ii);
  }

  /* Public calls to send/receive data via YAC */
  void receive_fields_from_yac();

  void send_fields_to_yac();
};

/* type satisfying CoupledDyanmics solver concept
specifically for thermodynamics and wind velocities
that are received from YAC */
struct YacCartesianDynamics {
 private:
  const unsigned int interval;
  const unsigned int end_time;
  std::shared_ptr<CartesianDynamics> dynvars;  // pointer to (thermo)dynamic variables

 public:
  YacCartesianDynamics(const Config& config, const unsigned int couplstep,
                       const std::array<size_t, 3> ndims, const unsigned int nsteps,
                       const CartesianDecomposition& decomp)
      : interval(couplstep),
        end_time(config.get_timesteps().T_END),
        dynvars(std::make_shared<CartesianDynamics>(config, ndims, nsteps, decomp)) {}

  auto get_couplstep() const { return interval; }

  void prepare_to_timestep() const {}

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void run_step(const unsigned int t_mdl, const unsigned int t_next) const {}

  int get_yac_coupling_flag() const { return dynvars->get_yac_coupling_flag(); }

  void receive_fields_from_yac() const { dynvars->receive_fields_from_yac(); }

  void send_fields_to_yac() const { dynvars->send_fields_to_yac(); }

  double get_press(const size_t ii) const { return dynvars->get_press(ii); }

  double get_temp(const size_t ii) const { return dynvars->get_temp(ii); }

  double get_qvap(const size_t ii) const { return dynvars->get_qvap(ii); }

  void set_temp_delta(const size_t ii, const double new_temp) const {
    dynvars->set_temp_delta(ii, new_temp);
  }

  void set_qvap_delta(const size_t ii, const double new_qvap) const {
    dynvars->set_qvap_delta(ii, new_qvap);
  }

  void set_qcloud_delta(const size_t ii, const double new_qcloud) const {
    dynvars->set_qcloud_delta(ii, new_qcloud);
  }

  void set_qrain_delta(const size_t ii, const double new_qrain) const {
    dynvars->set_qrain_delta(ii, new_qrain);
  }

  std::pair<double, double> get_wvel(const size_t ii) const { return dynvars->get_wvel(ii); }

  std::pair<double, double> get_uvel(const size_t ii) const { return dynvars->get_uvel(ii); }

  std::pair<double, double> get_vvel(const size_t ii) const { return dynvars->get_vvel(ii); }
};

// int YacCartesianDynamics::get_counter = 1;

#endif  // LIBS_COUPLDYN_YAC_YAC_CARTESIAN_DYNAMICS_HPP_
