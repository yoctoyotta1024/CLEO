/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: yac_cartesian_dynamics.cpp
 * Project: coupldyn_yac
 * Created Date: Friday 13th October 2023
 * Author: Wilton Loch (WL)
 * Additional Contributors: Clara Bayley (CB), Lakshmi Aparna Devulapalli (LAD)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for dynamics solver in CLEO
 * where coupling is one-way and dynamics
 * are read from a file via yac
 */

#include "coupldyn_yac/yac_cartesian_dynamics.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <cstring>
#include "cleoconstants.hpp"
extern "C" {
#include "yac.h"
}

enum { VERTICAL = 0, EASTWARD = 1, NORTHWARD = 2 };

namespace dlc = dimless_constants;
/* return (k,i,j) indicies from idx for a flattened 3D array
with ndims [nz, nx, ny]. kij is useful for then getting
position in of a variable in a flattened array defined on
the faces of the same grid. E.g for the w velocity defined
on z faces of the grid which therefore has dims [nz+1, nx, ny] */
std::array<size_t, 3> kijfromindex(const std::array<size_t, 3> &ndims, const size_t index) {
  const size_t j = index / (ndims[0] * ndims[1]);
  const size_t k = index % ndims[0];
  const size_t i = index / ndims[0] - ndims[1] * j;

  return std::array<size_t, 3>{k, i, j};
}

double get_dlon_from_metres(const double lower_latitude,
    double delta_east) {
  // double EarthRadiusMeters = 6378137.0;
  double EarthRadiusMeters = (100000)/ (2 * M_PI);   // semi‑major axis a
  double dLon = 0.0;

  double MetresToRadians = (2*M_PI)/(100000);

  // dLon = (delta_east*dlc::COORD0)*MetresToRadians;

  if (delta_east != 0.0) {
    double cosLat = std::cos(lower_latitude);
    dLon = (delta_east*dlc::COORD0) / (EarthRadiusMeters * cosLat);
    return dLon;
  }
  return 0.0;

  // return dLon;
}
double get_dlat_from_metres(const double lower_latitude,
                            double delta_north) {
  // double EarthRadiusMeters = 6378137.0;
  double EarthRadiusMeters = (100000)/ (2 * M_PI);   // semi‑major axis a
  double dLat = 0.0;
  double MetresToRadians = (2*M_PI)/(100000);
  // dLat = (delta_north*dlc::COORD0)*MetresToRadians;

  if (delta_north != 0.0) {
    dLat = (delta_north*dlc::COORD0) / EarthRadiusMeters;
    return dLat;
  }
  return 0.0;

  // return dLat;
}


void create_vertex_coordinates(const Config &config, const std::array<size_t, 3> partition_size,
                               std::vector<std::vector<double>> gridbox_bounds,
                               std::array<std::array<double, 3>, 2> domain_bounds,
                               std::vector<double> &vertex_longitudes,
                               std::vector<double> &vertex_latitudes) {
  const auto lower_longitude = config.get_yac_dynamics().lower_longitude;
  const auto upper_longitude = config.get_yac_dynamics().upper_longitude;
  const auto lower_latitude = config.get_yac_dynamics().lower_latitude;
  const auto upper_latitude = config.get_yac_dynamics().upper_latitude;

  auto vertex_longitudes_new = std::vector<double>(partition_size[EASTWARD] + 1, 0);
  auto vertex_latitudes_new = std::vector<double>(partition_size[NORTHWARD] + 1, 0);

  // Defines the vertex longitude and latitude values in radians for grid creation
  // The values are later permuted by YAC to generate all vertex coordinates

  // ********* 1 process ***************
  /*
  for (size_t i = 0; i < vertex_longitudes.size(); i++)
    vertex_longitudes[i] =
        lower_longitude + i * ((upper_longitude - lower_longitude) / partition_size[EASTWARD]);

  for (size_t i = 0; i < vertex_latitudes.size(); i++)
    vertex_latitudes[i] =
        lower_latitude + i * ((upper_latitude - lower_latitude) / partition_size[NORTHWARD]);
  */
  // ***** multiprocess **********
  for (size_t i = 0; i < vertex_longitudes.size(); i++) {
    auto dLon = get_dlon_from_metres(lower_latitude,
        (gridbox_bounds[EASTWARD][i] - domain_bounds[0][EASTWARD]));
    vertex_longitudes[i] = lower_longitude + dLon;
    std::cout << "dLon : " << dLon <<std::endl;
  }

  for (size_t i = 0; i < vertex_latitudes.size(); i++) {
    auto dLat = get_dlat_from_metres(lower_latitude,
         (gridbox_bounds[NORTHWARD][i] - domain_bounds[0][NORTHWARD]));
    vertex_latitudes[i] = lower_latitude + dLat;
    std::cout << "dLat : " << dLat <<std::endl;
    }

  for (size_t i = 0; i < vertex_longitudes.size(); i++)
    std::cout << "vertex_lon old vs new " << vertex_longitudes[i] << " vs "
      << vertex_longitudes_new[i] << std::endl;

  for (size_t i = 0; i < vertex_latitudes.size(); i++)
    std::cout << "vertex_lat old vs new " << vertex_latitudes[i] << " vs "
      << vertex_latitudes_new[i] << std::endl;
}

/* Creates the YAC grid and defines the cell and edge points based on ndims data */
void create_grid_and_points_definitions(const Config &config, const std::array<size_t, 3> ndims,
                                        const std::string grid_name, int &grid_id,
                                        int &cell_point_id, int &edge_point_id,
                                        std::array<size_t, 3> partition_size,
                                        std::array<size_t, 3> partition_origin,
                                        std::vector<std::vector<double>> gridbox_bounds,
                                        std::array<std::array<double, 3>, 2> domain_bounds) {
  int cyclic_dimension[2] = {0, 0};
  int total_cells[2] = {static_cast<int>(partition_size[EASTWARD]),
                        static_cast<int>(partition_size[NORTHWARD])};
  int total_vertices[2] = {static_cast<int>(partition_size[EASTWARD] + 1),
                           static_cast<int>(partition_size[NORTHWARD] + 1)};
  int total_edges[2] = {static_cast<int>(partition_size[EASTWARD] *
                        (partition_size[NORTHWARD] + 1)),
                        static_cast<int>(partition_size[NORTHWARD] *
                        (partition_size[EASTWARD] + 1))};

  auto vertex_longitudes = std::vector<double>(partition_size[EASTWARD] + 1, 0);
  auto vertex_latitudes = std::vector<double>(partition_size[NORTHWARD] + 1, 0);
  auto cell_center_longitudes = std::vector<double>(partition_size[EASTWARD]);
  auto cell_center_latitudes = std::vector<double>(partition_size[NORTHWARD]);

  std::vector<double> edge_centers_longitudes;
  std::vector<double> edge_centers_latitudes;
  int * global_cell_ids;
  create_vertex_coordinates(config, partition_size, gridbox_bounds, domain_bounds,
                            vertex_longitudes, vertex_latitudes);

  // Defines a regular 2D grid
  yac_cdef_grid_reg2d(grid_name.c_str(), total_vertices, cyclic_dimension,
      vertex_longitudes.data(), vertex_latitudes.data(), &grid_id);
  yac_cset_global_index(global_cell_ids,
                        YAC_LOCATION_CELL,
                        grid_id);
  // --- Point definitions ---
  // Defines the cell center longitude and latitude values in radians
  // The values are later permuted by YAC to generate all cell center coordinates
  for (size_t i = 0; i < cell_center_longitudes.size(); i++)
    cell_center_longitudes[i] =
        vertex_longitudes[i] + (vertex_longitudes[i + 1] - vertex_longitudes[i]) / 2;

  for (size_t i = 0; i < cell_center_latitudes.size(); i++)
    cell_center_latitudes[i] =
        vertex_latitudes[i] + (vertex_latitudes[i + 1] - vertex_latitudes[i]) / 2;

  // Defines the edge center longitude and latitude values in radians.
  // Since it is not possible to generate edge center coordinates with a single
  // permutation of longitude and latitude values, usage of the
  // yac_cdef_points_unstruct is required. The call then takes x and y arrays,
  // with the actual radian coordinates for each edge center point. Therefore,
  // these arrays will have a size equal to the number of edges.
  for (size_t lat_index = 0; lat_index < vertex_latitudes.size() * 2 - 1; lat_index++) {
    if (lat_index % 2 == 0) {
      edge_centers_longitudes.insert(edge_centers_longitudes.end(), cell_center_longitudes.rbegin(),
                                     cell_center_longitudes.rend());
      edge_centers_latitudes.insert(edge_centers_latitudes.end(), cell_center_longitudes.size(),
                                    vertex_latitudes[lat_index / 2]);
    } else {
      edge_centers_longitudes.insert(edge_centers_longitudes.end(), vertex_longitudes.rbegin(),
                                     vertex_longitudes.rend());
      edge_centers_latitudes.insert(edge_centers_latitudes.end(), vertex_longitudes.size(),
                                    cell_center_latitudes[(lat_index - 1) / 2]);
    }
  }

  yac_cdef_points_reg2d(grid_id, total_cells, YAC_LOCATION_CELL, cell_center_longitudes.data(),
                        cell_center_latitudes.data(), &cell_point_id);
  yac_cdef_points_unstruct(grid_id, total_edges[0] + total_edges[1], YAC_LOCATION_EDGE,
                           edge_centers_longitudes.data(), edge_centers_latitudes.data(),
                           &edge_point_id);
}

/*
fill's target_array with values from yac_raw_data at multiplied by their conversion factor
*/

void CartesianDynamics::receive_yac_field(unsigned int yac_field_id, double **yac_raw_data,
                                          std::vector<double> &target_array,
                                          const size_t ndims_north, const size_t ndims_east,
                                          const size_t vertical_levels,
                                          double conversion_factor = 1.0) const {
  int info, error;
  int action;
  yac_cget_action(yac_field_id, &action);

  if (action != YAC_ACTION_GET_FOR_RESTART) {
    std::cout << "Get Action: " << action << std::endl;
    yac_cget(yac_field_id, vertical_levels, yac_raw_data, &info, &error);
  } else {
    std::cout << "Last Get Action: " << action << std::endl;
  }

  for (size_t j = 0; j < ndims_north; j++) {
    for (size_t i = 0; i < ndims_east; i++) {
      for (size_t k = 0; k < vertical_levels; k++) {
        auto vertical_idx = k;
        auto source_idx = j * ndims_east + i;
        auto ii = (ndims_east * j + i) * vertical_levels + k;
        target_array[ii] = yac_raw_data[vertical_idx][source_idx] / conversion_factor;
      }
    }
  }
}

/* This subroutine is the main entry point for receiving data from YAC.
 * It checks the dimensionality of the simulation based on the config data. */
void CartesianDynamics::receive_fields_from_yac() {
  receive_yac_field(temp_yac_id, yac_raw_cell_data, temp, ndims[NORTHWARD], ndims[EASTWARD],
                    ndims[VERTICAL], dlc::TEMP0);
  receive_yac_field(pressure_yac_id, yac_raw_cell_data, press, ndims[NORTHWARD], ndims[EASTWARD],
                    ndims[VERTICAL], dlc::P0);
  receive_yac_field(qvap_yac_id, yac_raw_cell_data, qvap, ndims[NORTHWARD], ndims[EASTWARD],
                    ndims[VERTICAL]);
  receive_yac_field(qcond_yac_id, yac_raw_cell_data, qcond, ndims[NORTHWARD], ndims[EASTWARD],
                    ndims[VERTICAL]);

  receive_yac_field(vertical_wind_yac_id, yac_raw_vertical_wind_data, wvel, ndims[NORTHWARD],
                    ndims[EASTWARD], ndims[VERTICAL] + 1, dlc::W0);

  receive_yac_field(eastward_wind_yac_id, yac_raw_edge_data, uvel, ndims[NORTHWARD],
                    ndims[EASTWARD] + 1, ndims[VERTICAL], dlc::W0);

  receive_yac_field(northward_wind_yac_id, yac_raw_edge_data, vvel, ndims[NORTHWARD] + 1,
                    ndims[EASTWARD], ndims[VERTICAL], dlc::W0);
}

void CartesianDynamics::send_yac_field(int field_id, double* field_data,
    double conversion_factor = 1.0) {
  auto ndims_vertical = ndims[VERTICAL];
  auto ndims_north = ndims[NORTHWARD];
  auto ndims_east = ndims[EASTWARD];
  auto ncells = ndims_north * ndims_east;

  int info, ierror, action;

  send_buffer = new double **[ndims_vertical];

  for (int j = 0; j < ndims_vertical; ++j) {
      send_buffer[j] = new double*[1];
      send_buffer[j][0] = new double[ncells];
  }

  for (size_t j = 0; j < ndims_north; j++) {
    for (size_t i = 0; i < ndims_east; i++) {
      for (size_t k = 0; k < ndims_vertical; k++) {
        auto ii = (ndims_east * j + i) * ndims_vertical + k;
        auto vertical_idx = k;
        auto source_idx = j * ndims_east + i;
        send_buffer[vertical_idx][0][source_idx] = field_data[ii] * conversion_factor;
      }
    }
  }
  yac_cget_action(field_id, &action);

  if (action != YAC_ACTION_PUT_FOR_RESTART) {
    std::cout << "Put Action: " << action << std::endl;
    yac_cput(field_id, ndims_vertical, send_buffer, &info, &ierror);
  } else {
    std::cout << "Last Put Action: " << action << std::endl;
  }

  std::cout << "yac_cput completed with info as: " << info << std::endl;

  for (int j = 0; j < ndims_vertical; ++j) {
      delete send_buffer[j][0];
      delete send_buffer[j];
  }
  delete[] send_buffer;
}
void CartesianDynamics::send_fields_to_yac(double* h_temp,
                               double* h_qvap,
                               double* h_qcond) {
  // double temp_field = h_temp;
  send_yac_field(temp_yac_id2, h_temp, dlc::TEMP0);
  send_yac_field(qvap_yac_id2, h_qvap);
  send_yac_field(qcond_yac_id2, h_qcond);
}

CartesianDynamics::CartesianDynamics(const Config &config, const std::array<size_t, 3> i_ndims,
                              const unsigned int nsteps, const CartesianDecomposition& decomp)
    : ndims(i_ndims),
      config(config),
      get_wvel(nullwinds()),
      get_uvel(nullwinds()),
      get_vvel(nullwinds()) {
  std::cout << "\n--- coupled cartesian dynamics from file ---\n";

  // Get YAC Component id from the communicator init class
  int component_id = init_communicator::get_yac_comp_id();

  // --- Grid definition ---
  int grid_id = -1;
  int cell_point_id = -1;
  int edge_point_id = -1;
  std::string grid_name = "cleo_grid";

  partition_size = decomp.get_local_partition_size();
  partition_origin = decomp.get_local_partition_origin();
  gridbox_bounds = decomp.get_local_gridbox_bounds();
  domain_bounds = decomp.get_domain_bounds();

  std::cout << "Partition Size in z:" << partition_size[0] << std::endl;
  std::cout << "Partition Size in x:" << partition_size[1] << std::endl;
  std::cout << "Partition Size in y:" << partition_size[2] << std::endl;
  std::cout << "Partition Origin in z:" << partition_origin[0] << std::endl;
  std::cout << "Partition Origin in x:" << partition_origin[1] << std::endl;
  std::cout << "Partition Origin in y:" << partition_origin[2] << std::endl;

  std::cout << "ndims Size in z:" << ndims[0] << std::endl;
  std::cout << "ndims Size in x:" << ndims[1] << std::endl;
  std::cout << "ndims Size in y:" << ndims[2] << std::endl;


  std::cout << "gridbox_bounds in z direction: " << gridbox_bounds[0].front() << " ; "
            << gridbox_bounds[0].back() <<std::endl;
  std::cout << "gridbox_bounds in x direction: " << gridbox_bounds[1].front() << " ; "
            << gridbox_bounds[1].back() <<std::endl;
  std::cout << "gridbox_bounds in y direction: " << gridbox_bounds[2].front() << " ; "
            << gridbox_bounds[2].back() <<std::endl;

  std::cout << "Domain_bounds in z direction: " << domain_bounds[0][0] << " ; "
            << domain_bounds[1][0] <<std::endl;

  create_grid_and_points_definitions(config, ndims, grid_name, grid_id, cell_point_id,
                                     edge_point_id, partition_size, partition_origin,
                                     gridbox_bounds, domain_bounds);

  // --- Interpolation stack ---
  int interp_stack_id;
  // yac_cget_interp_stack_config(&interp_stack_id);
  // yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_AVG, 1, 0.0, 1.0);
  yac_cget_interp_stack_config_from_string_yaml(
    "[{\"nnn\":"
    "  {\"n\": 1,"
    "   \"weighted\": \"arithmetic_average\"}}]",
    &interp_stack_id);

  // --- Field definitions ---
  int num_point_sets = 1;
  int horizontal_fields_collection_size = ndims[VERTICAL];
  int vertical_winds_collection_size = ndims[VERTICAL] + 1;

  const char coupling_timestep[6] = "PT30S";  // TODO(CB): move these variables to config
  const char coupldyn_grid_name[16] = "icon_atmos_grid";

  // --- Field definitions for recieving data from ICON ---

  yac_cdef_field("pressure", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &pressure_yac_id);

  yac_cdef_field("temperature", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &temp_yac_id);

  yac_cdef_field("qvap", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &qvap_yac_id);

  yac_cdef_field("qcond", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &qcond_yac_id);

  yac_cdef_field("eastward_wind", component_id, &edge_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &eastward_wind_yac_id);

  yac_cdef_field("northward_wind", component_id, &edge_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &northward_wind_yac_id);

  yac_cdef_field("vertical_wind", component_id, &cell_point_id, num_point_sets,
                 vertical_winds_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &vertical_wind_yac_id);

  // --- Field definitions for sending data to ICON ---

  yac_cdef_field("temperature_send", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &temp_yac_id2);

  yac_cdef_field("qvap_send", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &qvap_yac_id2);

  yac_cdef_field("qcond_send", component_id, &cell_point_id, num_point_sets,
                 horizontal_fields_collection_size, coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                 &qcond_yac_id2);
  int icon_lag = 1;
  int cleo_lag = 0;
  // --- Field coupling definitions for receiving from ICON ---
  yac_cdef_couple("atm", coupldyn_grid_name, "pressure", "cleo", "cleo_grid", "pressure",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "temperature", "cleo", "cleo_grid", "temperature",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "qvap", "cleo", "cleo_grid", "qvap", coupling_timestep,
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "qcond", "cleo", "cleo_grid", "qcond",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "eastward_wind", "cleo", "cleo_grid", "eastward_wind",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "northward_wind", "cleo", "cleo_grid",
                  "northward_wind", coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                  YAC_REDUCTION_TIME_NONE, interp_stack_id, icon_lag, cleo_lag);

  yac_cdef_couple("atm", coupldyn_grid_name, "vertical_wind", "cleo", "cleo_grid", "vertical_wind",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, icon_lag, cleo_lag);

  // --- Field coupling definitions for sending to ICON ---

  yac_cdef_couple("cleo", "cleo_grid", "temperature_send", "atm", coupldyn_grid_name,
                  "temperature_receive", coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT,
                  YAC_REDUCTION_TIME_NONE, interp_stack_id, cleo_lag, icon_lag);
  yac_cdef_couple("cleo", "cleo_grid", "qvap_send", "atm", coupldyn_grid_name, "qvap_receive",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, cleo_lag, icon_lag);
  yac_cdef_couple("cleo", "cleo_grid", "qcond_send", "atm", coupldyn_grid_name, "qcond_receive",
                  coupling_timestep, YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, cleo_lag, icon_lag);

  // ---------------------------------------------------------

  yac_cset_config_output_file(
  "coupling_debug.yaml", YAC_CONFIG_OUTPUT_FORMAT_YAML,
  YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF, 1);
  // --- End of YAC definitions ---

  // yac_cset_grid_output_file(grid_name.c_str(), "cleo_grid_vis.nc" );
  yac_cenddef();

  int yac_instance_id = yac_cget_default_instance_id();
  std::cout << "YAC instance id: " << yac_instance_id << std::endl;
  auto comp_metadata = yac_cget_component_metadata_instance(yac_instance_id, "cleo");
  std::cout << "Component Metadata:" << comp_metadata << std::endl;

  const char* string_onewaycoupling = "oneway_coupling";
  const char* string_twowaycoupling = "twoway_coupling";

  if (std::strcmp(comp_metadata, string_onewaycoupling) == 0) {
    yac_coupling_flag = 1;
    std::cout << "coupling_flag = " << yac_coupling_flag << std::endl;
  }

  if (std::strcmp(comp_metadata, string_twowaycoupling) == 0) {
    yac_coupling_flag = 2;
    std::cout << "coupling_flag = " << yac_coupling_flag << std::endl;
  }


  size_t horizontal_cell_number = yac_cget_grid_size(YAC_LOCATION_CELL, grid_id);
  size_t horizontal_edge_number = yac_cget_grid_size(YAC_LOCATION_EDGE, grid_id);

  yac_raw_cell_data = new double *[ndims[VERTICAL]];
  yac_raw_edge_data = new double *[ndims[VERTICAL]];
  yac_raw_vertical_wind_data = new double *[ndims[VERTICAL] + 1];

  for (size_t i = 0; i < ndims[VERTICAL]; i++) {
    yac_raw_cell_data[i] = new double[horizontal_cell_number];
    yac_raw_edge_data[i] = new double[horizontal_edge_number];
  }

  for (size_t i = 0; i < ndims[VERTICAL] + 1; i++)
    yac_raw_vertical_wind_data[i] = new double[horizontal_cell_number];

  // Initialization of target containers for receiving data
  press = std::vector<double>(horizontal_cell_number * ndims[VERTICAL], 0);
  temp = std::vector<double>(horizontal_cell_number * ndims[VERTICAL], 0);
  qvap = std::vector<double>(horizontal_cell_number * ndims[VERTICAL], 0);
  qcond = std::vector<double>(horizontal_cell_number * ndims[VERTICAL], 0);
  uvel = std::vector<double>(ndims[NORTHWARD] * (ndims[EASTWARD] + 1) * ndims[VERTICAL], 0);
  vvel = std::vector<double>(ndims[EASTWARD] * (ndims[NORTHWARD] + 1) * ndims[VERTICAL], 0);
  wvel = std::vector<double>(horizontal_cell_number * (ndims[VERTICAL] + 1), 0);

  // Calls the first data retrieval from YAC to have thermodynamic data for first timestep
  // receive_fields_from_yac();

  std::cout << "Finished setting up YAC for receiving:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n";

  // Defines the functions that will be used to retrieve data from the containers
  // (Can probably be simplified)
  set_winds(config);
  std::cout << "--- cartesian dynamics from YAC: success ---\n";
}

CartesianDynamics::~CartesianDynamics() {
  for (size_t i = 0; i < ndims[VERTICAL]; i++) {
    delete yac_raw_cell_data[i];
    delete yac_raw_edge_data[i];
  }
  for (size_t i = 0; i < ndims[VERTICAL] + 1; i++) delete yac_raw_vertical_wind_data[i];

  delete[] yac_raw_cell_data;
  delete[] yac_raw_edge_data;
  delete[] yac_raw_vertical_wind_data;

  yac_cfinalize();
}

/* depending on nspacedims, read in data
for 1-D, 2-D or 3-D wind velocity components */
void CartesianDynamics::set_winds(const Config &config) {
  const auto nspacedims = config.get_nspacedims();

  switch (nspacedims) {
    case 0:
      std::cout << "0-D model has no wind data\n";
      break;

    case 1:
    case 2:
    case 3:  // 1-D, 2-D or 3-D model
    {
      const std::string windstr(set_winds_from_yac(nspacedims));
      std::cout << windstr;
    } break;

    default:
      throw std::invalid_argument("nspacedims for wind data is invalid");
  }
}

/* Read in data from binary files for wind
velocity components in 1D, 2D or 3D model
and check they have correct size */
std::string CartesianDynamics::set_winds_from_yac(const unsigned int nspacedims) {
  std::string infostart(std::to_string(nspacedims) + "-D model, wind velocity");

  std::string infoend;
  switch (nspacedims) {
    case 3:  // 3-D model
      get_vvel = get_vvel_from_yac();
      infoend = ", u";
      [[fallthrough]];
    case 2:  // 3-D or 2-D model
      get_uvel = get_uvel_from_yac();
      infoend = ", v" + infoend;
      [[fallthrough]];
    case 1:  // 3-D, 2-D or 1-D model
      get_wvel = get_wvel_from_yac();
      infoend = "w" + infoend;
  }

  return infostart + " = [" + infoend + "]\n";
}

/* nullwinds retuns an empty function 'func' that returns
{0.0, 0.0}. Useful for setting get_[X]vel[Y]faces functions
in case of non-existent wind component e.g. get_uvelyface
when setup is 2-D model (x and z only) */
CartesianDynamics::get_winds_func CartesianDynamics::nullwinds() const {
  const auto func = [](const unsigned int ii) { return std::pair<double, double>{0.0, 0.0}; };

  return func;
}

/* set function for retrieving wvel defined at zfaces of
a gridbox with index 'gbxindex' and return vector
containting wvel data from binary file */
CartesianDynamics::get_winds_func CartesianDynamics::get_wvel_from_yac() const {
  const auto func = [&](const unsigned int gbxindex) {
    const auto kij =
        kijfromindex(ndims, static_cast<size_t>(gbxindex));  // [k,i,j] of gridbox centre on 3D grid
    const size_t nzfaces(ndims[0] + 1);                      // no. z faces to same 3D grid

    size_t lpos(ndims[1] * nzfaces * kij[2] + nzfaces * kij[1] +
                kij[0]);           // position of z lower face in 1D wvel vector
    const size_t uppos(lpos + 1);  // position of z upper face

    return std::pair(wvel.at(lpos), wvel.at(uppos));
  };

  return func;
}

/* returns vector of yvel retrieved from binary
file called 'filename' where uvel is defined on
the x-faces (coord1) of gridboxes */
CartesianDynamics::get_winds_func CartesianDynamics::get_uvel_from_yac() const {
  const auto func = [&](const unsigned int gbxindex) {
    const auto kij =
        kijfromindex(ndims, static_cast<size_t>(gbxindex));  // [k,i,j] of gridbox centre on 3D grid
    const size_t nxfaces(ndims[1] + 1);                      // no. x faces to same 3D grid

    size_t lpos(nxfaces * ndims[0] * kij[2] + ndims[0] * kij[1] +
                kij[0]);                  // position of x lower face in 1D uvel vector
    const size_t uppos(lpos + ndims[0]);  // position of x upper face

    return std::pair(uvel.at(lpos), uvel.at(uppos));
  };

  return func;
}

CartesianDynamics::get_winds_func CartesianDynamics::get_vvel_from_yac() const {
  const auto func = [&](const unsigned int gbxindex) {
    const size_t lpos(static_cast<size_t>(gbxindex));  // position of y lower face in 1D vvel vector
    const size_t uppos(lpos + ndims[1] * ndims[0]);    // position of x upper face

    return std::pair(vvel.at(lpos), vvel.at(uppos));
  };

  return func;
}
