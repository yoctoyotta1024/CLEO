/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: yac_cartesian_dynamics.cpp
 * Project: coupldyn_yac
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 17th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for dynamics solver in CLEO
 * where coupling is one-way and dynamics
 * are read from file
 */

#include "coupldyn_yac/yac_cartesian_dynamics.hpp"

#include <iostream>

#include <mpi.h>
extern "C" {
#include "yac.h"
}

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

void create_vertex_coordinates(const Config &config,
                               const std::array<size_t, 3> ndims,
                               std::vector<double> & vertex_longitudes,
                               std::vector<double> & vertex_latitudes) {
  double lower_longitude = 0;
  double upper_longitude = ndims[0] * (2 * std::numbers::pi / (ndims[0] + 1));
  double lower_latitude = (-0.5 * std::numbers::pi * ndims[1]) / (ndims[1] + 2);
  double upper_latitude = (0.5 * std::numbers::pi * ndims[1]) / (ndims[1] + 2);

  if (config.get_yac_dynamics().lower_longitude != NaNVals::dbl() &&
      config.get_yac_dynamics().upper_longitude != NaNVals::dbl() &&
      config.get_yac_dynamics().lower_latitude != NaNVals::dbl() &&
      config.get_yac_dynamics().upper_latitude != NaNVals::dbl()) {
      lower_longitude = config.get_yac_dynamics().lower_longitude;
      upper_longitude = config.get_yac_dynamics().upper_longitude;
      lower_latitude = config.get_yac_dynamics().lower_latitude;
      upper_latitude = config.get_yac_dynamics().upper_latitude;
  }

  // Defines the vertex longitude and latitude values in radians for grid creation
  // The values are later permuted by YAC to generate all vertex coordinates
  for (size_t i = 0; i < vertex_longitudes.size(); i++)
    vertex_longitudes[i] = lower_longitude + i * ((upper_longitude - lower_longitude) / ndims[0]);

  for (size_t i = 0; i < vertex_latitudes.size(); i++)
    vertex_latitudes[i] = lower_latitude + i * ((upper_latitude - lower_latitude) / ndims[1]);
}

/* Creates the YAC grid and defines the cell and edge points based on ndims data */
void create_grid_and_points_definitions(const Config &config,
                                        const std::array<size_t, 3> ndims,
                                        const std::string grid_name,
                                        int & grid_id, int & cell_point_id, int & edge_point_id) {
  int cyclic_dimension[2] = {0, 0};
  int total_cells[2] = {static_cast<int>(ndims[0]), static_cast<int>(ndims[1])};
  int total_vertices[2] = {static_cast<int>(ndims[0] + 1), static_cast<int>(ndims[1] + 1)};
  int total_edges[2] = {static_cast<int>(ndims[0] * (ndims[1] + 1)),
                        static_cast<int>(ndims[1] * (ndims[0] + 1))};

  auto vertex_longitudes = std::vector<double>(ndims[0] + 1, 0);
  auto vertex_latitudes = std::vector<double>(ndims[1] + 1, 0);
  auto cell_center_longitudes = std::vector<double>(ndims[0]);
  auto cell_center_latitudes = std::vector<double>(ndims[1]);
  std::vector<double> edge_centers_longitudes;
  std::vector<double> edge_centers_latitudes;

  create_vertex_coordinates(config, ndims, vertex_longitudes, vertex_latitudes);

  // Defines a regular 2D grid
  yac_cdef_grid_reg2d(grid_name.c_str(), total_vertices, cyclic_dimension,
                      vertex_longitudes.data(), vertex_latitudes.data(), &grid_id);

  // --- Point definitions ---
  // Defines the cell center longitude and latitude values in radians
  // The values are later permuted by YAC to generate all cell center coordinates
  for (size_t i = 0; i < cell_center_longitudes.size(); i++)
    cell_center_longitudes[i] = vertex_longitudes[i]
                                + (vertex_longitudes[i + 1] - vertex_longitudes[i]) / 2;

  for (size_t i = 0; i < cell_center_latitudes.size(); i++)
    cell_center_latitudes[i] = vertex_latitudes[i]
                               + (vertex_latitudes[i + 1] - vertex_latitudes[i]) / 2;

  // Defines the edge center longitude and latitude values in radians.
  // Since it is not possible to generate edge center coordinates with a single
  // permutation of longitude and latitude values, usage of the
  // yac_cdef_points_unstruct is required. The call then takes x and y arrays,
  // with the actual radian coordinates for each edge center point. Therefore,
  // these arrays will have a size equal to the number of edges.
  for (size_t lat_index = 0; lat_index < vertex_latitudes.size() * 2 - 1; lat_index++) {
    if (lat_index % 2 == 0) {
      edge_centers_longitudes.insert(edge_centers_longitudes.end(), cell_center_longitudes.begin(),
                                     cell_center_longitudes.end());
      edge_centers_latitudes.insert(edge_centers_latitudes.end(), cell_center_longitudes.size(),
                                    vertex_latitudes[lat_index / 2]);
    } else {
      edge_centers_longitudes.insert(edge_centers_longitudes.end(), vertex_longitudes.begin(),
                                     vertex_longitudes.end());
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

void CartesianDynamics::receive_yac_field(unsigned int field_type,
                                          unsigned int yac_field_id,
                                          double ** yac_raw_data,
                                          std::vector<double> & target_array,
                                          size_t vertical_levels) {
  int info, error;
  unsigned int total_horizontal_cells = ndims[0] * ndims[1];
  bool edge_dimension = false;
  std::vector<double>::iterator target_it = target_array.begin();

  yac_cget(yac_field_id, vertical_levels, yac_raw_data, &info, &error);

  switch (field_type) {
    case 0:
      for (size_t i = 0; i < vertical_levels; i++)
        for (size_t j = 0; j < total_horizontal_cells; j++)
          target_array[i * total_horizontal_cells + j] = yac_raw_data[i][j];
      return;

    case 1:
      edge_dimension = false;
      break;

    case 2:
      edge_dimension = true;
      break;
  }

  for (size_t vertical_index = 0; vertical_index < ndims[2]; vertical_index++) {
    unsigned int source_index = 0;
    for (size_t lat_index = 0; lat_index < (ndims[1] + 1) * 2 - 1; lat_index++) {
      if (lat_index % 2 == edge_dimension) {
        for (size_t index = 0; index < ndims[0] + edge_dimension;
             index++, target_it++, source_index++)
          *target_it = yac_raw_data[vertical_index][source_index];
      } else
        source_index += ndims[0] + !edge_dimension;
    }
  }
}

/* This subroutine is the main entry point for receiving data from YAC.
 * It checks the dimensionality of the simulation based on the config data. */
void CartesianDynamics::receive_field_collections_from_yac() {
  enum field_types {
    CELL,
    U_EDGE,
    W_EDGE
  };

  receive_yac_field(CELL, temp_yac_id, yac_raw_cell_data, temp, ndims[2]);
  receive_yac_field(CELL, pressure_yac_id, yac_raw_cell_data, press, ndims[2]);
  receive_yac_field(CELL, qvap_yac_id, yac_raw_cell_data, qvap, ndims[2]);
  receive_yac_field(CELL, qcond_yac_id, yac_raw_cell_data, qcond, ndims[2]);

  receive_yac_field(CELL, vvel_yac_id, yac_raw_vertical_wind_data, vvel, ndims[2] + 1);

  receive_yac_field(U_EDGE, eastward_wind_yac_id, yac_raw_edge_data, uvel, ndims[2]);
  receive_yac_field(W_EDGE, northward_wind_yac_id, yac_raw_edge_data, wvel, ndims[2]);
}

CartesianDynamics::CartesianDynamics(const Config &config, const std::array<size_t, 3> i_ndims,
                                     const unsigned int nsteps)
    : ndims(i_ndims),
      config(config),
      get_wvel(nullwinds()),
      get_uvel(nullwinds()),
      get_vvel(nullwinds()) {
  std::cout << "\n--- coupled cartesian dynamics from file ---\n";

  // -- YAC initialization and calendar definitions ---
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("1850-01-01T00:00:00", "1850-12-31T00:00:00");

  // --- Component definition ---
  std::string component_name = "cleo";
  int component_id = -1;
  yac_cdef_comp(component_name.c_str(), &component_id);

  // --- Grid definition ---
  int grid_id = -1;
  int cell_point_id = -1;
  int edge_point_id = -1;
  std::string grid_name = "cleo_grid";

  create_grid_and_points_definitions(config, ndims, grid_name, grid_id,
                                     cell_point_id, edge_point_id);

  // --- Interpolation stack ---
  int interp_stack_id;
  yac_cget_interp_stack_config(&interp_stack_id);
  yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_AVG, 1, 1.0);

  // --- Field definitions ---
  int num_point_sets = 1;
  int horizontal_fields_collection_size = ndims[2];
  int vertical_winds_collection_size = ndims[2] + 1;

  yac_cdef_field("pressure", component_id, &cell_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &pressure_yac_id);

  yac_cdef_field("temperature", component_id, &cell_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &temp_yac_id);

  yac_cdef_field("qvap", component_id, &cell_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &qvap_yac_id);

  yac_cdef_field("qcond", component_id, &cell_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &qcond_yac_id);

  yac_cdef_field("eastward_wind", component_id, &edge_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &eastward_wind_yac_id);

  yac_cdef_field("northward_wind", component_id, &edge_point_id,
                 num_point_sets, horizontal_fields_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &northward_wind_yac_id);

  yac_cdef_field("vvel", component_id, &cell_point_id,
                 num_point_sets, vertical_winds_collection_size, "PT1M",
                 YAC_TIME_UNIT_ISO_FORMAT, &vvel_yac_id);

  // --- Field coupling definitions ---
  yac_cdef_couple("icon_data_reader", "bubble_grid", "pressure", "cleo", "cleo_grid", "pressure",
                  "PT1M", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "temperature", "cleo", "cleo_grid",
                  "temperature", "PT1M", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE,
                  interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "qvap", "cleo", "cleo_grid", "qvap", "PT1M",
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "qcond", "cleo", "cleo_grid", "qcond", "PT1M",
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "eastward_wind",
                  "cleo", "cleo_grid", "eastward_wind", "PT1M",
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "northward_wind",
                  "cleo", "cleo_grid", "northward_wind", "PT1M",
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  yac_cdef_couple("icon_data_reader", "bubble_grid", "vvel", "cleo", "cleo_grid", "vvel", "PT1M",
                  YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, interp_stack_id, 0, 0);

  // --- End of YAC definitions ---
  yac_cenddef();

  size_t horizontal_cell_number = yac_cget_grid_size(YAC_LOCATION_CELL, grid_id);
  size_t horizontal_edge_number = yac_cget_grid_size(YAC_LOCATION_EDGE, grid_id);

  yac_raw_cell_data = new double * [ndims[2]];
  yac_raw_edge_data = new double * [ndims[2]];
  yac_raw_vertical_wind_data = new double * [ndims[2] + 1];

  for (size_t i = 0; i < ndims[2]; i++) {
    yac_raw_cell_data[i] = new double[horizontal_cell_number];
    yac_raw_edge_data[i] = new double[horizontal_edge_number];
  }

  for (size_t i = 0; i < ndims[2] + 1; i++)
    yac_raw_vertical_wind_data[i] = new double[horizontal_cell_number];

  // Initialization of target containers for receiving data
  press = std::vector<double>(horizontal_cell_number * ndims[2], 0);
  temp = std::vector<double>(horizontal_cell_number * ndims[2], 0);
  qvap = std::vector<double>(horizontal_cell_number * ndims[2], 0);
  qcond = std::vector<double>(horizontal_cell_number * ndims[2], 0);
  uvel = std::vector<double>(ndims[0] * (ndims[1] + 1) * ndims[2], 0);
  wvel = std::vector<double>(ndims[1] * (ndims[0] + 1) * ndims[2], 0);
  vvel = std::vector<double>(horizontal_cell_number * (ndims[2] + 1), 0);

  // Calls the first data retrieval from YAC to have thermodynamic data for first timestep
  receive_field_collections_from_yac();

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
  for (size_t i = 0; i < ndims[2]; i++) {
    delete yac_raw_cell_data[i];
    delete yac_raw_edge_data[i];
  }

  for (size_t i = 0; i < ndims[2] + 1; i++)
    delete yac_raw_vertical_wind_data[i];

  delete [] yac_raw_cell_data;
  delete [] yac_raw_edge_data;
  delete [] yac_raw_vertical_wind_data;

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
    case 2:  // 3-D or 2-D model
      get_uvel = get_uvel_from_yac();
      infoend = ", v" + infoend;
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
