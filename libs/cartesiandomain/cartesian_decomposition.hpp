/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_decomposition.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 30th July 2023
 * Author: Wilton J. Loch
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class and suporting functions to perform domain decomposition
 * in a cartesian domain
 */

#ifndef LIBS_CARTESIANDOMAIN_CARTESIAN_DECOMPOSITION_HPP_
#define LIBS_CARTESIANDOMAIN_CARTESIAN_DECOMPOSITION_HPP_

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "../cleoconstants.hpp"
#include "cartesiandomain/domainboundaries.hpp"
#include "configuration/communicator.hpp"
#include "initialise/gbx_bounds_from_binary.hpp"

class CartesianDecomposition {
 private:
  int my_rank;
  // Number of dimensions of the global domain
  std::vector<size_t> ndims;
  // Global origins of all partitions
  std::vector<std::array<size_t, 3>> partition_origins;
  // Sizes of all partitions
  std::vector<std::array<size_t, 3>> partition_sizes;

  // Geometrical coordinates of the local partition begin and end
  std::array<double, 3> partition_begin_coordinates;
  std::array<double, 3> partition_end_coordinates;

  // Geometric bounds of a gridbox in z, x, y directions
  std::vector<std::vector<double>> gridbox_bounds;

  // Geometric bounds of the entire domain
  // First Index 0 : Lower bounds in z, x, y directions
  // First Index 1 : Upper bounds in z, x, y directions
  std::array<std::array<double, 3>, 2> domain_bounds;

  // Behavior of each dimension, being either periodic or finite
  std::array<size_t, 3> dimension_bound_behavior;

  // Which process neighbors the current one in each direction
  // (it can also be the same as the local one)
  std::map<std::array<int, 3>, int> neighboring_processes;

  // Domain decomposition factors for each dimension
  // (i.e. in how many parts each dimension is divided)
  std::array<size_t, 3> decomposition;
  // Number of local gridboxes
  size_t total_local_gridboxes;

  // Fill the partition_begin_coordinates and partition_end_coordinates arrays
  void calculate_partition_coordinates();

  // Fills the neighboring_processes array
  void calculate_neighboring_processes();

 public:
  CartesianDecomposition();
  ~CartesianDecomposition();

  // Creates the decomposition
  bool create(std::vector<size_t> ndims, GbxBoundsFromBinary gfb);

  // Local and global amount of gridboxes
  size_t get_total_local_gridboxes() const;
  size_t get_total_global_gridboxes() const;

  // Get the origin and size of local partition in terms of number of gridboxes
  std::array<size_t, 3> get_local_partition_origin() const;
  std::array<size_t, 3> get_local_partition_size() const;

  // Get partition index and partition coordinates
  int get_partition_index_from_slice(std::array<int, 3> slice_indices) const;
  std::array<int, 3> get_slice_indices_from_partition(int partition_index) const;
  // Checks whether a coordinate is bounded by one specific partition
  bool check_indices_inside_partition(std::array<size_t, 3> indices, int partition_index) const;

  // Gridbox related subroutines
  int local_to_global_gridbox_index(size_t local_gridbox_index, int process = -1) const;
  int global_to_local_gridbox_index(size_t global_gridbox_index) const;
  int get_gridbox_owner_process(size_t global_gridbox_index) const;
  unsigned int get_local_bounding_gridbox_index(std::array<double, 3> &coordinates) const;
  // void set_gridbox_size(double z_size, double x_size, double y_size);
  void set_gridbox_bounds(GbxBoundsFromBinary gfb);
  // Sets the behavior of all dimensions
  void set_dimensions_bound_behavior(std::array<size_t, 3> behaviors);
};

// Given the global domain, a global decomposition and a partition index,
// returns the partition origin and size
void construct_partition(const std::vector<size_t> ndims, std::vector<size_t> decomposition,
                         int partition_index, std::array<size_t, 3> &partition_origin,
                         std::array<size_t, 3> &partition_size);

// Adds all permutations of a particular decomposition and removes the ones that
// do not fit the global dimension sizes
void permute_and_trim_factorizations(std::vector<std::vector<size_t>> &factors,
                                     const std::vector<size_t> ndims);

// Finds the best decomposition given by the most even division of gridboxes among processes
int find_best_decomposition(std::vector<std::vector<size_t>> &factors,
                            const std::vector<size_t> ndims);

// Functions for getting indexes and coordinates in an arbitrary 3D gridbox domain
size_t get_index_from_coordinates(const std::vector<size_t> &ndims, const size_t k, const size_t i,
                                  const size_t j);
std::array<size_t, 3> get_coordinates_from_index(const std::vector<size_t> &ndims,
                                                 const size_t index);

// Support functions
std::vector<std::vector<size_t>> factorize(int n);
void factorize_helper(int n, int start, std::vector<size_t> &current,
                      std::vector<std::vector<size_t>> &result);
void heap_permutation(std::vector<std::vector<size_t>> &results, std::vector<size_t> arr, int size);
int get_multiplications_to_turn_int(double entry_value);
int binary_search(std::array<double, 3> &coordinates,
    int dimension, std::array<size_t, 3> partition_size,
    std::vector<std::vector<double>> gridbox_bounds);
#endif  // LIBS_CARTESIANDOMAIN_CARTESIAN_DECOMPOSITION_HPP_
