/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_decomposition.cpp
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

#include "cartesian_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "cartesiandomain/domainboundaries.hpp"
#include "domainboundaries.hpp"

CartesianDecomposition::CartesianDecomposition() {}
CartesianDecomposition::~CartesianDecomposition() {}

size_t CartesianDecomposition::get_total_global_gridboxes() const {
  return ndims[0] * ndims[1] * ndims[2];
};

size_t CartesianDecomposition::get_total_local_gridboxes() const { return total_local_gridboxes; };

std::array<size_t, 3> CartesianDecomposition::get_local_partition_origin() const {
  return partition_origins[my_rank];
};

std::array<size_t, 3> CartesianDecomposition::get_local_partition_size() const {
  return partition_sizes[my_rank];
};

void CartesianDecomposition::set_gridbox_size(double z_size, double x_size, double y_size) {
  gridbox_size[0] = z_size;
  gridbox_size[1] = x_size;
  gridbox_size[2] = y_size;
};

// Defines if a dimension is periodic or finite
void CartesianDecomposition::set_dimensions_bound_behavior(std::array<size_t, 3> behaviors) {
  dimension_bound_behavior = behaviors;
}

int CartesianDecomposition::get_partition_index_from_slice(std::array<int, 3> slice_indices) const {
  return slice_indices[0] * (decomposition[1] * decomposition[2]) +
         slice_indices[1] * decomposition[2] + slice_indices[2];
}

std::array<int, 3> CartesianDecomposition::get_slice_indices_from_partition(
    int partition_index) const {
  return {static_cast<int>(partition_index / (decomposition[1] * decomposition[2])),
          static_cast<int>((partition_index / decomposition[2]) % decomposition[1]),
          static_cast<int>(partition_index % decomposition[2])};
}

bool CartesianDecomposition::check_indices_inside_partition(std::array<size_t, 3> indices,
                                                            int partition_index) const {
  bool inside = true;
  // Checks for each dimension whether the gridbox coordinate indices are inside
  // the partition limits
  for (size_t dimension = 0; dimension < partition_origins[partition_index].size(); dimension++) {
    inside = inside && indices[dimension] >= partition_origins[partition_index][dimension] &&
             indices[dimension] < (partition_origins[partition_index][dimension] +
                                   partition_sizes[partition_index][dimension]);
  }
  return inside;
}

unsigned int CartesianDecomposition::get_local_bounding_gridbox(
    std::array<double, 3> &coordinates) const {
  std::array<size_t, 3> bounding_gridbox_coordinates;
  std::array<int, 3> external_direction = {0, 0, 0};
  bool local_coordinate = true;

  for (auto dimension : {0, 1, 2}) {
    // Tests whether the coordinate in that dimension is smaller than the
    // beginning of the partition
    if (coordinates[dimension] < partition_begin_coordinates[dimension]) {
      // If the dimension behavior is finite and the coordinate is smaller than
      // the beginning of the domain return the out_of_bounds value
      if (dimension_bound_behavior[dimension] == 0 && coordinates[dimension] < 0)
        return LIMITVALUES::oob_gbxindex;

      // The coordinate is inside of the domain but outside of the partition in that dimension
      external_direction[dimension] -= 1;
      local_coordinate = false;

      // Tests whether the coordinate for that dimension is larger than the end of the partition
    } else if (coordinates[dimension] > partition_end_coordinates[dimension]) {
      // If the dimension behavior is finite and the coordinate is larger than
      // the end of the domain return the out_of_bounds value
      if (dimension_bound_behavior[dimension] == 0 &&
          coordinates[dimension] > ndims[dimension] * gridbox_size[dimension])
        return LIMITVALUES::oob_gbxindex;

      // The coordinate is inside of the domain but outside of the partition in that dimension
      external_direction[dimension] += 1;
      local_coordinate = false;

      // If none of the tests above pass then the coordinate is inside of the
      // partition in that dimension
    } else if (local_coordinate) {
      int multiplications =
          std::max({get_multiplications_to_turn_int(coordinates[dimension]),
                    get_multiplications_to_turn_int(partition_begin_coordinates[dimension]),
                    get_multiplications_to_turn_int(gridbox_size[dimension])});
      int64_t integer_coordinate =
          std::round(coordinates[dimension] * std::pow(10, multiplications));
      int64_t integer_domain_begin =
          std::round(partition_begin_coordinates[dimension] * std::pow(10, multiplications));
      int64_t integer_gridbox_size =
          std::round(gridbox_size[dimension] * std::pow(10, multiplications));
      bounding_gridbox_coordinates[dimension] =
          (integer_coordinate - integer_domain_begin) / integer_gridbox_size;
    }
  }

  // If the coordinate is inside of the partition in all dimensions, returns the
  // index of its bounding gridbox
  if (local_coordinate) {
    return get_index_from_coordinates({get_local_partition_size()[0], get_local_partition_size()[1],
                                       get_local_partition_size()[2]},
                                      bounding_gridbox_coordinates[0],
                                      bounding_gridbox_coordinates[1],
                                      bounding_gridbox_coordinates[2]);

    // Otherwise, the coordinate is outside of the partition
  } else {
    bool corrected = false;
    for (auto dimension : {0, 1, 2}) {
      // Since the finite dimensions have already been checked before,
      // if the coordinate is outside of the domain for a dimension correct it to go around
      if (coordinates[dimension] < 0) {
        coordinates[dimension] += ndims[dimension] * gridbox_size[dimension];
        corrected = true;
      } else if (coordinates[dimension] > ndims[dimension] * gridbox_size[dimension]) {
        coordinates[dimension] -= ndims[dimension] * gridbox_size[dimension];
        corrected = true;
      }
    }

    // If the coordinates have been corrected but the target partition is local
    // then get the local bounding gridbox
    if (corrected && neighboring_processes.at(external_direction) == my_rank)
      return get_local_bounding_gridbox(coordinates);

    // If the coordinate is outside of the local partition, encode in the return
    // value which process contains it
    return (LIMITVALUES::oob_gbxindex - 1) - neighboring_processes.at(external_direction);
  }
}

// Given a global gridbox index returns which process ows it
int CartesianDecomposition::get_gridbox_owner_process(size_t global_gridbox_index) const {
  // Tests whether the gridbox index is out of the domain
  if (global_gridbox_index == LIMITVALUES::oob_gbxindex) return -1;

  // Get global gridbox coordinates
  auto gridbox_coordinates = get_coordinates_from_index(ndims, global_gridbox_index);

  // For each process tests whether its partition bounds the gridbox coordinates
  for (size_t process = 0; process < partition_origins.size(); process++)
    if (check_indices_inside_partition(gridbox_coordinates, process)) return process;

  return -1;
};

// Given a global gridbox index returns the corresponding local gridbox index
int CartesianDecomposition::global_to_local_gridbox_index(size_t global_gridbox_index) const {
  // Tests whether the gridbox index is out of the domain
  if (global_gridbox_index == LIMITVALUES::oob_gbxindex) return -1;

  // Tests whether the gridbox is owned by the local process
  if (my_rank != get_gridbox_owner_process(global_gridbox_index)) return -1;

  auto partition_origin = get_local_partition_origin();
  auto partition_size = get_local_partition_size();
  auto global_coordinates = get_coordinates_from_index(ndims, global_gridbox_index);
  std::array<size_t, 3> local_coordinates{global_coordinates[0] - partition_origin[0],
                                          global_coordinates[1] - partition_origin[1],
                                          global_coordinates[2] - partition_origin[2]};

  return get_index_from_coordinates({partition_size[0], partition_size[1], partition_size[2]},
                                    local_coordinates[0], local_coordinates[1],
                                    local_coordinates[2]);
};

// Given a local gridbox index returns the corresponding global gridbox index
int CartesianDecomposition::local_to_global_gridbox_index(size_t local_gridbox_index,
                                                          int process) const {
  // if no process is specified assumes the local one
  if (process == -1) process = my_rank;

  if (process == my_rank && local_gridbox_index > total_local_gridboxes) return -1;

  auto partition_size = partition_sizes[process];
  auto partition_origin = partition_origins[process];
  auto local_coordinates = get_coordinates_from_index(
      {partition_size[0], partition_size[1], partition_size[2]}, local_gridbox_index);

  return get_index_from_coordinates(ndims, local_coordinates[0] + partition_origin[0],
                                    local_coordinates[1] + partition_origin[1],
                                    local_coordinates[2] + partition_origin[2]);
};

// Main subroutine for the creation of the decomposition
bool CartesianDecomposition::create(std::vector<size_t> ndims, double gridbox_z_size,
                                    double gridbox_x_size, double gridbox_y_size) {
  this->ndims = ndims;
  set_gridbox_size(gridbox_z_size, gridbox_x_size, gridbox_y_size);

  int comm_size, decomposition_index;
  std::vector<std::vector<size_t>> factorizations;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // If the comm_size is equal to 1 there will be no suitable factorization,
  // so treat it as a special case
  if (comm_size == 1) {
    partition_origins.push_back({0, 0, 0});
    partition_sizes.push_back({ndims[0], ndims[1], ndims[2]});
    total_local_gridboxes =
        partition_sizes[my_rank][0] * partition_sizes[my_rank][1] * partition_sizes[my_rank][2];
    decomposition = {1, 1, 1};

    calculate_partition_coordinates();
    calculate_neighboring_processes();

    return true;
  }

  // Get all the possible factorizations of the total number of processes
  factorizations = factorize(comm_size);

  // conforms all factorizations to the number of dimensions by deleting factorizations larger
  // than ndims size and padding factorizaions smaller than ndims size
  for (size_t factorization = 0; factorization < factorizations.size();)
    if (factorizations[factorization].size() > ndims.size())
      factorizations.erase(factorizations.begin() + factorization);
    else {
      while (factorizations[factorization].size() < ndims.size())
        factorizations[factorization].push_back(1);
      factorization++;
    }

  // Gets all the permutations of the factorizations and removes the ones that
  // do not fit the global domain
  permute_and_trim_factorizations(factorizations, ndims);

  // Raise an error if there are no decompositions left after trimming
  assert(!factorizations.empty() &&
         "No domain decomposition found for the number of gridboxes and processes");

  // Finds the best (most even) decomposition of gridboxes among processes
  decomposition_index = find_best_decomposition(factorizations, ndims);

  decomposition = {factorizations[decomposition_index][0], factorizations[decomposition_index][1],
                   factorizations[decomposition_index][2]};

  // Saves the origin and sizes of the partitions of all processes
  for (int process = 0; process < comm_size; process++) {
    std::array<size_t, 3> partition_origin;
    std::array<size_t, 3> partition_size;
    construct_partition(ndims, factorizations[decomposition_index], process, partition_origin,
                        partition_size);
    partition_origins.push_back(partition_origin);
    partition_sizes.push_back(partition_size);
  }

  // Sets the number of local gridboxes for convenience
  total_local_gridboxes =
      partition_sizes[my_rank][0] * partition_sizes[my_rank][1] * partition_sizes[my_rank][2];

  calculate_partition_coordinates();
  calculate_neighboring_processes();

  return true;
}

// Calculates the neighboring processes in all directions
void CartesianDecomposition::calculate_neighboring_processes() {
  auto my_slice_indices = get_slice_indices_from_partition(my_rank);
  std::array<int, 3> target_slice = my_slice_indices;

  // cycles through all possible directions that a superdroplet can move
  for (auto k : {-1, 0, 1})
    for (auto i : {-1, 0, 1})
      for (auto j : {-1, 0, 1}) {
        if (i != 0 || j != 0 || k != 0) {
          // calculates indices of domain slices from the neighbors
          target_slice[0] = my_slice_indices[0] + k;
          target_slice[1] = my_slice_indices[1] + i;
          target_slice[2] = my_slice_indices[2] + j;

          // corrects neighbors that wrap (go around) the domain
          for (auto dimension : {0, 1, 2}) {
            if (target_slice[dimension] < 0) target_slice[dimension] = decomposition[dimension] - 1;
            if (target_slice[dimension] == static_cast<int>(decomposition[dimension]))
              target_slice[dimension] = 0;
          }

          // saves the neighboring process for the target domain slice
          neighboring_processes[{k, i, j}] = get_partition_index_from_slice(target_slice);
        }
      }
}

// Calculates the geometrical coordinates of the beginning and end of the local partition
void CartesianDecomposition::calculate_partition_coordinates() {
  int my_rank;
  auto partition_origin = get_local_partition_origin();
  auto partition_size = get_local_partition_size();

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  for (int dimension = 0; dimension < 3; dimension++) {
    int multiplications = get_multiplications_to_turn_int(gridbox_size[dimension]);
    int64_t integer_gridbox_size =
        std::round(gridbox_size[dimension] * std::pow(10, multiplications));
    partition_begin_coordinates[dimension] = partition_origin[dimension] * integer_gridbox_size;
    partition_begin_coordinates[dimension] /= std::pow(10, multiplications);
  }

  partition_end_coordinates[0] =
      partition_begin_coordinates[0] + partition_size[0] * gridbox_size[0];
  partition_end_coordinates[1] =
      partition_begin_coordinates[1] + partition_size[1] * gridbox_size[1];
  partition_end_coordinates[2] =
      partition_begin_coordinates[2] + partition_size[2] * gridbox_size[2];
}

// Returns how many multiplications by 10 are needed to turn a double to int
int get_multiplications_to_turn_int(double entry_value) {
  int total_multiplications = 0;
  while (std::round(entry_value) != entry_value) {
    entry_value *= 1e1;
    total_multiplications++;
  }
  return total_multiplications;
}

void permute_and_trim_factorizations(std::vector<std::vector<size_t>> &factorizations,
                                     const std::vector<size_t> ndims) {
  int original_factors_size = factorizations.size();

  // Find all permutations for each factorization
  for (int factorization = 0; factorization < original_factors_size; factorization++)
    heap_permutation(factorizations, factorizations[factorization], ndims.size());

  // Remove factorizations whose factors don't fit the dimension sizes
  for (size_t factorization = 0; factorization < factorizations.size();) {
    bool deleted = false;
    for (size_t dimension = 0; dimension < 3; dimension++)
      if (factorizations[factorization][dimension] > ndims[dimension]) {
        factorizations.erase(factorizations.begin() + factorization);
        deleted = true;
        break;
      }
    if (!deleted) factorization++;
  }
}

int find_best_decomposition(std::vector<std::vector<size_t>> &factors,
                            const std::vector<size_t> ndims) {
  std::array<size_t, 3> partition_origin, partition_size;
  int comm_size, best_factorization = -1;
  double vertical_split_penalization = 1.0;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // Calculates the ideal (most even possible) division and initializes the minimum error
  double ideal_division =
      static_cast<double>(ndims[0] * ndims[1] * ndims[2]) / static_cast<double>(comm_size);
  double smallest_mean_squared_error = 999999999;

  // For each factorization sums the difference of each process to the ideal division
  for (size_t factorization = 0; factorization < factors.size(); factorization++) {
    double mean_squared_error = 0;
    for (int process = 0; process < comm_size; process++) {
      construct_partition(ndims, factors[factorization], process, partition_origin, partition_size);
      mean_squared_error +=
          sqrt(pow(partition_size[0] * partition_size[1] * partition_size[2] - ideal_division, 2));
    }
    // Penalizes decompositions which split the vertical dimension
    mean_squared_error *= std::pow(factors[factorization][0], vertical_split_penalization);

    // Saves the smallest mean found so far and the index of the corresponding decomposition
    if (mean_squared_error < smallest_mean_squared_error) {
      smallest_mean_squared_error = mean_squared_error;
      best_factorization = factorization;
    }
  }

  return best_factorization;
}

void construct_partition(const std::vector<size_t> ndims, std::vector<size_t> decomposition,
                         int partition_index, std::array<size_t, 3> &partition_origin,
                         std::array<size_t, 3> &partition_size) {
  // Finds the slice index in each dimension for the partition index
  std::array<size_t, 3> slice_indices = {partition_index / (decomposition[1] * decomposition[2]),
                                         (partition_index / decomposition[2]) % decomposition[1],
                                         partition_index % decomposition[2]};

  // Uses the slice indices to calculate the origin and size in each dimension
  for (int dimension = 0; dimension < 3; dimension++) {
    partition_size[dimension] = ndims[dimension] / decomposition[dimension];
    partition_origin[dimension] = partition_size[dimension] * slice_indices[dimension];
    partition_origin[dimension] +=
        std::min(slice_indices[dimension], ndims[dimension] % decomposition[dimension]);

    // When the division is uneven, spreads the remainder through the first
    // slices (up to the remainder size)
    if (slice_indices[dimension] < ndims[dimension] % decomposition[dimension])
      partition_size[dimension]++;
  }
}

// Returns the index of a gridbox inside a domain of arbitrary size given the
// gridbox coordinate indices
size_t get_index_from_coordinates(const std::vector<size_t> &ndims, const size_t k, const size_t i,
                                  const size_t j) {
  return k + ndims[0] * (i + ndims[1] * j);
}

// Returns the coordinate indices of a gridbox inside a domain of arbitrary size
// given the gridbox index
std::array<size_t, 3> get_coordinates_from_index(const std::vector<size_t> &ndims,
                                                 const size_t index) {
  const size_t j = index / (ndims[0] * ndims[1]);
  const size_t k = index % ndims[0];
  const size_t i = index / ndims[0] - ndims[1] * j;

  return std::array<size_t, 3>{k, i, j};
}

// Returns all possible factorizations of an integer
std::vector<std::vector<size_t>> factorize(int n) {
  std::vector<std::vector<size_t>> result;
  std::vector<size_t> current;
  factorize_helper(n, 2, current, result);  // Start from 1 to include trivial factorizations
  return result;
}

// Helper function for the integer factorization
void factorize_helper(int n, int start, std::vector<size_t> &current,
                      std::vector<std::vector<size_t>> &result) {
  if (n == 1) {
    if (!current.empty()) {  // We want to consider all factorizations including those with 1
      result.push_back(current);
    }
    return;
  }

  for (int i = start; i <= n; ++i) {
    if (n % i == 0) {
      current.push_back(i);
      factorize_helper(n / i, i, current, result);
      current.pop_back();
    }
  }
}

// Implements the heap permutation algorithm (https://en.wikipedia.org/wiki/Heap%27s_algorithm)
// It is used to permute the factorizations through the dimensions
void heap_permutation(std::vector<std::vector<size_t>> &results, std::vector<size_t> arr,
                      int size) {
  if (size == 1) {
    for (auto i : results)
      if (arr == i) return;
    results.push_back(arr);
    return;
  }

  for (int i = 0; i < size; i++) {
    heap_permutation(results, arr, size - 1);

    // Swap logic
    if (size % 2 == 1) {
      std::swap(arr[0], arr[size - 1]);
    } else {
      std::swap(arr[i], arr[size - 1]);
    }
  }
}
