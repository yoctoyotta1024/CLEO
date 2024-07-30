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

#include <array>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <mpi.h>

#include "cartesiandomain/domainboundaries.hpp"
#include "domainboundaries.hpp"

CartesianDecomposition::CartesianDecomposition() {}
CartesianDecomposition::~CartesianDecomposition() {}

size_t CartesianDecomposition::get_total_global_gridboxes() const {
    return ndims[0] * ndims[1] * ndims[2];
};

size_t CartesianDecomposition::get_total_local_gridboxes() const {
    return total_local_gridboxes;
};

std::array<size_t, 3> CartesianDecomposition::get_local_partition_origin() const {
    return partition_origins[my_rank];
};

std::array<size_t, 3> CartesianDecomposition::get_local_partition_size() const {
    return partition_sizes[my_rank];
};

int CartesianDecomposition::get_gridbox_owner_process(size_t global_gridbox_index) const {
    // Tests whether the gridbox index is out of the domain
    if (global_gridbox_index == outofbounds_gbxindex())
        return -1;

    // Get global gridbox coordinates
    auto gridbox_coordinates = get_coordinates_from_index(ndims, global_gridbox_index);

    // For each process tests whether its partition bounds the gridbox coordinates
    for (size_t process = 0; process < partition_origins.size(); process++)
        if (check_coordinates_inside_partition(gridbox_coordinates, process))
            return process;

    return -1;
};

int CartesianDecomposition::global_to_local_gridbox_index(size_t global_gridbox_index) const {
    // Tests whether the gridbox index is out of the domain
    if (global_gridbox_index == outofbounds_gbxindex())
        return -1;

    // Tests whether the gridbox is owned by the local process
    if (my_rank != get_gridbox_owner_process(global_gridbox_index))
        return -2;

    auto partition_origin = get_local_partition_origin();
    auto partition_size = get_local_partition_size();
    auto global_coordinates = get_coordinates_from_index(ndims, global_gridbox_index);
    std::array<size_t, 3> local_coordinates{global_coordinates[0] - partition_origin[0],
                                            global_coordinates[1] - partition_origin[1],
                                            global_coordinates[2] - partition_origin[2]};


    return get_index_from_coordinates({partition_size[0], partition_size[1], partition_size[2]},
                                      local_coordinates[0],
                                      local_coordinates[1],
                                      local_coordinates[2]);
};

int CartesianDecomposition::local_to_global_gridbox_index(size_t local_gridbox_index) const {
    if (local_gridbox_index > total_local_gridboxes)
        return -1;

    auto partition_size = get_local_partition_size();
    auto partition_origin = get_local_partition_origin();
    auto local_coordinates = get_coordinates_from_index({ partition_size[0],
                                                          partition_size[1],
                                                          partition_size[2] },
                                                        local_gridbox_index);

    return get_index_from_coordinates(ndims,
                                      local_coordinates[0] + partition_origin[0],
                                      local_coordinates[1] + partition_origin[1],
                                      local_coordinates[2] + partition_origin[2]);
};


bool CartesianDecomposition::create(std::vector<size_t> ndims) {
    this->ndims = ndims;
    int comm_size, decomposition_index;
    std::vector<std::vector<size_t>> factorizations;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // If the comm_size is equal to 1 there will be no suitable factorization,
    // so treat it as a special case
    if (comm_size == 1) {
        partition_origins.push_back({0, 0, 0});
        partition_sizes.push_back({ndims[0], ndims[1], ndims[2]});
        total_local_gridboxes = partition_sizes[my_rank][0] *
                                partition_sizes[my_rank][1] *
                                partition_sizes[my_rank][2];
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

    // Saves the origin and sizes of the partitions of all processes
    for (int process = 0; process < comm_size; process++) {
      std::array<size_t, 3> partition_origin;
      std::array<size_t, 3> partition_size;
      construct_partition(ndims, factorizations[decomposition_index],
                          process, partition_origin, partition_size);
      partition_origins.push_back(partition_origin);
      partition_sizes.push_back(partition_size);
    }

    // Sets the number of local gridboxes for convenience
    total_local_gridboxes = partition_sizes[my_rank][0] *
                            partition_sizes[my_rank][1] *
                            partition_sizes[my_rank][2];
    return true;
}


void permute_and_trim_factorizations(std::vector<std::vector<size_t>> & factorizations,
                                     const std::vector<size_t> ndims) {
    int original_factors_size = factorizations.size();

    // Find all permutations for each factorization
    for (int factorization = 0; factorization < original_factors_size; factorization++)
        heapPermutation(factorizations, factorizations[factorization], ndims.size());

    // Remove factorizations whose factors don't fit the dimension sizes
    for (size_t factorization = 0; factorization < factorizations.size();) {
        bool deleted = false;
        for (size_t dimension = 0; dimension < 3; dimension++)
            if (factorizations[factorization][dimension] > ndims[dimension]) {
                factorizations.erase(factorizations.begin() + factorization);
                deleted = true;
                break;
            }
        if (!deleted)
           factorization++;
    }
}

int find_best_decomposition(std::vector<std::vector<size_t>> & factors,
                            const std::vector<size_t> ndims) {
    std::array<size_t, 3> partition_origin, partition_size;
    int comm_size, best_factorization = -1;
    double vertical_split_penalization = 0.5;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Calculates the ideal (most even possible) division and initializes the minimum error
    double ideal_division = static_cast<double>(ndims[0] * ndims[1] * ndims[2]) /
                            static_cast<double>(comm_size);
    double smallest_mean_squared_error = 999999999;

    // For each factorization sums the difference of each process to the ideal division
    for (size_t factorization = 0; factorization < factors.size(); factorization++) {
        double mean_squared_error = 0;
        for (int process = 0; process < comm_size; process++) {
           construct_partition(ndims, factors[factorization],
                               process, partition_origin, partition_size);
           mean_squared_error += sqrt(pow(partition_size[0] *
                                          partition_size[1] *
                                          partition_size[2] - ideal_division, 2));
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
                         int partition_index, std::array<size_t, 3> & partition_origin,
                         std::array<size_t, 3> & partition_size) {
  // Finds the slice index in each dimension for the partition index
  std::array<size_t, 3> slice_indices = {partition_index / (decomposition[1] * decomposition[2]),
                                         (partition_index / decomposition[2]) % decomposition[1],
                                         partition_index % decomposition[2]};

  // Uses the slice indices to calculate the origin and size in each dimension
  for (int dimension = 0; dimension < 3; dimension++) {
      partition_size[dimension] = ndims[dimension] / decomposition[dimension];
      partition_origin[dimension] = partition_size[dimension] * slice_indices[dimension];
      partition_origin[dimension] += std::min(slice_indices[dimension],
                                              ndims[dimension] % decomposition[dimension]);

      // When the division is uneven, spreads the remainder through the first
      // slices (up to the remainder size)
      if (slice_indices[dimension] < ndims[dimension] % decomposition[dimension])
        partition_size[dimension]++;
  }
}

bool CartesianDecomposition::check_coordinates_inside_partition(std::array<size_t, 3> coordinates,
                                                                int partition_index) const {
  bool inside = true;
  // Checks whether the coordinate for each dimension is inside the partition limits
  for (size_t dimension = 0; dimension < partition_origins[partition_index].size(); dimension++) {
    inside = inside &&
             coordinates[dimension] >= partition_origins[partition_index][dimension] &&
             coordinates[dimension] < (partition_origins[partition_index][dimension] +
                                       partition_sizes[partition_index][dimension]);
  }
  return inside;
}

size_t get_index_from_coordinates(const std::vector<size_t> & ndims,
                                  const size_t k, const size_t i, const size_t j) {
  return k + ndims[0] * (i + ndims[1] * j);
}

std::array<size_t, 3> get_coordinates_from_index(const std::vector<size_t> & ndims,
                                                 const size_t index) {
  const size_t j = index / (ndims[0] * ndims[1]);
  const size_t k = index % ndims[0];
  const size_t i = index / ndims[0] - ndims[1] * j;

  return std::array<size_t, 3>{k, i, j};
}

std::vector<std::vector<size_t>> factorize(int n) {
    std::vector<std::vector<size_t>> result;
    std::vector<size_t> current;
    factorizeHelper(n, 2, current, result);  // Start from 1 to include trivial factorizations
    return result;
}

void factorizeHelper(int n, int start, std::vector<size_t> & current,
                     std::vector<std::vector<size_t>> & result) {
    if (n == 1) {
        if (!current.empty()) {  // We want to consider all factorizations including those with 1
            result.push_back(current);
        }
        return;
    }

    for (int i = start; i <= n; ++i) {
        if (n % i == 0) {
            current.push_back(i);
            factorizeHelper(n / i, i, current, result);
            current.pop_back();
        }
    }
}

void heapPermutation(std::vector<std::vector<size_t>> & results,
                     std::vector<size_t> arr, int size) {
    if (size == 1) {
        for (auto i : results)
            if (arr == i)
                return;
        results.push_back(arr);
        return;
    }

    for (int i = 0; i < size; i++) {
        heapPermutation(results, arr, size - 1);

        // Swap logic
        if (size % 2 == 1) {
            std::swap(arr[0], arr[size - 1]);
        } else {
            std::swap(arr[i], arr[size - 1]);
        }
    }
}
