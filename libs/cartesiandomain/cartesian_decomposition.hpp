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

#include <cstddef>
#include <array>
#include <vector>

class CartesianDecomposition {
 private:
    int my_rank;
    // Number of dimensions of the global domain
    std::vector<size_t> ndims;
    // Global origins of all partitions
    std::vector<std::array<size_t, 3>> partition_origins;
    // Sizes of all partitions
    std::vector<std::array<size_t, 3>> partition_sizes;
    // Number of local gridboxes
    size_t total_local_gridboxes;

 public:
    CartesianDecomposition();
    ~CartesianDecomposition();

    // Creates the decomposition
    bool create(std::vector<size_t> ndims);

    size_t get_total_local_gridboxes() const;
    size_t get_total_global_gridboxes() const;
    std::array<size_t, 3> get_local_partition_origin() const;
    std::array<size_t, 3> get_local_partition_size() const;
    int get_gridbox_owner_process(size_t global_gridbox_index) const;
    int local_to_global_gridbox_index(size_t local_gridbox_index) const;
    int global_to_local_gridbox_index(size_t global_gridbox_index) const;

    // Checks whether a coordinate is bounded by one specific partition
    bool check_coordinates_inside_partition(std::array<size_t, 3> coordinates,
                                            int partition_index) const;
};

// Given the global domain, a global decomposition and a partition index,
// returns the partition origin and size
void construct_partition(const std::vector<size_t> ndims,
                         std::vector<size_t> decomposition, int partition_index,
                         std::array<size_t, 3> &partition_origin,
                         std::array<size_t, 3> &partition_size);

// Adds all permutations of a particular decomposition and removes the ones that
// do not fit the global dimension sizes
void permute_and_trim_factorizations(std::vector<std::vector<size_t>> &factors,
                                     const std::vector<size_t> ndims);

// Finds the best decomposition given by the most even division of gridboxes among processes
int find_best_decomposition(std::vector<std::vector<size_t>> &factors,
                            const std::vector<size_t> ndims);

size_t get_index_from_coordinates(const std::vector<size_t> &ndims, const size_t k, const size_t i,
                                  const size_t j);
std::array<size_t, 3> get_coordinates_from_index(const std::vector<size_t> &ndims,
                                                 const size_t index);
std::vector<std::vector<size_t>> factorize(int n);
void factorizeHelper(int n, int start, std::vector<size_t> &current,
                     std::vector<std::vector<size_t>> &result);
void heapPermutation(std::vector<std::vector<size_t>> &results, std::vector<size_t> arr,
                     int size);

#endif  // LIBS_CARTESIANDOMAIN_CARTESIAN_DECOMPOSITION_HPP_
