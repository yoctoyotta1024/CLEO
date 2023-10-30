/*
 * ----- CLEO -----
 * File: initsupers1.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 30th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for superdroplets' 
 * initial conditions for CLEO SDM
 * (e.g. superdroplet attributes)
 * which can be used by InitConds
 * struct as SuperdropInitConds type
 */

#ifndef INITSUPERS1_HPP
#define INITSUPERS1_HPP

#include <vector>

#include "./config.hpp"

// TODO 

struct InitSupers1
/* struct containing functions which return data
for the initial conditions needed to create
superdroplets e.g. via the CreateSupers struct */
{
private:
  size_t totnsupers; // total number of superdroplets (in kokkos view on device initially)
  unsigned int nspacedims; // number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D)

public:
  InitSupers1(const Config &config)
      : totnsupers(config.totnsupers), nspacedims(config.nspacedims) {}

  auto get_totnsupers() const { return totnsupers; }

  auto get_nspacedims() const { return nspacedims; }

  size_t get_size() const
  {
    return sdgbxindex().size();
  }

  std::vector<unsigned int> sdgbxindex() const
  {
    std::vector<unsigned int> sdgbxindex{0, 0, 0, 0, 0};

    // std::vector<unsigned int> sdgbxindex{1, 3, 5, 0, 0,
    //                                      5, 6, 6, 0, 4};

    return sdgbxindex;
  }

  std::vector<double> coord3() const
  {
    std::vector<double> coord3(totnsupers, 0.0);

    return coord3;
  }

  std::vector<double> coord1() const
  {
    std::vector<double> coord1(totnsupers, 0.0);

    return coord1;
  }

  std::vector<double> coord2() const
  {
    std::vector<double> coord2(totnsupers, 0.0);

    return coord2;
  }

  std::vector<double> radius() const
  {
    const double r(0.1e-6 / dlc::R0);
    std::vector<double> radius(totnsupers, r);
    return radius;
  }

  std::vector<double> msol() const
  {
    const double m(8.446695447951756e-18 / dlc::MASS0);
    std::vector<double> msol(totnsupers, m);

    return msol;
  }

  std::vector<unsigned long long> xi() const
  {
    const unsigned long long x(100000000);
    std::vector<unsigned long long> xi(totnsupers, x);

    return xi;
  }
};

#endif // INITSUPERS1_HPP