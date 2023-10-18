/*
 * ----- CLEO -----
 * File: initconds.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for initial conditions for CLEO SDM
 * e.g. for superdroplet attributes
 */

#ifndef INITCONDS_HPP
#define INITCONDS_HPP

#include "./config.hpp"

struct InitSupers
/* struct containing functions which return data
for the initial conditions needed to create
superdroplets e.g. via the CreateSupers struct */
{
private:
  int totnsupers; // total number of superdroplets (in kokkos view on device initially)
  int nspacedims; // number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D)

public:
  InitSupers(const int totnsupers, const int nspacedims)
      : totnsupers(totnsupers), nspacedims(nspacedims) {}

  int get_totnsupers() const { return totnsupers; }

  size_t get_size() const
  {
    return sdgbxindex().size();
  }

  std::vector<unsigned int> sdgbxindex() const
  {
    std::vector<unsigned int> sdgbxindex{1, 3, 5, 0, 0,
                                         5, 6, 6, 0, 4};

    return sdgbxindex;
  }

  std::vector<double> coord3() const
  {
    std::vector<double> coord3(totnsupers);

    return coord3;
  }

  std::vector<double> coord1() const
  {
    std::vector<double> coord1(totnsupers);

    return coord1;
  }

  std::vector<double> coord2() const
  {
    std::vector<double> coord2(totnsupers);

    return coord2;
  }

  std::vector<double> radius() const
  {
    std::vector<double> radius(totnsupers);

    return radius;
  }

  std::vector<double> msol() const
  {
    std::vector<double> msol(totnsupers);

    return msol;
  }

  std::vector<unsigned long long> xi() const
  {
    std::vector<unsigned long long> xi(totnsupers);

    return xi;
  }
};

struct InitConds
/* struct for functions to generate
intial conditions for CLEO */
{
  InitSupers initsupers; // initial conditions for creating superdroplets

  InitConds(const Config &config)
      : initsupers(config.totnsupers, config.nspacedims) {}
};

#endif // INITCONDS_HPP