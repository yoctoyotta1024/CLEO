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
  InitSupers(const Config &config)
      : totnsupers(config.totnsupers), nspacedims(config.nspacedims) {}

  auto get_totnsupers() const { return totnsupers; }
  
  auto get_nspacedims() const { return nspacedims; }

  size_t get_size() const
  {
    return sdgbxindex().size();
  }

  std::vector<unsigned int> sdgbxindex() const
  {
    std::vector<unsigned int> sdgbxindex{1, 3, 5, 0, 0,
                                         5, 6, 6, 0, 4};

    // std::vector<unsigned int> sdgbxindex{0, 1, 2, 3, 4,
    //                                      5, 6, 7, 8, 9};

    return sdgbxindex;
  }

  std::vector<double> coord3() const
  {
    std::vector<double> coord3(totnsupers, 0.3);

    return coord3;
  }

  std::vector<double> coord1() const
  {
    std::vector<double> coord1(totnsupers, 0.1);

    return coord1;
  }

  std::vector<double> coord2() const
  {
    std::vector<double> coord2(totnsupers, 0.2);

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

struct InitGbxs
/* struct containing functions which return data
for the initial conditions needed to create
gridboxes e.g. via the CreateGbxs struct */
{
private:
  int ngbxs;
  int getsizetmp; //TODO: erradicate

public:
  InitGbxs(const Config &config)
      : ngbxs(config.ngbxs) {}

  auto get_ngbxs() const { return ngbxs; }
  size_t get_size() const { return getsizetmp; }
};

struct InitConds
/* struct for functions to generate
intial conditions for CLEO */
{
  InitSupers initsupers; // initial conditions for creating superdroplets
  InitGbxs initgbxs; // initial conditions for creating gridboxes

  InitConds(const Config &config)
      : initsupers(config),
        initgbxs(config) {}
};

#endif // INITCONDS_HPP