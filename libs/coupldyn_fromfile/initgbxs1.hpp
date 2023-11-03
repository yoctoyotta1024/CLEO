/*
 * ----- CLEO -----
 * File: initgbxs1.hpp
 * Project: coupldyn_fromfile
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for griboxes'
 * initial conditions for CLEO SDM
 * (e.g. thermodynamics) which can be used 
 * by InitConds struct as GbxInitConds type
 */

#ifndef INITGBXS1_HPP
#define INITGBXS1_HPP

#include <vector>
#include <utility>

#include "initialise/config.hpp"

// TODO 

struct InitGbxs1
/* struct containing functions which return data
for the initial conditions needed to create
gridboxes e.g. via the create_gbxs function */
{
private:
  size_t ngbxs;

public:
  InitGbxs1(const Config &config)
      : ngbxs(config.ngbxs) {}

  auto get_ngbxs() const { return ngbxs; }

  // size_t get_size() const
  // {
  //   return volume().size();
  // }

  std::vector<double> temp() const
  {
    std::vector<double> temp(ngbxs, 1000);

    return temp;
  }

  std::vector<double> press() const
  {
    std::vector<double> press(ngbxs, 1000);

    return press;
  }

  std::vector<double> qvap() const
  {
    std::vector<double> qvap(ngbxs, 1000);

    return qvap;
  }

  std::vector<double> qcond() const
  {
    std::vector<double> qcond(ngbxs, 1000);

    return qcond;
  }

  std::vector<std::pair<double, double>> wvel() const
  {
    auto w = std::make_pair<double, double>(0.0, 0.0);
    std::vector<std::pair<double, double>> wvel(ngbxs, w);

    return wvel;
  }

  std::vector<std::pair<double, double>> uvel() const
  {
    auto u = std::make_pair<double, double>(0.0, 0.0);
    std::vector<std::pair<double, double>> uvel(ngbxs, u);

    return uvel;
  }

  std::vector<std::pair<double, double>> vvel() const
  {
    auto v = std::make_pair<double, double>(0.0, 0.0);
    std::vector<std::pair<double, double>> vvel(ngbxs, v);

    return vvel;
  }
};

#endif // INITGBXS1_HPP