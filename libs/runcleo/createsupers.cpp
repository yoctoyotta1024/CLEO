/*
 * ----- CLEO -----
 * File: createsupers.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
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
 * file for functions to create a view of
 * superdroplets (on device) from some
 * initial conditions
 */

#include "./createsupers.hpp"

void CreateSupers::ensure_initialisation_complete(viewd_constsupers supers,
                                                  const size_t size) const
/* ensure the number of superdrops in the view matches the
size according to the initial conditions */
{
  if (supers.extent(0) < size)
  {
    const std::string err("Fewer superdroplets were created than were"
                          " given by initialisation data ie. " +
                          std::to_string(supers.extent(0)) + " < " +
                          std::to_string(size));
    throw std::invalid_argument(err);
  }

  if (is_sorted(supers) == 0)
  {
    const std::string err("supers ordered incorrectly "
                          "(ie. not sorted by asceding sdgbxindex");
    throw std::invalid_argument(err);
  }
}

std::array<double, 3>
CreateSupers::GenSuperdrop::coords_at(const unsigned int kk) const
/* returns superdroplet spatial coordinates. A coordinate is
only copied from the corresponding coords vector if that
coordinate is consistent with number of spatial dimensions of
model. Otherwise coordinate = 0. E.g. if model is 1-D,
only coord3 obtained from vectorr (coord1 = coord2 = 0.0) */
{
  std::array<double, 3> coords312{0.0, 0.0, 0.0};

  if (nspacedims > 0)
  {
    coords312[0] = coord3s.at(kk);

    if (nspacedims > 1)
    {
      coords312[1] = coord1s.at(kk);

      if (nspacedims == 3)
      {
        coords312[2] = coord2s.at(kk);
      }
    }
  }

  return coords312;
}