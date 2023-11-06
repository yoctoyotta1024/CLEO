/*
 * ----- CLEO -----
 * File: createsupers.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 6th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * some functionality for structure(s) to
 * create a view of superdroplets (on device)
 * using some initial conditions
 */

#include "./createsupers.hpp"

void is_sdsinit_complete(const viewd_constsupers supers,
                        const size_t size)
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

void print_supers(const viewd_constsupers supers)
/* print superdroplet information */
{
  auto h_supers = Kokkos::create_mirror_view(supers); // mirror of supers in case view is on device memory
  Kokkos::deep_copy(h_supers, supers);

  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    std::cout << "SD: " << h_supers(kk).sdId.value
              << " [gbx, (coords), (attrs)]: [ "
              << h_supers(kk).get_sdgbxindex() << ", ("
              << h_supers(kk).get_coord3() << ", "
              << h_supers(kk).get_coord1() << ", "
              << h_supers(kk).get_coord2() << "), ("
              << h_supers(kk).is_solute() << ", "
              << h_supers(kk).get_radius() << ", "
              << h_supers(kk).get_msol() << ", "
              << h_supers(kk).get_xi() << ") ] \n";
  }
}

std::array<double, 3>
GenSuperdrop::coords_at(const unsigned int kk) const
/* returns superdroplet spatial coordinates. A coordinate is
only copied from the corresponding coords vector if that
coordinate is consistent with number of spatial dimensions of
model. Otherwise coordinate = 0. E.g. if model is 1-D,
only coord3 obtained from vectorr (coord1 = coord2 = 0.0) */
{
  std::array<double, 3> coords312{0.0, 0.0, 0.0};

  switch (nspacedims)
  {
  case 3: // 3-D model
    coords312[2] = initdata.coord2s.at(kk);
  case 2: // 3-D or 2-D model
    coords312[1] = initdata.coord1s.at(kk);
  case 1: // 3-D, 2-D or 1-D model
    coords312[0] = initdata.coord3s.at(kk);
  }

  return coords312;
}

SuperdropAttrs GenSuperdrop::attrs_at(const unsigned int kk) const
/* helper function to return a superdroplet's attributes
at position kk in the initial conditions data. All
superdroplets created with same solute properties */
{
  const double radius(initdata.radii.at(kk));
  const double msol(initdata.msols.at(kk));
  const unsigned long long xi(initdata.xis.at(kk));
  const SoluteProperties solute(initdata.solutes.at(0));

  return SuperdropAttrs(solute, xi, radius, msol);
}

Superdrop GenSuperdrop::operator()(const unsigned int kk) const
{
  const unsigned int sdgbxindex(initdata.sdgbxindexes.at(kk));
  const std::array<double, 3> coords312(coords_at(kk));
  const SuperdropAttrs attrs(attrs_at(kk));
  const auto sdId(sdIdGen->next());

  return Superdrop(sdgbxindex, coords312[0],
                   coords312[1], coords312[2],
                   attrs, sdId);
}