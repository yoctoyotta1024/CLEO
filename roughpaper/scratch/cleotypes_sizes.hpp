/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cleotypes_sizes.hpp
 * Project: scratch
 * Created Date: Friday 21st June 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#ifndef ROUGHPAPER_SCRATCH_CLEOTYPES_SIZES_HPP_
#define ROUGHPAPER_SCRATCH_CLEOTYPES_SIZES_HPP_

#include <Kokkos_Core.hpp>
#include <cstdint>
#include <iostream>

#include "gridboxes/gbxindex.hpp"
#include "gridboxes/gridbox.hpp"
#include "kokkosaliases.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

struct Gridbox2 {
  /* gbx with different ordering of members to see if padding is reduced */
  Gbxindex gbxindex;
  SupersInGbx supersingbx;
  State state;
};

struct SupersInGbx2 {
  /* supersingbx with different ordering of members to see if padding is reduced */
  unsigned int idx;
  kkpair_size_t refs;
  viewd_supers totsupers;
};

struct Superdrop2 {
  /* superdrop with different ordering of members to see if padding is reduced */
  SuperdropAttrs attrs;
  double coord3;
  double coord1;
  double coord2;
  unsigned int sdgbxindex;
  using IDType = IntID;
  // using IDType = EmptyID;
  [[no_unique_address]] IDType sdId;
};

struct SuperdropAttrs2 {
  SoluteProperties solute;
  uint64_t xi;
  double radius;
  double msol;
};

void print_type_sizes(int argc, char *argv[]) {
  Kokkos::initialize();
  {
    std::cout << "GBx: " << sizeof(Gridbox) << "\n";
    std::cout << "re-ordered GBx: " << sizeof(Gridbox2) << "\n";
    std::cout << "no padding: " << sizeof(State) + sizeof(SupersInGbx) + sizeof(Gbxindex) << "\n";
    std::cout << "  State: " << sizeof(State) << "\n";
    std::cout << "  SupersInGBx: " << sizeof(SupersInGbx) << "\n";
    std::cout << "  gbxindex: " << sizeof(Gbxindex) << "\n";

    std::cout << "\nSupersInGBx: " << sizeof(SupersInGbx) << "\n";
    std::cout << "re-ordered SupersInGBx: " << sizeof(SupersInGbx2) << "\n";
    std::cout << "no padding: "
              << sizeof(viewd_supers) + sizeof(kkpair_size_t) + sizeof(unsigned int) << "\n";
    std::cout << "  View: " << sizeof(viewd_supers) << "\n";
    std::cout << "  refs: " << sizeof(kkpair_size_t) << "\n";
    std::cout << "  idx: " << sizeof(unsigned int) << "\n";

    std::cout << "\nSD: " << sizeof(Superdrop) << "\n";
    std::cout << "re-ordered SD: " << sizeof(Superdrop2) << "\n";
    std::cout << "no padding: "
              << sizeof(unsigned int) + sizeof(double) * 3 + sizeof(SuperdropAttrs) + sizeof(IntID);
    std::cout << " or "
              << sizeof(unsigned int) + sizeof(double) * 3 + sizeof(SuperdropAttrs) +
                     sizeof(EmptyID)
              << "\n";
    std::cout << "  sdgbxindex: " << sizeof(unsigned int) << "\n";
    std::cout << "  coords: " << sizeof(double) * 3 << "\n";
    std::cout << "  attrs: " << sizeof(SuperdropAttrs) << "\n";
    std::cout << "  id: " << sizeof(IntID) << " or " << sizeof(EmptyID) << "\n";

    std::cout << "\nSDAttrs: " << sizeof(SuperdropAttrs) << "\n";
    std::cout << "re-ordered SDAttrs: " << sizeof(SuperdropAttrs2) << "\n";
    std::cout << "no padding: " << sizeof(uint64_t) + 2 * sizeof(double) + sizeof(SoluteProperties)
              << "\n";
    std::cout << "  xi: " << sizeof(uint64_t) << "\n";
    std::cout << "  radius: " << sizeof(double) << "\n";
    std::cout << "  msol: " << sizeof(double) << "\n";
    std::cout << "  SoluteProperties: " << sizeof(SoluteProperties) << "\n";
  }
  Kokkos::finalize();
}

#endif  // ROUGHPAPER_SCRATCH_CLEOTYPES_SIZES_HPP_
