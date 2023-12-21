/*
 * ----- CLEO -----
 * File: cartesianmotion_withreset.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 21st December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodic boundary
 * conditions and reset of superdroplets that leave
 * the domain through the lower domain boundary
 */

#ifndef CARTESIANMOTION_WITHRESET_HPP
#define CARTESIANMOTION_WITHRESET_HPP

#include <functional>
#include <random>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../kokkosaliases.hpp"
#include "./cartesianboundaryconds.hpp"
#include "./cartesianmaps.hpp"
#include "./cartesianmotion.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "superdrops/urbg.hpp"
#include "gridboxes/predcorrmotion.hpp"

struct ResetSuperdrop
{
  Kokkos::View<GenRandomPool[1]> genpool4reset;
  uint64_t nbins;
  Kokkos::pair<unsigned int, unsigned int> gbxidxs;

  ResetSuperdrop(const uint64_t i_nbins,
                 const unsigned int ngbxs,
                 const unsigned int ngbxs4reset)
      : genpool4reset("genpool4reset"),
        nbins(i_nbins),
        gbxidxs({ngbxs - ngbxs4reset, ngbxs})
  { 
    auto h_genpool4reset = Kokkos::create_mirror_view(genpool4reset);
    h_genpool4reset(0) = GenRandomPool(std::random_device{}());
    Kokkos::deep_copy(genpool4reset, h_genpool4reset);
  }

  KOKKOS_FUNCTION unsigned int
  reset_position(const CartesianMaps &gbxmaps,
                        URBG<ExecSpace> &urbg,
                        Superdrop &drop) const
  /* randomly update position of superdroplet by 
  randomly selecting a gbxindex from gbxidxs and then
  randomly selecting a coord3 with that gbx's bounds */
  {
    const auto sdgbxindex = urbg(gbxidxs.first,
                                 gbxidxs.second); // randomly selected gbxindex in range {incl., excl.} 
   
    const auto bounds = gbxmaps.coord3bounds(sdgbxindex);
    const auto coord = urbg.drand(bounds.first, bounds.second); // random coord within gbx bounds

    std::cout << "pos: "<< sdgbxindex << ", " << coord <<" \n";

    drop.set_sdgbxindex(sdgbxindex);
    drop.set_coord3(coord);

    return sdgbxindex;
  }

  KOKKOS_FUNCTION void
  reset_attributes(URBG<ExecSpace> &urbg,
                   Superdrop &drop) const
  /* reset radius and multiplicity of superdroplet
  by randomly sampling from binned distributions */
  {
    const auto bin = urbg(0, nbins);                    // index of randomly selected bin
    std::cout << "bin: "<< bin <<" \n";
  }

  KOKKOS_FUNCTION unsigned int
  operator()(const CartesianMaps &gbxmaps,
             Superdrop &drop) const
  {
    URBG<ExecSpace> urbg{genpool4reset(0).get_state()}; // thread safe random number generator
    
    const auto sdgbxindex = reset_position(gbxmaps, urbg, drop);
    reset_attributes(urbg, drop);

    genpool4reset(0).free_state(urbg.gen);

    return sdgbxindex;
  }
};

KOKKOS_FUNCTION unsigned int
change_if_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                const CartesianMaps &gbxmaps,
                                unsigned int idx,
                                Superdrop &drop);

struct CartesianChangeIfNghbrWithReset
/* wrapper of functions for use in PredCorrMotion's
ChangeToNghbr type for deciding if a superdroplet should move
to a neighbouring gbx in a cartesian domain and then updating the
superdroplet appropriately. Struct has three functions, one
for each direction (coord3 = z, coord1 = x, coord2 = y). For each,
the superdrop's coord is compared to gridbox bounds given by gbxmaps
for the current gbxindex 'idx'. If superdrop coord lies outside
bounds, forward or backward neighbour functions are called to
update sdgbxindex (and possibly other superdrop attributes).
Struct is same as CartesianChangeIfNghbr except for in 
coord3(...){...} function */
{
  ResetSuperdrop reset_superdrop;

  CartesianChangeIfNghbrWithReset(const uint64_t nbins,
                                  const unsigned int ngbxs,
                                  const unsigned int ngbxs4reset)
      : reset_superdrop(ResetSuperdrop(nbins, ngbxs, ngbxs4reset)) {}

  KOKKOS_INLINE_FUNCTION unsigned int
  coord3(const CartesianMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  {
    return change_if_coord3nghbr_withreset(reset_superdrop,
                                           gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int
  coord1(const CartesianMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  {
    return change_if_coord1nghbr(gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int
  coord2(const CartesianMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  {
    return change_if_coord2nghbr(gbxmaps, idx, drop);
  }
};

template <VelocityFormula TV>
inline PredCorrMotion<CartesianMaps, TV,
                      CartesianChangeIfNghbrWithReset,
                      CartesianCheckBounds>
CartesianMotionWithReset(const unsigned int motionstep,
                         const std::function<double(unsigned int)> int2time,
                         const TV terminalv,
                         const uint64_t nbins,
                         const unsigned int ngbxs,
                         const unsigned int ngbxs4reset)
/* returned type satisfies motion concept for motion of a
superdroplet using a predictor-corrector method to update
a superdroplet's coordinates and then updating it's
sdgbxindex as appropriate for a cartesian domain */
{
  const auto cin = CartesianChangeIfNghbrWithReset(nbins, ngbxs, ngbxs4reset);
  return PredCorrMotion<CartesianMaps, TV,
                        CartesianChangeIfNghbrWithReset,
                        CartesianCheckBounds>(motionstep,
                                              int2time,
                                              terminalv,
                                              cin,
                                              CartesianCheckBounds{});
}

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION unsigned int
change_to_backwards_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                          const unsigned int idx,
                                          const CartesianMaps &gbxmaps,
                                          Superdrop &superdrop);

KOKKOS_FUNCTION unsigned int
change_to_forwards_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                         const unsigned int idx,
                                         const CartesianMaps &gbxmaps,
                                         Superdrop &superdrop);

KOKKOS_FUNCTION unsigned int
change_if_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                const CartesianMaps &gbxmaps,
                                unsigned int idx,
                                Superdrop &drop)
/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord3 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
Note: backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
{
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx),
                                    drop.get_coord3()); // if value != 0 idx needs to change
  switch (flag)
  {
  case 1:
    idx = change_to_backwards_coord3nghbr_withreset(reset_superdrop,
                                                    idx, gbxmaps, drop);
    break;
  case 2:
    idx = change_to_forwards_coord3nghbr_withreset(reset_superdrop,
                                                   idx, gbxmaps, drop);
    break;
  }
  return idx;
}

KOKKOS_FUNCTION unsigned int
change_to_backwards_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                          const unsigned int idx,
                                          const CartesianMaps &gbxmaps,
                                          Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
if its coord3 has exceeded the z lower domain boundary */
{
  auto nghbr = (unsigned int)gbxmaps.coord3backward(idx);

  const auto incre = (unsigned int)1;                         // increment
  if (beyond_domainboundary(idx, incre, gbxmaps.get_ndim(0))) // drop was at lower z edge of domain (now moving below it)
  {
    nghbr = reset_superdrop(gbxmaps, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr; // gbxindex of z backwards (down) neighbour
};

KOKKOS_FUNCTION unsigned int
change_to_forwards_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                         const unsigned int idx,
                                         const CartesianMaps &gbxmaps,
                                         Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if superdrop has exceeded the z upper domain boundary */
{
  auto nghbr = (unsigned int)gbxmaps.coord3forward(idx);

  const auto incre = (unsigned int)1;                                 // increment
  if (beyond_domainboundary(idx + incre, incre, gbxmaps.get_ndim(0))) // drop was upper z edge of domain (now moving above it)
  {
    nghbr = reset_superdrop(gbxmaps, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr; // gbxindex of z forwards (up) neighbour
};

#endif // CARTESIANMOTION_WITHRESET_HPP

