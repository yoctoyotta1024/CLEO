/*
 * Copyright (c) 2023 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: yac_comms.hpp
 * Project: coupldyn_yac
 * Created Date: Tuesday 31st October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Wilton Loch (WL); Aparna Dev (AD)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * send and receive dynamics communication functions for SDM when coupled to the yac
 * dynamics solver
 */

#ifndef LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_
#define LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "coupldyn_yac/yac_cartesian_dynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/supersindomain.hpp"
#include "superdrops/state.hpp"
#include "superdrops/thermodynamic_equations.hpp"

namespace KCS = KokkosCleoSettings;

/* 1-way coupling from coupldyn to CLEO's gridboxes where
coupldyn sends information to h_gbxs but doesn't
receive any back. Struct obeys coupling comms concept */
struct YacComms {
 private:
  /*
   * calculate mixing ratios for "cloud" water and "rain" water
   * calculation could potentially be on device instead of host
   */
  void calculate_cloud_rain_masses(const viewd_constgbx d_gbxs, const SupersInDomain& allsupers,
                                   Kokkos::DualView<double*> mass_cloud,
                                   Kokkos::DualView<double*> mass_rain) const;

  /*
   * updates the state of a gridbox using information received from YacCartesianDynamics solver
   * for 1-way or 2-way coupling to from Dynamics to CLEO SDM
   */
  void update_gridbox_state(const YacCartesianDynamics& ffdyn, const size_t ii, Gridbox& gbx) const;

  /*
   * fills YacComms send buffers for a 2-way coupling to CLEO SDM with the change in
   * state of a gridbox compared to data in dynamics
   */
  void update_send_buffers(const YacCartesianDynamics& ffdyn, const size_t ii, const double temp,
                           const double qvap, const double qcloud, const double qrain) const;

 public:
  /*
   * functor for calculate mass of "cloud" and "rain" water in each gridbox
   * "cloud" water is defined as droplets with radius < 0.38 mm (equivalent to ICON 1-mom scheme)
   */
  struct CalculateCloudRainMassesFunctor {
    subviewd_constsupers domainsupers;
    viewd_constgbx d_gbxs;
    Kokkos::DualView<double*>::t_dev d_mass_cloud;
    Kokkos::DualView<double*>::t_dev d_mass_rain;

    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember& team_member) const {
      constexpr double cloud_radius = 0.38e-3 / dlc::R0;  // 0.38 mm (as in ICON 1-mom scheme)
      const auto ii = team_member.league_rank();
      const auto supers = d_gbxs(ii).supersingbx.readonly(domainsupers);
      const size_t nsupers(supers.extent(0));

      double mass_cloud = 0.0;
      double mass_rain = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team_member, nsupers),
          [=](const size_t kk, double& m_cloud, double& m_rain) {
            const auto& drop = supers(kk);
            if (drop.get_radius() < cloud_radius) {
              // "cloud" water is defined as droplets with radius < 0.38 mm (ICON 1-mom)
              m_cloud += drop.condensate_mass() * drop.get_xi();
            } else {
              // "rain" water is any water in droplets > "cloud" radius 0.38 mm
              m_rain += drop.condensate_mass() * drop.get_xi();
            }
          },
          mass_cloud, mass_rain);

      Kokkos::single(Kokkos::PerTeam(team_member), [=, this]() {
        d_mass_cloud(ii) = mass_cloud;
        d_mass_rain(ii) = mass_rain;
      });
    }
  };

  /* send information from Gridboxes' states
  to coupldyn is null for YacCartesianDynamics*/
  template <typename GbxMaps, typename CD = YacCartesianDynamics>
  void send_dynamics(const GbxMaps& gbxmaps, const dualview_constgbx gbxs,
                     const SupersInDomain& allsupers, const YacCartesianDynamics& ffdyn) const;

  /* update Gridboxes' states using information
  received from YacCartesianDynamics solver for
  1-way coupling to CLEO SDM */
  template <typename GbxMaps, typename CD = YacCartesianDynamics>
  void receive_dynamics(const GbxMaps& gbxmaps, const YacCartesianDynamics& ffdyn,
                        const viewh_gbx h_gbxs) const;
};

#endif  // LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_
