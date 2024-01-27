/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: runcleo.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 27th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Generic functions for timestepping CLEO coupled model
 * (CLEO SDM coupled one-way/two-ways to a Dynamics Solver)
 */

#ifndef LIBS_RUNCLEO_RUNCLEO_HPP_
#define LIBS_RUNCLEO_RUNCLEO_HPP_

#include <concepts>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Random.hpp>

#include "../kokkosaliases.hpp"
#include "./coupleddynamics.hpp"
#include "./couplingcomms.hpp"
#include "./creategbxs.hpp"
#include "./createsupers.hpp"
#include "./initialconditions.hpp"
#include "./sdmmethods.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

template <CoupledDynamics CD, GridboxMaps GbxMaps, MicrophysicalProcess Microphys,
          Motion<GbxMaps> M, Observer Obs, CouplingComms<CD> Comms>
class RunCLEO {
 private:
  const SDMMethods<GbxMaps, Microphys, M, Obs> &sdm;
  CD &coupldyn;
  const Comms &comms;

  /* prepare CLEO SDM and Coupled Dyanmics for timestepping */
  int prepare_to_timestep(const dualview_constgbx gbxs) const {
    std::cout << "\n--- prepare timestepping ---\n";

    coupldyn.prepare_to_timestep();
    sdm.prepare_to_timestep(gbxs.view_host());

    std::cout << "--- prepare timestepping: success ---\n";
    return 0;
  }

  /* check  coupling of CLEO SDM and Coupled Dyanmics is correct.
  For example ensuring they have the same timestep for coupling */
  void check_coupling() const {
    if (sdm.get_couplstep() != coupldyn.get_couplstep()) {
      const std::string err(
          "coupling timestep of dynamics "
          "solver and CLEO SDM are not equal");
      throw std::invalid_argument(err);
    }
  }

  /* timestep CLEO from t=0 to t=t_end */
  int timestep_cleo(const unsigned int t_end, const dualview_gbx gbxs, const viewd_supers totsupers,
                    GenRandomPool genpool) const {
    std::cout << "\n--- timestepping ---\n";

    unsigned int t_mdl(0);
    while (t_mdl <= t_end) {
      /* start step (in general involves coupling) */
      const auto t_next = (unsigned int)start_step(t_mdl, gbxs);

      /* advance dynamics solver (optionally concurrent to SDM) */
      coupldyn_step(t_mdl, t_next);

      /* advance SDM (optionally concurrent to dynamics solver) */
      sdm_step(t_mdl, t_next, gbxs, totsupers, genpool);

      /* proceed to next step (in general involves coupling) */
      t_mdl = proceed_to_next_step(t_mdl, t_next, gbxs);
    }

    std::cout << "--- timestepping: success ---\n";
    return 0;
  }

  /** Start of every timestep: 1) communication of thermodynamic
  state from dynamics solver to CLEO's Gridboxes. 2) call
  sdm at_start_step function (e.g. to make observations).
  3) Return size of step to take given current timestep, t_mdl */
  unsigned int start_step(const unsigned int t_mdl, dualview_gbx gbxs) const {
    if (t_mdl % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.receive_dynamics(coupldyn, gbxs.view_host());
      gbxs.modify_host();
    }

    gbxs.sync_host();
    sdm.at_start_step(t_mdl, gbxs.view_host());

    return get_next_step(t_mdl);
  }

  /**
  * returns size of next step to take given current
  * timestep, t_mdl, such that next timestep is
  * sooner out of next timestep for obs or coupl
  */
  unsigned int get_next_step(const unsigned int t_mdl) const {
    const auto next_couplstep = [&, t_mdl]() {
      const auto interval = (unsigned int)sdm.get_couplstep();
      return ((t_mdl / interval) + 1) * interval;
    };

    /* t_next is sooner out of time for next coupl or obs */
    const auto next_coupl = (unsigned int)next_couplstep();
    const auto next_obs = (unsigned int)sdm.obs.next_obs(t_mdl);
    const auto t_next(!(next_coupl < next_obs)
                          ? next_obs
                          : next_coupl);  // return smaller of two unsigned ints (see std::min)

    return t_next;  // stepsize = t_next - t_mdl
  }

  /* run CLEO SDM (on host and device) */
  void sdm_step(const unsigned int t_mdl, unsigned int t_next, dualview_gbx gbxs,
                const viewd_supers totsupers, GenRandomPool genpool) const {
    gbxs.sync_device();  // get device up to date with host
    sdm.run_step(t_mdl, t_next, gbxs.view_device(), totsupers, genpool);
    gbxs.modify_device();  // mark device view of gbxs as modified
  }

  /* run coupled dynamics solver (on host) */
  void coupldyn_step(const unsigned int t_mdl, const unsigned int t_next) const {
    coupldyn.run_step(t_mdl, t_next);
  }

  /* returns incremented timestep 't_mdl' to 't_next'.
  Also this is where communication from CLEO SDM Gridbox
  States to coupled dynamics solver may occur */
  unsigned int proceed_to_next_step(unsigned int t_mdl, unsigned int t_next,
                                    dualview_gbx gbxs) const {
    if (t_mdl % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.send_dynamics(gbxs.view_host(), coupldyn);
    }

    return t_next;
  }

 public:
  RunCLEO(const SDMMethods<GbxMaps, Microphys, M, Obs> &sdm, CD &coupldyn, const Comms &comms)
      : sdm(sdm), coupldyn(coupldyn), comms(comms) {
    check_coupling();
  }

  ~RunCLEO() { sdm.obs.after_timestepping(); }

  /* create gridboxes and superdrops using initial
  conditions, then prepare and do timestepping. */
  int operator()(const InitialConditions auto &initconds, const unsigned int t_end) const {
    // create runtime objects
    viewd_supers totsupers(create_supers(initconds.initsupers));  // all the superdrops in domain
    dualview_gbx gbxs(create_gbxs(sdm.gbxmaps, initconds.initgbxs, totsupers));
    GenRandomPool genpool(std::random_device {}());

    // prepare CLEO for timestepping
    prepare_to_timestep(gbxs);

    // do timestepping from t=0 to t=t_end
    timestep_cleo(t_end, gbxs, totsupers, genpool);

    return 0;
  }
};

#endif  // LIBS_RUNCLEO_RUNCLEO_HPP_
