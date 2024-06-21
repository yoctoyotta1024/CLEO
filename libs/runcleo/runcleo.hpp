/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: runcleo.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TK)
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Generic templated class for timestepping CLEO SDM coupled (one-way/two-way)
 * a Dynamics Solver
 */

#ifndef LIBS_RUNCLEO_RUNCLEO_HPP_
#define LIBS_RUNCLEO_RUNCLEO_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Random.hpp>
#include <concepts>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "initialise/initialconditions.hpp"
#include "observers/observers.hpp"
#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/creategbxs.hpp"
#include "runcleo/createsupers.hpp"
#include "runcleo/sdmmethods.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @class RunCLEO
 * @brief Generic templated class for timestepping CLEO SDM coupled
 * (one-way/two-way) a Dynamics Solver
 *
 * This class orchestrates the timestepping of CLEO coupled model,
 * which consists of `SDM Methods` coupled one-way or two-way
 * to `Coupled Dynamics` with communication handled by `Coupling Comms`.
 *
 * @tparam CD Type of CoupledDynamics.
 * @tparam GbxMaps Type of GridboxMaps.
 * @tparam Microphys Type of MicrophysicalProcess.
 * @tparam M Type of Motion.
 * @tparam Obs Type of Observer.
 * @tparam Comms Type of CouplingComms.
 */
template <CoupledDynamics CD, GridboxMaps GbxMaps, MicrophysicalProcess Microphys,
          Motion<GbxMaps> M, typename BoundaryConditions, Observer Obs, CouplingComms<CD> Comms>
class RunCLEO {
 private:
  const SDMMethods<GbxMaps, Microphys, M, BoundaryConditions, Obs> &sdm; /**< SDMMethods object. */
  CD &coupldyn;       /**< CoupledDynamics object.  */
  const Comms &comms; /**< CouplingComms object. */

  /**
   * @brief Prepare SDM and Coupled Dynamics for timestepping.
   *
   * This function prepares the SDM and Coupled Dynamics for timestepping. It
   * calls the `prepare_to_timestep` function of both the Coupled Dynamics and
   * SDMMethods objects.
   *
   * @param gbxs DualView of gridboxes.
   * @return 0 on success.
   */
  int prepare_to_timestep(const dualview_constgbx gbxs) const {
    std::cout << "\n--- prepare timestepping ---\n";

    coupldyn.prepare_to_timestep();
    sdm.prepare_to_timestep(gbxs.view_device());

    std::cout << "--- prepare timestepping: success ---\n";
    return 0;
  }

  /**
   * @brief Check if coupling between SDM and Coupled Dynamics is correct.
   *
   * This function checks if the coupling timestep of the Dynamics Solver and
   * SDM are equal. Throws an exception if they are not.
   */
  void check_coupling() const {
    if (sdm.get_couplstep() != coupldyn.get_couplstep()) {
      const std::string err(
          "coupling timestep of dynamics "
          "solver and CLEO SDM are not equal");
      throw std::invalid_argument(err);
    }
  }

  /**
   * @brief Timestep CLEO from t=0 to t=t_end.
   *
   * This function performs the main timestepping loop for CLEO from the initial
   * time (t_mdl=0) to the specified end time (t_mdl=t_end). It calls RunCLEO's
   * `start_step`, `coupldyn_step`, `sdm_step`, and `proceed_to_next_step`
   * functions in a loop until timestepping is complete.
   *
   * @param t_end End time for timestepping.
   * @param gbxs DualView of gridboxes.
   * @param totsupers View of all superdroplets (both in and out of bounds of domain).
   * @return 0 on success.
   */
  int timestep_cleo(const unsigned int t_end, const dualview_gbx gbxs,
                    const viewd_supers totsupers) const {
    std::cout << "\n--- timestepping ---\n";

    unsigned int t_mdl(0);
    while (t_mdl <= t_end) {
      /* start step (in general involves coupling) */
      const auto t_next = (unsigned int)start_step(t_mdl, gbxs);

      /* advance dynamics solver (optionally concurrent to SDM) */
      coupldyn_step(t_mdl, t_next);

      /* advance SDM (optionally concurrent to dynamics solver) */
      sdm_step(t_mdl, t_next, gbxs, totsupers);

      /* proceed to next step (in general involves coupling) */
      t_mdl = proceed_to_next_step(t_next, gbxs);
    }

    std::cout << "--- timestepping: success ---\n";
    return 0;
  }

  /**
   * @brief Start of every timestep.
   *
   * This function is called at the start of every timestep. It includes
   * 1) communication of dynamics fields from the Dynamics Solver to the States
   * of CLEO's Gridboxes, 2) calling the `at_start_step` function of SDMMethods
   * (e.g. to make observations), and 3) returning the size of the timestep to
   * take now given the current timestep `t_mdl`.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   * @return Size of the next timestep.
   */
  unsigned int start_step(const unsigned int t_mdl, dualview_gbx gbxs) const {
    if (t_mdl % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.receive_dynamics(coupldyn, gbxs.view_host());
      gbxs.modify_host();
    }

    gbxs.sync_device();
    sdm.at_start_step(t_mdl, gbxs);

    return get_next_step(t_mdl);
  }

  /**
   * @brief Get the size of the next timestep.
   *
   * This function calculates and returns the next step size to take based on the
   * current model time, `t_mdl` and the coupling and obs times
   * obtained from the `sdm` object; `t_coupl` and `t_obs` respectively.
   *
   * @param t_mdl The current timestep of the model.
   * @return The size of the next timestep.
   *
   * @details
   * The size of the next timestep is determined by finding the smaller out of the
   * step to the next coupling time and the next observation time. The next coupling
   * time is calculated after receiving the size of the coupling timestep (a constant)
   * using the `sdm.get_couplstep()` function. The time of the next observation
   * is obtained from the `sdm.obs.next_obs()` function.
   *
   * The size of the next timestep is then calculated as `t_next - t_mdl`,
   * where `t_next` is the time closer to `t_mdl` out of `next_coupl`
   * and `next_obs`. The function uses Kokkos' version of C++ standard
   * library's `std::min` to find `t_next` (also GPU compatible).
   *
   * @see SDMMethods::get_couplstep()
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

  /**
   * @brief Run timestep of CLEO's Super-Droplet Model (SDM).
   *
   * This function runs SDM on both host and device from `t_mdl` to `t_next`.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param t_next Next timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   * @param totsupers View of all superdrops (both in and out of bounds of domain).
   */
  void sdm_step(const unsigned int t_mdl, unsigned int t_next, dualview_gbx gbxs,
                const viewd_supers totsupers) const {
    gbxs.sync_device();  // get device up to date with host
    sdm.run_step(t_mdl, t_next, gbxs.view_device(), totsupers);
    gbxs.modify_device();  // mark device view of gbxs as modified
  }

  /**
   * @brief Run timestep of Coupled Dynamics.
   *
   * This function runs the Coupled Dynamics on host from t_mdl to t_next.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param t_next Next timestep of the coupled model.
   */
  void coupldyn_step(const unsigned int t_mdl, const unsigned int t_next) const {
    coupldyn.run_step(t_mdl, t_next);
  }

  /**
   * @brief Proceed to the next timestep.
   *
   * This function returns the incremented timestep (`t_mdl` to `t_next`).
   * It is also where communication from the States of CLEO's
   * Gridboxes to the Coupled Dynamics may occur.
   *
   * @param t_next Next timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   * @return Incremented timestep.
   */
  unsigned int proceed_to_next_step(unsigned int t_next, dualview_gbx gbxs) const {
    if (t_next % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.send_dynamics(gbxs.view_host(), coupldyn);
    }

    return t_next;
  }

 public:
  /**
   * @brief Constructor for RunCLEO.
   *
   * Initializes the RunCLEO object with the provided SDMMethods, CoupledDynamics,
   * and CouplingComms objects. Checks if coupling between SDM and Dynamics is correct.
   *
   * @param sdm SDMMethods object.
   * @param coupldyn CoupledDynamics object.
   * @param comms CouplingComms object.
   */
  RunCLEO(const SDMMethods<GbxMaps, Microphys, M, BoundaryConditions, Obs> &sdm, CD &coupldyn,
          const Comms &comms)
      : sdm(sdm), coupldyn(coupldyn), comms(comms) {
    check_coupling();
  }

  /**
   * @brief Destructor for RunCLEO.
   *
   * Calls the after_timestepping function of the SDM observer.
   */
  ~RunCLEO() { sdm.obs.after_timestepping(); }

  /* create  then prepare and do timestepping. */
  /**
   * @brief Operator () for RunCLEO.
   *
   * Creates runtime objects, gridboxes, superdrops and random number generators
   * using initial conditions, then prepares and performs CLEO timestepping.
   *
   * @param initconds InitialConditions object containing initial conditions.
   * @param t_end End time for timestepping.
   * @return 0 on success.
   */
  int operator()(const InitialConditions auto &initconds, const unsigned int t_end) const {
    // create runtime objects
    viewd_supers totsupers(create_supers(initconds.initsupers));
    dualview_gbx gbxs(create_gbxs(sdm.gbxmaps, initconds.initgbxs, totsupers));

    // prepare CLEO for timestepping
    prepare_to_timestep(gbxs);

    // do timestepping from t=0 to t=t_end
    timestep_cleo(t_end, gbxs, totsupers);

    return 0;
  }
};

#endif  // LIBS_RUNCLEO_RUNCLEO_HPP_
