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
 * Last Modified: Wednesday 26th June 2024
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
 * @tparam BoundaryConditions Type of boundary conditions for super-droplet movement.
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
   * @brief Prepare SDM and Coupled Dynamics for timestepping.
   *
   * This function prepares the SDM and Coupled Dynamics for timestepping. It
   * calls the `prepare_to_timestep` function of both the Coupled Dynamics and
   * SDMMethods objects.
   *
   * @param gbxs DualView of gridboxes.
   */
  void prepare_to_timestep(const dualview_constgbx gbxs) const {
    std::cout << "\n--- prepare timestepping ---\n";

    coupldyn.prepare_to_timestep();
    sdm.prepare_to_timestep(gbxs.view_device());

    std::cout << "--- prepare timestepping: success ---\n";
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
   * and `next_obs`. The function uses explicit implementation of C++ standard
   * library's `std::min` to find `t_next` to make it GPU compatible.
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
   * @brief Timestep CLEO from t=0 to t=t_end.
   *
   * This function performs the main timestepping loop for CLEO from the initial
   * time (t_mdl=0) to the specified end time (t_mdl=t_end). It calls RunCLEO's
   * `get_next_step`, `coupldyn_step`, `start_sdm_step`, `sdm_step`, and
   * `end_sdm_step` functions in a loop until timestepping is complete.
   *
   * @param t_end End time for timestepping.
   * @param gbxs DualView of gridboxes.
   * @param totsupers View of all superdroplets (both in and out of bounds of domain).
   */
  void timestep_cleo(const unsigned int t_end, const dualview_gbx gbxs,
                     const viewd_supers totsupers) const {
    std::cout << "\n--- timestepping ---\n";

    unsigned int t_mdl = 0;
    at_initial_conditions(t_mdl, gbxs);
    while (t_mdl < t_end) {
      const auto t_next = get_next_step(t_mdl);

      /* advance Dynamics Solver (optionally asynchronous to SDM) */
      coupldyn_step(t_mdl, t_next);

      /* start step (in general involves coupling: Dynamics -> SDM) */
      start_sdm_step(t_next, gbxs);

      /* advance SDM (optionally asynchronous to Dynamics Solver) */
      sdm_step(t_mdl, t_next, gbxs, totsupers);

      /* end SDM step (in general involves coupling SDM -> Dynamics) */
      end_sdm_step(t_next, gbxs);

      /* proceed to next step */
      t_mdl = t_next;
    }

    std::cout << "--- timestepping: success ---\n";
  }

  /**
   * @brief Execute at initial conditions before timestepping routine begins.
   *
   * This function is called only once at the start of the simulation before
   * any actions have occured in the timestepping routine. It calls the
   * `at_initial_conditions` function of the SDMMethods (e.g. to make an observation).
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   */
  void at_initial_conditions(const unsigned int t_mdl, dualview_gbx gbxs) const {
    sdm.at_initial_conditions(t_mdl, gbxs);
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
   * @brief Start of every SDM timestep.
   *
   * This function is called at the start of every timestep. It includes
   * communication of dynamics fields from the Dynamics Solver to the States
   * of CLEO's Gridboxes.
   *
   * @param t_next Next timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   */
  void start_sdm_step(const unsigned int t_next, dualview_gbx gbxs) const {
    if (t_next % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.receive_dynamics(coupldyn, gbxs.view_host());
      gbxs.modify_host();
    }
  }

  /**
   * @brief End of every SDM timestep.
   *
   * This function is called at the end of every timestep, i.e. when t_mdl = t_next.
   * It includes 1) calling the `at_end_step` function of SDMMethods (e.g. to make observations),
   * and 2) possible communication from the States of CLEO's Gridboxes to the Coupled Dynamics.
   *
   * @param t_mdl Current timestep of the coupled model.
   * @param gbxs DualView of gridboxes.
   */
  void end_sdm_step(const unsigned int t_mdl, dualview_gbx gbxs) const {
    if (t_mdl % sdm.get_couplstep() == 0) {
      gbxs.sync_host();
      comms.send_dynamics(gbxs.view_host(), coupldyn);
    }

    gbxs.sync_device();
    sdm.at_end_step(t_mdl, gbxs);
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
   */
  void operator()(const InitialConditions auto &initconds, const unsigned int t_end) const {
    // create runtime objects
    viewd_supers totsupers(create_supers(initconds.initsupers));
    dualview_gbx gbxs(create_gbxs(sdm.gbxmaps, initconds.initgbxs, totsupers));

    // prepare CLEO for timestepping
    prepare_to_timestep(gbxs);

    // do timestepping from t=0 to t=t_end
    timestep_cleo(t_end, gbxs, totsupers);
  }
};

#endif  // LIBS_RUNCLEO_RUNCLEO_HPP_
