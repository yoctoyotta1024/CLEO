/*
 * ----- CLEO -----
 * File: runcleo.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Generic functions for timestepping CLEO coupled model
 * (CLEO SDM coupled one-way/two-ways to a Dynamics Solver)
 */

#ifndef RUNCLEO_HPP
#define RUNCLEO_HPP

#include <string>
#include <stdexcept>
#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "./coupleddynamics.hpp"
#include "./createsupers.hpp"
#include "./creategbxs.hpp"
#include "./runtimestats.hpp"
#include "./sdmmethods.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"

template <CoupledDynamics CD, GridboxMaps GbxMaps,
          MicrophysicalProcess Microphys,
          Motion M, Observer Obs>
class RunCLEO
{
private:
  const CD &coupldyn;
  const SDMMethods<CD, GbxMaps, Microphys, M, Obs> &sdm;

  int prepare_timestepping() const
  /* prepare CLEO SDM and Coupled Dyanmics for timestepping */
  {
    coupldyn.prepare_to_timestep();
    sdm.prepare_to_timestep();

    return 0;
  }

  void check_coupling() const
  /* check  coupling of CLEO SDM and Coupled Dyanmics is correct.
  For example ensuring they have the same timestep for coupling */
  {
    if (sdm.get_couplstep() != coupldyn.get_couplstep())
    {
      const std::string err("coupling timestep of dyanmics "
                            "solver and CLEO SDM are not equal");
      throw std::invalid_argument(err);
    }
  }

  int timestep_cleo(const unsigned int t_end,
                    RunStats &stats,
                    const dualview_gbx gbxs,
                    const viewd_supers supers) const
  /* timestep CLEO from t=0 to t=t_end */
  {
    unsigned int t_mdl(0);
    while (t_mdl <= t_end)
    {
      /* start step (in general involves coupling) */
      const unsigned int t_next(start_step(t_mdl, gbxs));

      /* advance dynamics solver (optionally concurrent to SDM) */
      coupldyn_step(t_mdl, t_next);
  
      /* advance SDM (optionally concurrent to dynamics solver) */
      sdm_step(t_mdl, t_next, gbxs, supers);

      /* proceed to next step (in general involves coupling) */
      t_mdl = proceed_to_next_step(t_mdl, t_next, gbxs);
    }

    return 0;
  }

  unsigned int start_step(const unsigned int t_mdl,
                          dualview_gbx gbxs) const
  /* Start of every timestep: 1) communication of thermodynamic state
  from dynamics solver to CLEO's Gridboxes. 2) Make observation.
  3) Return size of step to take given current timestep, t_mdl */
  {
    if (t_mdl % sdm.get_couplstep() == 0)
    {
      gbxs.sync_host();
      sdm.receive_dynamics(coupldyn, gbxs.view_host());
      gbxs.modify_host();
    }

    gbxs.sync_host();
    sdm.obs.at_start_step(t_mdl, gbxs.view_host());

    return get_next_step(t_mdl);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int get_next_step(const unsigned int t_mdl) const
  /* returns size of next step to take given current
  timestep, t_mdl, such that next timestep is
  sooner out of next timestep for obs or coupl */
  {
    const auto next_couplstep = [&, t_mdl]()
    {
      const unsigned int interval(sdm.get_couplstep());
      return ((t_mdl / interval) + 1) * interval;
    };

    /* t_next is sooner out of time for next coupl or obs */
    const unsigned int next_coupl(next_couplstep());
    const unsigned int next_obs(sdm.obs.next_obs(t_mdl));
    const auto t_next(!(next_coupl < next_obs) ? next_obs : next_coupl); // return smaller of two unsigned ints (see std::min)
    
    return t_next;                                                       // stepsize = t_next - t_mdl
  }

  void sdm_step(const unsigned int t_mdl,
                unsigned int t_next,
                dualview_gbx gbxs,
                const viewd_supers supers) const
  /* run CLEO SDM (on device) */
  {
    gbxs.sync_device();
    sdm.run_step(t_mdl, t_next, gbxs.view_device(), supers);
    gbxs.modify_device();
  }

  void coupldyn_step(const unsigned int t_mdl,
                     const unsigned int t_next) const
  /* run coupled dynamics solver (on host) */
  {
    coupldyn.run_step(t_mdl, t_next);
  }

  unsigned int proceed_to_next_step(unsigned int t_mdl,
                                    unsigned int t_next,
                                    dualview_gbx gbxs) const

  /* returns incremented timestep 't_mdl' to 't_next'.
  Also this is where communication from CLEO SDM Gridbox
  States to coupled dynamics solver may occur */
  {
    if (t_mdl % sdm.get_couplstep() == 0)
    {
      gbxs.sync_host();
      sdm.send_dynamics(coupldyn, gbxs.view_host());
    }

    return t_next;
  }

public:
  RunCLEO(const CD &coupldyn,
          const SDMMethods<CD, GbxMaps, Microphys, M, Obs> &sdm)
      : coupldyn(coupldyn), sdm(sdm)
  {
    check_coupling(); 
  }

  template <typename IC>
  int operator()(const IC &initconds,
                 const unsigned int t_end) const
  /* create gridboxes and superdrops using initial conditions,
  then prepare and do timestepping. Meanwhile there is the
  option to record some runtime statistics using "stats" */
  {
    // create runtime objects
    RunStats stats;
    viewd_supers supers(CreateSupers{}(initconds.initsupers));
    dualview_gbx gbxs(CreateGbxs{}(initconds.initgbxs, supers));

    // prepare CLEO for timestepping
    prepare_timestepping();
    stats.before_timestepping();

    // do timestepping from t=0 to t=t_end
    timestep_cleo(t_end, stats, gbxs, supers);
    stats.after_timestepping();
    return 0;
  }
};

#endif // RUNCLEO_HPP