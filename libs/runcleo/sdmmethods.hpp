/*
 * ----- CLEO -----
 * File: sdmmethods.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct wrapping the core ingredients of the Super-droplet Model
 * (SDM) part of CLEO to enact on super-droplets and gridboxes
 */

#ifndef SDMMETHODS_HPP
#define SDMMETHODS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_Random.hpp>

#include "./kokkosaliases.hpp"
#include "./coupleddynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "observers/observers.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/motion.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/urbg.hpp"

template <CoupledDynamics CD, GridboxMaps GbxMaps,
          MicrophysicalProcess Microphys,
          Motion M, Observer Obs>
struct SDMMethods
{
// private:
public: // private except that GPU compatible Kokkos requries public access during calls for thread parallelisation

  unsigned int couplstep;           // coupled timestep
  
  KOKKOS_INLINE_FUNCTION
  unsigned int next_sdmstep(const unsigned int t_sdm,
                            const unsigned int next_mdl) const
  /* given current timestep, t_sdm, work out which event
  (motion or one complete step) is next to occur and return
  the time of the sooner event, (ie. next t_move or t_mdl) */
  {
    const unsigned int next_move(movesupers.next_step(t_sdm));        // t of next sdm movement
    const unsigned int t_next(!(next_mdl <
                                next_move)
                                  ? next_move
                                  : next_mdl);

    return t_next; // return smaller of two unsigned ints (see std::min)
  }

  KOKKOS_INLINE_FUNCTION
  void superdrops_movement(const unsigned int t_sdm,
                           const viewd_gbx d_gbxs,
                           const viewd_supers supers) const
  /* move superdroplets (including movement between
  gridboxes) according to movesupers struct */
  {
    movesupers.run_step(t_sdm, gbxmaps, d_gbxs, supers);
  }

  void sdm_microphysics(const unsigned int t_sdm,
                        const unsigned int t_next,
                        const viewd_gbx d_gbxs,
                        Kokkos::Random_XorShift64_Pool<ExecSpace> genpool) const
  /* enact SDM microphysics for each gridbox
  (using sub-timestepping routine) */
  {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "sdm_microphysics",
        Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
        KOKKOS_CLASS_LAMBDA(const size_t ii) 
    {
      // URBG<ExecSpace> urbg(genpool.get_state()); // thread safe random number generator
      // URBG<ExecSpace> urbg; // thread safe random number generator

      auto gbx = d_gbxs(ii);
      for (unsigned int subt = t_sdm; subt < t_next;
           subt = microphys.next_step(subt))
      {
        // microphys.run_step(subt, gbx.supersingbx(), gbx.state, urbg); // TODO use returned supers here
      }

      // genpool.free_state(urbg.gen);
    });
  }

public:
  GbxMaps gbxmaps;                  // maps from gridbox indexes to domain coordinates
  Microphys microphys;              // microphysical process
  MoveSupersInDomain<M> movesupers; // super-droplets' motion in domain
  Obs obs;                          // observer

  SDMMethods(const CD &coupldyn,
             const GbxMaps gbxmaps,
             const Microphys microphys,
             const MoveSupersInDomain<M> movesupers,
             const Obs obs)
      : couplstep(coupldyn.get_couplstep()),
        gbxmaps(gbxmaps),
        microphys(microphys),
        movesupers(movesupers),
        obs(obs) {}

  KOKKOS_INLINE_FUNCTION
  auto get_couplstep() const { return couplstep; }

  void prepare_to_timestep(const viewh_constgbx h_gbxs) const
  /* prepare CLEO SDM for timestepping */
  {
    obs.before_timestepping(h_gbxs);
  }

  void receive_dynamics(const CD &coupldyn,
                        const viewh_gbx h_gbxs) const {}
  /* update Gridboxes' states using information
  received from coupldyn */

  void send_dynamics(const CD &coupldyn,
                     const viewh_constgbx h_gbxs) const {}
  /* send information from Gridboxes' states to coupldyn */

  void run_step(const unsigned int t_mdl,
                const unsigned int t_mdl_next,
                const viewd_gbx d_gbxs,
                const viewd_supers supers,
                Kokkos::Random_XorShift64_Pool<ExecSpace> genpool) const
  /* run CLEO SDM (on device) from time t_mdl to t_mdl_next
  with sub-timestepping routine for super-droplets'
  movement and microphysics */
  {
    unsigned int t_sdm(t_mdl);
    while (t_sdm < t_mdl_next)
    {
      unsigned int t_sdm_next(next_sdmstep(t_sdm, t_mdl_next));

      superdrops_movement(t_sdm, d_gbxs, supers);
      sdm_microphysics(t_sdm, t_sdm_next, d_gbxs, genpool);

      t_sdm = t_sdm_next;
    }
  }
};

#endif // SDMMETHODS_HPP


// /home/m/m300950/testCLEOfire/libs/observers/printobs.hpp(61): warning #20011-D: calling a __host__ function("std::function<double (int)> ::function(const ::std::function<double (int)> &)") from a __host__ __device__ function("PrintObserver::PrintObserver") is not allowed

// Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

// /home/m/m300950/testCLEOfire/libs/observers/./observers.hpp(173): warning #20011-D: calling a __host__ function("std::function<double (int)> ::function(const ::std::function<double (int)> &)") from a __host__ __device__ function("DoTimeObs::DoTimeObs") is not allowed

// /home/m/m300950/testCLEOfire/libs/runcleo/./sdmmethods.hpp(91): warning #20011-D: calling a __host__ function("URBG< ::Kokkos::Cuda> ::URBG( ::Kokkos::Random_XorShift64< ::Kokkos::Cuda> )") from a __host__ __device__ function("SDMMethods<    ::FromFileDynamics,     ::CartesianMaps,     ::ConstTstepMicrophysics<    ::DoCondensation> ,     ::PredCorrMotion,     ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver< ::CombinedObserver<    ::PrintObserver,     ::ConstTstepObserver<    ::DoTimeObs> > ,     ::GbxindexObserver> ,  ::ConstTstepObserver<    ::DoNsupersObs> > ,  ::ConstTstepObserver<    ::DoNrainsupersObs> > ,  ::ConstTstepObserver<    ::DoTotNsupersObs> > ,  ::ConstTstepObserver<    ::DoMassMomentsObs> > ,  ::ConstTstepObserver<    ::DoRainMassMomentsObs> > ,  ::ConstTstepObserver<    ::DoStateObs> > ,  ::ConstTstepObserver<    ::DoSupersAttrsObs<    ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers< ::CombinedSuperdropsBuffers<    ::SdIdBuffer,     ::XiBuffer> ,     ::MsolBuffer> ,     ::RadiusBuffer> ,     ::Coord3Buffer> ,     ::Coord1Buffer> ,     ::Coord2Buffer> ,     ::SdgbxindexBuffer> > > > > ::sdm_microphysics(unsigned int, unsigned int,  ::Kokkos::View<    ::Gridbox *, void, void, void > ,  ::Kokkos::Random_XorShift64_Pool< ::Kokkos::Cuda> ) const::[lambda(unsigned long) (instance 1)]::operator () const") is not allowed

// /home/m/m300950/testCLEOfire/libs/runcleo/runcleo.hpp(132): warning #20011-D: calling a __host__ function("std::_Function_base::~_Function_base()") from a __host__ __device__ function("std::_Function_base::~_Function_base [subobject]") is not allowed
