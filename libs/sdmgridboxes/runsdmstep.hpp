// Author: Clara Bayley
// File: runsdmstep.hpp
/* Header file for class to
run 1 timestep of SDM */

#ifndef RUNSDMSTEP_HPP
#define RUNSDMSTEP_HPP

#include <concepts>
#include <random>
#include <vector>
#include <algorithm>

#include "./gridbox.hpp"
#include "./sdmotion.hpp"
#include "./maps4gridboxes.hpp"
#include "./movesuperdropsindomain.hpp"
#include "./superdropwithgbxindex.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "observers/observers.hpp"

template <SdMotion M, SdmProcess P, Observer O>
class RunSDMStep
{
private:
  int onestep_or_motion(const int t_sdm, const int onestep,
                      const MoveSuperdropsInDomain<M> &sdmmotion) const
  /* given current timestep, t_sdm, work out which event
  (motion or one complete step) is next to occur and return the
  time of the sooner event */
  {
    const int next_one = ((t_sdm / onestep) + 1) * onestep; // t of next output
    const int next_motion = sdmmotion.next_step(t_sdm);           // t of next sdmmotion

    return std::min(next_one, next_motion);
  }

public:
  const MoveSuperdropsInDomain<M> sdmmotion;
  const Maps4GridBoxes &gbxmaps;
  const P &sdmprocess;
  const O &observer;
  const size_t ngridboxes;

  RunSDMStep(const Maps4GridBoxes &gbxmaps, const M &movesd,
             const P &sdmprocess, const O &observer)
      : sdmmotion(movesd), gbxmaps(gbxmaps), sdmprocess(sdmprocess),
      observer(observer), ngridboxes(gbxmaps.gbxidxs.size()) {}

  void run_sdmstep(const int t_mdl, const int onestep,
                   std::mt19937 &gen,
                   std::vector<GridBox> &gridboxes,
                   std::vector<SuperdropWithGbxindex> &SDsInGBxs) const
  /* run SDM for each gridbox from time t_mdl to t_mdl+onestep
  with subtimestepping such that each step (onestep) can be subdivided
  to allow the movement of superdroplets between gridboxes and the
  SDM process to occur at smaller time intervals */
  {

    int t_sdm(t_mdl);

    /* sdm model time is incremented until >= t_mdl+onestep
    allowing for motion and process subtimestepping*/
    while (t_sdm < t_mdl + onestep)
    {
      int nextt = onestep_or_motion(t_sdm, onestep, sdmmotion);

      sdmmotion.run_step(t_sdm, gbxmaps, SDsInGBxs, gridboxes);

      /* run SDM process for each gridbox
      using sdmprocess subttimestepping routine */
      for (auto &gbx : gridboxes)
      {
        for (int subt = t_sdm; subt < nextt;
             subt = sdmprocess.next_step(subt))
        {
          sdmprocess.run_step(subt, gbx.span4SDsinGBx,
                              gbx.state, gen);
        }
      }

      t_sdm = nextt;
    }
  }
};

#endif // RUNSDMSTEP_HPP