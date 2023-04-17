// Author: Clara Bayley
// File: run_sdmstep.hpp
/* Header file for class to run
1 timestep of (uncoupled) SDM */

#ifndef RUN_SDMSTEP_HPP
#define RUN_SDMSTEP_HPP

#include <concepts>
#include <random>
#include <vector>
#include <algorithm>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./movesuperdropsindomain.hpp"
#include "./superdropwithgbxindex.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmotion.hpp"

template <SdMotion M>
class RunSDMStep
{
private:
  int coupl_or_motion(const int t_sdm, const int couplstep,
                      const MoveSuperdropsInDomain<M> &sdmmotion)
  /* given current timestep, t_sdm, work out which event
  (motion or coupling) is next to occur and return the time
  of the sooner event */
  {
    const int next_coupl = ((t_sdm / couplstep) + 1) * couplstep; // t of next output
    const int next_motion = sdmmotion.next_step(t_sdm);           // t of next sdmmotion

    return std::min(next_coupl, next_motion);
  }

public:
  RunSDMStep() {}

  void run_sdmstep(const int t_mdl, const int couplstep,
                   const Maps4GridBoxes &gbxmaps,
                   const MoveSuperdropsInDomain<M> &sdmmotion,
                   const SdmProcess auto &sdmprocess,
                   std::mt19937 &gen,
                   std::vector<GridBox> &gridboxes,
                   std::vector<SuperdropWithGbxindex> &SDsInGBxs)
  /* run SDM for each gridbox from time t_mdl to t_mdl+couplstep
  with subtimestepping such that each coupled timestep (couplstep)
  can be subdivided to allow the movement of superdroplets between
  gridboxes and the SDM process to occur at smaller time intervals */
  {

    int t_sdm(t_mdl);

    /* sdm model time is incremented until >= t_mdl+couplstep
    allowing for motion and process subtimestepping*/
    while (t_sdm < t_mdl + couplstep)
    {
      int nextt = coupl_or_motion(t_sdm, couplstep, sdmmotion);

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

#endif // RUN_SDMSTEP_HPP