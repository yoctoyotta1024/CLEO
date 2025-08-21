/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_testikstraubcolls.cpp
 * Project: src
 * Created Date: Sunday 16th June 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for 0-D box model with coalescence, rebound and
 * breakup with flag decided based on section 4 of Testik et al. 2011 (figure 12) as well
 * as coalescence efficiency from Straub et al. 2010 and Schlottke et al. 2010.
 * After make/compiling, execute for example via:
 * ./src/testikstraubcolls ../src/config/config.yaml
 */

#include "./main_impl.hpp"
#include "superdrops/collisions/breakup.hpp"
#include "superdrops/collisions/breakup_nfrags.hpp"
#include "superdrops/collisions/coalbure.hpp"
#include "superdrops/collisions/coalbure_flag.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/longhydroprob.hpp"

struct TestikStraubCreateMicrophysics {
  MicrophysicalProcess auto operator()(const Config &config, const Timesteps &tsteps) const {
    const PairProbability auto collprob = LongHydroProb();
    const NFragments auto nfrags = CollisionKineticEnergyNFrags{};
    const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag{};
    const MicrophysicalProcess auto colls =
        CoalBuRe(tsteps.get_collstep(), &step2realtime, collprob, nfrags, coalbure_flag);
    return colls;
  }
};

int main(int argc, char *argv[]) {
  return generic_microphysics_main(argc, argv, TestikStraubCreateMicrophysics{});
}
