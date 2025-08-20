/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_lowlistcolls.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for 0-D box model with Low and List's kernel.
 * After make/compiling, execute for example via:
 * ./src/lowlistcolls ../src/config/config.yaml
 */

#include "./main_impl.hpp"
#include "superdrops/collisions/breakup.hpp"
#include "superdrops/collisions/breakup_nfrags.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/lowlistprob.hpp"

struct LowListCreateMicrophysics {
  MicrophysicalProcess auto operator()(const Config &config, const Timesteps &tsteps) const {
    const auto c = config.get_breakup();
    const PairProbability auto buprob = LowListBuProb();
    const NFragments auto nfrags = ConstNFrags(c.constnfrags.nfrags);
    const MicrophysicalProcess auto bu =
        CollBu(tsteps.get_collstep(), &step2realtime, buprob, nfrags);

    const PairProbability auto coalprob = LowListCoalProb();
    const MicrophysicalProcess auto coal =
        CollCoal(tsteps.get_collstep(), &step2realtime, coalprob);

    return coal >> bu;
  }
};

int main(int argc, char *argv[]) {
  return generic_microphysics_main(argc, argv, LowListCreateMicrophysics{});
}
