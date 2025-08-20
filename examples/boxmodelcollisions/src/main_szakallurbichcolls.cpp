/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_szakallurbichcolls.cpp
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
 * breakup with flag decided based on section 2.2 of Szakáll and Urbich 2018.
 * After make/compiling, execute for example via:
 * ./src/szakallurbichcolls ../src/config/config.yaml
 */

#include "./main_impl.hpp"
#include "superdrops/collisions/breakup.hpp"
#include "superdrops/collisions/breakup_nfrags.hpp"
#include "superdrops/collisions/coalbure.hpp"
#include "superdrops/collisions/coalbure_flag.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/longhydroprob.hpp"

struct SzakallUrbichCreateMicrophysics {
  MicrophysicalProcess auto operator()(const Config &config, const Timesteps &tsteps) const {
    const auto c = config.get_breakup();

    const PairProbability auto collprob = LongHydroProb();
    const NFragments auto nfrags = ConstNFrags(c.constnfrags.nfrags);
    const CoalBuReFlag auto coalbure_flag = SUCoalBuReFlag{};
    const MicrophysicalProcess auto colls =
        CoalBuRe(tsteps.get_collstep(), &step2realtime, collprob, nfrags, coalbure_flag);
    return colls;
  }
};

int main(int argc, char *argv[]) {
  return generic_microphysics_main(argc, argv, SzakallUrbichCreateMicrophysics{});
}
