/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: main_golcolls.cpp
 * Project: src
 * Created Date: Thursday 12th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * runs the CLEO super-droplet model (SDM) for 0-D box model with Golovin's kernel.
 * After make/compiling, execute for example via:
 * ./src/golcolls ../src/config/config.yaml
 */

#include "./main_impl.hpp"
#include "superdrops/collisions/coalescence.hpp"
#include "superdrops/collisions/golovinprob.hpp"

struct GolovinCreateMicrophysics {
  MicrophysicalProcess auto operator()(const Config &config, const Timesteps &tsteps) const {
    const PairProbability auto prob = GolovinProb();
    const MicrophysicalProcess auto colls = CollCoal(tsteps.get_collstep(), &step2realtime, prob);
    return colls;
  }
};

int main(int argc, char *argv[]) {
  return generic_microphysics_main(argc, argv, GolovinCreateMicrophysics{});
}
