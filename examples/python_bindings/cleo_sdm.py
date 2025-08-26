"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: cleo_sdm.py
Project: python_bindings
Created Date: Thursday 12th June 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
class and functions for handing setup and running of CLEO via python bindings
"""


def create_null_sdm(cleo, tsteps, gbxmaps):
    print("CLEO STATUS: creating Null Observer")
    obs = cleo.NullObserver()

    print("CLEO STATUS: creating Null Microphysical Process")
    micro = cleo.NullMicrophysicalProcess()  # no microphysics

    print("CLEO STATUS: creating Null Superdroplet Movement")
    motion = cleo.NullMotion()
    transport = cleo.CartesianTransportAcrossDomain()
    boundary_conditions = cleo.NullBoundaryConditions()
    move = cleo.CartesianNullMoveSupersInDomain(motion, transport, boundary_conditions)

    print("CLEO STATUS: creating Null SDM Methods")
    sdm = cleo.CartesianNullSDMMethods(
        tsteps.get_couplstep(), gbxmaps, micro, move, obs
    )  # no microphysics

    return sdm


def create_sdm(cleo, cleo_config, tsteps, is_sdm_null):
    print("CLEO STATUS: creating GridboxMaps")
    gbxmaps = cleo.create_cartesian_maps(
        cleo_config.get_ngbxs(),
        cleo_config.get_nspacedims(),
        cleo_config.get_grid_filename(),
    )

    if is_sdm_null:
        sdm = create_null_sdm(cleo, tsteps, gbxmaps)
        store, dataset = None, None
    else:
        print("CLEO STATUS: creating Observer")
        store = cleo.FSStore(cleo_config.get_zarrbasedir())
        dataset = cleo.SimpleDataset(store)
        obs = cleo.pycreate_observer(cleo_config, tsteps, dataset, store)

        print("CLEO STATUS: creating Microphysical Process")
        micro = cleo.pycreate_microphysical_process(
            cleo_config, tsteps
        )  # config gives microphysics

        print("CLEO STATUS: creating Superdroplet Movement")
        motion = cleo.create_cartesian_predcorr_motion(
            cleo_config, tsteps.get_motionstep()
        )
        # motion = cleo.create_cartesian_predcorr_motion(cleo_config, False)
        transport = cleo.CartesianTransportAcrossDomain()
        boundary_conditions = cleo.NullBoundaryConditions()
        move = cleo.CartesianMoveSupersInDomain(motion, transport, boundary_conditions)

        print("CLEO STATUS: creating SDM Methods")
        sdm = cleo.CartesianSDMMethods(
            tsteps.get_couplstep(), gbxmaps, micro, move, obs
        )  # microphysics determined by settings for microphysics given in config

    print(f"CLEO STATUS: SDM created with couplstep = {sdm.get_couplstep()}")
    return sdm, dataset, store


def prepare_to_timestep_sdm(cleo, cleo_config, sdm):
    print("CLEO STATUS: creating superdroplets")
    initsupers = cleo.InitSupersFromBinary(
        cleo_config.get_initsupersfrombinary(), sdm.gbxmaps
    )
    allsupers = cleo.create_supers_from_binary(
        initsupers, sdm.gbxmaps.get_local_ngridboxes_hostcopy()
    )

    print("CLEO STATUS: creating gridboxes")
    initgbxs = cleo.InitGbxsNull(sdm.gbxmaps.get_local_ngridboxes_hostcopy())
    gbxs = cleo.create_gbxs_cartesian_null(sdm.gbxmaps, initgbxs, allsupers)

    print("CLEO STATUS: preparing sdm")
    sdm.prepare_to_timestep(gbxs, allsupers)

    print("CLEO STATUS: preparation complete")
    return sdm, gbxs, allsupers


class CleoSDM:
    def __init__(
        self,
        cleo,  # CLEO python bindings module
        cleo_config,
        t_start,
        timestep,
        press,
        temp,
        qvap,
        qcond,
        wvel,
        uvel,
        vvel,
        is_sdm_null=False,
    ):
        self.cleo = cleo

        tsteps = cleo.pycreate_timesteps(cleo_config)
        assert (
            cleo.realtime2step(timestep) == tsteps.get_couplstep()
        ), "timestep and SDM coupling not equal"

        self.t_sdm = cleo.realtime2step(
            t_start
        )  # convert from seconds to model timesteps (!)

        self.coupldyn = cleo.coupldyn_numpy.NumpyDynamics(
            tsteps.get_couplstep(),
            press,
            temp,
            qvap,
            qcond,
            wvel,
            uvel,
            vvel,
        )
        self.comms = cleo.coupldyn_numpy.NumpyComms()

        self.sdm, self.dataset, self.store = create_sdm(
            cleo, cleo_config, tsteps, is_sdm_null
        )
        self.sdm, self.gbxs, self.allsupers = prepare_to_timestep_sdm(
            cleo, cleo_config, self.sdm
        )

    def do_step(self, timestep, thermo):
        timestep = self.cleo.realtime2step(
            timestep
        )  # convert from seconds to model timesteps (!)
        t_mdl_next = self.t_sdm + timestep
        assert t_mdl_next == self.sdm.next_couplstep(
            self.t_sdm
        ), "SDM out of sync with coupling"

        print(f"CLEO STATUS: start t_sdm = {self.t_sdm} [model timesteps]")
        while self.t_sdm < t_mdl_next:
            t_sdm_next = min(
                self.sdm.next_couplstep(self.t_sdm), self.sdm.obs.next_obs(self.t_sdm)
            )

            if self.t_sdm % self.sdm.get_couplstep() == 0:
                self.comms.receive_dynamics(self.sdm.gbxmaps, self.coupldyn, self.gbxs)

            self.sdm.at_start_step(self.t_sdm, self.gbxs, self.allsupers)

            self.coupldyn.run_step(self.t_sdm, t_sdm_next)

            self.sdm.run_step(self.t_sdm, t_sdm_next, self.gbxs, self.allsupers)

            if self.t_sdm % self.sdm.get_couplstep() == 0:
                self.comms.send_dynamics(self.sdm.gbxmaps, self.gbxs, self.coupldyn)

            self.t_sdm = t_sdm_next
        print(f"CLEO STATUS: end t_sdm = {self.t_sdm} [model timesteps]")
        return thermo
