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


def create_null_sdm(pycleo, tsteps, gbxmaps):
    print("PYCLEO STATUS: creating Null Observer")
    obs = pycleo.NullObserver()

    print("PYCLEO STATUS: creating Null Microphysical Process")
    micro = pycleo.NullMicrophysicalProcess()  # no microphysics

    print("PYCLEO STATUS: creating Null Superdroplet Movement")
    motion = pycleo.NullMotion()
    transport = pycleo.CartesianTransportAcrossDomain()
    boundary_conditions = pycleo.NullBoundaryConditions()
    move = pycleo.CartesianNullMoveSupersInDomain(
        motion, transport, boundary_conditions
    )

    print("PYCLEO STATUS: creating Null SDM Methods")
    sdm = pycleo.CartesianNullSDMMethods(
        tsteps.get_couplstep(), gbxmaps, micro, move, obs
    )  # no microphysics

    return sdm


def create_sdm(pycleo, cleo_config, tsteps, is_sdm_null):
    print("PYCLEO STATUS: creating GridboxMaps")
    gbxmaps = pycleo.create_cartesian_maps(
        cleo_config.get_ngbxs(),
        cleo_config.get_nspacedims(),
        cleo_config.get_grid_filename(),
    )

    if is_sdm_null:
        sdm = create_null_sdm(pycleo, tsteps, gbxmaps)
        store, dataset = None, None
    else:
        print("PYCLEO STATUS: creating Observer")
        store = pycleo.FSStore(cleo_config.get_zarrbasedir())
        dataset = pycleo.SimpleDataset(store)
        obs = pycleo.pycreate_observer(cleo_config, tsteps, dataset, store)

        print("PYCLEO STATUS: creating Microphysical Process")
        micro = pycleo.pycreate_microphysical_process(
            cleo_config, tsteps
        )  # config gives microphysics

        print("PYCLEO STATUS: creating Superdroplet Movement")
        motion = pycleo.create_cartesian_predcorr_motion(
            cleo_config, tsteps.get_motionstep()
        )
        # motion = pycleo.create_cartesian_predcorr_motion(cleo_config, False)
        transport = pycleo.CartesianTransportAcrossDomain()
        boundary_conditions = pycleo.NullBoundaryConditions()
        move = pycleo.CartesianMoveSupersInDomain(
            motion, transport, boundary_conditions
        )

        print("PYCLEO STATUS: creating SDM Methods")
        sdm = pycleo.CartesianSDMMethods(
            tsteps.get_couplstep(), gbxmaps, micro, move, obs
        )  # microphysics determined by settings for microphysics given in config

    print(f"PYCLEO STATUS: SDM created with couplstep = {sdm.get_couplstep()}")
    return sdm, dataset, store


def prepare_to_timestep_sdm(pycleo, cleo_config, sdm):
    print("PYCLEO STATUS: creating superdroplets")
    initsupers = pycleo.InitSupersFromBinary(
        cleo_config.get_initsupersfrombinary(), sdm.gbxmaps
    )
    allsupers = pycleo.create_supers_from_binary(
        initsupers, sdm.gbxmaps.get_local_ngridboxes_hostcopy()
    )

    print("PYCLEO STATUS: creating gridboxes")
    initgbxs = pycleo.InitGbxsNull(sdm.gbxmaps.get_local_ngridboxes_hostcopy())
    gbxs = pycleo.create_gbxs_cartesian_null(sdm.gbxmaps, initgbxs, allsupers)

    print("PYCLEO STATUS: preparing sdm")
    sdm.prepare_to_timestep(gbxs, allsupers)

    print("PYCLEO STATUS: preparation complete")
    return sdm, gbxs, allsupers


class CleoSDM:
    def __init__(
        self,
        pycleo,  # CLEO python bindings module
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
        self.pycleo = pycleo

        tsteps = pycleo.pycreate_timesteps(cleo_config)
        assert (
            pycleo.realtime2step(timestep) == tsteps.get_couplstep()
        ), "timestep and SDM coupling not equal"

        self.t_sdm = pycleo.realtime2step(
            t_start
        )  # convert from seconds to model timesteps (!)

        self.coupldyn = pycleo.coupldyn_numpy.NumpyDynamics(
            tsteps.get_couplstep(),
            press,
            temp,
            qvap,
            qcond,
            wvel,
            uvel,
            vvel,
        )
        self.comms = pycleo.coupldyn_numpy.NumpyComms()

        self.sdm, self.dataset, self.store = create_sdm(
            pycleo, cleo_config, tsteps, is_sdm_null
        )
        self.sdm, self.gbxs, self.allsupers = prepare_to_timestep_sdm(
            pycleo, cleo_config, self.sdm
        )

    def do_step(self, timestep, thermo):
        timestep = self.pycleo.realtime2step(
            timestep
        )  # convert from seconds to model timesteps (!)
        t_mdl_next = self.t_sdm + timestep
        assert t_mdl_next == self.sdm.next_couplstep(
            self.t_sdm
        ), "SDM out of sync with coupling"

        print(f"PYCLEO STATUS: start t_sdm = {self.t_sdm} [model timesteps]")
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
        print(f"PYCLEO STATUS: end t_sdm = {self.t_sdm} [model timesteps]")
        return thermo
