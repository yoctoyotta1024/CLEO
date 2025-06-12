"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: python_bindings.py
Project: python_bindings
Created Date: Thursday 5th June 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 10th June 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import argparse
import shutil
import sys
from mpi4py import MPI

from pathlib import Path
from ruamel.yaml import YAML

import python_bindings_inputfiles

parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for PySD)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "config_filename", type=Path, help="Absolute path to configuration YAML file"
)
parser.add_argument(
    "--do_inputfiles",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Generate initial condition binary files",
)
parser.add_argument(
    "--do_run_executable",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Run fromfile executable",
)
parser.add_argument(
    "--do_plot_results",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Plot results of fromfile example",
)
args = parser.parse_args()

path2CLEO = args.path2CLEO
path2build = args.path2build
src_config_filename = args.config_filename
isfigures = [False, False]

do_inputfiles = True
if args.do_inputfiles == "FALSE":
    do_inputfiles = False
do_run_executable = True
if args.do_run_executable == "FALSE":
    do_run_executable = False
do_plot_results = True
if args.do_plot_results == "FALSE":
    do_plot_results = False

sys.path.append(str(path2build / "pycleo"))
import pycleo

yaml = YAML()
with open(src_config_filename, "r") as file:
    config = yaml.load(file)
config_filename = (
    Path(config["python_initconds"]["paths"]["tmppath"]) / src_config_filename.name
)

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
if do_inputfiles:
    ### --- ensure build, share and bin directories exist --- ###
    pyconfig = config["python_initconds"]
    tmppath = Path(pyconfig["paths"]["tmppath"])
    sharepath = Path(pyconfig["paths"]["sharepath"])
    binpath = Path(pyconfig["paths"]["binpath"])
    savefigpath = Path(pyconfig["paths"]["savefigpath"])
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
    thermofiles = Path(pyconfig["thermo"]["thermofiles"])

    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")
    else:
        path2build.mkdir(exist_ok=True)
        tmppath.mkdir(exist_ok=True)
        sharepath.mkdir(exist_ok=True)
        binpath.mkdir(exist_ok=True)
        savefigpath.mkdir(exist_ok=True)
        shutil.copy(src_config_filename, config_filename)

        ### --- delete any existing initial conditions --- ###
        shutil.rmtree(grid_filename, ignore_errors=True)
        shutil.rmtree(initsupers_filename, ignore_errors=True)
        all_thermofiles = thermofiles.parent / Path(
            f"{thermofiles.stem}*{thermofiles.suffix}"
        )
        shutil.rmtree(all_thermofiles, ignore_errors=True)

        python_bindings_inputfiles.main(
            path2CLEO,
            path2build,
            config_filename,
            grid_filename,
            initsupers_filename,
            thermofiles,
            isfigures=isfigures,
        )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
def mpi_info(comm):
    print("\n--- PYCLEO STATUS: MPI INFORMATION ---")
    print(f"MPI version: {MPI.Get_version()}")
    print(f"Processor name: {MPI.Get_processor_name()}")
    print(f"Total processes: {comm.Get_size()}")
    print(f"Process rank: {comm.Get_rank()}")
    print("-----------------------")


def create_sdm(config, tsteps):
    print("PYCLEO STATUS: creating GridboxMaps")
    gbxmaps = pycleo.create_cartesian_maps(
        ngbxs=ngbxs, nspacedims=nspacedims, grid_filename=grid_filename
    )

    print("PYCLEO STATUS: creating Observer")
    obs = pycleo.NullObserver()

    print("PYCLEO STATUS: creating MicrophysicalProcess")
    micro = pycleo.NullMicrophysicalProcess()

    print("PYCLEO STATUS: creating Superdroplet Movement")
    motion = pycleo.NullMotion()
    transport = pycleo.CartesianTransportAcrossDomain()
    boundary_conditions = pycleo.NullBoundaryConditions()
    move = pycleo.CartesianNullMoveSupersInDomain(
        motion, transport, boundary_conditions
    )

    print("PYCLEO STATUS: creating SDM Methods")
    sdm = pycleo.CartesianNullSDMMethods(
        tsteps.get_couplstep(), gbxmaps, micro, move, obs
    )

    print(f"PYCLEO STATUS: SDM created with couplstep = {sdm.get_couplstep()}")
    return sdm


def prepare_to_timestep_sdm(config, sdm):
    print("PYCLEO STATUS: creating superdroplets")
    initsupers = pycleo.InitSupersFromBinary(
        config.get_initsupersfrombinary(), sdm.gbxmaps
    )
    allsupers = pycleo.create_supers_from_binary(
        initsupers, sdm.gbxmaps.get_local_ngridboxes_hostcopy()
    )

    print("PYCLEO STATUS: creating gridboxes")
    initgbxs = pycleo.InitGbxsNull(sdm.gbxmaps.get_local_ngridboxes_hostcopy())
    gbxs = pycleo.create_gbxs_cartesian_null(sdm.gbxmaps, initgbxs, allsupers)

    print("PYCLEO STATUS: preparing sdm")
    sdm.prepare_to_timestep(gbxs)

    print("PYCLEO STATUS: preparation complete")
    return sdm, gbxs, allsupers


def timestep_sdm(tsteps, sdm, gbxs, allsupers):
    t_mdl, t_end = 0, tsteps.get_t_end()

    print(f"PYCLEO STATUS: timestepping SDM from {t_mdl} to {t_end} [model timesteps]")
    while t_mdl <= t_end:
        print(f"PYCLEO STATUS: t = {t_mdl}")
        t_mdl_next = min(sdm.next_couplstep(t_mdl), sdm.obs.next_obs(t_mdl))

        # dynamics -> SDM: if (t_mdl % sdm.get_couplstep() == 0):
        # ``comms.receive_dynamics(sdm.gbxmaps, coupldyn, gbxs.view_host());``

        sdm.at_start_step(t_mdl, gbxs, allsupers)

        # here would be ``dynamics.run_step();``

        sdm.run_step(t_mdl, t_mdl_next, gbxs, allsupers)

        # SDM -> dynamics: if (t_mdl % sdm.get_couplstep() == 0):
        # ``comms.send_dynamics(sdm.gbxmaps, gbxs.view_host(), coupldyn);``

        t_mdl = t_mdl_next


def run_cleo_sdm(config):
    tsteps = pycleo.pycreate_timesteps(config)
    sdm = create_sdm(config, tsteps)
    sdm, gbxs, allsupers = prepare_to_timestep_sdm(config, sdm)
    timestep_sdm(tsteps, sdm, gbxs, allsupers)


def run_exec():
    mpi_info(MPI.COMM_WORLD)
    config = pycleo.Config(config_filename)
    pycleo.pycleo_initialize(config)
    run_cleo_sdm(config)
    pycleo.pycleo_finalize()


if do_run_executable:
    nspacedims = config["domain"]["nspacedims"]
    ngbxs = config["domain"]["ngbxs"]
    grid_filename = Path(config["inputfiles"]["grid_filename"])

    print(f"PYCLEO STATUS: i+j={pycleo.test_python_bindings(i=1, j=2)}")

    run_exec()


### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
if do_plot_results:
    print("no plotting script for python bindings example")
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
