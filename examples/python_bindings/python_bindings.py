"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: python_bindings.py
Project: python_bindings
Created Date: Thursday 5th June 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

# %%
### -------------------------------- IMPORTS ------------------------------- ###
import argparse
import numpy as np
import shutil
import subprocess
import sys
from mpi4py import MPI
from pathlib import Path
from ruamel.yaml import YAML

from cleo_sdm import CleoSDM
from thermodynamics import Thermodynamics

# %%
### --------------------------- PARSE ARGUMENTS ---------------------------- ###
parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "src_config_filename",
    type=Path,
    help="Absolute path to source configuration YAML file",
)
parser.add_argument(
    "--do_inputfiles",
    action="store_true",  # default is False
    help="Generate initial condition binary files",
)
parser.add_argument(
    "--do_run_executable",
    action="store_true",  # default is False
    help="Run executable",
)
parser.add_argument(
    "--do_plot_results",
    action="store_true",  # default is False
    help="Plot results of example",
)
args = parser.parse_args()

# %%
### -------------------------- INPUT PARAMETERS ---------------------------- ###
### --- command line parsed arguments --- ###
path2CLEO = args.path2CLEO
path2build = args.path2build
src_config_filename = args.src_config_filename

### --- additional/derived arguments --- ###
tmppath = path2build / "tmp"
sharepath = path2build / "share"
binpath = path2build / "bin"
savefigpath = binpath

config_filename = path2build / "tmp" / "pybind_config.yaml"
config_params = {
    "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
    "grid_filename": str(sharepath / "pybind_dimlessGBxboundaries.dat"),
    "initsupers_filename": str(sharepath / "pybind_dimlessSDsinit.dat"),
    "setup_filename": str(binpath / "pybind_setup.txt"),
    "zarrbasedir": str(binpath / "pybind_sol.zarr"),
}

isfigures = [False, False]  # booleans for [showing, saving] initialisation figures


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def inputfiles(
    path2CLEO,
    path2build,
    tmppath,
    sharepath,
    binpath,
    savefigpath,
    src_config_filename,
    config_filename,
    config_params,
    isfigures,
):
    from cleopy import editconfigfile

    ### --- ensure build, share and bin directories exist --- ###
    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")
    else:
        path2build.mkdir(exist_ok=True)
        tmppath.mkdir(exist_ok=True)
        sharepath.mkdir(exist_ok=True)
        binpath.mkdir(exist_ok=True)
        savefigpath.mkdir(exist_ok=True)

    ### --- copy src_config_filename into tmp and edit parameters --- ###
    config_filename.unlink(missing_ok=True)  # delete any existing config
    shutil.copy(src_config_filename, config_filename)
    editconfigfile.edit_config_params(config_filename, config_params)

    ### --- delete any existing initial conditions --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["inputfiles"]["grid_filename"]).unlink(missing_ok=True)
    Path(config["initsupers"]["initsupers_filename"]).unlink(missing_ok=True)

    ### --- input binary files generation --- ###
    # equivalent to ``import python_bindings_inputfiles`` followed by
    # ``python_bindings_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "python_bindings" / "python_bindings_inputfiles.py"
    )
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
    ]
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


# %%
### --------------- FUNCTION DEFINITIONS FOR RUN_EXECUTABLE ---------------- ###
def mpi_info(comm):
    print("\n--- CLEO STATUS: MPI INFORMATION ---")
    print(f"MPI version: {MPI.Get_version()}")
    print(f"Processor name: {MPI.Get_processor_name()}")
    print(f"Total processes: {comm.Get_size()}")
    print(f"Process rank: {comm.Get_rank()}")
    print("--------------------------------------")


def create_thermodynamics(python_config):
    ngbxs = python_config["domain"]["ngbxs"]

    temp = np.repeat(288.15, ngbxs)
    rho = np.repeat(1.225, ngbxs)
    press = np.repeat(101325, ngbxs)
    qvap = np.repeat(0.01, ngbxs)
    qcond = np.repeat(0.002, ngbxs)
    qice = np.repeat(0.003, ngbxs)
    qrain = np.repeat(0.004, ngbxs)
    qsnow = np.repeat(0.005, ngbxs)
    qgrau = np.repeat(0.006, ngbxs)

    return Thermodynamics(temp, rho, press, qvap, qcond, qice, qrain, qsnow, qgrau)


def create_winds(python_config):
    nfaces = (
        2 * python_config["domain"]["ngbxs"]
    )  # each gridbox has 2 faces in each direction

    wvel = np.repeat(1.0, nfaces)
    uvel = np.repeat(0.0, nfaces)
    vvel = np.repeat(0.0, nfaces)

    return wvel, uvel, vvel


def timestep_example(t_mdl, t_end, timestep, thermo, cleo_sdm):
    print(
        f"CLEO STATUS: timestepping SDM from {t_mdl}s to {t_end}s (timestep = {timestep}s)"
    )

    print("\n--- CLEO STATUS: THERMO INFORMATION ---")
    thermo.print_state()
    print("\n-----------------------------------------")

    while t_mdl <= t_end:
        print(f"CLEO STATUS: t = {t_mdl}s")

        thermo.temp[0] += 10  # mock example of changing dynamics

        thermo = cleo_sdm.do_step(timestep, thermo)
        print("temp[0:2]:", thermo.temp[0:2])
        print("qvap[0:2]:", thermo.massmix_ratios[0][0:2])

        t_mdl += timestep


def cleo_sdm_example(cleo, python_config, cleo_config):
    t_mdl, t_end = 0, python_config["timesteps"]["T_END"]  # [s]
    timestep = python_config["timesteps"]["COUPLTSTEP"]  # [s]

    thermo = create_thermodynamics(python_config)
    wvel, uvel, vvel = create_winds(python_config)
    cleo_sdm = CleoSDM(
        cleo,
        cleo_config,
        t_mdl,
        timestep,
        thermo.press,
        thermo.temp,
        thermo.massmix_ratios[0],
        thermo.massmix_ratios[1],
        wvel,
        uvel,
        vvel,
        is_sdm_null=False,
    )

    timestep_example(t_mdl, t_end, timestep, thermo, cleo_sdm)


def run_sdm_example(cleo, python_config, config_filename):
    cleo_config = cleo.Config(config_filename)
    cleo.cleo_initialize(cleo_config)
    cleo_sdm_example(cleo, python_config, cleo_config)


def run_exectuable(path2build, config_filename):
    sys.path.append(str(path2build / "cleo_python_bindings"))
    import cleo_python_bindings as cleo
    from cleo_python_bindings import coupldyn_numpy

    yaml = YAML()
    with open(config_filename, "r") as file:
        python_config = yaml.load(file)

    print(f"CLEO STATUS: 2+3={cleo.test_cleo_python_bindings(i=3, j=2)}")
    print(f"COUPLDYN_NUMPY STATUS: 2*3={coupldyn_numpy.test_coupldyn_numpy(i=3, j=2)}")

    mpi_info(MPI.COMM_WORLD)
    run_sdm_example(cleo, python_config, config_filename)


# %%
### ----------------------------- RUN EXAMPLE ------------------------------ ###
if args.do_inputfiles:
    inputfiles(
        path2CLEO,
        path2build,
        tmppath,
        sharepath,
        binpath,
        savefigpath,
        src_config_filename,
        config_filename,
        config_params,
        isfigures,
    )

if args.do_run_executable:
    run_exectuable(path2build, config_filename)

if args.do_plot_results:
    print("\nno plotting script for python bindings example")
