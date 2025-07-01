"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: python_bindings.py
Project: python_bindings
Created Date: Thursday 5th June 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 1st July 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import argparse
import numpy as np
import shutil
import sys
from mpi4py import MPI
from pathlib import Path
from ruamel.yaml import YAML

import python_bindings_inputfiles
from cleo_sdm import CleoSDM
from thermodynamics import Thermodynamics

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
from pycleo import coupldyn_numpy

yaml = YAML()
with open(src_config_filename, "r") as file:
    python_config = yaml.load(file)
config_filename = (
    Path(python_config["python_initconds"]["paths"]["tmppath"])
    / src_config_filename.name
)

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
if do_inputfiles:
    ### --- ensure build, share and bin directories exist --- ###
    pyinit = python_config["python_initconds"]
    tmppath = Path(pyinit["paths"]["tmppath"])
    sharepath = Path(pyinit["paths"]["sharepath"])
    binpath = Path(pyinit["paths"]["binpath"])
    savefigpath = Path(pyinit["paths"]["savefigpath"])
    grid_filename = Path(python_config["inputfiles"]["grid_filename"])
    initsupers_filename = Path(python_config["initsupers"]["initsupers_filename"])

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

        python_bindings_inputfiles.main(
            path2CLEO,
            path2build,
            config_filename,
            grid_filename,
            initsupers_filename,
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


def timestep_example(t_mdl, t_end, timestep, thermo, cleo_sdm):
    print(
        f"PYCLEO STATUS: timestepping SDM from {t_mdl}s to {t_end}s (timestep = {timestep}s)"
    )

    print("\n--- PYCLEO STATUS: THERMO INFORMATION ---")
    thermo.print_state()
    print("\n-----------------------------------------")

    while t_mdl <= t_end:
        print(f"PYCLEO STATUS: t = {t_mdl}s")

        thermo.temp[0] += 10  # mock example of changing dynamics

        thermo = cleo_sdm.do_step(timestep, thermo)
        print("temp[0:2]:", thermo.temp[0:2])
        print("qvap[0:2]:", thermo.massmix_ratios[0][0:2])

        t_mdl += timestep


def cleo_sdm_example(python_config, cleo_config):
    t_mdl, t_end = 0, python_config["timesteps"]["T_END"]  # [s]
    timestep = python_config["timesteps"]["COUPLTSTEP"]  # [s]

    thermo = create_thermodynamics(python_config)
    cleo_sdm = CleoSDM(
        pycleo,
        cleo_config,
        t_mdl,
        timestep,
        thermo.press,
        thermo.temp,
        thermo.massmix_ratios[0],
        thermo.massmix_ratios[1],
        is_sdm_null=False,
    )

    timestep_example(t_mdl, t_end, timestep, thermo, cleo_sdm)


def run_exec(python_config, config_filename):
    cleo_config = pycleo.Config(config_filename)
    pycleo.pycleo_initialize(cleo_config)
    cleo_sdm_example(python_config, cleo_config)


if do_run_executable:
    print(f"PYCLEO STATUS: 2+3={pycleo.test_pycleo(i=3, j=2)}")
    print(f"COUPLDYN_NUMPY STATUS: 2*3={coupldyn_numpy.test_coupldyn_numpy(i=3, j=2)}")

    mpi_info(MPI.COMM_WORLD)
    run_exec(python_config, config_filename)


### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
if do_plot_results:
    print("no plotting script for python bindings example")
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
