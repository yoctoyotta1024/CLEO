"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: kokkostools.py
Project: kokkostools
Created Date: Thursday 21st August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Script generates input files, runs program with "spdtest" executable, and post-processes
kokkos tools kernel timer profiling outful to test performance of CLEO using Kokkos tools
for a particular buildtype
"""


# %%
### -------------------------------- IMPORTS ------------------------------- ###
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from ruamel.yaml import YAML

from kp_kernel_timer import KpKernelTimer

# %%
### --------------------------- PARSE ARGUMENTS ---------------------------- ###
parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "path2kokkostools",
    type=Path,
    help="Absolute path to kokkos tools installation libkp_[XXX]",
)
parser.add_argument(
    "src_config_filename",
    type=Path,
    help="Absolute path to source configuration YAML file",
)
parser.add_argument(
    "postproc_filename",
    type=Path,
    help="Absolute path for profiler data .txt file(s)",
)
parser.add_argument(
    "--nruns",
    type=int,
    default="2",
    help="Number of times to run program",
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
path2kokkostools = args.path2kokkostools
path2build = args.path2build
src_config_filename = args.src_config_filename
postproc_filename = args.postproc_filename
nruns = args.nruns

### --- additional/derived arguments --- ###
tmppath = path2build / "tmp"
sharepath = path2build / "share"
binpath = path2build / "bin"
savefigpath = binpath

run_configs = {}  # run number: [config_filename, config_params]
thermofiles = sharepath / "kokkostools_dimlessthermo.dat"

for run in range(nruns):
    cf = path2build / "tmp" / f"kokkostools_config_run{run}.yaml"
    cp = {
        "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
        "grid_filename": str(sharepath / "kokkostools_dimlessGBxboundaries.dat"),
        "initsupers_filename": str(
            sharepath / f"kokkostools_dimlessSDsinit_run{run}.dat"
        ),
        "setup_filename": str(binpath / f"kokkostools_setup_run{run}.txt"),
        "zarrbasedir": str(binpath / f"kokkostools_sol_run{run}.zarr"),
    }
    run_configs[run] = [cf, cp]

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
    thermofiles,
    postproc_filename,
    isfigures,
):
    from cleopy import editconfigfile

    ### --- ensure build, share and bin directories exist --- ###
    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")
    path2build.mkdir(exist_ok=True)
    tmppath.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    postproc_filename.parent.mkdir(exist_ok=True)
    if savefigpath is not None:
        savefigpath.mkdir(exist_ok=True)

    ### --- add names of thermofiles to config_params --- ###
    for var in ["press", "temp", "qvap", "qcond", "wvel", "uvel", "vvel"]:
        config_params[var] = str(
            thermofiles.parent / Path(f"{thermofiles.stem}_{var}{thermofiles.suffix}")
        )

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
    all_thermofiles = thermofiles.parent.glob(
        f"{thermofiles.stem}*{thermofiles.suffix}"
    )
    for file in all_thermofiles:
        file.unlink(missing_ok=True)

    ### --- input binary files generation --- ###
    # equivalent to ``import kokkostools_inputfiles`` followed by
    # ``kokkostools_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "kokkostools" / "kokkostools_inputfiles.py"
    )
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
        thermofiles,
    ]
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def run_exectuable(path2kokkostools, path2build, config_filename, postproc_filename):
    ### --- delete any existing output dataset and setup files --- ###
    ### --- Note: profiler and post-processes data is not deleted --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["outputdata"]["setup_filename"]).unlink(missing_ok=True)
    shutil.rmtree(Path(config["outputdata"]["zarrbasedir"]), ignore_errors=True)
    all_postproc_filenames = postproc_filename.parent.glob(
        f"{postproc_filename.stem}*{postproc_filename.suffix}"
    )
    for file in all_postproc_filenames:
        file.unlink(missing_ok=True)

    ### --- run exectuable with given config file --- ###
    os.chdir(path2build / "bin")
    profiler = KpKernelTimer(path2kokkostools)
    executable = path2build / "examples" / "kokkostools" / "src" / "spdtest"
    cmd = [executable, config_filename]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)
    profiler.postprocess(Path.cwd(), postproc_filename)


# %%
### ----------------------------- RUN EXAMPLE ------------------------------ ###
for run, [config_filename, config_params] in run_configs.items():
    print(f"---------- RUN NUMBER: {run} ----------")
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
            thermofiles,
            postproc_filename,
            isfigures,
        )

    if args.do_run_executable:
        run_exectuable(path2kokkostools, path2build, config_filename, postproc_filename)

    if args.do_plot_results:
        print("\nno plotting script for kokkostools example")
    print("-----------------------------------")
