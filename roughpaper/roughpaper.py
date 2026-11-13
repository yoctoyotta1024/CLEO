"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- CLEO -----
File: roughpaper.py
Project: roughpaper
Created Date: Wednesday 7th October 2026
Author: Clara Bayley (CB)
Additional Contributors: Harshada Balasubramanian (HB)
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files and then runs the roughpaper CLEO executable
"cleocoupledsdm" (see roughpaper/src/main.cpp) with its configuration file. Use it as a
starting point for running your own setup of CLEO.
"""

# %%
### -------------------------------- IMPORTS ------------------------------- ###
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from ruamel.yaml import YAML

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

config_filename = path2build / "tmp" / "roughpaper_config.yaml"
thermofiles = sharepath / "dimlessthermo.dat"
config_params = {
    "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
    "grid_filename": str(sharepath / "dimlessGBxboundaries.dat"),
    "initsupers_filename": str(sharepath / "dimlessSDsinit.dat"),
    "setup_filename": str(binpath / "setup.txt"),
    "zarrbasedir": str(binpath / "SDMdata.zarr"),
}

isfigures = [False, True]  # booleans for [showing, saving] initialisation figures


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
    gen_config,
    gen_gbxs,
    gen_supers,
    gen_thermo,
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
    if savefigpath is not None:
        savefigpath.mkdir(exist_ok=True)

    ### --- add names of thermofiles to config_params --- ###
    for var in ["press", "temp", "qvap", "qcond", "wvel", "uvel", "vvel"]:
        config_params[var] = str(
            thermofiles.parent / Path(f"{thermofiles.stem}_{var}{thermofiles.suffix}")
        )

    ### --- copy src_config_filename into tmp and edit parameters --- ###
    if gen_config:
        config_filename.unlink(missing_ok=True)  # delete any existing config
        shutil.copy(src_config_filename, config_filename)
        editconfigfile.edit_config_params(config_filename, config_params)

    ### --- delete any existing initial conditions --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    if gen_gbxs:
        Path(config["inputfiles"]["grid_filename"]).unlink(missing_ok=True)
    if gen_supers:
        Path(config["initsupers"]["initsupers_filename"]).unlink(missing_ok=True)
    if gen_thermo:
        all_thermofiles = thermofiles.parent.glob(
            f"{thermofiles.stem}*{thermofiles.suffix}"
        )
        for file in all_thermofiles:
            file.unlink(missing_ok=True)

    ### --- input binary files generation --- ###
    # equivalent to ``import roughpaper_inputfiles`` followed by
    # ``roughpaper_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = path2CLEO / "roughpaper" / "roughpaper_inputfiles.py"
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
        thermofiles,
    ]
    if gen_gbxs:
        cmd.append("--gen_gbxs")
    if gen_supers:
        cmd.append("--gen_supers")
    if gen_thermo:
        cmd.append("--gen_thermo")
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def run_exectuable(path2build, config_filename):
    ### --- delete any existing output dataset and setup files --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["outputdata"]["setup_filename"]).unlink(missing_ok=True)
    shutil.rmtree(Path(config["outputdata"]["zarrbasedir"]), ignore_errors=True)

    ### --- run exectuable with given config file --- ###
    executable = path2build / "roughpaper" / "src" / "cleocoupledsdm"
    cmd = [executable, config_filename]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def plot_results():
    print("\nno plotting script for roughpaper")


# %%
### ----------------------------- RUN ROUGHPAPER --------------------------- ###
if args.do_inputfiles:
    gen_config = True
    gen_gbxs = True
    gen_supers = True
    gen_thermo = True
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
        gen_config,
        gen_gbxs,
        gen_supers,
        gen_thermo,
        isfigures,
    )

if args.do_run_executable:
    run_exectuable(path2build, config_filename)

if args.do_plot_results:
    plot_results()
