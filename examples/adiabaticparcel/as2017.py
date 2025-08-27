"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: as2017.py
Project: adiabaticparcel
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generate input files, runs CLEO adia0d executable to create data and then
creates plots for adiabatic parcel example similar to Figure 5 of
"On the CCN (de)activation nonlinearities" S. Arabas and S. Shima 2017 to show
example of adaibatic parcel expansion and contraction.
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

isfigures = [False, True]  # booleans for [showing, saving] initialisation figures

### configuration parameters for yaml file for each of three runs of three
### different initial radii and number concentration runs (nine runs overall)
run_configs = {}  # run number: [config_filename, config_params]

params1 = {
    "W_avg": 1,
    "TAU_half": 150,
    "T_END": 300,
    "COUPLTSTEP": 1,
    "OBSTSTEP": 2,
}
params2 = {
    "W_avg": 0.5,
    "TAU_half": 300,
    "T_END": 600,
    "COUPLTSTEP": 1,
    "OBSTSTEP": 2,
}
params3 = {
    "W_avg": 0.002,
    "TAU_half": 75000,
    "T_END": 150000,
    "COUPLTSTEP": 3,
    "OBSTSTEP": 750,
}

runnum = 0
for icond in range(3):  # 3 diff initial conditions
    for params in [params1, params2, params3]:  # 3 different run parameters
        cf = path2build / "tmp" / f"as2017_config_run{runnum}.yaml"
        cp = {
            "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
            "grid_filename": str(sharepath / "as2017_dimlessGBxboundaries.dat"),
            "initsupers_filename": str(
                sharepath / f"as2017_dimlessSDsinit_icond{icond}.dat"
            ),
            "setup_filename": str(binpath / f"as2017_setup_run{runnum}.txt"),
            "zarrbasedir": str(binpath / f"as2017_sol_run{runnum}.zarr"),
        }
        for key, value in params.items():
            cp[key] = value

        run_configs[runnum] = [cf, cp]
        runnum += 1


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
    gen_gbxs,
    gen_supers,
    icond,
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

    ### --- copy src_config_filename into tmp and edit parameters --- ###
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

    ### --- input binary files generation --- ###
    # equivalent to ``import as2017_inputfiles`` followed by
    # ``as2017_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "adiabaticparcel" / "as2017_inputfiles.py"
    )
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
    ]
    if gen_gbxs:
        cmd.append("--gen_gbxs")
    if gen_supers:
        cmd.append("--gen_supers")
        cmd.append(f"--icond={icond}")
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
    executable = path2build / "examples" / "adiabaticparcel" / "src" / "adia0d"
    cmd = [executable, config_filename]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def plot_results(path2CLEO, savefigpath, config_filenames, runnums):
    plotting_script = path2CLEO / "examples" / "adiabaticparcel" / "as2017_plotting.py"
    python = sys.executable

    setupfiles, datasets = [], []
    for config_filename in config_filenames:
        yaml = YAML()
        with open(config_filename, "r") as file:
            config = yaml.load(file)
        setupfiles.append(Path(config["outputdata"]["setup_filename"]))
        datasets.append(Path(config["outputdata"]["zarrbasedir"]))
    grid_filename = Path(
        config["inputfiles"]["grid_filename"]
    )  # same grid for all datasets

    # equivalent to ``import as2017_plotting`` followed by
    # ``as2017_plotting.main(path2CLEO, savefigpath, ...)``
    cmd = [
        python,
        plotting_script,
        f"--path2CLEO={path2CLEO}",
        f"--savefigpath={savefigpath}",
        f"--grid_filename={grid_filename}",
    ]
    cmd += ["--setupfiles"] + setupfiles
    cmd += ["--datasets"] + datasets
    cmd += ["--runnums"] + [str(r) for r in runnums]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


# %%
### ----------------------------- RUN EXAMPLE ------------------------------ ###
if args.do_inputfiles:
    for runnum, [config_filename, config_params] in run_configs.items():
        if runnum == 0:
            gen_gbxs = True
        else:
            gen_gbxs = False

        if runnum % 3 == 0:  # 3 = number of different params dicts
            gen_supers, icond = (
                True,
                runnum // 3,
            )  # 3 = number of different params dicts
        else:
            gen_supers, icond = False, None

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
            gen_gbxs,
            gen_supers,
            icond,
            isfigures,
        )

if args.do_run_executable:
    cfs = list([cfgs[0] for cfgs in run_configs.values()])
    for config_filename in cfs:
        run_exectuable(path2build, config_filename)

if args.do_plot_results:
    runnums = list(run_configs.keys())
    cfs = list([cfgs[0] for cfgs in run_configs.values()])
    plot_results(path2CLEO, savefigpath, cfs, runnums)
