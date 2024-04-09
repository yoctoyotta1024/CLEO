# %% [markdown]
# # Zarr files to netCDF
#
# With this script, the zarr files from CLEOS output can be transformed into a single netcdf file which has the dimensions
# - time
# - sd_id (super droplet id)

# %%
import sys
import os
import numpy as np
from pathlib import Path
import yaml
import xarray as xr

# from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *
from pySD.sdmout_src import sdtracing
from pySD.sdmout_src.supersdata import SupersData
from pySD.initsuperdropsbinary_src import *


from sdm_eurec4a.visulization import set_custom_rcParams

set_custom_rcParams()

#%%
path2CLEO = Path("/home/m/m301096/CLEO")
path2sdm_eurec4a = Path("/home/m/m301096/repositories/sdm-eurec4a")
rawdirectory = path2CLEO / "data/output/raw/rain/"
processeddirectory = path2CLEO / "data/output/processed/rain/"
yamldirectory = path2sdm_eurec4a / "data/model/input/all_rain_clusters_shallow"

# use paths to files
path2build = path2CLEO / "build"
configfile = path2CLEO / "eurec4a/experiment_03/src/config/rain1d_config.txt"
# yaml_config_file = path2sdm_eurec4a / "data/model/input/example_input_18.yaml"


# %%
def main(yaml_config_file):
    path2dataoutput = path2CLEO / "data/output"
    with open(yaml_config_file, 'r') as f:
        config_yaml = yaml.safe_load(f)

    identification_type = config_yaml['cloud']["identification_type"]
    cloud_id = config_yaml['cloud']["cloud_id"]

    OUTPUT_DIR = processeddirectory / f"{identification_type}_{cloud_id}"
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    OUTPUT_FILEPATH = OUTPUT_DIR / "full_dataset.nc"
    ### ----------------------- INPUT PARAMETERS ----------------------- ###
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constsfile    = path2CLEO / "libs/cleoconstants.hpp"
    # path and file names for plotting results
    setupfile     = rawdirectory / f"{identification_type}_{cloud_id}/rain1d_setup.txt"
    dataset       = rawdirectory / f"{identification_type}_{cloud_id}/rain1d_sol.zarr"


    # Get the config, constants and also an initial superdroplet dataset



    if os.path.exists(OUTPUT_FILEPATH):
        print(f"The file {OUTPUT_FILEPATH} exists. Skip processing.")
    else:
        print(f"The file {OUTPUT_FILEPATH} does not exist. Create it.")

        # read in constants and intial setup from setup .txt file
        config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
        consts = pysetuptxt.get_consts(setupfile, isprint=False)
        # Create a first simple dataset to have the coordinates for later netcdf creation
        sddata = pyzarr.get_supers(str(dataset), consts)
        lagrange = sddata.to_Dataset()
        lagrange.to_netcdf(OUTPUT_FILEPATH)


import glob
yaml_files = glob.glob(str(yamldirectory / "*.yaml"))
yaml_files.sort()
for yaml_file in yaml_files:
    try:
        main(yaml_file)
    except Exception as e:
        print(f"Error in {yaml_file}")
        print(e)

# %%
