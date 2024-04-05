# %% [markdown]
# # Zarr files to netCDF
# 
# With this script, the zarr files from CLEOS output can be transformed into a single netcdf file which has the dimensions
# - time
# - sd_id (super droplet id)

# %%
import sys
import numpy as np
from pathlib import Path
import yaml
import xarray as xr
import tqdm

# from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *
from pySD.sdmout_src import sdtracing
from pySD.initsuperdropsbinary_src import *


from sdm_eurec4a.visulization import set_custom_rcParams

set_custom_rcParams()

path2CLEO = Path("/home/m/m301096/CLEO")
path2sdm_eurec4a = Path("/home/m/m301096/repositories/sdm-eurec4a")
sys.path.append(path2CLEO)  # for imports from pySD package
# sys.path.append(path2CLEO / "examples/exampleplotting") # for imports from example plotting package

# use paths to files
path2build = path2CLEO / "build"
configfile = path2CLEO / "eurec4a/experiment_02/src/config/rain1d_config.txt"
yaml_config_file = sys.argv[1]

# yaml_config_file = path2sdm_eurec4a / "data/model/input/example_input_18.yaml"



# %%
with open(yaml_config_file, 'r') as f:
    config_yaml = yaml.safe_load(f)

identification_type = config_yaml['cloud']["identification_type"]
cloud_id = config_yaml['cloud']["cloud_id"]

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile    = path2CLEO / "libs/cleoconstants.hpp"
# path and file names for plotting results
setupfile     = path2CLEO / f"data/output/raw/long_run/{identification_type}_{cloud_id}/rain1d_setup.txt"
dataset       = path2CLEO / f"data/output/raw/long_run/{identification_type}_{cloud_id}/rain1d_sol.zarr"

# get cloud imformation
cloud_id = config_yaml['cloud']['cloud_id']
identification_type = config_yaml['cloud']['identification_type']
savefigpath = path2CLEO / "results/experiment_02" / f"{identification_type}_{cloud_id}" # directory for saving figures
savefigpath.mkdir(exist_ok=True, parents=True)

# %% [markdown]
# Get the config, constants and also an initial superdroplet dataset

# %%
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)

# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(str(dataset), consts)
simple_ds = xr.open_dataset(dataset, engine="zarr",
                                consolidated=False);

# %% [markdown]
# ## Create a xarray dataset from the zarr file
# 
# 
# For this, the ``sdtracing.attribute_for_superdroplets_sample`` function will be used for a subset of superdroplet ids.
# Using a size of 100 makes sense, and is not too slow.
# 
# The subdatasets will be stored in a temporary folder.

# %% [markdown]
# To do so, we create a list of temporaray filenames

# %%
path2dataoutput = path2CLEO / "data/output"
TEMPORARY_DIR = path2dataoutput / "temporary"
OUTPUT_DIR = path2dataoutput / "processed/long_run/" f"{identification_type}_{cloud_id}"
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# %%
all_ids = np.arange(config['totnsupers'], dtype = int)
len(all_ids)
# all_ids_reshaped = all_ids.reshape(100)
# use each gridcell on its own
n = config['totnsupers'] / 256
all_ids = all_ids.reshape(int(n), -1)
# temporary file names based on seed

seed = 2
np.random.seed(seed)
hashs = np.random.choice(int(1e6), size = all_ids.shape[0])
str_func = np.vectorize(lambda x : TEMPORARY_DIR / f"{hex(x)}.nc")
temp_filepaths = str_func(hashs)

assert len(temp_filepaths) == all_ids.shape[0]


print(f"Temporary files will be stored in {TEMPORARY_DIR}")

# %%
attributes = ["xi", "radius", "coord3", "sdgbxindex", "msol"]

print(f"Processing {all_ids.shape[0]} files")

for ids, temp_filepath in tqdm.tqdm(zip(all_ids, temp_filepaths)):
    print(f"Processing {ids.min()}-{ids.max()} to \t{temp_filepath}")
    ids = list(ids)
    list_dataarrays = []
    for attr in attributes:
        data = sdtracing.attribute_for_superdroplets_sample(
            sddata,
            attr,
            ids=ids,
        )
        da = xr.DataArray(
            data=data,
            dims=["time", "sd_id"],
            coords={"time": simple_ds.time, "sd_id": ids}
        )
        da.name = attr
        list_dataarrays.append(da)

    ds = xr.merge(list_dataarrays)
    ds.to_netcdf(temp_filepath)

# %% [markdown]
# ## Combine datasets and store in folder

# %%
print(f"Combine datasets and store in folder {OUTPUT_DIR}")

full_dataset = xr.open_mfdataset(temp_filepaths, parallel = True)
full_dataset.to_netcdf(OUTPUT_DIR / "full_dataset.nc")


