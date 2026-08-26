# First Steps

### Clone ICON repository
``
mkdir /home/m/m300950/icon
git clone --recursive git@git.gitlab.dkrz.de:icon/icon.git
``

or to update ICON repo:
- pull ``main`` branch and rebase your relevant branch(es)
- ``git submodule update --init --recursive``

### Switch to CLEO two-way coupling branch
```
git switch cleo-twoway-coupling
```

### Build ICON in build directory
```
mkdir -p /work/mh0731/m300950/icon/build_icon
cd /work/mh0731/m300950/icon/build_icon
/home/m/m300950/icon/icon/config/dkrz/levante.gcc-11.2.0 --enable-openmp
make -j 16
```


# Steps to run default ICON bubble with 1-moment microphysics

### Install mkexp and set path to build

([see ICON docs](https://docs.icon-model.org/documentation/buildrun/buildrun_running.html#ref-buildrun-running))

```
cd /home/m/m300950/icon/icon/
mamba activate clouds
python -m pip install utils/mkexp/

export MKEXP_PATH=.:/work/mh0731/m300950/icon/build_icon/run/
```

### Generate Default ICON bubble run script
```
cd /home/m/m300950/icon/icon/run
cp ./examples/bubble.config ./bubble_1mom.config
```
then add ACCOUNT into .config file ``ACCOUNT = mh0731`` beneath ``EXP_TYPE = torus``.

then generate run script(s):
```
/home/m/m300950/icon/icon/utils/mkexp/mkexp bubble_1mom.config
```

### run ICON
```
cd /home/m/m300950/icon/icon/experiments/bubble_1mom/scripts
sbatch bubble_1mom.run_start
```

### Important Directories including Output
The output will be in the "work directory" given by mkexp after ``[...]/mkexp bubble_1mom.config``,
e.g.
```
Script directory: '/home/m/m300950/icon/icon/experiments/bubble_1mom/scripts'
Data directory: '/work/mh0731/m300950/icon/icon/experiments/bubble_1mom/outdata'
Work directory: '/scratch/m/m300950/icon/icon/experiments/bubble_1mom/work'
Log directory: '/work/mh0731/m300950/icon/icon/experiments/bubble_1mom/log'
```

e.g. probably you then want to move the bubble_1mom_atm .nc files into outdata e.g.
```
mv /scratch/m/m300950/icon/icon/experiments/bubble_1mom/work/run_20080801T000000-20080801T015930/bubble_1mom_atm_* \
  /work/mh0731/m300950/icon/icon/experiments/bubble_1mom/outdata/

```


# Steps to run ICON bubble with CLEO microphysics (default one-way coupling)

### Build CLEO bubble3d executable and create CLEO input files
```
git switch testing-mpi_yac_step3_twoway_step2_withmpi-bubble-helping
```

edit the ``bubble.sh`` bash script settings: ``path2build``, ``path2experiment``, ``path2iconfiles``,
and set ``script_args="${src_config_filename} ${src_yac_config_filename} ${path2experiment} ${path2iconfiles --do_inputfiles"``:
```
vim /home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

*NOTE*: to run ICON-YAC-CLEO you will need to have the same version of YAC for CLEO as that used by ICON (in ``icon/externals/yac``).
At the time of writing this is YAC v3.20.2 ([see here](https://gitlab.dkrz.de/icon/icon/-/tree/main/externals?ref_type=heads)). The easiest way
is to build YAC and YAXT from the icon source code with the ``install_yac.sh`` bash script edited as follows:
```
# yaxt_tag=0.12.1
# yaxt_version=yaxt-${yaxt_tag}
# yaxt_release_tag=release-${yaxt_tag}
# yaxt_source=https://gitlab.dkrz.de/dkrz-sw/yaxt/-/archive/$yaxt_release_tag/$yaxt_version.tar.gz
export yaxt_src=/home/m/m300950/icon/icon/externals/yaxt

# yac_tag=v3.20.2
# yac_version=yac_$yac_tag
# yac_source=https://gitlab.dkrz.de/dkrz-sw/yac/-/archive/$yac_tag/$yac_version.tar.gz
export yac_src=/home/m/m300950/icon/icon/externals/yac

[...]

  ### --------------------- install YAXT ------------------- ###
  # mkdir ${root4YAC}/${yaxt_version}
  # cd ${root4YAC}/${yaxt_version} &&  pwd
  # curl -s -L ${yaxt_source} | tar xvz --strip-components=1
  # mkdir build && cd build
  mkdir -p ${root4YAC}/yaxt/build && cd ${root4YAC}/yaxt/build && pwd
  ${yaxt_src}/configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    --without-regard-for-quality \
    --without-example-programs \
    --without-perf-programs \
    --with-pic \
    --prefix=${root4YAC}/yaxt
  make -j 8
  make install
  # cd ${root4YAC} && rm -rf ${yaxt_version}
  ### ------------------------------------------------------ ###

  ## --------------------- install YAC -------------------- ###
  # python bindings made in yac_version directory (note this is not yac directory!)
  # mkdir ${root4YAC}/${yac_version}
  # cd ${root4YAC}/${yac_version} && pwd
  # curl -s -L ${yac_source} | tar xvz --strip-components=1
  # mkdir build && cd build
  mkdir -p ${root4YAC}/yac/build && cd ${root4YAC}/yac/build && pwd
  ${yac_src}/configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    LDFLAGS="-lm" \
    PYTHON=${python} \
    --disable-mpi-checks \
    --with-yaxt-root=${root4YAC}/yaxt \
    --with-netcdf-root=${netcdf_root} \
    --with-fyaml-root=${fyaml_root} \
    --enable-python-bindings \
    --enable-rpaths \
    --with-pic \
    --prefix=${root4YAC}/yac
  make -j 8
  make install

  # mv ${root4YAC}/${yac_version}/build/python ${root4YAC}/yac/
  # cd ${root4YAC} && rm -rf ${yac_version}
  ### ------------------------------------------------------ ###
```

Then make sure Cleo's build for the bubble points to these YAC and YAXT installations
by setting the ``CLEO_YACYAXTROOT`` accordingly in
``/home/m/m300950/CLEO/scripts/levante/examples/build_compile_run_plot.sh``.

Make the relevant directories for the bubble input/output (must match ICON run_start, see below)
```
cd /work/mh0731/m300950/icon/icon/experiments/ && mkdir bubble_cleo
cd bubble_cleo && mkdir bin log outdata scripts share tmp work
```

Finally, run the ``bubble3d`` example's script,
```
cd /home/m/m300950/CLEO
/home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

### Copy ICON CLEO-bubble run_start script (and create empty logfiles)

``/home/m/m300950/CLEO/bubble_1mom.run_start`` and ``/home/m/m300950/CLEO/bubble_cleo.run_start`` are
possible drafts you could use for running ICON with its 1 moment scheme and/or Cleo one-way/two-way coupled.
Otherwise you can start from your own ``bubble_1mom.run_start``:

```
mkdir -p /home/m/m300950/icon/icon/experiments/bubble_cleo/scripts
cp /home/m/m300950/icon/icon/experiments/bubble_1mom/scripts/bubble_1mom.run_start /home/m/m300950/icon/icon/experiments/bubble_cleo/scripts/bubble_cleo.run_start
```


### Change CLEO params ``bubble_cleo.run_start``

First change any instances of ``bubble_1mom`` to ``bubble_cleo``. Then adapt the run to turn on the cleo coupling:
```
#SBATCH --account=mh0731
#SBATCH --nodes=2
#SBATCH --time=00:20:00
[...]

# Environment variables for the target system
[...]
# export paths for CLEO microphysics
# export PYTHON="/home/m/m300950/CLEO/.venv/bin/python3"
# export PYTHONPATH="/work/mh0731/m300950/yacyaxt/gcc/yac/python:${PYTHONPATH}"
export LD_LIBRARY_PATH="/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib:${LD_LIBRARY_PATH}"
export ICON_MODEL="/work/mh0731/m300950/icon/build/bin/icon"
export CLEO_MODEL="/work/mh0731/m300950/icon/build_cleo/examples/bubble3d/src/bubble3d"
export CLEO_CONFIGFILE="/work/mh0731/m300950/icon/icon/experiments/bubble_cleo/tmp/bubble3d_config.yaml"

[...]

&coupling_mode_nml
    coupled_to_cleo = .TRUE.
    coupled_to_ocean = .false.
/

[...]

# Call processes
srun -l --kill-on-bad-exit=1 --cpu-bind=quiet,cores --distribution=block:block --propagate=STACK,CORE -c 4 -n 2 \
    $ICON_MODEL : -n 1 $CLEO_MODEL $CLEO_CONFIGFILE
```

### run ICON with CLEO
```
cd /home/m/m300950/icon/icon/experiments/bubble_cleo/scripts
sbatch bubble_cleo.run_start
```

### Important Directories including Output
The output will be in the "work directory" of ``bubble_cleo`` analagously to
that given by mkexp after ``[...]/mkexp bubble_1mom.config``,
e.g.
```
Script directory: '/home/m/m300950/icon/icon/experiments/bubble_cleo/scripts'
Data directory: '/work/mh0731/m300950/icon/icon/experiments/bubble_cleo/outdata'
Work directory: '/scratch/m/m300950/icon/icon/experiments/bubble_cleo/work/run_[XXX]'
Log directory: '/work/mh0731/m300950/icon/icon/experiments/bubble_cleo/log'
```

e.g. probably you then want to move the bubble_1mom_atm .nc files into outdata e.g.
```
mv /scratch/m/m300950/icon/icon/experiments/bubble_cleo/work/run_20080801T000000-20080801T015930/bubble_cleo_atm_* \
  /work/mh0731/m300950/icon/icon/experiments/bubble_cleo/outdata
```


# Steps to run ICON bubble with CLEO two-way coupling

Edit Cleo's YAC coupling YAML file to activate two-way coupling; in
in ``examples/bubble3d/src/config/yac_icon_cleo_coupling_config.yaml`` simply un-comment
the two-way coupling settings:

```
# CLEO -> ATM
  - <<: [*cleo2atm, *nnn_interp_stack]
    field:
      - src: liquid_water_mixing_ratio_out
        tgt: liquid_water_mixing_ratio_from_cleo
      - src: humidity_mixing_ratio_out
        tgt: humidity_mixing_ratio_from_cleo
  - <<: [*cleo2atm, *nnn_interp_stack]
    field:
      - src: air_temperature_out
        tgt: air_temperature_from_cleo
```


# Steps to run ICON bubble with CLEO > 1 MPI process

in your ``bubble_cleo.run_start`` run script, simply edit ``#SBATCH --nodes=X+1`` (one extra for ICON)
and ``srun [...] : -n X $CLEO_MODEL $CLEO_CONFIGFILE`` to the number of desired MPI processes for CLEO, ``X``:

*Note:* you may also need to comment out the MassMomentsObservers (``obs4`` and ``obs5``)
in ``main_bubble3d.cpp``

# Debugging

By adding to ``yac_icon_cleo_coupling_config.yaml`` the following lines:
```
debug:
  coordinates_mismatch_is_fatal: false   # warns instead of aborting
  global_config:
    enddef: /work/mh0731/m300950/icon/icon/experiments/bubble_cleo/bin/coupling_config_debug.yaml
  output_grids:
    - grid_name: cleo_cartesian_grid
      file_name: /work/mh0731/m300950/icon/icon/experiments/bubble_cleo/bin/all_grids_debug.nc
    - grid_name: icon_atmos_grid
      file_name: /work/mh0731/m300950/icon/icon/experiments/bubble_cleo/bin/all_grids_debug.nc
```
and removing from ``bubble3d_config.yaml``:
```
  yac_debug_config_file: [...]
  yac_debug_grid_file: [...]
```

You can output the full grid and coupling configuration of ICON-Cleo. You can then plot the grids, e.g.
with ``python ~/icon/icon/externals/yac/tools/yac_plot_grids.py all_grids_debug_gid.nc --projection platecarree --radius 0``.
*Note:* you will probably first need to post-process the Cleo grid output with ``grids_gid_fix.py``.

Similarly you can visualise the coupling configuration with:
``python ~/icon/icon/externals/yac/tools/yac_plot_coupling_config.py coupling_config_debug.yaml``
which produces ``coupling.svg``.

You can output the weights for a particular coupling too with
``python ~/icon/icon/externals/yac/tools/yac_plot_weights.py weights_temp_debug.nc all_grids_debug_gid.nc --projection platecarree --radius 0``
after you have added the weights file to ``yac_icon_cleo_coupling_config.yaml`` for a particle field, e.g.
```yaml
coupling:
# ATM -> CLEO
  - <<: [*atm2cleo, *nnn_interp_stack]
    field:
      - src: eastward_wind_to_cleo
        tgt: eastward_wind_in
      - src: northward_wind_to_cleo
        tgt: northward_wind_in
      - src: upward_air_velocity_to_cleo
        tgt: upward_air_velocity_in
    weight_file:
      name: /work/mh0731/m300950/icon/icon/experiments/bubble_cleo/bin/weights_winds_debug.nc
      on_existing: overwrite
```
