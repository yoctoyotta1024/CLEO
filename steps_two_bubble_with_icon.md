# First Steps

### Clone ICON-MPIM repository
``
mkdir /home/m/m300950/icon-mpim
git clone --recursive git@gitlab.dkrz.de:icon/icon-mpim.git
``

or to update ICON-MPIM repo:
- pull ``master`` branch and rebase your relevant branch(es)
- ``git submodule update --init --recursive``

### Switch to CLEO two-way coupling branch
```
git switch cleo-twoway-branch2
```

### Build ICON in build directory
```
mkdir /work/mh0731/m300950/icon-mpim/build
cd /work/mh0731/m300950/icon-mpim/build
/home/m/m300950/icon-mpim/config/dkrz/levante.gcc-11.2.0 --enable-openmp
make -j 16
```


# Steps to run default ICON bubble with 1-moment microphysics

### Install mkexp and set path to build

([see ICON docs](https://docs.icon-model.org/documentation/buildrun/buildrun_running.html#ref-buildrun-running))

```
cd /home/m/m300950/icon-mpim/
mamba activate clouds
python -m pip install utils/mkexp/

export MKEXP_PATH=.:/work/mh0731/m300950/icon-mpim/build/run/
```

### Generate Default ICON bubble run script
```
cd /home/m/m300950/icon-mpim/run
cp ./examples/bubble.config ./bubble_1mom.config
```
then add ACCOUNT into .config file ``ACCOUNT = mh0731`` beneath ``EXP_TYPE = torus``.

then generate run script(s):
```
/home/m/m300950/icon-mpim//utils/mkexp/mkexp bubble_1mom.config
```

### run ICON
```
cd /home/m/m300950/icon-mpim/experiments/bubble_1mom/scripts
sbatch bubble_1mom.run_start
```

### Important Directories including Output
The output will be in the "work directory" given by mkexp after ``[...]/mkexp bubble_1mom.config``,
e.g.
```
Script directory: '/home/m/m300950/icon-mpim/experiments/bubble_1mom/scripts'
Data directory: '/work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/outdata'
Work directory: '/scratch/m/m300950/icon-mpim/experiments/bubble_1mom/work/run_[XXX]'
Log directory: '/work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/log'
```

e.g. probably you then want to move it into outdata e.g.
```
mv /scratch/m/m300950/icon-mpim/experiments/bubble_1mom/work/run_20080801T000000-20080801T015930/bubble_1mom_atm_cgrid_ml.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/outdata/
mv /scratch/m/m300950/icon-mpim/experiments/bubble_1mom/work/run_20080801T000000-20080801T015930/bubble_1mom_atm_2d_ml_20080801T000000Z.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/outdata/
mv /scratch/m/m300950/icon-mpim/experiments/bubble_1mom/work/run_20080801T000000-20080801T015930/bubble_1mom_atm_3d_ml_20080801T000000Z.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/outdata/
```


# Steps to run ICON bubble with CLEO microphysics (default one-way coupling)

### Build CLEO bubble3d executable and create CLEO input files
```
git switch testing-mpi_yac_step3_twoway_step2_withmpi-bubble-helping
```

edit ``path2build=${HOME}/CLEO/build_bubble3d/`` and
``script_args="${src_config_filename} --do_inputfiles"`` in ``bubble3d.sh``,
```
vim /home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

*NOTE*: to run ICON-YAC-CLEO you will need to have the same version of YAC for CLEO as that used by ICON (in ``icon-mpim/externals/yac``).
At the time of writing this is YAC v3.16.0 ([see here](https://gitlab.dkrz.de/dkrz-sw/yac/-/blob/master/CMakeLists.txt?ref_type=heads)). The easiest way
to ensure this is to build Cleo with ICON's YAC and YAXT by setting the ``CLEO_YACYAXTROOT`` to point to the icon-mpim externals directory,
e.g. in ``/path/to/cleo/scripts/levante/examples/build_compile_run_plot.sh``, set ``yacyaxtroot=/home/m/m300950/icon-mpim/externals``.

run script bubble3d script,
```
/home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

### Copy ICON CLEO-bubble run_start script (and create empty logfiles)

``bubble_1mom.run_start`` and ``bubble_cleo_iconopenmp.run_start`` are
possible drafts you could use for running ICON with its 1 moment scheme +- Cleo coupled.
Otherwise you can start from your own ``bubble_1mom.run_start``:

```
cp /home/m/m300950/icon-mpim/experiments/bubble_1mom/scripts/bubble_1mom.run_start /home/m/m300950/icon-mpim/experiments/bubble_cleo/scripts/bubble_cleo.run_start
```

create empty logfiles in icon-mpim experiment source location,
e.g. ``bubble_cleo.dump  bubble_cleo.log  bubble_cleo.run.log`` in ``/home/m/m300950/icon-mpim/experiments/bubble_cleo``


### Change CLEO params ``bubble_cleo.run_start``

First change any instances of ``bubble_1mom`` to ``bubble_cleo``. Then adapt the run to turn on the cleo coupling:
```
#SBATCH --account=mh0731
#SBATCH --nodes=2

[...]

# Environment variables for the target system
[...]
# export paths for CLEO microphysics
# export PYTHON="/home/m/m300950/CLEO/.venv/bin/python3"
# export PYTHONPATH="/work/mh0731/m300950/yacyaxt/gcc/yac/python:${PYTHONPATH}"
export LD_LIBRARY_PATH="/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib:${LD_LIBRARY_PATH}"
export ICON_MODEL="/work/mh0731/m300950/icon-mpim/build/bin/icon"
export CLEO_MODEL="/work/mh0731/m300950/icon-mpim/build_cleo/examples/bubble3d/src/bubble3d"
export CLEO_CONFIGFILE="/work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/tmp/bubble3d_config.yaml"

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
cd /home/m/m300950/icon-mpim/experiments/bubble_cleo/scripts
sbatch bubble_cleo.run_start
```

### Important Directories including Output
The output will be in the "work directory" of ``bubble_cleo`` analagously to
that given by mkexp after ``[...]/mkexp bubble_1mom.config``,
e.g.
```
Script directory: '/home/m/m300950/icon-mpim/experiments/bubble_cleo/scripts'
Data directory: '/work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/outdata'
Work directory: '/scratch/m/m300950/icon-mpim/experiments/bubble_cleo/work/run_[XXX]'
Log directory: '/work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/log'
```

e.g. probably you then want to move it into outdata e.g.
```
mkdir /work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/outdata
mv /scratch/m/m300950/icon-mpim/experiments/bubble_cleo/work/run_20080801T000000-20080801T015930/bubble_cleo_atm_cgrid_ml.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/outdata/
mv /scratch/m/m300950/icon-mpim/experiments/bubble_cleo/work/run_20080801T000000-20080801T015930/bubble_cleo_atm_2d_ml_20080801T000000Z.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/outdata/
mv /scratch/m/m300950/icon-mpim/experiments/bubble_cleo/work/run_20080801T000000-20080801T015930/bubble_cleo_atm_3d_ml_20080801T000000Z.nc \
  /work/mh0731/m300950/icon-mpim/experiments/bubble_cleo/outdata/
```


# Steps to run ICON bubble with CLEO microphysics, two-way coupling

### Edit ICON source code to activate two-way coupling

Edit ``oneway_coupling`` --> ``twoway_coupling`` at TWO locations of ICON source code

In ``icon-mpim/src/coupling/mo_atmo_coupling_frame.f90``:

```
CALL construct_atmo_cleo_coupling_post_sync(p_patch(jg),
  comp_id, cell_point_id(1), patch_horz%n_patch_cells, timestepstring,
  twoway_coupling)
```

In ``icon-mpim/src/atm_phy_aes/mo_aes_phy_main.f90``:
```
if(is_coupled_to_cleo()) CALL interface_cleo(jg, twoway_coupling)
```

Then (Re-)Compile and follow steps to run ICON bubble with CLEO microphysics as above.

Note: This test of coupling uses CLEO by default with NullMicrophysicalProcess. Two-way coupling
with CLEO's microphysical procsses enabled (i.e. not Null) would need to first
de-activate one-moment scheme...


# Steps to run ICON bubble with CLEO microphysics, one/two-way coupling > 1 MPI process

in your ``bubble_cleo.run_start`` run script, simply edit ``#SBATCH --nodes=X+1`` (one extra for ICON)
and ``srun [...] : -n X $CLEO_MODEL $CLEO_CONFIGFILE`` to the number of desired MPI processes for CLEO, ``X``:

Note: you may also need to comment out the MassMomentsObservers (``obs4`` and ``obs5``)
in ``main_bubble3d.cpp``
