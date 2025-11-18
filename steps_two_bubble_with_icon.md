# First Steps

### Clone ICON repository
mkdir /home/m/m300950/icon-mpim
git clone --recursive git@gitlab.dkrz.de:icon/icon-mpim.git

### Switch to CLEO two-way coupling branch
```
git switch cleo-twoway-branch2
```

### Build ICON in build directory
```
mkdir /work/bm1183/m300950/icon-mpim/build
cd /work/bm1183/m300950/icon-mpim/build
/home/m/m300950/icon-mpim/config/dkrz/levante.gcc-11.2.0 --enable-openmp
make -j 16
```


# Steps to run default ICON bubble with 1-moment microphysics

### Generate Default ICON bubble run script
```
cd /work/bm1183/m300950/icon-mpim/build
./make_runscripts aes_bubble
```

### Change account in ``exp.aes_bubble.run`` and maybe output directory
```
#SBATCH --account=bm1183
#SBATCH --output=/work/bm1183/m300950/icon-mpim/build/run/LOG.exp.aes_bubble.run.%j.o
```

### run ICON
```
cd /work/bm1183/m300950/icon-mpim/build/run
sbatch ./exp.aes_bubble.run
```


# Steps to run ICON bubble with CLEO microphysics (default one-way coupling)

### Build CLEO bubble3d executable and create CLEO input files
```
git switch mpi_yac_step3_twoway_step2-withmpi-bubble_helping
```

edit ``path2build=${HOME}/CLEO/build_bubble3d/`` and
``script_args="${src_config_filename} --do_inputfiles"`` in ``bubble3d.sh``,
```
vim /home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

run script bubble3d script,
```
/home/m/m300950/CLEO/scripts/levante/examples/bubble3d.sh
```

### Copy ICON CLEO-bubble run script

_hint_: ``exp.aes_bubble_cleo_iconopenmp.run`` and ``exp.aes_bubble_cleo_iconopenmp.run`` are
possible drafts you could use for ``exp.aes_bubble_cleo.run``.

```
cp /path/to/your/copy/exp.aes_bubble_cleo.run /work/bm1183/m300950/icon-mpim/build/run/exp.aes_bubble_cleo.run
```

### Change CLEO params ``exp.aes_bubble_cleo.run``
```
#SBATCH --account=bm1183
#SBATCH --partition=compute
#SBATCH --nodes=2

module load openmpi/4.1.2-gcc-11.2.0 # for CLEO

[...]

# export the python interpreter
export PYTHON="/home/m/m300950/CLEO/.venv/bin/python3"
export PYTHONPATH=":/work/bm1183/m300950/icon-mpim/build/externals/mtime/build/python"
export PYTHONPATH="/work/bm1183/m300950/yacyaxt/gcc/yac/python:${PYTHONPATH}"

# export paths for CLEO microphysics
export LD_LIBRARY_PATH="/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib:${LD_LIBRARY_PATH}"
export CLEO_PATH2BUILD="/work/bm1183/m300950/icon-mpim/build_cleo"
export CLEO_MODEL="${CLEO_PATH2BUILD}/examples/bubble3d/src/bubble3d"
export CLEO_CONFIGFILE="${CLEO_PATH2BUILD}/tmp/bubble3d_config.yaml"

[...]

&coupling_mode_nml
  coupled_to_cleo = .TRUE.
/

[...]

START_MODEL_CLEO="-n 1 $CLEO_MODEL $CLEO_CONFIGFILE"
START_MODEL="${START_MODEL:=$START $MODEL : $START_MODEL_CLEO}"
```

### run ICON with CLEO
```
cd /work/bm1183/m300950/icon-mpim/build/run
sbatch ./exp.aes_bubble_cleo.run
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

in your ``exp.aes_bubble_cloe.run`` run script, simply edit ``--ntasks=X+1`` and
``START_MODEL_CLEO="-n X [...]`` to the number of desired MPI processes for CLEO, ``X``:

```
#SBATCH --ntasks=X+1 (one extra for ICON)

[...]

START_MODEL_CLEO="-n 1 $CLEO_MODEL $CLEO_CONFIGFILE"
```

Note: you may also need to comment out the MassMomentsObservers (``obs4`` and ``obs5``)
in ``main_bubble3d.cpp``
