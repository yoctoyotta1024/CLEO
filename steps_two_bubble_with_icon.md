# First Steps

### Clone ICON repository
mkdir /home/m/m300950/icon-mpim
git clone --recursive git@gitlab.dkrz.de:icon/icon-mpim.git

### Switch to CLEO two-way coupling branch
```
git switch cleo-twoway-coupling
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

### Change account in ``exp.aes_bubble.run``
```
#SBATCH --account=bm1183
```

### run ICON
```
cd /work/bm1183/m300950/icon-mpim/build/run
sbatch ./exp.aes_bubble.run
```


# Steps to run ICON bubble with CLEO microphysics

### Build CLEO bubble3d executable and create CLEO input files
``
git switch mpi_yac_step3_twoway_step1
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
