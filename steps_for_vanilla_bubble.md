# BUILD ICON
Do this before icon-mpim generic configure wrapper:

1. brew install cmake eccodes gcc hdf5 libfyaml libxml2 netcdf netcdf-fortran open-mpi
2. export OMPI_CC=gcc-16 export OMPI_CXX=g++-16 export OMPI_FC=gfortran-16
3. which gcc-16 g++-16 gfortran-16 mpicc mpicxx mpif90
4. cd /Users/yoctoyotta1024/Documents/d2_springsummer2026/debug_bubble/build
5. /Users/yoctoyotta1024/Documents/icon-mpim/config/generic/gcc --enable-openmp --enable-yaxt
6. make -j4 2>&1 | tee make.log

#### Notes:
- Maybe you need ``export ICON_SW_PREFIX=/opt/homebrew``
* Since you're on OpenMPI (which has known macOS quirks per the project notes), two things worth having ready:
    * Oversubscription: OpenMPI refuses to run more ranks than physical cores by default, which bites on a laptop. export OMPI_MCA_rmaps_base_oversubscribe=1 (relevant at runtime, not configure time).
    * If you hit MPI startup hangs later, the fix is export TMPDIR=/tmp at runtime (and possibly at configure time via BUILD_ENV="export TMPDIR='/tmp';").


# RUN 1mom BUBBLE
1.  mamba activate cleoenv
2.  python -m pip install jinja2 six
3.  cd run
4.  cp examples/bubble.config bubble_1mom.config
5.  /Users/yoctoyotta1024/Documents/icon-mpim/utils/mkexp/mkexp bubble_1mom.config
6.  cd /Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_1mom/scripts
7.  bash bubble_1mom.run_start > ../log/bubble_1mom.run.log


Output:
Script directory: '/Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_1mom/scripts'
Data directory: '/Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_1mom/outdata'
Work directory: '/Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_1mom/work'
Log directory: '/Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_1mom/log'

Note Settings:
```
export OMP_NUM_THREADS='1'
[...]
mpi_icon_procs=4
```
=> 4 MPI ranks * 1 OpenMP Threads => 4 logical cores for ICON (check compatible with ``sysctl -n hw.ncpu`` (FYI: can also see sysctl -n hw.physicalcpu)


# BUILD CLEO

#### Compile and Setup Cleo (one-way or two-way coupled)

Compile cleo on ``bubble_vanilla`` Cleo branch with ``vanilla_inputfiles_compile_bubble.patch`` applied,
and optionally the ``vanilla_40min_debugging.patch`` or ``vanilla_twoway.patch`` patches.

#### Install ICON YAXT and YAC versions

*NOTE*: you may get the CMake error:
```
CMake Error at /opt/homebrew/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:290 (message):
  Could NOT find YAXT (missing: YAXT_Fortran_LIBRARY YAXT_C_LIBRARY)
Call Stack (most recent call first):
  /opt/homebrew/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:654 (_FPHSA_FAILURE_MESSAGE)
  libs/coupldyn_yac/cmake/FindYAXT.cmake:24 (find_package_handle_standard_args)
  libs/coupldyn_yac/cmake/FindYAC.cmake:1 (find_package)
  libs/configuration/CMakeLists.txt:26 (find_package)
```
and/or
```
CMake Error at /opt/homebrew/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:290 (message):
  Could NOT find YAC (missing: YAC_C_LIBRARY YAC_C_INCLUDE_DIR
  YAC_C_MTIME_LIBRARY)
Call Stack (most recent call first):
  /opt/homebrew/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:654 (_FPHSA_FAILURE_MESSAGE)
  libs/coupldyn_yac/cmake/FindYAC.cmake:30 (find_package_handle_standard_args)
  libs/configuration/CMakeLists.txt:26 (find_package)
```

Then you will have to build YAXT and YAC from ICON as follows:
``` bash
# compilers
export CC=gcc-16
export CXX=g++-16
export OMPI_CC=gcc-16
export OMPI_CXX=g++-16
export OMPI_FC=gfortran-16

# src directories
export yaxt_src=/Users/yoctoyotta1024/Documents/icon-mpim/externals/yaxt
export yac_src=/Users/yoctoyotta1024/Documents/icon-mpim/externals/yac

# build directory
export root4YAC=/Users/yoctoyotta1024/Documents/d2_springsummer2026/debug_bubble/build/yacyaxt
mkdir ${root4YAC}

# checks:
echo ${yaxt_src} && ls ${yaxt_src}
echo ${yac_src} && ls ${yac_src}
echo ${root4YAC} && ls ${root4YAC}

# install YAXT
mkdir -p ${root4YAC}/yaxt/build && cd ${root4YAC}/yaxt/build && pwd
${yaxt_src}/configure \
  CC=mpicc FC=mpif90 \
  CFLAGS="-O0 -g -Wall" \
  FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
  --without-regard-for-quality \
  --without-example-programs \
  --without-perf-programs \
  --with-pic \
  --prefix=${root4YAC}/yaxt
make -j 8
make install

# install YAC
mkdir -p ${root4YAC}/yac/build && cd ${root4YAC}/yac/build && pwd
${yac_src}/configure \
    CC=mpicc FC=mpif90 \
    CFLAGS="-O0 -g -Wall" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    LDFLAGS="-lm" \
    PYTHON=/Users/yoctoyotta1024/Documents/CLEO/.venv/bin/python \
    --disable-mpi-checks \
    --with-yaxt-root=${root4YAC}/yaxt \
    --with-netcdf-root=/opt/homebrew \
    --with-fyaml-root=/opt/homebrew \
    --enable-python-bindings \
    --enable-rpaths \
    --with-pic \
    --prefix=${root4YAC}/yac
make -j 8
make install
```


# RUN AND PLOT ICON-CLEO (one-way or two-way coupled)

- run with ``bash bubble_cleo.run_start > ../log/bubble_cleo.run.log`` in ``/Users/yoctoyotta1024/Documents/icon-mpim/experiments/bubble_cleo/scripts``
- move relevatn *.nc files from bubble_cleo/work to bubble_cleo/outdata, like for bubble_1mom
- plot results on ``testing-mpi_yac_bubble_example-bubble-helping`` Cleo branch with ``vanilla_plotting.patch`` applied.
