Requirements
============

The following requirements ensure CLEO's build, compilation and execution on DKRZ's Levante HPC.
If they do not work, please :ref:`contact us <contact>` or `open a new
issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository.

Of course other architectures, other compilers, versions etc. are possible, but we leave this for
you to discover.

CMake
-----
CMake minimum version 3.18.0.

On Levante it's best to use version 3.23.1 which can be loaded e.g. via the command

.. code-block:: console

  $ spack load cmake@3.23.1%gcc

Compilers
---------
A C++ compiler with the C++20 standard library is the absolute minimum.

On Levante you can use the latest gcc compilers. At the time of writing this is gcc 11.2.0, e.g.

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0

To compile with CUDA, use Levante's latest nvhpc compilers, e.g.

.. code-block:: console

  $ module load nvhpc/23.9-gcc-11.2.0

Python
------
To use PySD you need Python minimum version 3.10.4. We advise you to :ref:`create an
environment<environment>` using our envirnoment.yml file. This environment should automatically
include all the additional packages you may require. If not, please :ref:`contact us <contact>` or
`open a new issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository.

On Levante it's a good idea to load the python3 module, e.g.

.. code-block:: console

  $ module load python3/2022.01-gcc-11.2.0

To use PySD and to run CLEO's examples, particular Python packages are needed. These are included in
our environment.yml file and are the following: ``matplotlib``, ``numpy``, ``scipy``, ``xarray``,
``zarr``, and ``awkward``. If there are other dependencies not listed here, you will have to install
them too. We kindly ask that you also :ref:`contact us <contact>` or `open a new
issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository to notify us.

You can install Python packages to an existing Conda (or Mamba) environment via:

.. code-block:: console

  $ conda activate [your conda environment]
  $ python -m pip install [package name(s)]

YAC
---

YAC is one of the :doc:`external libraries<extern>` which CLEO may require in order to
couple to dynamics and/or have MPI domain decomposition.

Note that YAC (and its YAXT dependency) need to be installed manually before you can build
CLEO with them. You can find instructions on how to do install YAC (and YAXT) in the
external libraries section.

YAC also requires some additional MPI, NetCDF and yaml libraries alongside the compatible gcc
compiler which you can load on Levante via:

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0 netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
  $ spack load openblas@0.3.18%gcc@=11.2.0

When you want to run CLEO with YAC, you will also need to export some additional paths:

.. code-block:: console

  $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
  $ export PYTHONPATH=/your/path/to/yac/python/
  $ spack load py-numpy
