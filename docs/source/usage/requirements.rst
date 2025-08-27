Requirements
============

The following requirements ensure CLEO's build, compilation and execution on DKRZ's Levante HPC.
If they do not work, please :ref:`contact us <contact>` or `open a new
issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository.

Of course other architectures, other compilers, versions etc. are possible, but we leave this for
you to discover.


The full list of packages CLEO uses for builds using the intel or gcc compiler
on Levante can be found in our
`levante packages script <https://github.com/yoctoyotta1024/CLEO/blob/main/scripts/levante/bash/src/levante_packages.sh>`_
and similarly for JUWELS in our
`juwels packages script <https://github.com/yoctoyotta1024/CLEO/blob/main/scripts/juwels/bash/src/juwels_packages.sh>`_

Compilers
---------
A C++ compiler with the C++20 standard library and MPI is the absolute minimum.

On Levante you can use the latest MPI compiler wrappers for the gcc or intel compilers.
At the time of writing for gcc 11.2.0, these are:

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0

To compile with C++ and CUDA, use Levante's nvhpc compilers too, e.g.

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0
  $ spack load cuda@12.2.0%gcc@=11.2.0

(Note you still need to use the gcc openmpi library as opposed to the cuda (nvhpc)
ones because CLEO uses gcc not nvhpc where this is relevant).

CMake
-----
CMake minimum version 3.18.0.

On Levante it's best to use version 3.26.3 which can be loaded e.g. for gcc 11.2.0 via:

.. code-block:: console

  $ cmake@3.26.3%gcc@=11.2.0/fuvwuhz


YAC
---

YAC is one of the :doc:`external libraries<extern>` which CLEO requires for its configuration
library (for MPI domain decomposition with/without YAC) and in order to couple to dynamics via YAC.

YAC (and its YAXT dependency) need to be installed manually before you can build CLEO with them.
Please refer to the instructions on how to do install YAC (and YAXT) in our
:doc:`external libraries<extern>` page.

YAC also requires some additional MPI, NetCDF and yaml libraries alongside the compatible gcc/intel
compiler. For example for gcc 11.2.0, you can load these on Levante via:

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0 openmpi@4.1.2%gcc@11.2.0 netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
  $ spack load openblas@0.3.18%gcc@=11.2.0

When you want to run CLEO, you will also need to export the fyaml library and python bindings paths,
e.g.

.. code-block:: console

  $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
  $ export PYTHONPATH=${PYTHONPATH}:/your/path/to/yac/python/


Python
------
To use ``cleopy`` you need Python minimum version 3.13.0. We advise you to use the
``pyproject.toml`` file via `uv <https://docs.astral.sh/uv/>`_ to manage the python version and
packages you need for ``cleopy`` and other parts of CLEO, e.g to run the examples.
If using the ``pyproject.toml`` does not automatically include all the version control and
dependencies you require, please :ref:`contact us <contact>` or
`open a new issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository.

To install python with all of CLEO's dependencies, please do
``uv sync --extra examples --extra yac``. If you only require only the dependencies for ``cleopy``,
you can instead do ``uv sync --no-dev``. If there are any missing version controls or dependencies,
you will have to install them too. We kindly ask that you also :ref:`contact us <contact>`
or `open a new issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub
repository to notify us.

To run your python for CLEO or find it's location you can do  ``uv run python`` and
``uv python find`` respectively.\
