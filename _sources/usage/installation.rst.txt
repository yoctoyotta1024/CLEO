.. _installation:

Installation
============

First clone CLEO's GitHub repository. Everything will be much easier for you if you clone CLEO in
your home directory. It’s not essential, but if you choose to do otherwise you may have to change
some extra paths in the bash and Python scripts.

.. code-block:: console

  $ git clone https://github.com/yoctoyotta1024/CLEO.git


Everything will also be much easier for you if you setup and use CLEO with a python version and
all the packages as listed in our ``pyproject.toml`` file. You can acheive this quickly and
simply using `uv <https://docs.astral.sh/uv/>`_ (*hint*: optionally, you can install ``uv`` within
a Conda/Micromamba environment using our ``environment.yml`` file,
e.g. ``mamba create --file=environment.yml``). Using ``uv`` to setup python, simply do:

.. code-block:: console

  $ uv sync --extra examples --extra yac

Alternatively, to only install the python dependencies required by CLEO's python package,
``cleopy``, you can just do ``uv sync --no-dev``.

*Note*: on Levante/HPCs you may need to set the paths to your mpi wrapper/libraries before
installing mpi4py in order to be able to run MPI via mpi4py. On Levante if you want to use openMPI
from ``module load openmpi/4.1.2-gcc-11.2.0``, you will need to uninstall the default
``mpi4py`` installation from ``uv sync [...]`` and re-install with the correct paths, e.g.

.. code-block:: console

  $ uv pip uninstall mpi4py
  $ export MPI4PY_BUILD_MPICC=/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc  # mpicc wrapper
  $ export MPI4PY_BUILD_MPILD=/sw/spack-levante/openmpi-4.1.2-mnmady/lib  # path to mpi libraries (libmpi.so)
  $ uv pip install --no-cache-dir --no-binary=mpi4py mpi4py
  $ uv run python -c 'from mpi4py import MPI; print(f"MPI version: {MPI.Get_version()}")'  # check installation

See the `mpi4py documentation <https://mpi4py.readthedocs.io/en/stable/install.html#build-backends>`_
for more information

For advanced users/developers, even if you already have ``uv`` installed elsewhere you may
want to install CLEO's non-python dependencies sourced from conda-forge via our ``environment.yml``,
(e.g. ``doxygen``), for examples with Micromamba: ``mamba create --file=environment.yml``.

Finally we suggest you use pre-commit. You can install our hooks via:

.. code-block:: console

  $ pre-commit install


That's it, you're done!

Now maybe you want to run some of :doc:`CLEO's examples <examples>` or try out
the :doc:`quickstart<quickstart>`...
