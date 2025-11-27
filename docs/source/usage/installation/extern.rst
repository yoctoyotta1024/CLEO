External Libraries
==================

CLEO depends upon ``Kokkos`` and ``YAC`` and may depend upon some additional external libraries
such as ``yaml-cpp`` or ``CVODE`` depending on your setup.

YAC
---
YAC is required if CLEO couples to dynamics using YAC and/or uses MPI domain decompoisiton. You can
find more information about it from `its documentation <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_.

To build CLEO with dependency on YAC, you will first need to install YAXT and YAC manually.
(YAXT is a dependency of YAC.)

How to Install YAC (and YAXT)
#############################

To help install YAXT and YAC, we have created bash scripts for installation on
DKRZ's Levante HPC or on a generic/arbitrary, so-called "vanilla", computer in
``scripts/levante/bash/`` and ``scripts/vanilla/bash/`` respectively. First you need to ensure
your python version has pip installed first, e.g. with ``uv`` via

.. code-block:: console

  $ uv run python -m ensurepip --upgrade

Then to use ``install_yac.sh`` you will need to provide the path to the directory where you want
to put the installations, as well as the compiler type (e.g. ``intel`` or ``gcc``), and the python
interpreter you want to use for YAC's python bindings.

.. admonition:: E.g. On Levante

  for intel:

  .. code-block:: console

    $ scripts/levante/bash/install_yac.sh /work/bm1183/m300950/yacyaxt/intel/ intel /home/m/m300950/CLEO/.venv/bin/python3

  for gcc:

  .. code-block:: console

    $ scripts/levante/bash/install_yac.sh /work/bm1183/m300950/yacyaxt/gcc/ gcc /home/m/m300950/CLEO/.venv/bin/python3


.. admonition:: E.g. On "vanilla" computer

  (With a gcc compiler:) First you may need to edit the paths to certain packages
  in ``scripts/vanilla/bash/vanilla_packages.sh`` to match your system setup. You may also need to install YAC/YAXT's
  dependencies, `as described in YAC's documentation <https://dkrz-sw.gitlab-pages.dkrz.de/yac/d1/d9f/installing_yac.html>`_
  Then you should be able to run:

  .. code-block:: console

    $ scripts/vanilla/bash/install_yac.sh ${CLEO_YACYAXTROOT} gcc ${CLEO_PYTHON}


Alternatively you can download `YAXT <https://swprojects.dkrz.de/redmine/>`_ and
`YAC <https://gitlab.dkrz.de/dkrz-sw/yac/>`_ as compressed files and then configure and compile
them yourself.

Kokkos
------
All builds of CLEO require Kokkos in order to implement thread parallelism. You can read more about
how we use Kokkos on :doc:`our page about Kokkos<../../intro/kokkos>`.

The Kokkos libaries for CLEO are automatically built using CMAKE and compiled if required.

yaml-cpp
--------
CLEO's ``configuration`` library depends on the ```yaml-cpp``` package to read and write YAML files. You
can find more information about it from `the yaml-cpp repository: <https://github.com/jbeder/yaml-cpp>`_.

The yaml-cpp library for CLEO is automatically built using CMAKE and compiled if required.

CVODE
-----
CLEO's ``coupldyn_cvode`` library requires the SUNDIALS CVODE package. You can find more information
about it from `its webpage: <https://computing.llnl.gov/projects/sundials/cvode>`_.

The CVODE libraries for CLEO are automatically built using CMAKE and compiled if required.

pybind11
--------
CLEO's ``cleo_python_bindings`` library requires pybind11 to create python binding for selected parts of CLEO's
libraries. You can find more information about it from `the pybind11 repository: <https://github.com/pybind/pybind11>`_.

The pybind11 library for CLEO is automatically built using CMAKE and compiled if required.
*Note*: you can avoid making CLEO's python bindings, and thereby avoid the pybind11
dependency on your build, by passing a non-empty "CLEO_NO_PYBINDINGS" flag to
cmake. E.g. ``cmake [...] -DCLEO_NO_PYBINDINGS=true``
