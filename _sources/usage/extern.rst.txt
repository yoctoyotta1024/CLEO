External Libraries
==================

CLEO depends upon Kokkos and may depend upon some additional external libraries such as ``YAC`` and
``yaml-cpp`` depending on your setup.

Kokkos
------
All builds of CLEO require Kokkos in order to implement thread parallelism. You can read more about
how we use Kokkos on :doc:`our page about Kokkos<../intro/kokkos>`.

The Kokkos libaries for CLEO are automatically built using CMAKE and compiled if required.

YAC
---
YAC is required if CLEO couples to dynamics using YAC and/or uses MPI domain decompoisiton. You can
find more information about it from `its documentation: <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_.

To build CLEO with dependency on YAC, you will first need to install YAXT and YAC manually.
(YAXT is a dependency of YAC.)

.. note::
  The installation of YAC for CLEO is currently in active development and so may not be exactly as written here.

How to Install YAC (and YAXT)
#############################

The easiest way to install YAXT and YAC is to run the ``install_yac.sh`` bash script found in
``scripts/bash/``. Note you will need to provide the path to the directory where you want
to put the installations.

Alternatively you can download `YAXT <https://swprojects.dkrz.de/redmine/>`_ and
`YAC <https://gitlab.dkrz.de/dkrz-sw/yac/>`_ as compressed files and then configure and compile
them yourself.

yaml-cpp
--------
CLEO's ``initialise`` library depends on the ```yaml-cpp``` package to read and write YAML files. You
can find more information about it from `its repository: <https://github.com/jbeder/yaml-cpp>`_.

The yaml-cpp library for CLEO is automatically built using CMAKE and compiled if required.

CVODE
-----
CLEO's ``coupldyn_cvode`` library requires the SUNDIALS CVODE package. You can find more information
about it from `its webpage: <https://computing.llnl.gov/projects/sundials/cvode>`_.

The CVODE libraries for CLEO are automatically built using CMAKE and compiled if required.
