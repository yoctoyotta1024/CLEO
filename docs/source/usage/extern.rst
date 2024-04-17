External Libraries
==================

CLEO depends upon Kokkos and may depend upon some additional external libraries such as YAC and
yaml-cpp depending on your setup. These are automatically built using CMAKE and compiled if
required.

.. note::
  The installation of YAC for CLEO is currently in development and may rquire some manual installation.

Kokkos
------
All builds of CLEO require Kokkos in order to implement thread parallelism. You can read more about
how we use Kokkos on :doc:`our page about Kokkos<intro/kokkos>`.

YAC
---
YAC is required if CLEO couples to dynamics using YAC and/or uses MPI domain decompoisiton. You can
find more information about it from `its documentation: <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_.

TODO(all): Detail how to install YAC for CLEO.

yaml-cpp
--------
The ```initialise``` library depends on the ```yaml-cpp``` package to read and write YAML files. You
can find more information about it from `its repository: <https://github.com/jbeder/yaml-cpp>`_.

CVODE
-----
The ```coupldyn_cvode``` library requires the SUNDIALS CVODE package. You can find more information
about it from `its webpage: <https://computing.llnl.gov/projects/sundials/cvode>`_.
