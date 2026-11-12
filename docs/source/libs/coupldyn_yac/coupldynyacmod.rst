Coupldyn YAC Module
=====================

This module of Cleo implements a ``CoupledDynamics``/``CouplingComms`` pair (see
:doc:`../runcleo/coupleddynamics`, :doc:`../runcleo/couplingcomms` and :doc:`../../intro/coupling`)
that couples Cleo to an external dynamical core (e.g. ICON) via the
`YAC <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_ coupler, including interpolation between
different grids and MPI communication between independent domain decompositions.


Contents:
---------

.. toctree::
   :maxdepth: 1

   yac_cartesian_dynamics
   yac_comms
