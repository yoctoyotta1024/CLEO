YAC Cartesian Dynamics
========================

Header file: ``<libs/coupldyn_yac/yac_cartesian_dynamics.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/coupldyn_yac/yac_cartesian_dynamics.hpp>`_

A ``CoupledDynamics`` (see :doc:`../runcleo/coupleddynamics`) for a Cartesian domain whose
thermodynamic and wind fields are received from an external dynamical core via YAC.

.. note::
   This header also defines its own internal helper struct called ``CartesianDynamics`` that
   holds the thermodynamics/wind arrays received via YAC (constructed internally by
   ``YacCartesianDynamics`` below) — unrelated to, but identically named to, the internal struct
   in :doc:`../coupldyn_fromfile/fromfile_cartesian_dynamics`. Because Sphinx's C++ domain cannot
   distinguish two identically-named, non-namespaced structs, it is described here in prose only
   (see that page if you want the generated API entry for the other one).

.. doxygenstruct:: YacCartesianDynamics
   :project: coupldyn_yac
   :private-members:
   :protected-members:
   :members:
   :undoc-members:
