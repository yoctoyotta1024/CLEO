FromFile Cartesian Dynamics
=============================

Header file: ``<libs/coupldyn_fromfile/fromfile_cartesian_dynamics.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/coupldyn_fromfile/fromfile_cartesian_dynamics.hpp>`_

A ``CoupledDynamics`` (see :doc:`../runcleo/coupleddynamics`) reading thermodynamics and wind
fields for a Cartesian domain from binary files (see
:doc:`../../cleopy/thermobinary_src/thermobinary_src`), one timestep at a time.

.. note::
   This header also defines an internal helper struct called ``CartesianDynamics`` that holds
   the thermodynamics/wind arrays read from file (constructed internally by
   ``FromFileDynamics`` below). :doc:`../coupldyn_yac/yac_cartesian_dynamics` defines an
   unrelated struct of the same name for the YAC coupling option. Because Sphinx's C++ domain
   cannot distinguish two identically-named, non-namespaced structs, we have given neither struct
   a documented API entry; please see the source code itself for the ``CartesianDynamics`` API for
   the struct in ``fromfile_cartesian_dynamics.hpp``.

.. doxygenstruct:: FromFileDynamics
   :project: coupldyn_fromfile
   :private-members:
   :protected-members:
   :members:
   :undoc-members:
