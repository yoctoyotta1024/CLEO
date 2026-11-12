Create Super-droplets Sub-Module
================================

Header file: ``<libs/runcleo/createsupers.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/runcleo/createsupers.hpp>`_

This sub-module defines functions and classes called by RunCLEO to create and initialise
super-droplets in Kokkos view (on device). Generating a single super-droplet's initial
attributes is delegated to :doc:`gensuperdrop`.

.. doxygenfunction:: create_supers
   :project: runcleo

.. doxygenfunction:: initialise_supers
   :project: runcleo

.. doxygenfunction:: initialise_supers_on_host
   :project: runcleo

.. doxygenfunction:: is_sdsinit_complete
   :project: runcleo

.. doxygenfunction:: print_supers
   :project: runcleo
