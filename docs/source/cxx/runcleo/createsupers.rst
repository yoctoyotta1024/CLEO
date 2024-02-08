CreateSupers Sub-Module
=======================

Header file: ``<libs/runcleo/createsupers.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/runcleo/createsupers.hpp>`_

.. cpp:module:: createsupers
   :project: runcleo

   This module defines functions and classes called by RunCLEO to create and initialise
   super-droplets in Kokkos view (on device).

.. doxygenclass:: GenSuperdrop
   :project: runcleo
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. cpp:function:: create_supers
   :project: runcleo

.. cpp:function:: initialise_supers
   :project: runcleo

.. cpp:function:: initialise_supers_on_host
   :project: runcleo

.. cpp:function:: is_sdsinit_complete
   :project: runcleo

.. cpp:function:: print_supers
   :project: runcleo
