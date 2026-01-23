Create Gridboxes Sub-Module
===========================

Header file: ``<libs/runcleo/creategbxs.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/runcleo/creategbxs.hpp>`_

This sub-module defines functions and classes called by RunCLEO to create and
initialise Gridboxes (i.e. GBxs).

.. doxygenclass:: GenGridbox
   :project: runcleo
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. doxygenfunction:: create_gbxs
   :project: runcleo

.. doxygenfunction:: initialise_gbxs
   :project: runcleo

.. doxygenfunction:: initialise_gbxs_on_host
   :project: runcleo

.. doxygenfunction:: is_gbxinit_complete
   :project: runcleo

.. doxygenfunction:: print_gbxs
   :project: runcleo
