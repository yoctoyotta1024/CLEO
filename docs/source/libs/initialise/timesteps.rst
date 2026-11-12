Timesteps
===========

Header file: ``<libs/initialise/timesteps.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/initialise/timesteps.hpp>`_

Converts between real (physical) time and Cleo's dimensionless, integer model timestep, and
stores the various timestep intervals (e.g. coupling, observer, motion) a run configures.

.. doxygenclass:: Timesteps
   :project: initialise
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. doxygenfunction:: realtime2step
   :project: initialise

.. doxygenfunction:: realtime2dimless
   :project: initialise

.. doxygenfunction:: step2realtime
   :project: initialise

.. doxygenfunction:: step2dimlesstime
   :project: initialise
