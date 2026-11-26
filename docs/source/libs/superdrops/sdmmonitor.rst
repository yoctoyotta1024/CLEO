SDMMonitor
============

Header file: ``<libs/superdrops/sdmmonitor.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/superdrops/sdmmonitor.hpp>`_

The C++ concept for a monitor that a microphysical process can write diagnostic data into
during a timestep (e.g. condensation rate, precipitation), independently of a gridbox's own
state. See :doc:`../observers/sdmmonitor/sdmmonitormod` for the observers that read a monitor
out to a dataset.

.. doxygenconcept:: SDMMonitor
   :project: superdrops

.. doxygenstruct:: CombinedSDMMonitor
   :project: superdrops
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. doxygenstruct:: NullSDMMonitor
   :project: superdrops
   :private-members:
   :protected-members:
   :members:
   :undoc-members:
