SDMMonitor Observers Module
=============================

This submodule of ``observers`` contains observers that read out data from an ``SDMMonitor``
(see :doc:`../../superdrops/sdmmonitor`) — a side-channel that microphysical processes can write
diagnostic data into during a timestep (e.g. the condensation/evaporation rate, or precipitation
flux), separate from a gridbox's own state. These observers write that monitored data to a
dataset.


Contents:
---------

.. toctree::
   :maxdepth: 1

   do_sdmmonitor_obs
   monitor_condensation_observer
   monitor_massmoments_change
   monitor_massmoments_change_observer
   monitor_precipitation_observer
