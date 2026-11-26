CartesianDomain Module
=======================

This module of Cleo contains the implementation of a Cartesian (0-D, 1-D, 2-D or 3-D) domain:
its ``GridboxMaps`` (see :doc:`../gridboxes/gridboxmaps`), its MPI domain decomposition, and
the super-droplet motion and transport-across-domain logic specific to a Cartesian grid.


Contents:
---------

.. toctree::
   :maxdepth: 1

   movement/movementmod
   cartesian_decomposition
   cartesianmaps
   createcartesianmaps
   domainboundaries
   doubly_periodic_domain
   collect_data_for_collective_dataset
