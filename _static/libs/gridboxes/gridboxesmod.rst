Gridboxes Module
==================

This module of Cleo contains the domain-independent implementation of a gridbox (the
``GridboxMaps`` concept, the ``Gridbox`` struct itself, and everything to do with finding,
sorting and moving the super-droplets within and between gridboxes). The concrete, Cartesian
implementation of ``GridboxMaps`` lives in :doc:`../cartesiandomain/cartesiandomainmod`.


Contents:
---------

.. toctree::
   :maxdepth: 1

   gridboxmaps
   gridbox
   gbxindex
   supersingbx
   supersindomain
   findrefs
   sortsupers
   movesupersindomain
   predcorr
   predcorrmotion
   boundary_conditions
   transport_across_domain
   cfl_criteria
