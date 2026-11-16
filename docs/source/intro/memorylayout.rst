.. _memlayout:

Memory Layout
=============

This page describes Cleo's memory layout, which is also explained more thoroughly in
:cite:`bayley2026a` (Sect. 2).

The key structures which compose Cleo and their memory layout on each MPI process are
shown schematically here:

.. figure:: ../figs/cleomemlayout.png
   :alt: Schematic for Cleo's memory layout
   :align: center

   A schematic for the memory layout of Cleo's key structures.

Each MPI process simulates a sub-portion of the total domain composed of superdroplets and gridboxes.
The gridboxes exist in an ordered, contigous block of memory (managed by a Kokkos view-like type)
and each gridbox contains a "view" to a seperate contigous block of memory of superdroplets, also
ordered by gridbox-index. The gridboxes exchange information about their dynamical state
(thermodynamics and wind velocity) via a "Coupler" to/from a coupled "host" dynamical model. The
coupled dynamics and coupler are abstract ideas (i.e. they are C++ concepts for constraining
templated types; see :doc:`coupling`), but get instantiated at compilation, for example as
`ICON <https://www.icon-model.org>`_ and `YAC <https://dkrz-sw.gitlab-pages.dkrz.de/yac/index.html>`_ respectively.

Each gridbox is composed of certain elements:

.. figure:: ../figs/gridbox.png
   :alt: Schematic showing the elements of a gridbox in Cleo
   :align: center

   Schematic showing the elements of a gridbox in Cleo

For reasons of computational efficiency, superdroplets are not stored within gridboxes. Instead all
the superdroplets whose positions lie within the portion of the domain on that node at a given timestep
are stored in another distinct contigous block of memory (also managed by a
Kokkos view-like type) and that block is ordered by the gridbox each superdroplet is associated with
(i.e. the superdroplets are ordered based on their location in the domain). Each gridbox then
contains a pointer-like object which identifies another contigous block of memory,
seperate from the gridboxes, that contains the superdroplets associated with that particular gridbox
at a given timestep.

Every superdroplet has spatial coordinates and an index which identifies which gridbox those
coordinates lie within. Furthermore the superdroplet has attributes which defines what "real"
droplet it represents (e.g. the superdroplet may have the attributes to define its radius and
multiplicity). Additionally, to aid analysis involving Lagrangian particle tracking, each superdroplet
can optionally be given a unique identifier. This is illustrated in the schematic below.

.. figure:: ../figs/superdroplet.png
   :alt: Schematic showing the elements of a superdroplet in Cleo
   :align: center

   Schematic showing the elements of a superdroplet in Cleo
