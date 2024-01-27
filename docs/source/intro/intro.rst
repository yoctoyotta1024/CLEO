.. CLEO Introduction documentation master file, created by
   sphinx-quickstart on Mon Nov 20 12:27:54 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Programming Guide
=================

.. note::
   Please consider that this project is under active development
   and our documentation is a work in progress.

It seems apparent that a new implementation of SDM is required;
capable of modelling warm rain in LES with realistic boundary
conditions and large scale forcings, and capable of application
in large regional domains, with horizontal extents O(100km). CLEO is
an attempt to build such a SDM. It strives to be a library for SDM
to model warm clouds with exceptional computational performance.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   background
   motivation
   memorylayout
   monoids
   coupling
   timestepping
   kokkos

Questions?
----------
Yes please! Simply :ref:`contact us! <contact>`
