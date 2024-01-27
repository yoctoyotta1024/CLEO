.. CLEO Introduction documentation master file, created by
   sphinx-quickstart on Mon Nov 20 12:27:54 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Programming Guide
=================

.. note::
   Please consider that this project is under active development
   and our documentation is a work in progress.

It seems apparent that a new implementation of SDM is required; capable of modelling warm rain
in LES with realistic boundary conditions and large scale forcings, and capable of application
in large regional domains, with horizontal extents O(100km). CLEO is an attempt to build such
a SDM. It strives to be a library for SDM to model warm clouds with exceptional computational
performance.

The fundamental basis for computational performance in CLEO is through efficient memory access
patterns. Primarily this is acheived by ensuring super-droplets which occupy the same gridbox are
always located contiguously in memory. We also try to avoid new memory allocation and we use
parallelised loops over gridboxes and super-droplets and simplistic microphysics for low cost
at run-time.

A key novel feature of CLEO is the construction of monoids. We use C++20 concepts to constrain
templated types for microphysics and observers. This ensures they satisfy monoid set properties
and can be combined in well-defined ways. The purpose is to allow for several microphysics
processes (and observers) to be combined easily and with adaptive-timestepping whilst avoiding
the use of conditional branches in the code. This enables massive model flexibility without
additional run-time cost. It also helps with maintaining readable and modifyable code.

For portable performance portable thread parallelism we embrace Kokkos. As a consequence,
Kokkos' macros and functions are littered throughout our code and many of our key data structures,
for example Gridboxes and super-droplets, are contained within Kokkos Views. For those seeking
advanced understanding, we refer to `Kokkos' github repositories <https://github.com/kokkos>`_
and documentation therein.

TODO: a bit of schpiel on coupling and dynamics solver

TODO: a bit of schpiel on the final time-stepping of the model

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
