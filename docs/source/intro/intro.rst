.. CLEO Introduction documentation master file, created by
   sphinx-quickstart on Mon Nov 20 12:27:54 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Programming Guide
=================

.. note::
   Please consider that this project is under active development
   and our documentation is a work in progress.

Motivation
----------
It seems apparent that a new implementation of SDM is required; capable of modelling warm rain
in LES with realistic boundary conditions and large scale forcings, and capable of application
in large regional domains, with horizontal extents O(100km). CLEO is an attempt to build such
a SDM. It strives to be a library for SDM to model warm clouds with exceptional computational
performance.

Memory Layout
-------------
The fundamental basis for computational performance in CLEO is through efficient memory access
patterns. Primarily this is acheived by ensuring super-droplets which occupy the same gridbox are
always located contiguously in memory. We also try to avoid new memory allocation and cache misses
through our organisation of gridboxes and super-droplets and we use simplistic microphysics for
low cost at run-time.

Timestepping
------------
CLEO's monoidal structures (see below) are designed to allow adaptive-timestepping,
meaning different microphysical processes and observers may have arbitrary time-steps which bare
no relation to one another and can in general change during run-time. This flexibility is contained
in a sub-timestepping routine which is part of CLEO's larger timestepping routine to run SDM
coupled to dynamics.

Coupling to Dynamics
--------------------
Whilst CLEO handles the transport of super-droplets throughout the domain, it cannot perform
advection of Gridboxes' dyanmic variables itself (temperature, pressure, winds etc.). Instead,
CLEO can be one-way or two-way coupled to a dyanmical core capable of advection.
Nevertheless, such a dynamical core is not a necessity. There are many ways for CLEO to receive
dynamics, which may even just involve reading data from binary files.

Monoids
-------
A key novel feature of CLEO is the construction of monoids. We use C++20 concepts to constrain
templated types for microphysics and observers. This ensures they satisfy monoid set properties
and can be combined in well-defined ways. The purpose is to allow for several microphysics
processes (and likewise observers) to be combined simply and flexibly whilst ensuring
adaptive-timestepping and avoiding the use of conditional branches in the code. This enables
extra-ordinary model flexibility without additional run-time cost. It also helps with maintaining
readable and modifyable code.

Kokkos Thread Parallelism
-------------------------
For performance portable thread parallelism we embrace Kokkos. As a consequence,
Kokkos' macros and functions are littered throughout our code and many of our key data structures,
for example Gridboxes and super-droplets, are contained within Kokkos Views. For those seeking
advanced understanding, we defer to Kokkos' GitHub repositories and documentation therein.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   background
   motivation
   memorylayout
   timestepping
   coupling
   monoids
   kokkos

Questions?
----------
Yes please! Simply :ref:`contact us! <contact>`
