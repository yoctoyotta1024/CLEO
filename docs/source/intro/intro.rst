Programming Guide
=================

This programming guide gives a short summary of what is explained more thoroughly in Cleo's model
description papers :cite:`bayley2026a, bayley2026b`, and mainly focusses on Cleo's fundamental
design and code structure rather than its microphysics capabilities.

Motivation
----------
It seems apparent that a new implementation of SDM is required; capable of modelling warm rain
in LES with realistic boundary conditions and large scale forcings, and capable of application
in large regional domains, with horizontal extents O(100km). Cleo is an attempt to build such
a SDM. It strives to be a library for SDM to model warm clouds with exceptional computational
performance.

Read more about Cleo's :ref:`background <background>` and :ref:`motivation <motivation>`.

Memory Layout
-------------
The fundamental basis for computational performance in Cleo is through efficient memory access
patterns. Primarily this is acheived by ensuring super-droplets which occupy the same gridbox are
always located contiguously in memory. We also try to avoid new memory allocation and cache misses
through our organisation of gridboxes and super-droplets and we use simplistic microphysics for
low cost at run-time.

Read more about Cleo's :ref:`memory layout <memlayout>`.

Timestepping
------------
Cleo's monoidal structures (see below) are designed to allow adaptive-timestepping,
meaning different microphysical processes and observers may have arbitrary time-steps which bare
no relation to one another and can in general change during run-time. This flexibility is contained
in the logic for how different microphysical processes and observers are combined. Cleo's overall
timestepping routine therefore calls a single (combined) observer and a single (combined)
microphysical process. Overall, the outermost timestepping routine calls the observer, the dynamics
Cleo is coupled to, and SDM. Within the call to SDM there is a sub-timestepping routine which
calls superdroplet motion and the microphysical process.

Read more about Cleo's :ref:`timestepping <timestepping>`.

Coupling to Dynamics
--------------------
Cleo is pretty much agnostic to the thermodynamics model it is coupled to. The coupled thermodynamics
can be thought of as a "black box" which provides thermodynamics data to Cleo at each of its coupling
time-steps. This "black box" can range from simply a structure which reads data from binary
files, all the way to a fully-fledged dynamical core. The bridge between Cleo and any specific
thermodynamics model (a specific "black box"), is made by a certain coupler. This coupler controls
if the coupling is one-way, meaning Cleo receives but doesn't send back data, or two-way, meaning Cleo both
receives data from and sends data to the coupled thermodynamics model. Note that whilst Cleo
handles the transport of super-droplets throughout the domain, it cannot perform
advection of Gridboxes' thermodynamical variables itself (temperature, pressure, winds etc.). For
thermodynamic advection, Cleo must be coupled to a thermodynamics model capable of advection.
Read more about Cleo's :ref:`coupling to dynamics <dynamics_coupling>`.

Monoids
-------
A key novel feature of Cleo is the construction of monoids. We use C++20 concepts to constrain
templated types for microphysics and observers. This ensures they satisfy monoid set properties
and can be combined in well-defined ways. The purpose is to allow for several microphysics
processes (and likewise observers) to be combined simply and flexibly whilst ensuring
adaptive-timestepping and avoiding the use of conditional branches in the code. This enables
extra-ordinary model flexibility without additional run-time cost. It also helps with maintaining
readable and modifyable/extendable code.

Read more about Cleo's :ref:`monoids <monoids>`.

Kokkos Thread Parallelism
-------------------------
For performance portable thread parallelism we embrace Kokkos. As a consequence,
Kokkos' macros and functions are littered throughout our code and many of our key data structures,
for example Gridboxes and super-droplets, are contained within Kokkos Views. For those seeking
advanced understanding, we defer to Kokkos' GitHub repositories and documentation therein.

Read more about Cleo's :ref:`thread parallelism <kokkos>`.

In More Detail:
---------------
.. toctree::
   :maxdepth: 1

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
