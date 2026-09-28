.. _dynamics_coupling:

Coupling to Dynamics
====================

This page explains the concept of Cleo's coupling to a dynamical model, as explained in more detail
in :cite:`bayley2026a` (Sect. 3.3). See the :doc:`coupled dynamics deep dive<../deep_dive/coupled_dynamics>`
to understand more are the code level (the actual C++ concepts and options).

Cleo for SDM Microphyiscs, Dynamical Model for Dynamics
-------------------------------------------------------

Cleo is a stand-alone Super-Droplet Model (SDM): it models cloud microphysics and moves super-droplets
through a domain based on the thermodynamics of each grid-box. It is *not* a fully-fledged climate
model — it has no representation of how pressure, temperature, moisture and
wind fields evolve under advection, nor any parameterisation of radiation or subgrid-scale turbulence.
For those, Cleo is designed to run concurrently with a separate piece of software, some kind of
"Dynamics Solver" which is responsible for the fluid dynamics and any other physical processes
beyond SDM microphysics. This Dynamics solver, computationally defined by the ``CoupledDynamics``
concept (see below), is essentially a black box; it may contain a dynamical core and all sorts of
complicated parameterisations, or it may simply be a set of pre-computed fields
read from a file. Cleo's own job is confined to using whatever thermodynamic and wind fields it gets
from the Dynamics Solver in order to enact SDM microphysics and moves super-droplets around the domain
(recall superdroplet transport is Cleo's responsibility even though field advection is not),
and if two-way coupled, return updated thermodynamics to the Dynamics Solver.

Because the two run concurrently as separate programs (in general with their own, independent
grids, MPI domain decompositions, and no shared memory — more on this below), something has to move
data between Cleo and a certain Dynamics Solver at each coupling timestep. This is the
"Dynamics Coupler", defined computationally by the ``CouplingComms`` in Cleo (see below).
The coupler is also responsible for controlling if Cleo is one-way or two-way coupled to a
certain Dynamics Solver.

One-Way or Two-Way Coupling
---------------------------

The Dynamics Coupler can move information in one or both directions:

* **One-way coupling**: Cleo *receives* the state of each gridbox (wind velocity and
  thermodynamics — pressure, temperature, vapour and liquid mixing ratios) from the
  Dynamics Solver via the Coupler, but nothing flows back. The simplest example is what Cleo's
  ``coupldyn_fromfile`` option implements: thermodynamics for a given timestep are read from
  arrays stored in a file, standing in for a Dynamics Solver that Cleo has no ability to
  influence. This mode is sometimes called *piggybacking*: Cleo rides along on a pre-existing
  dynamical solution, useful when you want realistic-looking fields without a two-way
  feedback, or when no live dynamical core is required.
* **Two-way coupling**: Cleo not only receives, but also *sends* the gridbox state back to the
  Dynamics Solver, so that the effect of microphysics (for example, the latent heating from
  condensation or moisture change from precipitation) can feed back into the evolving dynamical
  fields. This is the physically more consistent since it allows microphysics to shape the
  circulation and thermodynamics, rather than merely respond to them.

Given the state of a gridbox, Cleo can enact microphysical processes on the super-droplets within
it, and it can move super-droplets according to the (possibly newly-received) wind field —
updating their coordinates and which gridbox they belong to, and re-ordering the super-droplet
arrays accordingly (see :doc:`memorylayout`). The physics (numerical methods) for exactly how
condensation, collisions and motion are computed from a gridbox's state is the subject of Cleo's
second model description paper :cite:`bayley2026b`.

Why the Coupler?
----------------

* Because Cleo and the Dynamics Solver are, in general, separate MPI programs with no shared
  memory, they need not share the same grid, nor MPI domain decomposition** — i.e. how the domain is
  split geometrically and across compute node. Coupling to a certain Dynamics Solver with `YAC
  <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_ (Yet Another Coupler) specifically supports this:
  YAC handles the MPI communication and interpolates variables between two different grids.
* Allowing different grids and is a deliberate design choice so that the Dynamics Solver is free to
  lay out its grid however suits its fluid dynamics, while Cleo is free to lay out gridboxes however
  suits SDM — for example a nested grid, or gridbox boundaries chosen to simplify or reduce the
  complexity of super-droplet motion/microphysics calculations.
* The trade-off is communication cost: with no shared memory, exchanging state between two
  independent domain decompositions is costlier than it would be if Cleo and the Dynamics Solver
  were tightly integrated into one program with one decomposition. This is however, worthwhile,
  because it maximises the freedom to load-balance and allocate computational
  resources independently for the two very different workloads (SDM's is superdroplet-based;
  a dynamical core's is grid-based) — see :doc:`../deep_dive/code_structure` for how this
  independence is reflected in Cleo's own MPI domain decomposition.

From Theory to Code
-------------------

Cleo's software mirrors this design directly. The Dynamics Solver is what Cleo's code calls
a type satisfying the ``CoupledDynamics`` C++ concept; the Dynamics Coupler is what Cleo's code
calls a type satisfying the ``CouplingComms`` concept. One-way coupling is simply the case where
a ``CouplingComms`` type's "send" half does nothing; two-way coupling is when both directions do real
work. Cleo currently ships four concrete ``CoupledDynamics`` options (no dynamics at all,
reading from file, a 0-D CVODE solver, and YAC coupling to an external dynamical core such as ICON)
plus a fifth route via Python/numpy bindings. For the concept definitions themselves, a table of
all the options and when to use each, and a worked, annotated code example switching between two of
them, see :doc:`../deep_dive/coupled_dynamics`.

----

.. note::
   This page was expanded by Claude (Anthropic), at the request of Clara Bayley, based on
   Sect. 3.3 ("Coupling Cleo to a Host Dynamical Driver") of :cite:`bayley2026a`.
   Please be critical of its accuracy and open an issue or get in touch if you spot anything wrong
   or out of date.
