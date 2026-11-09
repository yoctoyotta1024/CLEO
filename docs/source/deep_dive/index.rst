.. _deep_dive:

Cleo Deep Dive
================================

This section is aimed at climate researchers or software engineers who are comfortable with
Python and scientific computing, but have little or no C++ experience and have not yet worked inside
Cleo's source code. It assumes you have already followed
:doc:`installation <../usage/installation/installation>` and at least run one of
:doc:`the examples <../usage/examples/examples>`.

What this section covers (and what it doesn't)
------------------------------------------------

The rest of the documentation tells you *what Cleo is* (:doc:`../intro/intro`), *how to install
it* (:doc:`../usage/installation/installation`), and *how to run the examples that already
exist* (:doc:`../usage/examples/examples`). It does not yet explain how the C++ source code is
organised, how to read or modify ``main.cpp``, or how to build a brand new example from scratch.
That is the gap this section fills.

Concretely, this section explains:

1. (This page) Scope and a glossary of terms used throughout.
2. :doc:`code_structure` — a map of the repository and the C++ ideas you'll meet in it.
3. :doc:`microphysics_and_observers` — how microphysics and observers are combined in code.
4. :doc:`coupled_dynamics` — what "coupled dynamics" means and the available options.
5. :doc:`main_cpp` — a walkthrough of ``main.cpp``/``main_impl.hpp`` and how to adapt them.
6. :doc:`tutorial_new_example` — a start-to-finish tutorial building a new example.
7. :doc:`external_project` — fetching and building Cleo as a dependency of another project.

It deliberately does not re-explain the physics of warm-cloud microphysics in depth — for that,
see the two papers describing Cleo:

* Bayley et al. (2026), *Cleo: the fundamental design of a new superdroplet model for high
  computational performance (v0.39.0)*, Geosci. Model Dev., 19, 6099-6120 — referred to
  throughout as **the design paper**.
* Bayley et al. (2026), *Cleo: the numerical methods of a new superdroplet model including a
  droplet breakup algorithm (v0.52.0)*, Geosci. Model Dev., 19, 6121-6142 — referred to
  throughout as **the numerics paper**.

Nor does it duplicate CLEO's Doxygen/Breathe-generated :doc:`API reference <../libs/libs>` —
instead it cross-references specific classes there so you can look up full signatures.

.. toctree::
   :maxdepth: 1
   :caption: In This Section:

   code_structure
   microphysics_and_observers
   coupled_dynamics
   main_cpp
   tutorial_new_example
   external_project

.. _glossary:

Glossary
--------

These terms are used freely from here on. Cloud-microphysics terms first, then Cleo-specific
software terms.

Cloud microphysics
~~~~~~~~~~~~~~~~~~~

.. glossary::

   Super-droplet
      A Lagrangian particle that statistically represents a large number of real, physical
      droplets ("real droplets") sharing similar properties. See its :term:`Multiplicity`.

   Multiplicity
      Denoted :math:`\xi`, the number of real droplets a single super-droplet represents.
      Multiplicity lets one computational particle stand in for a whole population, which is
      what makes SDM computationally tractable.

   Super-Droplet Method (SDM)
      The Lagrangian, Monte-Carlo method Cleo implements for warm-cloud microphysics, in which
      condensation/evaporation is applied deterministically to each super-droplet and collisions
      between super-droplets are resolved probabilistically. See :doc:`../intro/background`.

   Condensation / evaporation
      The growth or shrinkage of a droplet's water mass by phase change with the surrounding
      water vapour. Modelled deterministically for every super-droplet.

   Collision-coalescence
      Two droplets colliding and merging into one. In SDM, whether a collision is enacted is
      decided probabilistically for a randomly-sampled subset of super-droplet pairs each
      timestep (see the numerics paper, Sect. 3).

   Warm rain
      Precipitation formed entirely through liquid-phase processes (condensation and
      collision-coalescence), without any ice-phase microphysics.

   Rebound / breakup
      Alternative outcomes of a droplet collision to coalescence: the two droplets bounce apart
      unchanged (rebound), or merge briefly and then fragment into several smaller droplets
      (breakup). Cleo can model both; see the numerics paper, Sect. 4, and
      :doc:`microphysics_and_observers`.

   Gridbox
      One cell of the spatial domain over which thermodynamic quantities (pressure, temperature,
      vapour/liquid mixing ratio, wind) are defined, and which owns a subset of the
      super-droplets currently located within its boundaries.

Cleo-specific software terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. glossary::

   Coupled dynamics
      The external (or internal) source of a gridbox's thermodynamic and wind fields at every
      coupling timestep. Cleo itself only models microphysics and super-droplet transport; it
      always needs *something* to supply/advect the dynamical fields. See
      :doc:`coupled_dynamics`.

   Coupling
      The exchange of information between Cleo's Super-Droplet Model and the coupled dynamics:
      "receiving" dynamics fields into gridboxes, and (for two-way coupling) "sending" the
      SDM's effect on those fields back out.

   Observer
      A Cleo software abstraction for anything that reads out data during a run (writing to a
      Zarr dataset, printing to the screen, feeding a monitor). See
      :doc:`microphysics_and_observers`.

   Microphysical process
      A Cleo software abstraction for anything that changes super-droplets' attributes
      during a timestep (e.g. condensation, collision-coalescence). See
      :doc:`microphysics_and_observers`.

   Monoid
      The C++ pattern Cleo uses to let multiple microphysical processes (or multiple observers)
      be combined, in any combination, into one object of the same kind — see
      :doc:`../intro/monoids` for the conceptual introduction and
      :doc:`microphysics_and_observers` for worked examples.

   Concept (C++20)
      A named, compiler-checked set of requirements a type must satisfy (which member functions
      it has, and their signatures). Cleo uses concepts (``GridboxMaps``, ``MicrophysicalProcess``,
      ``Observer``, ``CoupledDynamics``, ``Motion``, ...) instead of e.g. Python's duck typing or
      abstract base classes, to select and combine behaviour entirely at compile time. See
      :doc:`code_structure`.

   GridboxMaps
      The C++ concept (and its concrete implementation, ``CartesianMaps``) providing the
      geometry of the domain: converting between a gridbox's index and its physical coordinate
      bounds, and finding neighbouring gridboxes. See :doc:`code_structure`.

   SDMMethods
      The class that owns one complete, ready-to-run instance of Cleo's SDM: its ``GridboxMaps``,
      its combined microphysics, its super-droplet motion, and its combined observer. See
      :doc:`main_cpp`.

   RunCLEO
      The class that owns the outer timestepping loop coupling an ``SDMMethods`` instance to a
      ``CoupledDynamics`` instance via a ``CouplingComms`` instance. See :doc:`main_cpp`.

   cleopy
      Cleo's Python package (``cleopy``), used for generating YAML configuration files and
      binary input files for a run, and for reading Cleo's Zarr output back into ``xarray``.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
