.. _deep_dive:

Cleo Deep Dive
================

This section is aimed at climate researchers or software engineers who are comfortable with
Python and scientific computing, but have little or no C++ experience and have not yet worked inside
Cleo's source code. It assumes you have already followed
:doc:`installation <../usage/installation/installation>` and at least run one of
:doc:`the examples <../usage/examples/examples>`.

What this section covers
-------------------------

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

It deliberately does not re-explain the physics of warm-cloud microphysics nor Cleo's computational
performance in depth — for that, see the two papers describing Cleo:

* **the design paper:** Bayley et al. (2026), *Cleo: the fundamental design of a new superdroplet model for high
  computational performance (v0.39.0)*, Geosci. Model Dev., 19, 6099-6120 — :cite:`bayley2026a`
  referred to throughout as "B2026a" .
* **the numerics paper:** Bayley et al. (2026), *Cleo: the numerical methods of a new superdroplet model including a
  droplet breakup algorithm (v0.52.0)*, Geosci. Model Dev., 19, 6121-6142 — :cite:`bayley2026b`
  referred to throughout as "B2026b" .

C++ Ideas You'll Meet
----------------------

Cleo is written in C++20 and makes heavy use of a small number of language features. None of them
have an exact Python equivalent, but each has a reasonable analogy.

Header (``.hpp``) and source (``.cpp``) files
    In Python, a module's interface and its implementation live in the same ``.py`` file. In C++
    they are conventionally split in two: a **header file** (``.hpp``) *declares* what a class or
    function looks like (its name, its parameters, its return type) so that other files can use it
    without seeing how it works; a **source file** (``.cpp``) *defines* how it actually works.
    Anywhere Cleo does a substantial piece of work you'll find a ``foo.hpp``/``foo.cpp`` pair
    (e.g. ``superdrops/condensation.hpp`` and ``superdrops/condensation.cpp``). The exception is
    templated code (see below), which is usually written entirely in the header, because the
    compiler needs the full definition available at every place it's used.

Compilation and templates
    Python is interpreted: the same function can be called with an ``int`` or a ``str`` and the
    interpreter figures out what to do at call-time. C++ is compiled ahead of time, and a
    **template** is C++'s way of writing one piece of code that the compiler turns into many
    concrete versions, one per type it's used with, before the program ever runs. You will see
    template code everywhere in Cleo, for example::

        template <MicrophysicalProcess Microphys>
        class SDMMethods { ... };

    reads as "``SDMMethods`` is a class parameterised by some type ``Microphys``, and that type
    must satisfy the ``MicrophysicalProcess`` concept". Nothing here is decided at runtime — when
    you write your ``main.cpp`` and pick, say, ``Condensation`` as your microphysics, the compiler
    generates a version of ``SDMMethods`` specifically for the microphysical process ``Condensation``.
    This is why swapping microphysics or observers costs nothing at runtime (see :doc:`../intro/kokkos`
    and B2026a :cite:`bayley2026a`, Sect. 2): there is no dynamic dispatch, no virtual function calls, no ``if``/``else``
    branching on "which kind of microphysics is this" anywhere in the timestepping loop.

Concepts
    A **concept** names a set of requirements — usually, "has a member function with this
    signature" — that a type must satisfy. It is checked by the compiler, not at runtime. It's the
    closest C++ analogue to Python's duck typing ("if it walks like a duck and quacks like a
    duck..."), except the check happens when you compile, not when you run. For example, from
    ``libs/superdrops/microphysicalprocess.hpp``::

        template <typename P>
        concept MicrophysicalProcess =
            requires(P p, const TeamMember& tm, const unsigned int t, subviewd_supers supers,
                     State& state, const NullSDMMonitor mo) {
              { p.next_step(t) } -> std::convertible_to<unsigned int>;
              { p.on_step(t) } -> std::same_as<bool>;
              { p.run_step(tm, t, supers, state, mo) } -> std::same_as<subviewd_supers>;
            };

    Any type — ``Condensation``, ``CollCoal``, or one you write yourself — that has these three
    member functions with these signatures automatically *is* a ``MicrophysicalProcess``, with no
    inheritance or registration needed. Cleo defines one such concept for each kind of
    "pluggable" component: ``GridboxMaps``, ``MicrophysicalProcess``, ``Observer``,
    ``CoupledDynamics``, ``CouplingComms``, ``Motion``, and more. Wherever you see
    ``SomeConcept auto x = ...;`` in Cleo's code, read it as "``x`` is some type that satisfies
    ``SomeConcept`` — I don't need to know or write out its exact (often very long,
    template-instantiated) type name".

Monoids
    Building on concepts, Cleo defines **monoids**: sets of types that can (a) be created, (b) be
    combined pairwise into another member of the same set, and (c) have a "do nothing" null
    member. Microphysical processes and observers are both monoids in Cleo. This is covered in
    depth, with worked examples, in :doc:`microphysics_and_observers`; the :doc:`conceptual
    introduction to monoids<../intro/monoids>` is found in the programming guide.

Namespaces
    A C++ ``namespace`` groups related names to avoid clashes, similar to how a Python package
    groups modules. ``KokkosCleoSettings`` (aliased ``KCS``), ``dimless_constants`` (``dlc``),
    ``dimmed_constants`` (``DC``), and ``LIMITVALUES`` are namespaces you'll frequently see used in Cleo.

Kokkos views and functions
    Cleo uses the :doc:`kokkos<../intro/kokkos>` libraries for performance-portable thread
    parallelism; a ``Kokkos::View`` is roughly an array (or better said, a pointer to an array)
    that can live on CPU or GPU memory, and ``KOKKOS_INLINE_FUNCTION`` or ``KOKKOS_FUNCTION`` marks
    a function so it can run inside a parallel region on either. You don't need to understand Kokkos
    internals to read ``main.cpp``, but you will see its types (``viewd_gbx``, ``subviewd_supers``,
    ...) passed around — these are Cleo's own aliases defined in ``libs/kokkosaliases.hpp`` and
    ``libs/superdrops/kokkosaliases_sd.hpp``.


.. _glossary:

Glossary
--------

These cloud-microphysics and Cleo-specific terms are used repeatedly from here on and you should
become familiar with them.

Cloud microphysics
~~~~~~~~~~~~~~~~~~~

.. glossary::

   Super-droplet
      (Also "Superdroplet" or "Superdrop".) A Lagrangian particle that statistically represents
      a large number of real, physical droplets ("real droplets") sharing similar properties.

   Multiplicity
      Denoted ":math:`\xi`" or "``xi``", the number of real droplets a single super-droplet represents.
      Multiplicity lets one computational particle stand in for a similar-looking subset of the
      total population of real droplets, which is what makes SDM computationally tractable. The higher
      the multiplicity, the coarser the representation of the real droplet population.

   Super-Droplet Method (SDM)
      (Also "Superdroplet Model".) The Lagrangian, particle-in-cell, Cleo
      implements for warm-cloud microphysics, in which condensation/evaporation is applied
      deterministically to each super-droplet and collisions between super-droplets are resolved
      probabilistically via an "all or nothing"  Monte-Carlo algorithm.
      See :doc:`../intro/background`.

   Condensation / evaporation
      The growth or shrinkage of a droplet's water mass by phase change with the surrounding
      water vapour. Modelled deterministically for every super-droplet.

   Collision-coalescence
      Two droplets colliding and merging into one. In SDM, whether a collision is enacted is
      decided probabilistically for a randomly-sampled subset of super-droplet pairs each
      timestep (see B2026b :cite:`bayley2026b`, Sect. 3).

   Warm rain / clouds
      Clouds or precipitation formed entirely through liquid-phase processes (condensation and
      collision-coalescence), i.e. without any ice-phase microphysics. Warm clouds are often
      described interchangeably as "shallow clouds" because they are thin and low-lying; they rise
      only up-to around the top of the boundary layer, where the inversion seperating the
      free toposphere above acts like a cap on them. This is in contrast to deep convective clouds,
      which can reach from within the boundary layer all the way up to the tropopause and, above the
      freezing level, contain ice particles (including not just ice crystals, but graupel and
      snow too).

   Coalescence / rebound / breakup
      Alternative outcomes of a droplet collision. In coalescence two real droplets merge into one
      (the number of superdroplets remains unchanged, but the multiplicity and mass of the
      superdroplets change); in rebound the two droplets bounce apart unchanged; in breakup
      real droplets may merge briefly but then fragment into several smaller droplets (again, the
      number of superdroplets remains unchanged in Cleo, but the multiplicity and mass of the
      superdroplets change). Cleo can model all three processes for warm-clouds, see B2026b
      :cite:`bayley2026b`, Sect. 4, and :doc:`microphysics_and_observers`. Collision-coalescence
      is usually assumed rather than model breakup and rebound too because the unceratinties
      surrounding our knowledge of droplet breakup are very large and rebound is either
      neglected or accounted for in the collision efficiency parameterisation whem modelling
      collision-coalescence.

Cleo-specific software terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. glossary::

   Gridbox
      One cell of the spatial domain over which thermodynamic quantities (pressure, temperature,
      vapour/liquid mixing ratio, wind) are defined, and which owns a subset of the
      super-droplets currently located within its boundaries. The gridbox size defines the volume
      and the superdroplets over which collision-coalescence is enacted probabilistically, and
      over which droplet motion is resolved; as explained in B2026b :cite:`bayley2026b`.

   Coupled dynamics
      A Cleo software abstraction for the "thing" Cleo is coupled to. The coupled dynamics "thing"
      is the source of a gridbox's thermodynamic and wind fields at every coupling timestep.
      Cleo itself only models microphysics and super-droplet transport; it always needs *something*
      to supply/advect the dynamical fields. See :doc:`coupled_dynamics`. The definition in code
      of the coupled dynamics is the ``CoupledDynamics`` concept
      (``libs/runcleo/coupleddynamics.hpp``).

   Coupling
      A Cleo software abstraction for the the communication of information between Cleo's
      Super-Droplet Model and the coupled dynamics. The coupling is responsible for "receiving"
      dynamics fields into gridboxes, and, for a  two-way coupling, "sending" Cleo's effect on
      those fields back to the coupled dynamics. The definition in code of the coupling is
      the ``CouplingComms`` concept (``libs/runcleo/couplingcomms.hpp``).

   Observer and Monitor
      An observer is a Cleo software abstraction for anything that reads out data during a run
      (writing to a Zarr dataset, printing to the screen / terminal etc.). A monitor is a
      sub-component of an observer which is specifically for outputing data from within Cleo's SDM
      sub-timestepping routines (instead of at the start/end of an observer's timestep).
      See :doc:`microphysics_and_observers`.

   Microphysical process
      A Cleo software abstraction for anything that changes super-droplets' attributes
      during a timestep (e.g. condensation, collision-coalescence). See
      :doc:`microphysics_and_observers`.

   Monoid
      An idea taken from set-theory in mathematics which we've made computationally in Cleo
      using C++ concepts and templating. In mathematics, a monoid is a group (or set) of members
      which produce another member of the group when combined together in a certain way,
      e.g. positive integers when added or multiplied. The monoid group also contains a member
      which does nothing under that combination (0 for positive integer addition, 1 for
      multiplication). In the same way, in Cleo, we have defined types/objects for microphysical
      processses which act as monoids, meaning each type/object can be combined with other
      microphysical processes in any order to produce another microphysical process which is acts as
      their combination. We have also defined a "Null" microphysical process which does nothing,
      and can be combined with other microphysical processes and leave them unchanged.
      Likewise, we have defined monoids for observers and monitors in Cleo too — see
      :doc:`../intro/monoids` for the conceptual introduction and :doc:`microphysics_and_observers`
      for worked examples.

   Concept (C++20)
      A named set of requirements a type must satisfy during compilation (e.g. which member
      functions a type has, and their signatures). Cleo uses concepts (``GridboxMaps``,
      ``MicrophysicalProcess``, ``Observer``, ``CoupledDynamics``, ``Motion``, etc. etc.) so that
      when a templated type is instantiated at compile time we can be more specific about what
      that type must be able to do (or, more precicely, what the minimum requirements are for how
      a type's interface must look). This is especially useful for Cleo's flexible components
      (microphysics, observers, coupled dynamics, motion), and its definition of monoids. Concepts
      in C++ are somewhat similar to Python's duck typing or abstract base classes (in case that
      helps you understand).

   GridboxMaps
      The C++ concept providing the geometry of the domain: converting between a gridbox's index
      and its physical coordinate bounds, and finding neighbouring gridboxes. See ``CartesianMaps``
      for its concrete implementation for a Cartesian C-grid/

   SDMMethods
      The class that owns one complete, ready-to-run instance of Cleo's SDM: its ``GridboxMaps``,
      its microphysical process (which may be a combination of several microphysical processes),
      its super-droplet motion, and its (combined) observer. See :doc:`main_cpp`.

   RunCLEO
      The class that owns the outer timestepping routine coupling a ``SDMMethods`` instance to a
      ``CoupledDynamics`` instance via a ``CouplingComms`` instance. See :doc:`main_cpp`.

   cleopy
      Cleo's Python package (``cleopy``), used for generating the required YAML configuration files
      and binary input files for a simulation, and for help handling the reading of Cleo's Zarr
      output via ``xarray`` and ``awkward`` in Python (in combination with the input grid ``.dat``
      binary file and setup ``.txt`` file).

   plotcleo
      Cleo's additional Python package (``plotcleo``) in ``CLEO/examples/exampleplotting`` used
      by the examples to help plot Cleo's Zarr output.


In The Following Sections...
-----------------------------

.. toctree::
   :maxdepth: 1

   code_structure
   main_cpp
   tutorial_new_example
   external_project
   microphysics_and_observers
   coupled_dynamics

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
