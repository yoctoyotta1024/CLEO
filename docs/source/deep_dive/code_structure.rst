.. _code_structure:

Cleo's Code Structure
======================

This page maps the repository and introduces the handful of C++ ideas you need in order to read
that map. If you already know what a header file is and what a C++ concept is, skip ahead to
:ref:`repo_map`.

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
    (e.g. ``superdrops/condensation.hpp`` and ``superdrops/condensation.cpp``). Templated code
    (see below) is usually written entirely in the header, because the compiler needs the full
    definition available at every place it's used.

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
    generates a version of ``SDMMethods`` specifically for ``Condensation``. This is why swapping
    microphysics or observers costs nothing at runtime (see :doc:`../intro/kokkos` and the design
    paper, Sect. 2): there is no dynamic dispatch, no virtual function calls, no ``if``/``else``
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
    depth, with worked examples, in :doc:`microphysics_and_observers`; the conceptual
    introduction (with the "wet oil paint" analogy) is in :doc:`../intro/monoids`.

Namespaces
    A C++ ``namespace`` groups related names to avoid clashes, similar to how a Python package
    groups modules. ``KokkosCleoSettings`` (aliased ``KCS``) and ``LIMITVALUES`` are examples
    you'll see used in Cleo.

Kokkos views and ``KOKKOS_INLINE_FUNCTION``
    Cleo uses the `Kokkos <https://github.com/kokkos>`_ library for performance-portable thread
    parallelism (see :doc:`../intro/kokkos`); a ``Kokkos::View`` is roughly an array that can live
    on CPU or GPU memory, and ``KOKKOS_INLINE_FUNCTION`` marks a function so it can run inside a
    parallel region on either. You don't need to understand Kokkos internals to read
    ``main.cpp``, but you will see its types (``viewd_gbx``, ``subviewd_supers``, ...) passed
    around — these are Cleo's own aliases defined in ``libs/kokkosaliases.hpp`` and
    ``libs/superdrops/kokkosaliases_sd.hpp``.

.. _repo_map:

Top-Level Repository Map
--------------------------

::

    CLEO/
    ├── libs/                    # Cleo's C++ library code (the model itself)
    │   ├── configuration/       # reading YAML config files into typed C++ structs
    │   ├── initialise/          # turning config + binary input files into initial conditions
    │   ├── superdrops/          # super-droplets themselves and microphysical processes
    │   │   └── collisions/      # collision-coalescence, breakup and rebound schemes
    │   ├── gridboxes/           # gridboxes, super-droplet motion between/within them
    │   ├── cartesiandomain/     # a concrete GridboxMaps: a Cartesian (C-grid) domain
    │   ├── runcleo/             # the concepts/classes that tie SDM to coupled dynamics
    │   ├── coupldyn_null/       # coupled dynamics: no dynamics at all
    │   ├── coupldyn_fromfile/   # coupled dynamics: read from binary files (one-way)
    │   ├── coupldyn_cvode/      # coupled dynamics: 0-D CVODE ODE solver (two-way)
    │   ├── coupldyn_yac/        # coupled dynamics: YAC-coupled external dynamical core
    │   ├── cleo_python_bindings/# pybind11 bindings exposing parts of Cleo to Python/numpy
    │   ├── observers/           # observers: writing/streaming data out during a run
    │   └── zarr/                # a minimal Zarr writer used by observers for output
    ├── examples/                # one directory per example, each with its own main*.cpp,
    │   │                        #   CMakeLists.txt, and Python driver script(s)
    │   └── exampleplotting/     # the `plotcleo` Python package used to plot example output
    ├── roughpaper/src/          # main.cpp / main_impl.hpp used as a template/scratchpad —
    │                            #   see :doc:`main_cpp`
    ├── cleopy/                  # Cleo's installable Python package: config editing, binary
    │                            #   input-file generation, output reading
    ├── scripts/                 # bash scripts to build/compile/run on Levante, JUWELS or a
    │                            #   generic ("vanilla") computer
    ├── docs/                    # this documentation (Sphinx + Doxygen/Breathe)
    ├── extern/                  # vendored third-party dependencies (Kokkos, yaml-cpp, ...)
    └── CMakeLists.txt           # top-level build configuration, see :doc:`external_project`

Two directory names are easy to misread: ``libs/gridboxes`` contains the *generic* gridbox and
super-droplet-transport machinery (the ``GridboxMaps`` concept, motion between gridboxes,
boundary conditions), whereas ``libs/cartesiandomain`` contains one specific, concrete
implementation of that machinery for a Cartesian C-grid domain (``CartesianMaps``). If Cleo grows
another domain type in future (e.g. unstructured), it would live alongside
``cartesiandomain``, implementing the same ``GridboxMaps`` concept.

.. note::
   Ignore any ``.patch`` files and any ``build*`` directories (e.g. ``build_bubble3d``) you see
   at the repository root — these are local build artefacts and patches, not part of Cleo's
   source.

Why It's Structured This Way
-------------------------------

The split mirrors the physical/software separation the design paper draws out (Sect. 2): the
things that are *always* needed regardless of setup (``configuration``, ``initialise``,
``gridboxes``, ``runcleo``, ``observers``, ``zarr``, generic ``superdrops`` machinery) live in
``libs/`` and are built into every executable, per ``libs/CMakeLists.txt``. The things that vary
between setups — which coupled dynamics you use, which concrete domain/gridbox geometry — are
each their own subdirectory, built only when selected, controlled by two CMake cache variables
read in ``libs/CMakeLists.txt``:

.. code-block:: cmake

    # from libs/CMakeLists.txt
    set(valid_coupled_dynamics "null" "fromfile" "cvode" "yac" "numpy")
    # ...
    set(valid_domain "cartesian")

You choose between them when you invoke CMake, e.g. ``-DCLEO_COUPLED_DYNAMICS=fromfile
-DCLEO_DOMAIN=cartesian`` (the default, with no flags, builds both ``fromfile`` and ``cvode``).
This is also why an executable's ``CMakeLists.txt`` (e.g.
``examples/rainshaft1d/src/CMakeLists.txt``) explicitly lists which of these optional libraries
it links against:

.. code-block:: cmake

    set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)
    add_executable(rshaft1d EXCLUDE_FROM_ALL "main_rshaft1d.cpp")
    target_link_libraries(rshaft1d PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
    target_link_libraries(rshaft1d PUBLIC Kokkos::kokkos)

Subsystem-by-Subsystem
------------------------

``superdrops/`` — super-droplets and microphysics
    Defines the ``Superdrop`` type and its attributes (radius, multiplicity, solute mass,
    position, optional id), the ``MicrophysicalProcess`` concept and its monoid combinators
    (``microphysicalprocess.hpp``), and concrete microphysical processes: ``Condensation``
    (``condensation.hpp``/``.cpp``, using an implicit Euler solver in ``impliciteuler.hpp``) and,
    under ``superdrops/collisions/``, collision-coalescence, breakup and rebound schemes. This is
    where the physics described in the numerics paper is implemented. See
    :doc:`microphysics_and_observers`.

``gridboxes/`` — domain decomposition and super-droplet transport
    Defines ``Gridbox`` (one cell's ``State`` — pressure, temperature, vapour/liquid mixing
    ratio, wind — plus the super-droplets currently inside it), the ``GridboxMaps`` concept for
    converting between gridbox index and physical coordinates, and ``MoveSupersInDomain`` /
    ``movesupersindomain.hpp`` for moving super-droplets both within and between gridboxes each
    sub-timestep (including across the boundary of the whole domain via ``BoundaryConditions``).
    This is where the memory layout described in :doc:`../intro/memorylayout` is realised: all
    gridboxes are one contiguous ``Kokkos::View``, all super-droplets are another, ordered by
    gridbox, and each gridbox holds a reference/pointer-like object (``SupersInGbx``) into the
    slice of the super-droplet array that belongs to it.

``cartesiandomain/`` — a concrete Cartesian ``GridboxMaps``
    Implements ``CartesianMaps``, satisfying the ``GridboxMaps`` concept for a structured C-grid
    domain, plus ``CartesianMotion`` (satisfying the ``Motion`` concept) for advecting
    super-droplets through that grid, and MPI domain-decomposition logic
    (``cartesian_decomposition.hpp``) for splitting the domain across nodes.

``observers/`` — getting data out
    Defines the ``Observer`` concept and its monoid combinators (``observers.hpp``), plus
    concrete observers for writing time, gridbox index, super-droplet attributes, mass moments,
    wind and thermodynamic state to a Zarr dataset (via ``zarr/``), or for printing progress to
    the screen (``StreamOutObserver``). ``observers/sdmmonitor/`` holds special observers that
    "monitor" (record diagnostics on) the microphysics itself as it runs. See
    :doc:`microphysics_and_observers`.

``runcleo/`` — tying SDM to coupled dynamics
    The smallest but most structurally important directory. Defines the ``CoupledDynamics`` and
    ``CouplingComms`` concepts (see :doc:`coupled_dynamics`), the ``SDMMethods`` class (one
    complete, self-contained SDM instance), and the ``RunCLEO`` class (the outer timestepping
    loop). See :doc:`main_cpp`.

``configuration/`` — reading YAML config into typed structs
    Reads a YAML configuration file into a ``Config`` object built from two structs,
    ``RequiredConfigParams`` (no sensible default — must be in the YAML) and
    ``OptionalConfigParams`` (has a default). See :doc:`tutorial_new_example`.

``initialise/`` — building the first Gridboxes and super-droplets
    Reads Cleo's custom binary input-file format (produced by ``cleopy``, see
    :doc:`tutorial_new_example`) to construct the initial ``Gridbox`` array and initial
    ``Superdrop`` array before timestepping 0.

``coupldyn_*/`` — coupled dynamics implementations
    One directory per option for where the thermodynamic/wind fields come from. See
    :doc:`coupled_dynamics` for what each one does and when to use it.

``zarr/`` — output backend
    A minimal, dependency-light implementation of the `Zarr <https://zarr.dev>`_ chunked-array
    storage format, used by ``observers/`` to write output that ``cleopy`` (or ``xarray`` with
    the ``zarr`` engine) can read back.

``cleo_python_bindings/`` — pybind11 bindings
    Exposes selected pieces of Cleo (enough to drive an SDM run) to Python via
    `pybind11 <https://github.com/pybind/pybind11>`_, so Cleo can be driven and its state
    inspected from a numpy/Python-based coupled model. This is what the ``"numpy"`` coupled
    dynamics option (mentioned above) refers to.

The Simulation Loop, One Level Up
-----------------------------------

Reading these subsystems in isolation can obscure how they fit together at runtime. In one
sentence: ``RunCLEO`` repeatedly calls ``CoupledDynamics::run_step`` to advance the external
dynamics, calls ``SDMMethods::run_step`` to advance super-droplet motion and microphysics, and
uses ``CouplingComms`` to exchange gridbox state between the two — with the whole loop's step
size chosen adaptively based on the coupling and observation timesteps. The full walkthrough,
with the relevant code from ``libs/runcleo/runcleo.hpp`` and ``sdmmethods.hpp``, is in
:doc:`main_cpp`; the adaptive-timestepping mechanism specifically is described in
:doc:`../intro/timestepping`.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
