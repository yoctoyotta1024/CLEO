.. _code_structure:

Cleo's Code Structure
======================

This page introduces Cleo's code/library structure.

.. _repo_map:

Top-Level Repository Map
--------------------------

::

    CLEO/
    ├── libs/                      # Cleo's C++ library code (the model itself)
    │   ├── configuration/         # reading YAML config files into typed C++ structs
    │   ├── initialise/            # turning config + input files into initial conditions
    │   ├── runcleo/               # the concepts/classes that run Cleo SDM with CoupledDynamics
    │   ├── superdrops/            # super-droplets definitions and microphysical processes
    │   │   └── collisions/        # collision-coalescence, breakup and rebound schemes
    │   ├── gridboxes/             # gridboxes (as in volumes of space) and their super-droplet
    |   |                             handling
    │   ├── cartesiandomain/       # an implementatation of the GridboxMaps concept: a Cartesian
    |   |                             (C-grid) domain, including droplet motion and MPI domain
    |   |                             decomposition
    │   ├── coupldyn_null/         # an implementatation of the CoupledDynamics concept:
    |   |                             no dynamics at all
    │   ├── coupldyn_fromfile/     # an implementatation of the CoupledDynamics concept:
    |   |                             read from binary files (one-way)
    │   ├── coupldyn_cvode/        # an implementatation of the CoupledDynamics concept:
    |   |                             CVODE ODE solver (two-way; for a 0-D parcel model only)
    │   ├── coupldyn_yac/          # an implementatation of the CoupledDynamics concept:
    |   |                             YAC-coupled external dynamical core (one-way or two-way; ICON)
    │   ├── observers/             # observers: for writing/streaming data output
    │   ├── zarr/                  # a minimal Zarr dataset writer used by some observers for output
    |   └── cleo_python_bindings/  # pybind11 bindings exposing parts of Cleo to Python and an
    |                                 (one-way or two-way) implementatation of the CoupledDynamics
    |                                 concept: for numpy arrays
    ├── examples/                  # one directory per example, each with its own main*.cpp,
    │   │                          #   CMakeLists.txt, and Python driver script(s)
    │   └── exampleplotting/       # the `plotcleo` Python package used to plot example output
    ├── roughpaper/src/            # main.cpp / main_impl.hpp used as a template/scratchpad
    ├── cleopy/                    # Cleo's installable Python package: config editing, binary
    │                              #   input-file generation, and reading output
    ├── scripts/                   # bash scripts to build/compile/run on Levante, Jupiter or a
    │                              #   generic "vanilla" computer
    ├── docs/                      # this documentation (Sphinx + Doxygen/Breathe)
    ├── extern/                    # third-party dependencies (Kokkos, yaml-cpp, ...)
    └── CMakeLists.txt             # top-level build configuration

.. _why_code_structure:

Why It's Structured This Way
-------------------------------

The directory structure and their dependencies are conceptually important. Superdroplets and their
microphysical processes exist and act on a certain (thermodynamic) state, however, since
superdroplets, microphysics, and a thermodynamic state are somewhat abstract numerical/physical
descriptions, they are independent of the grid-boxes and domain geometries which we use in a
simulation. Therefore they live, self-contained, in ``libs/superdrops``. The gridboxes depend on
``libs/superdrops`` since they contain a state and (a view to) super-droplets, but they are
independent of a particular domain geometry. Therefore, gridboxes they live in ``libs/gridboxes``,
alongside the ``GridboxMaps`` concept which defines how gridboxes can be mapped to an arbitrary
domain geometry (without specifying any specific geometry). A concrete domain geometry and
its associated methods for super-droplet motion are dependent on the particular implementation of the
GridboxMaps concept, they thus depend on ``libs/gridboxes`` and are stored in a seperate library,
e.g. ``libs/cartesiandomain`` for a cartesian domain. If Cleo grows another domain type in the
future (e.g. unstructured), it would live alongside ``cartesiandomain``, implementing its own
version of the ``GridboxMaps`` concept. The directory for coupled dynamics
implementations (``libs/coupldyn_*``) are separate from all the SDM libraries (superdroplets,
gridboxes, observers etc.), because the coupled dynamics is independent of SDM. It is the
``libs/runcleo`` library which ties everything together: the model configuration, initialisation,
coupled dynamics (and coupler), superdroplets (including their motion and microphysics), gridboxes,
domain, and observers (including Zarr dataset output).

Occasionally this seperation between ``superdroplets``, ``gridboxes`` and a particular
``GridboxMaps`` implementation (i.e. domain geometry) requires a component of the model
to be divided into seperate parts across different libraries. This is the case for super-droplet
motion. The abstract concept for super-droplet ``Motion`` is defined in ``libs/superdroplets``
along with the ``VelocityFormula`` concept and several formulations of a superdroplet's terminal
velocity which satisfy that concept. The structure for moving superdroplets through a domain,
``MoveSupersInDomain``, is then defined in ``libs/gridboxes``, since this library handles how
superdroplets relate to gridboxes, however, this is a templated structure which is agnostic to a
particular domain geometry. Therefore, the actual implementation of the ``MoveSupersInDomain``
structure only comes together in a particular domain geometry's library
(e.g. in ``libs/cartesiandomain/movement``). This is where ``MoveSupersInDomain`` is filled-in with
an implmentation of the ``Motion`` concept, as well as implementations for how to transport
superdroplets across the domain and how to enact boundary conditions for that particular geometry.

Some libraries are *always* needed regardless of setup: ``configuration``, ``initialise``,
``gridboxes``, ``runcleo``, ``observers``, ``zarr``, and ``superdrops``, and are built into every
executable, as per ``libs/CMakeLists.txt``. The libraries that can vary between setups —
which coupled dynamics and which domain geometry you use — are built only when selected,
controlled by two CMake cache variables read in ``libs/CMakeLists.txt``:

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

Library Overview
----------------

Reading each library in isolation can obscure how they fit together at runtime. In one
sentence: ``RunCLEO`` repeatedly calls ``CoupledDynamics::run_step`` to advance the external
dynamics, calls ``SDMMethods::run_step`` to advance super-droplet motion and microphysics, and
uses ``CouplingComms`` to exchange gridbox state between the two — with the whole loop's step
size chosen adaptively based on the coupling and observation timesteps. The full walkthrough,
with the relevant code from ``libs/runcleo/runcleo.hpp`` and ``sdmmethods.hpp``, is in
:doc:`main_cpp`; the adaptive-timestepping mechanism specifically is described in
:doc:`../intro/timestepping`.

Nevertheless, here is each library in a bit more detail...

``superdrops/`` — super-droplets and microphysics
    Defines the ``Superdrop`` type and its attributes (radius, multiplicity, solute mass,
    position, optional identifier, etc.), the ``State`` (pressure, temperature, vapour  mixing
    ratio, wind, etc.), and the ``MicrophysicalProcess`` concept. The monoid set definition for
    microphysical processes is also defined here, as well as concrete microphysical
    processes: ``Condensation`` (``condensation.hpp``/``.cpp``, using an implicit Euler solver
    in ``impliciteuler.hpp``) and, under ``superdrops/collisions/``, collision-coalescence
    (often denoted``CollCoal``), and  collision-coalescence,-breakup,-rebound (``CoalBuRe``).
    This is where the physics described in the numerics paper (B2026b :cite:`bayley2026b`)
    is implemented. Refer to :doc:`microphysics_and_observers` for how to combine different
    microphysical processes monoidally.

``gridboxes/`` — domain decomposition and super-droplet transport
    Defines ``Gridbox`` (one cell's ``State`` and the super-droplets currently inside it),
    the ``GridboxMaps`` concept for converting between gridbox indexes and physical coordinates,
    and the ``MoveSupersInDomain`` framework for moving super-droplets both within and between
    gridboxes each sub-timestep (including the ``BoundaryConditions`` concept for creating or
    removing superdroplets at the boundaries of the domain). This is where
    Cleo's :doc:`../intro/memorylayout` is realised: all gridboxes are one contiguous
    ``Kokkos::View``, all super-droplets are another, ordered by gridbox, and each gridbox holds a
    reference/pointer-like object (``SupersInGbx``) into the slice of the super-droplet array that
    belongs to it.

``cartesiandomain/`` — a concrete implementatation of the ``GridboxMaps`` concept
    Implements ``CartesianMaps`` which satisfies the ``GridboxMaps`` concept for a structured
    cartesian C-grid domain including MPI domain-decomposition logic for splitting the domain
    across MPI processes. This library also includes a concrete implementation of
    ``MoveSupersInDomain`` for this cartisian domain which pertains of ``CartesianMotion``
    (satisfying the ``Motion`` concept) for advecting super-droplets through the cartesian grid,
    ``CartesianTransportAcrossDomain`` for updating super-droplet gridbox-indexes when they move
    across the domain (in x, y, or z), and additionally, an optional implementation of
    the ``BoundaryConditions`` concept, ``AddSupersToDomain``, which removes and then adds
    superdroplets to the top or bottom of the domain.

``observers/`` — data output
    Defines the ``Observer`` concept and its monoid combinators (``observers.hpp``), plus
    concrete observers for writing output to a Zarr dataset (via ``zarr/``), e.g. for time,
    gridbox indexes, super-droplet attributes, moments of the droplet mass distribution,
    thermodynamic state or wind fields. ``observers/sdmmonitor/`` holds special observers that
    "monitor" within the microphysics/motion sub-timestepping routine itself, rather than at the
    start or end of each timestep, e.g. to record diagnostics of precipitation or
    condenstaion/evaporation rates, and then output them to a Zarr dataset.
    See :doc:`microphysics_and_observers` for more information. This library also contains methods
    to output progress to the screen/terminal during a run (``StreamOutObserver``).

``runcleo/`` — tying SDM to coupled dynamics
    The smallest but most structurally important directory. Defines the ``CoupledDynamics`` and
    ``CouplingComms`` concepts (see :doc:`coupled_dynamics`), the ``SDMMethods`` class (one
    complete, self-contained SDM instance), and the ``RunCLEO`` class (the outer timestepping
    loop). See :doc:`main_cpp`. Also defines how initial conditions are read in to create gridboxes
    and super-droplets.

``configuration/`` — reading YAML config into typed structs
    Reads a YAML configuration file into a ``Config`` object built from two structs,
    ``RequiredConfigParams`` (no sensible default — must be in the YAML) and
    ``OptionalConfigParams`` (has a default). See :doc:`tutorial_new_example`.

``initialise/`` — building the first Gridboxes and super-droplets
    Reads Cleo's custom binary input-files (produced by ``cleopy``, see
    :doc:`tutorial_new_example`) to construct the initial conditions (``InitConds``) comprised of
    Gridbox initial-data arrays (``GbxInitConds``) and Superdroplet initial-data arrays
    (``SuperdropInitConds``). These arrays are then used by ``RunCleo`` to create gridboxes
    and superdroplets.

``coupldyn_*/`` — coupled dynamics implementations
    One directory per option for where the thermodynamic/wind fields come from. See
    :doc:`coupled_dynamics` for what each one does and when to use it.

``zarr/`` — output backend
    A minimal, dependency-light implementation of the
    `Zarr Storage Specification Version 2 <https://zarr-specs.readthedocs.io/en/v3.1.0/v2/v2.0.html>`_
    chunked-array storage format, used by ``observers/`` to write output, e.g. that ``xarray`` with
    the ``zarr`` engine (through ``cleopy`` helper modules) can read back.

``cleo_python_bindings/`` — pybind11 bindings
    Exposes selected pieces of Cleo to Python (enough to run CLEO-SDM coupled to numpy arrays) via
    `pybind11 <https://github.com/pybind/pybind11>`_, so Cleo can be driven and its state
    inspected from a numpy/Python-based coupled model. This is what the ``"numpy"`` coupled
    dynamics option (mentioned above) refers to - see ``CLEO/examples/python_bindings/`` for a
    working example of how these bindings can be used.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
