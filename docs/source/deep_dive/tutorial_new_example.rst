.. _tutorial_new_example:

Tutorial: Building a New Cleo Example From Scratch
=====================================================

This tutorial walks through everything needed to create a new Cleo example, using one running
example throughout: a minimal 1-D "rainshaft" — a column of gridboxes stacked vertically, with
thermodynamics and wind prescribed from file (one-way coupled), condensation and
collision-coalescence microphysics, and super-droplets falling under gravity/terminal velocity.
This mirrors (and is a simplified version of) the real ``examples/rainshaft1d`` example already
in the repository — cross-referenced throughout so you can compare against working code.

.. note::
   This tutorial assumes you have already built Cleo at least once — see
   :doc:`../usage/installation/installation` — and installed ``cleopy``
   (``uv sync --extra examples``, per :doc:`../usage/installation/installation`).

Step 1: Plan the Simulation
------------------------------

Before writing any code, decide:

* **Domain**: how many spatial dimensions, how many gridboxes, what geometry. *Our example*: 1-D
  (vertical only), a stack of ``ngbxs`` gridboxes with evenly-spaced ``z`` boundaries.
* **Coupled dynamics**: where do pressure/temperature/moisture/wind come from? See
  :doc:`coupled_dynamics`. *Our example*: ``fromfile`` — prescribed hydrostatic profiles with a
  sinusoidal updraught, generated once by a Python script rather than coming from a real
  dynamical core.
* **Microphysics**: which processes, and with what parameters? See
  :doc:`microphysics_and_observers`. *Our example*: condensation/evaporation plus
  collision-coalescence with a hydrodynamic kernel.
* **Initial super-droplets**: how many, what distribution of radii/positions? *Our example*: a
  fixed number of super-droplets per gridbox, radii sampled log-uniformly, multiplicities from a
  lognormal aerosol-size distribution.
* **Observers**: what output do you need to answer your scientific question? See
  :doc:`microphysics_and_observers`. *Our example*: time, super-droplet attributes, and mass
  moments, at a minimum.

Step 2: Create the Configuration File
-----------------------------------------

Cleo reads a single YAML file at startup into a ``Config`` object
(``libs/configuration/config.hpp``), itself built from two structs:

* ``RequiredConfigParams`` (``libs/configuration/required_config_params.hpp``) — fields with *no*
  default; the YAML must set them, or the program throws.
* ``OptionalConfigParams`` (``libs/configuration/optional_config_params.hpp``) — fields with a
  default value, only overridden if present in the YAML (e.g. under a
  ``coupled_dynamics: fromfile`` block).

The required fields, with units, and their effect:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Field (YAML path)
     - Units
     - Effect
   * - ``inputfiles.constants_filename``
     - --
     - Path to ``libs/cleoconstants.hpp`` (physical constants used to non-dimensionalise
       quantities; see :doc:`../intro/memorylayout`).
   * - ``inputfiles.grid_filename``
     - --
     - Path to the binary file of gridbox boundaries (see Step 3).
   * - ``outputdata.setup_filename``
     - --
     - Where a copy of the config + constants used for this run is written, for reproducibility.
   * - ``outputdata.zarrbasedir``
     - --
     - Base directory of the Zarr store observers write into.
   * - ``outputdata.maxchunk``
     - elements
     - Maximum number of elements per Zarr array chunk.
   * - ``domain.nspacedims``
     - --
     - Number of spatial dimensions modelled (``1`` for our example).
   * - ``domain.ngbxs``
     - --
     - Total number of gridboxes.
   * - ``domain.maxnsupers``
     - --
     - Maximum number of super-droplets Cleo will allocate memory for (must be large enough to
       accommodate any super-droplets entering the domain during the run, not just the initial
       count).
   * - ``timesteps.CONDTSTEP``, ``COLLTSTEP``, ``MOTIONTSTEP``, ``COUPLTSTEP``, ``OBSTSTEP``
     - seconds
     - Sub-timesteps for condensation, collisions, super-droplet motion, coupling, and
       observations respectively — see :doc:`../intro/timestepping`.
   * - ``timesteps.T_END``
     - seconds
     - Total simulated time, from ``t=0``.

Selected optional fields relevant to our example (see ``optional_config_params.hpp`` for the
complete set, including Kokkos settings, breakup parameters, and every coupled-dynamics block):

.. code-block:: yaml

    condensation:
      do_alter_thermo: true      # cond/evap feeds back onto gridbox thermodynamic state
      maxniters: 50              # max. Newton-Raphson iterations
      rtol: 1e-6                 # relative tolerance
      atol: 1e-6                 # absolute tolerance
      MINSUBTSTEP: 0.001         # [s] minimum condensation sub-timestep

    coalescence:
      constcoaleff:
        coaleff: 1.0             # constant coalescence efficiency used by LongHydroProb

    fromfiledynamics:
      nspacedims: 1
      press: [path to pressure binary]
      temp: [path to temperature binary]
      qvap: [path to vapour mixing ratio binary]
      qcond: [path to liquid mixing ratio binary]
      wvel: [path to vertical wind binary]

    initsupersfrombinary:
      initsupers_filename: [path to super-droplets binary]
      initnsupers: [initial number of super-droplets]

.. tip::
   You never have to write the YAML file entirely by hand every time: ``cleopy.editconfigfile``
   (used by every example's driver script, see Step 5) copies a source config template and
   programmatically edits the paths within it — see ``rainshaft1d.py``'s
   ``editconfigfile.edit_config_params(config_filename, config_params)`` call.

Step 3: Create the Input Files
----------------------------------

Cleo does not read NetCDF or YAML for its bulk numerical input — gridbox boundaries, initial
super-droplet attributes, and (for ``fromfile`` coupling) the dynamics fields are all read from
Cleo's own simple flat binary format, generated by ``cleopy``. You do this from a Python script,
not by hand. The pattern (mirroring ``examples/rainshaft1d/rainshaft1d_inputfiles.py``):

.. code-block:: python

    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import crdgens, rgens, dryrgens, probdists, attrsgen
    from cleopy.thermobinary_src import thermogen, windsgen, thermodyngen

    # 1. gridbox boundaries: evenly-spaced z, single-cell x and y (1-D column)
    zgrid = [0, zgrid_max, zgrid_max / ngbxs]     # [zmin, zmax, delta] in metres
    xgrid = [0, xgrid_max]
    ygrid = [0, ygrid_max]
    geninitconds.generate_gridbox_boundaries(
        grid_filename, zgrid, xgrid, ygrid, constants_filename)

    # 2. thermodynamics: hydrostatic pressure/temperature profile with prescribed lapse rates,
    #    plus a sinusoidal updraught for the vertical wind
    thermog = thermogen.HydrostaticLapseRates(
        config_filename, constants_filename, PRESSz0, TEMPz0, qvapz0, Zbase,
        TEMPlapses, qvaplapses, qcond=0.0)
    windsg = windsgen.SinusoidalUpdraught(WVEL, None, None, Wlength)
    geninitconds.generate_thermodynamics_conditions_fromfile(
        thermofiles, thermodyngen.ThermodynamicsGenerator(thermog, windsg),
        config_filename, constants_filename, grid_filename)

    # 3. initial super-droplets: radii sampled log-uniformly, multiplicities from a lognormal
    #    aerosol-size distribution, positions sampled randomly within the column
    initattrsgen = attrsgen.AttrsGenerator(
        radiigen=rgens.SampleLog10RadiiGen(rspan),
        dryradiigen=dryrgens.ScaledRadiiGen(dryr_sf),
        xiprobdist=probdists.LnNormal(geomeans, geosigs, scalefacs),
        coord3gen=crdgens.SampleCoordGen(True), coord1gen=None, coord2gen=None)
    geninitconds.generate_initial_superdroplet_conditions(
        initattrsgen, initsupers_filename, config_filename, constants_filename,
        grid_filename, nsupers, numconc)

Each ``generate_*`` call writes one (or a small family of) binary file(s) whose paths must match
what you put in your YAML config (Step 2). ``cleopy``'s generators are composable in the same
spirit as Cleo's C++ monoids — you pick a distribution/generator object for each attribute and
pass it in, rather than writing bespoke code per example.

Step 4: Wire Up Microphysics, Observers and Coupled Dynamics in Code
-------------------------------------------------------------------------

Copy ``roughpaper/src/main.cpp`` and ``main_impl.hpp`` (or ``examples/rainshaft1d/src/
main_rshaft1d.cpp`` as a single-file starting point) into your new example's ``src/`` directory,
then, per :doc:`main_cpp`:

1. Adjust ``create_gbxmaps`` if your domain differs (it likely won't, for a 1-D column).
2. Adjust ``create_microphysics``/``config_condensation``/``config_collisions`` for your chosen
   processes and their config-driven parameters (:doc:`microphysics_and_observers`).
3. Adjust ``create_coupldyn`` and the ``CouplingComms<...>`` type in ``main()`` for your chosen
   coupled dynamics (:doc:`coupled_dynamics`) — ``fromfile``/``FromFileDynamics``/
   ``FromFileComms`` for our running example.
4. Adjust ``create_observer`` for what you need written out.
5. Adjust ``create_initconds`` only if your initial-conditions binary format/generator differs
   from the default ``InitSupersFromBinary``/``InitGbxsNull`` pair.

Then write a ``CMakeLists.txt`` for your example, modelled on
``examples/rainshaft1d/src/CMakeLists.txt``:

.. code-block:: cmake

    project("myexample")
    set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)
    add_executable(myexample EXCLUDE_FROM_ALL "main_myexample.cpp")
    target_link_libraries(myexample PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
    target_link_libraries(myexample PUBLIC Kokkos::kokkos)
    target_include_directories(myexample PRIVATE "${CLEO_SOURCE_DIR}/libs")
    set_target_properties(myexample PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON
      CMAKE_CXX_EXTENSIONS ON CXX_STANDARD 20)

...and add ``add_subdirectory(myexample)`` to ``examples/CMakeLists.txt`` so the top-level build
picks it up.

Step 5: Build and Run
--------------------------

1. Configure and build (see :doc:`../usage/installation/installation` for full compiler/Kokkos
   flags):

   .. code-block:: console

     $ cmake -S . -B build -DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian
     $ cmake --build build --target myexample

2. Write a Python driver script mirroring ``rainshaft1d.py``'s three-stage pattern
   (``--do_inputfiles`` / ``--do_run_executable`` / ``--do_plot_results``): generate the config
   and input files (Step 2/3), run the compiled executable with the config file as its one
   argument, then plot. Cleo executables are always invoked the same way:

   .. code-block:: console

     $ ./build/examples/myexample/src/myexample /path/to/config.yaml

3. Optionally add a bash script under ``scripts/vanilla/examples/`` following the pattern of
   ``rainshaft1d.sh``, so your example fits into the existing
   :doc:`../usage/examples/examples_vanilla` workflow.

Step 6: Verify / Sanity-Check the Output
--------------------------------------------

Cleo writes output as a Zarr store, readable with ``cleopy`` or plain ``xarray``
(``xr.open_zarr(...)``). At minimum, check:

* The run completed and produced the number of output timesteps you expect
  (``T_END / OBSTSTEP + 1``).
* Total super-droplet multiplicity-weighted liquid water mass is conserved up to
  expected physical changes (precipitation falling out of the domain, if applicable) — a sudden,
  unphysical jump usually indicates a units or config error rather than a real result.
* Individual super-droplet radii stay within physically sensible bounds (not NaN, not negative,
  not orders of magnitude larger than the domain).
* If comparing to a known solution (e.g. an analytic collision-kernel test), plot against it — see
  ``examples/exampleplotting`` for utilities.

Common Errors at Each Step
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Step
     - Likely cause and fix
   * - Config file
     - ``std::invalid_argument``/YAML parse errors usually mean a required field
       (``required_config_params.hpp``) is missing, or a path is wrong/relative when an absolute
       path was expected — check every path in the YAML resolves from wherever you run the
       executable, not from the config file's own directory.
   * - Input files
     - A "wrong file size" or similar exception when Cleo reads a binary file almost always means
       a mismatch between what ``cleopy`` generated and what ``domain.ngbxs``/``nspacedims`` in
       the YAML say to expect — regenerate both from the same source of truth (a single Python
       config dict, as ``rainshaft1d.py`` does) rather than editing either independently.
   * - Wiring in ``main_impl.hpp``
     - A very long compiler error, usually naming a ``concept`` and ``requires`` clause, means
       some type you passed doesn't satisfy the concept expected at that point (e.g. a
       ``CouplingComms`` that doesn't match your chosen ``CoupledDynamics`` type — see the
       pitfall in :doc:`coupled_dynamics`). Read the *first* concept-related error in the output;
       later ones are often knock-on effects of the first.
   * - Build/CMake
     - "undefined reference" at link time usually means a needed ``coupldyn_*``/``cartesiandomain``
       library isn't listed in your executable's ``target_link_libraries``, even though the
       header compiled fine (see the pitfall in :doc:`main_cpp`).
   * - Runtime
     - A silent hang or a Kokkos-related abort often means the coupling timestep of your
       coupled-dynamics object and your SDM disagree — ``RunCLEO``'s constructor checks this
       (``check_coupling()`` in ``libs/runcleo/runcleo.hpp``) and throws
       ``"coupling timestep of dynamics solver and CLEO SDM are not equal"`` if so; make sure
       ``COUPLTSTEP`` in your YAML and however your ``create_coupldyn`` computes its coupling
       step agree.
   * - Output
     - Missing or empty output variables usually mean the corresponding observer wasn't included
       in your ``create_observer`` combination (:doc:`microphysics_and_observers`) — check the
       ``>>``-chain includes it.

See also :doc:`../usage/examples/examples_vanilla` for the general mechanics of running any
example's bash script once your example is wired up.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
