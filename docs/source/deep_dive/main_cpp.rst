.. _main_cpp_page:

Understanding and Adapting ``main.cpp``
==========================================

Every Cleo executable is, in essence, one ``main.cpp`` file (optionally split into a thin
``main.cpp`` plus a ``main_impl.hpp`` holding the "configuration" functions, as
``roughpaper/src/`` does). This page walks through ``roughpaper/src/main.cpp`` and
``main_impl.hpp`` block by block. Treat ``roughpaper/src/`` as a template/scratchpad: it is built
with ``-DCLEO_NO_ROUGHPAPER=false`` (the default) and is meant to be copied from and experimented
with, rather than being "the" canonical example — see :doc:`../usage/quickstart`.

Overall Control Flow
-----------------------

At the top level, every Cleo executable's ``main()`` does the same six things, in the same order:

1. Parse the path to a YAML configuration file from the command line and build a ``Config``.
2. Initialise MPI (via a ``Communicator``) and Kokkos.
3. Build one ``SDMMethods`` object: everything needed for the SDM part of the model
   (``GridboxMaps``, combined microphysics, super-droplet motion, combined observer).
4. Build one object satisfying ``CoupledDynamics``, and one satisfying
   ``CouplingComms<YourGbxMaps, YourCoupledDynamics>``.
5. Build the initial conditions, construct a ``RunCLEO`` from the ``SDMMethods``, the coupled
   dynamics, and the comms, then call it (this runs the entire simulation).
6. Finalise Kokkos and exit.

``main.cpp`` and ``main_impl.hpp`` Together
----------------------------------------------

``main.cpp`` itself is deliberately short — it only ``#include``s ``main_impl.hpp`` and defines
``main()``:

.. code-block:: cpp

    #include "./main_impl.hpp"

    int main(int argc, char* argv[]) {
      if (argc < 2) {
        throw std::invalid_argument("configuration file(s) not specified");
      }

      Kokkos::Timer kokkostimer;

      /* Read input parameters from configuration file(s) */
      const std::filesystem::path config_filename(argv[1]);
      const Config config(config_filename);

      /* Initialize Communicator here */
      init_communicator init_comm(argc, argv, config);

      /* Initialise Kokkos parallel environment */
      Kokkos::initialize(config.get_kokkos_initialization_settings());
      {
        Kokkos::print_configuration(std::cout);

        const Timesteps tsteps(config.get_timesteps());

        auto store = FSStore(config.get_zarrbasedir());
        auto dataset = SimpleDataset(store);

        /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
        const SDMMethods sdm(create_sdm(config, tsteps, dataset, store));

        /* Solver of dynamics coupled to CLEO SDM */
        CoupledDynamics auto coupldyn =
            create_coupldyn(config, sdm.gbxmaps, tsteps.get_couplstep(), tsteps.get_t_end());

        /* coupling between coupldyn and SDM */
        const CouplingComms<CartesianMaps, FromFileDynamics> auto comms = FromFileComms{};

        /* Initial conditions for CLEO run */
        const InitialConditions auto initconds = create_initconds(config, sdm.gbxmaps);

        /* Run CLEO (SDM coupled to dynamics solver) */
        const RunCLEO runcleo(sdm, coupldyn, comms);
        runcleo(initconds, tsteps.get_t_end());
      }
      Kokkos::finalize();

      const auto ttot = double{kokkostimer.seconds()};
      std::cout << "CLEO Total Program Duration: " << ttot << "s\n";
      return 0;
    }

Everything Kokkos-related must happen *between* ``Kokkos::initialize`` and ``Kokkos::finalize``
(hence the extra ``{ }`` block) — this is a Kokkos requirement, not a Cleo-specific one. Notice
that every "ingredient" (``sdm``, ``coupldyn``, ``comms``, ``initconds``) is declared using its
*concept* name with ``auto`` (``CoupledDynamics auto coupldyn = ...``), never a concrete class
name — this is what lets you swap out, e.g., ``FromFileDynamics`` for ``NullDynamics`` by editing
only ``main_impl.hpp``'s factory functions, and never having to touch ``main()`` itself.

``main_impl.hpp`` holds the ``create_*`` factory functions that ``main()`` calls. Reading it
top-to-bottom mirrors the order things are needed in ``main()``:

.. code-block:: cpp

    inline GridboxMaps auto create_gbxmaps(const Config& config) {
      return create_cartesian_maps(config.get_ngbxs(), config.get_nspacedims(),
                                   config.get_grid_filename());
    }

Builds the domain geometry from the config's ``ngbxs``, number of spatial dimensions, and the
gridbox-boundaries binary file (see :doc:`tutorial_new_example`).

.. code-block:: cpp

    inline Motion<CartesianMaps> auto create_motion(const unsigned int motionstep) {
      const auto terminalv = RogersGKTerminalVelocity{};
      return CartesianMotion(motionstep, &step2dimlesstime, terminalv);
    }

    inline BoundaryConditions<CartesianMaps> auto create_boundary_conditions(const Config& config) {
      return NullBoundaryConditions{};
    }

    template <GridboxMaps GbxMaps>
    inline auto create_movement(const Config& config, const Timesteps& tsteps, const GbxMaps& gbxmaps) {
      const Motion<CartesianMaps> auto motion = create_motion(tsteps.get_motionstep());
      const BoundaryConditions<CartesianMaps> auto boundary_conditions =
          create_boundary_conditions(config);
      return cartesian_movement(gbxmaps, motion, boundary_conditions);
    }

Builds how super-droplets move: a terminal-velocity formula (here, Rogers & Yau / "GK" —
commented-out alternatives such as ``NullTerminalVelocity`` or ``SimmelTerminalVelocity`` show
where you'd swap it), wrapped in ``CartesianMotion`` (satisfying the ``Motion<CartesianMaps>``
concept), plus what happens to super-droplets that reach the domain boundary
(``NullBoundaryConditions``, or ``AddSupersToDomain`` to replenish super-droplets leaving the
top/bottom, commented out immediately above).

.. code-block:: cpp

    inline MicrophysicalProcess auto config_condensation(const Config& config, const Timesteps& tsteps) {
      const auto c = config.get_condensation();
      return Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                          c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);
    }

    inline MicrophysicalProcess auto config_collisions(const Config& config, const Timesteps& tsteps) {
      const PairProbability auto collprob = LongHydroProb();
      const NFragments auto nfrags = ConstNFrags(config.get_breakup().constnfrags.nfrags);
      const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag(RogersGKTerminalVelocity{});
      return CoalBuRe(tsteps.get_collstep(), &step2realtime, collprob, nfrags, coalbure_flag);
    }

    inline MicrophysicalProcess auto create_microphysics(const Config& config, const Timesteps& tsteps) {
      const MicrophysicalProcess auto cond = config_condensation(config, tsteps);
      const MicrophysicalProcess auto colls = config_collisions(config, tsteps);
      return colls >> cond;
    }

This is where the microphysics monoid combination from :doc:`microphysics_and_observers` is
assembled — here, collision/coalescence/rebound/breakup (via ``CoalBuRe``, using a
collision-outcome flag ``TSCoalBuReFlag`` based on Testik & Straub, and a number-of-fragments
model ``ConstNFrags``) combined with condensation. Note ``main_impl.hpp`` keeps several
alternative configurations as commented-out code directly above ``config_collisions`` — this is
the pattern the file uses throughout to suggest what else you could plug in.

.. code-block:: cpp

    inline CoupledDynamics auto create_coupldyn(const Config& config, const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep, const unsigned int t_end) {
      const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
      const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});
      const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);
      return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
    }

    template <GridboxMaps GbxMaps>
    inline InitialConditions auto create_initconds(const Config& config, const GbxMaps& gbxmaps) {
      const auto initsupers = InitSupersFromBinary(config.get_initsupersfrombinary(), gbxmaps);
      const auto initgbxs = InitGbxsNull(gbxmaps.get_local_ngridboxes_hostcopy());
      return InitConds(initsupers, initgbxs);
    }

Builds the coupled dynamics (see :doc:`coupled_dynamics` for the alternatives) and the initial
conditions — the super-droplets' starting attributes/positions from a binary file
(``InitSupersFromBinary``), and the gridboxes' starting thermodynamic state (here
``InitGbxsNull``, meaning: don't set it here, it will come from the coupled dynamics'
``receive_dynamics`` on the very first step).

.. code-block:: cpp

    template <typename Dataset, typename Store>
    inline Observer auto create_observer(const Config& config, const Timesteps& tsteps,
                                         Dataset& dataset, Store& store) {
      // ... builds and combines ~9 observers, see :doc:`microphysics_and_observers`
    }

    template <typename Dataset, typename Store>
    inline auto create_sdm(const Config& config, const Timesteps& tsteps, Dataset& dataset, Store& store) {
      const auto couplstep = (unsigned int)tsteps.get_couplstep();
      const GridboxMaps auto gbxmaps = create_gbxmaps(config);
      const MicrophysicalProcess auto microphys = create_microphysics(config, tsteps);
      const MoveSupersInDomain movesupers = create_movement(config, tsteps, gbxmaps);
      const Observer auto obs = create_observer(config, tsteps, dataset, store);
      return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
    }

``create_sdm`` is the function that ties the previous four factory functions together into one
``SDMMethods`` object — this is exactly the object ``main()`` assigns to ``sdm``.

The Timestepping Loop Itself
--------------------------------

``main()`` hands ``sdm``, ``coupldyn`` and ``comms`` to a ``RunCLEO`` and calls it. Inside
``libs/runcleo/runcleo.hpp``, ``RunCLEO::operator()`` creates the initial gridboxes/super-droplets
and then runs:

.. code-block:: cpp

    int timestep_cleo(const unsigned int t_end, const dualview_gbx gbxs,
                      SupersInDomain& allsupers) const {
      unsigned int t_mdl(0);
      while (t_mdl <= t_end) {
        const auto t_next = start_step(t_mdl, gbxs, allsupers);   // 1. receive dynamics, observe
        coupldyn_step(t_mdl, t_next);                              // 2. advance coupled dynamics
        sdm_step(t_mdl, t_next, gbxs, allsupers);                  // 3. advance SDM (motion+microphysics)
        t_mdl = proceed_to_next_step(t_next, gbxs);                // 4. send dynamics back out
      }
    }

with ``t_next`` chosen adaptively each iteration as the sooner of the next coupling time and the
next observation time (``get_next_step``, mirroring the general pattern in
:doc:`../intro/timestepping`), and with ``sdm_step`` internally running its own nested
sub-timestepping loop (in ``SDMMethods::run_step``, ``libs/runcleo/sdmmethods.hpp``) so that
super-droplet motion and each microphysical process can each use their own, independently
adaptive sub-timestep, down to the granularity of ``t_next - t_mdl``.

.. figure:: ../intro/figs/timestepping.png
   :alt: Diagram for Cleo's timestepping routine
   :align: center

   The same routine-calling diagram from :doc:`../intro/timestepping`, now with each labelled
   step traceable to a specific function above: "start step" → ``start_step``, "run coupled
   dynamics" → ``coupldyn_step``, "run SDM" → ``sdm_step`` (itself nesting motion and
   microphysics sub-steps), "proceed" → ``proceed_to_next_step``.

Where Things Plug In — Summary
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - To change...
     - ...edit this factory function in ``main_impl.hpp``
   * - Domain geometry / resolution
     - ``create_gbxmaps`` (and the gridbox-boundaries input file, see :doc:`tutorial_new_example`)
   * - Super-droplet motion / terminal velocity / boundary conditions
     - ``create_motion``, ``create_boundary_conditions``, ``create_movement``
   * - Microphysics (condensation, collisions, breakup, rebound)
     - ``config_condensation``, ``config_collisions``, ``create_microphysics`` — see
       :doc:`microphysics_and_observers`
   * - Coupled dynamics (null / from-file / CVODE / YAC)
     - ``create_coupldyn``, and the ``CouplingComms<...>`` type used in ``main()`` — see
       :doc:`coupled_dynamics`
   * - Initial conditions
     - ``create_initconds``
   * - What gets written out (observers)
     - ``create_observer`` and its helpers (``create_superdrops_observer``,
       ``create_gridboxes_observer``, ...) — see :doc:`microphysics_and_observers`

What Typically Should *Not* Change
--------------------------------------

* The overall structure of ``main()`` (the six steps listed at the top of this page) — Cleo's
  timestepping and coupling logic lives in ``RunCLEO``/``SDMMethods``, not in ``main()``, and
  should not be reimplemented per-example.
* The ``Kokkos::initialize``/``{ ... }``/``Kokkos::finalize`` structure and its ordering relative
  to everything else.
* Signatures of the ``create_*`` functions, unless you are also changing what calls them — many
  of them are templated on ``Dataset``/``Store``/``GbxMaps`` specifically so they can be reused
  unchanged across examples with different concrete types.

Before/After: Swapping the Coupled Dynamics
----------------------------------------------

As a concrete "before/after" (condensed from :doc:`coupled_dynamics`'s worked example): switching
an example from reading dynamics from file to having no dynamics at all only requires changing
two things in ``main_impl.hpp``/``main.cpp``, everything else is untouched:

.. code-block:: cpp

    // before: one-way coupled to binary files
    inline CoupledDynamics auto create_coupldyn(const Config& config, const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep, const unsigned int t_end) {
      const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
      const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});
      const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);
      return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
    }
    // in main(): const CouplingComms<CartesianMaps, FromFileDynamics> auto comms = FromFileComms{};

.. code-block:: cpp

    // after: no coupled dynamics at all
    inline CoupledDynamics auto create_coupldyn(const Config& config, const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep, const unsigned int t_end) {
      return NullDynamics(couplstep);
    }
    // in main(): const CouplingComms<CartesianMaps, NullDynamics> auto comms = NullDynComms{};

You would also update your executable's ``CMakeLists.txt`` to link ``coupldyn_null`` instead of
``coupldyn_fromfile``, and (if you no longer need it) remove the ``fromfiledynamics`` block from
your YAML configuration file (see :doc:`tutorial_new_example`).

.. admonition:: Common pitfall

   Forgetting to update the ``CMakeLists.txt`` link line after changing which ``coupldyn_*``
   header you ``#include`` produces a *linker* error (undefined reference), not a compiler error —
   if you see undefined references to symbols that plainly exist in a header you're including,
   check which libraries your executable's ``target_link_libraries`` lists.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
