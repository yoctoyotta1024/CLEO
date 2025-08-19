Quickstart
==========

Start by following the steps to :doc:`install CLEO <installation>`.

To use CLEO as an SDM coupled to a dynamics solver, essentially your task is to contruct a
``main.cpp`` file.

Have a look at ``roughpaper/src/main.cpp`` and ``roughpaper/src/main_impl.hpp`` for inspiration.
Depending on the setup you desire, you need to include different instantiations of the core concepts
of CLEO in your main.cpp. For example for the SDM part of CLEO you need to specify the coupling
timestep, the gridbox maps for the domain, and the type of super-droplet motion, microphysics and
observer.

Once you’ve done that, your should have constructed all the necesary ingredients for CLEO so that
the end of your ``main.cpp`` contains a main function that looks something like this:

.. code-block:: c++

  int main(int argc, char *argv[]) {
    if (argc < 2) {
      throw std::invalid_argument("configuration file(s) not specified");
    }

    /* Read input parameters from configuration file(s) */
    const auto config_filename = std::filesystem::path(argv[1]);  // path to configuration file
    const auto config = Config(config_filename);

    /* Initialize Communicator here */
    init_communicator init_comm(argc, argv, config);

    /* Initialise Kokkos parallel environment */
    Kokkos::initialize(config.get_kokkos_initialization_settings());
    {
      Kokkos::print_configuration(std::cout);

      /* Create timestepping parameters from configuration */
      const auto tsteps = Timesteps(config.get_timesteps());

      /* Create Xarray dataset wit Zarr backend for writing output data to a store */
      auto store = FSStore(config.get_zarrbasedir());
      auto dataset = SimpleDataset(store);

      /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
      const SDMMethods sdm = create_sdm(config, tsteps, dataset);

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

    return 0;
  }

where for example the SDM is set-up as:

.. code-block:: c++

  template <typename Dataset, typename Store>
  inline auto create_sdm(const Config &config, const Timesteps &tsteps,
                         Dataset &dataset, Store &store) {
    const auto couplstep = (unsigned int)tsteps.get_couplstep();
    const GridboxMaps auto gbxmaps = create_gbxmaps(config);
    const MicrophysicalProcess auto microphys = create_microphysics(config, tsteps);
    const MoveSupersInDomain movesupers = create_movement(config, tsteps, gbxmaps);
    const Observer auto obs = create_observer(config, tsteps, dataset);

    return SDMMethods(couplstep, gbxmaps, microphys, movesupers, obs);
  }

and for example the ``create_microphysics`` function sets up the microphysics as:

.. code-block:: c++

  inline MicrophysicalProcess auto create_microphysics(const Config &config,
                                                      const Timesteps &tsteps) {
    const auto c = config.get_condensation();
    const MicrophysicalProcess auto cond =
      Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                  c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);

    const PairProbability auto coalprob = LongHydroProb(1.0);
    const MicrophysicalProcess auto colls = CollCoal(tsteps.get_collstep(), &step2realtime, coalprob);

    return colls >> cond;
  }

Have a look at ``roughpaper/src/main.cpp`` and ``roughpaper/src/main_impl.hpp`` for inspiration.
