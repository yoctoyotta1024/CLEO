Quickstart
==========

Start by following the steps to :doc:`install CLEO <installation>`.

To use CLEO as an SDM coupled to a dynamics solver, essentially your
task is to contruct a ``main.cpp`` file.

Have a look at ``src/main.cpp`` for inspiration. Depending on the
setup you desire, you need to include different instantiations of
the core concepts of CLEO in your ``main.cpp``. For example you
need to specify the type of Super-Droplet motion, microphysics,
dynamics solver and coupling.

Once you've done that, your should have constructed all the
necesary ingredients for CLEO so that the end of
your ``main.cpp`` can look something like this:

.. code-block:: c++

  /* Initialise Kokkos parallel environment */
  Kokkos::initialize(argc, argv);
  {
    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, fsstore));

    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, t_end);
  }
  Kokkos::finalize();
