.. _tutorial_new_example:

Tutorial: Building a New Cleo Experiment From Scratch
======================================================

This tutorial walks through everything needed to create a new Cleo example, using one running
example throughout: a minimal 1-D "rainshaft" — a column of gridboxes stacked vertically, with
thermodynamics and wind prescribed from file (one-way coupled). In terms of SDM microphysics, this
example models condensation and collision-coalescence, and super-droplets fall according to their
terminal velocity. This setup mirrors (and is a simplified version of) the real
``examples/rainshaft1d`` example already in the repository.

This tutorial assumes you have already built Cleo at least once — see
:doc:`../usage/installation/installation` — and installed ``cleopy``
(``uv sync --extra examples``, per :doc:`../usage/installation/installation`).


Step 1: Plan the Simulation
------------------------------

Before writing any code, decide:

* **Domain**: how many spatial dimensions, how many gridboxes, what geometry do you want?
  This tutorial uses a 1-D (vertical only) stack of ``ngbxs`` gridboxes with evenly-spaced
  ``z`` boundaries.
* **Initial super-droplet conditions**: how many superdroplets do you want, and
  what initial distribution of their radii and positions? This tutorial initialises a fixed
  number of super-droplets per gridbox as dry aerosol throughout the entire domain, with radii
  sampled log-uniformly and multiplicities given according to a lognormal aerosol-size distribution.
* **Coupled dynamics**: how do you want to define the pressure/temperature/moisture/wind fields? See
  :ref:`coupled_dynamics_available_options`. This tutorial uses ``fromfile`` dynamics to prescribe
  a hydrostatic profiles with a sinusoidal updraught, generated before the simulation by a
  Python script rather than coming from a real dynamical core.
* **Microphysics**: which processes, and with what parameters do you want? See
  :doc:`microphysics_and_observers`. This tutorial models condensation/evaporation and
  collision-coalescence with a hydrodynamic kernel.
* **Observers**: what output do you need? See :doc:`microphysics_and_observers`. This tutorial
  outputs time, super-droplet attributes, and mass moments, amongst other things...

Now create a new directory under ``CLEO/experiments/[your_new_experiment]``, and within there
create ``src/config`` directories. Your folder structure should now mirror that of Cleo's examples,
like so:

::

    CLEO/
    ├── cleopy/
    ├── docs/
    ├── examples/
    ├── experiments/
    │   └── your_new_experiment/
    │       └── src/
    │           └── config/
    ├── extern/
    ├── libs/
    [...]

Now you can proceeed to create a ``config.yaml`` file in
``CLEO/experiments/[your_new_experiment]/src/config``...

.. _create_config:

Step 2: Create the Configuration File
-----------------------------------------

Cleo reads a single YAML file at startup into a ``Config`` object, itself built from two parts:

* ``RequiredConfigParams`` — fields with *no* default; your the YAML must set them.
* ``OptionalConfigParams`` — fields with a default value, only overridden if present in your YAML.
  Note that optional parameters are often only optional because not all setups require them,
  however, they may be required for your particular setup. For example, if you use ``fromfile``
  dynamics, you will need the ``coupled_dynamics`` block with ``type : fromfile`` and the subsequent
  list of binary files you intend to supply for the thermodynamics. Don't worry though, Cleo is
  expected to throw an error if you do not supply optional config parameters which your
  simulation setup requires - so your simulation should fail rather than run nonsense!

Your goal is to create a YAML file, e.g.
``CLEO/experiments/[your_new_experiment]/src/config/config.yaml``, which contains all the required
and optional parameters you need for your setup, as explained next.

Here is the list of required configuration parameters every Cleo config file *must* have:

.. code-block:: yaml

    ### SDM Runtime Parameters ###
    domain:
      nspacedims : X                                        # no. of spatial dimensions to model
      ngbxs : X                                             # total number of Gbxs
      maxnsupers: X                                         # maximum number of SDs

    timesteps:
      CONDTSTEP : X                                         # time between SD condensation [s]
      COLLTSTEP : X                                         # time between SD collision [s]
      MOTIONTSTEP : X                                       # time between SDM motion [s]
      COUPLTSTEP : X                                        # time between dynamic couplings [s]
      OBSTSTEP : X                                          # time between SDM observations [s]
      T_END : X                                             # time span of integration from 0s to T_END [s]

    ### Initialisation Parameters ###
    inputfiles:
      constants_filename : /path/to/libs/cleoconstants.hpp  # name of file for values of physical constants
      grid_filename : /path/to/dimlessGBxboundaries.dat     # binary filename for initialisation of GbxMaps

    ### Output Parameters ###
    outputdata:
      setup_filename : /path/for/output/setup.txt           # .txt filename to copy configuration to
      zarrbasedir : /path/for/output/sol.zarr               # zarr store base directory
      maxchunk : X (~2500000 is generally a good choice)    # maximum no. of elements in chunks of zarr store array


Here are the selected optional fields relevant to our 1-D column example (*hint:* have a look at
other examples or ``optional_config_params.hpp`` for the complete set of optional parameters, e.g.
breakup parameters, other coupled-dynamics settings, etc.):

.. code-block:: yaml

    ### Kokkos Initialization Parameters ###
    kokkos_settings: # (best to delete block if you don't want to control computational performance)
      num_threads: X                                      # number of threads for host parallel backend

    ### Initialisation Parameters ###
    inputfiles:
      [...] # (required parameters, shown above)
    initsupers:
      type: frombinary                                    # type of initialisation of super-droplets
      initsupers_filename : /path/to/dimlessSDsinit.dat   # binary filename for initialisation of SDs

    ### Microphysics Parameters ###
    microphysics:
      condensation:
        do_alter_thermo : true/false                      # true = cond/evap alters the thermodynamic state
        maxniters : X                                     # maximum no. iterations of Newton Raphson Method
        MINSUBTSTEP : X                                   # minimum subtimestep in cases of substepping [s]
        rtol : X                                          # relative tolerance for implicit Euler integration
        atol : X                                          # absolute tolerance for implicit Euler integration
      coalescence:
        constcoaleff:
          coaleff: X                                       # constant coalescence efficiency

    ### Coupled Dynamics Parameters ###
    coupled_dynamics:
      type : fromfile                                      # type of coupled dynamics to configure
      press : /path/to/dimlessthermo_press.dat             # binary filename for pressure
      temp : /path/to/dimlessthermo_temp.dat               # binary filename for temperature
      qvap : /path/to/dimlessthermo_qvap.dat               # binary filename for vapour mixing ratio
      qcond : /path/to/dimlessthermo_qcond.dat             # binary filename for liquid mixing ratio
      wvel : /path/to/dimlessthermo_wvel.dat               # binary filename for vertical (coord3) velocity

.. tip::
   Occasionally a simulation will fail due to the Newton Raphson method for condensation /
   evaporation, with an error such as ``No root converged upon within max number of iterations``.
   To fix this, in your settings for ``condensation`` you need to increase ``maxniters``, and/or
   ``rtol`` and/or ``atol``, and/or you need to reduce ``MINSUBTSTEP``. Be cautious
   though, too large tolerances will lead to inaccuracies in the solver, and too large ``maxniters``
   or too small  ``MINSUBTSTEP`` will be computationally expensive (slow!).

If you are running many simulations with similar config files, you don't have to write the YAML
entirely by hand every time. You can instead copy a draft config file (in Python:
``shutil.copy(src_config_filename, dest_config_filename)``) and then use ``cleopy``'s Python
module, ``editconfigfile``, to re-write the config parameters you want to change. Simply pass
a dictionary into ``edit_config_params`` with the ``name: value`` items for the parameters you
want to re-write
(i.e. ``editconfigfile.edit_config_params(dest_config_filename, config_params_dict)``). In
fact, this is what all of Cleo's example scripts are already doing - see e.g. the call
to ``edit_config_params`` in ``rainshaft1d.py``.

In Cleo's examples you may also notice an additional block in the YAML file

.. code-block:: yaml

    ### Initialization Parameters for rshaft1d_inputfiles.py ###
    python_inputfiles:
      [...]

This block is not used by the actual Cleo model itself, (i.e. not by any C++ code). It is
specific to each example to make it easier to create the input files for Cleo in Python,
as explained in step 4.


Step 3: Create your Cleo C++ executable
-----------------------------------------

Into your ``CLEO/experiments/your_new_experiment/src`` directory, you now need to create a
``main.cpp`` file and associated ``CMakeLists.txt`` which will compile the executable you want
for your experiment. The easiest thing to do is copy the ``main.cpp`` and ``CMakeLists.txt`` from
a Cleo example which is the most similar to the setup you have in mind
(e.g. ``examples/rainshaft1d/src/main_rshaft1d.cpp`` for a 1-D rainshaft).

Your folder strucutre should now look as follows:

::

    CLEO/
    ├── cleopy/
    ├── docs/
    ├── examples/
    ├── experiments/
    │   └── your_new_experiment/
    │       └── src/
    │           ├── config/
    │           |   └── config.yaml
    |           ├── CMakeLists.txt
    |           └── main.cpp
    ├── extern/
    ├── libs/
    [...]


Editing your ``main.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, as per the :ref:`main_cpp_summary`, adjust the following C++ functions to create the setup
you have in mind:

1. Adjust ``create_gbxmaps`` if your domain geometry differs (e.g. you don't want a cartesian domain).
2. Adjust ``create_movement`` if you want a different type of movement (e.g. for a non-cartesian
   domain), a different terminal velocity parameterisation, or different boundary conditions.
3. Adjust ``create_microphysics`` to model the microphysical processes you require
   (see :doc:`microphysics_and_observers`).
4. Adjust ``create_observer`` for what you want to output.
5. Adjust ``create_coupldyn`` and the ``CouplingComms<...>`` type in ``main()`` for your chosen
   coupled dynamics (see :doc:`coupled_dynamics`).
6. Adjust ``create_initconds`` if your initial conditions differ from the conventional setup with
   ``InitSupersFromBinary`` and ``InitGbxsNull``.

Editing your ``CMakeLists.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you're finished designing your ``main.cpp``, correct your ``CMakeLists.txt`` to match your
experiment by adapting the following lines:

.. code-block:: cmake

    # [...]
    project("[your_new_experiment]")
    # [...]
    set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)
    # [...]
    add_executable([your_executable] EXCLUDE_FROM_ALL "main.cpp")
    # [...]
    target_link_libraries([your_executable] PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
    target_link_libraries([your_executable] PUBLIC Kokkos::kokkos)
    target_include_directories([your_executable] PRIVATE "${CLEO_SOURCE_DIR}/libs")
    set_target_properties([your_executable] PROPERTIES CMAKE_CXX_STANDARD_REQUIRED ON
                                                       CMAKE_CXX_EXTENSIONS ON CXX_STANDARD 20)

Note: if your ``main.cpp`` includes a different domain geometry or coupled dynamics, you
will need to link to those libraries instead of ``coupldyn_fromfile`` and ``cartesiandomain``.

Finally, add ``add_subdirectory(experiments/your_new_experiment)`` to ``CLEO/CMakeLists.txt``
(e.g. after ``add_subdirectory(libs)``) so Cleo's top-level build picks up your new executable.

.. _create_input:

Step 4: Create the Input Files
----------------------------------

Your executable and config file are now ready to compile and run, however you need to provide Cleo
input data for your experiment. The input data Cleo requires is for:

I. your setup's grid,
II. your initial super-droplet conditions,
III. your thermodynamical fields, if you're using ``fromfile`` dynamics.

Cleo does not read NetCDF or YAML for this input data, they are all read
from Cleo's own simple flat binary ``.dat`` files. These files can be generated with ``cleopy``
in Python.

Your next step therefore it to create a Python script e.g.
``CLEO/experiments/your_new_experiment/create_inputfiles.py``, which uses the ``cleopy`` modules to
generate the input data you want for your experiment.

Cleo's examples show how different grids / initial superdroplet conditions / ``fromfile``
thermodynamics can be generated, have a look at the ``XXX_inputfiles.py`` scripts for inspiration.
Alternatively, all options for creating grids are shown in ``create_gbxboundariesbinary_script.py``,
all options for creating initial superdroplet conditions are shown in
``create_initsuperdropsbinary_script.py``, and all options for creating thermodynamic conditions
(for setups using ``fromfile`` coupled dynamics) are shown in ``create_thermobinaries_script.py``.

I. Creating Cleo's Grid
~~~~~~~~~~~~~~~~~~~~~~~~~

Cleo's grid is created by reading a "gridbox boundaries" ``.dat`` binary file.
Note that despite its name, the gridbox boundaries file does not initialise gridboxes themselves,
rather it initialises the GridboxMaps of your simulation, to define Cleo's grid.

To create the ``.dat`` binary grid file, you goal is to run:

.. code-block:: python

    from cleopy import geninitconds

    [...]

    geninitconds.generate_gridbox_boundaries(
            /path/to/dimlessGBxboundaries.dat,
            zgrid,
            xgrid,
            ygrid,
            /path/to/libs/cleoconstants.hpp,
            isprintinfo=True/False,
            isfigures=[True/False, True/False],
            savefigpath=/path/to/folder/to/save/figures/,
            savelabel="optional_label_for_saved_figures",
        )

You thus need to create 3 python arrays or lists for ``zgrid``, ``xgrid``, and ``ygrid``.
See Cleo's examples or the script ``create_gbxboundariesbinary_script.py`` for how this can be done.

Note the call to ``generate_gridbox_boundaries`` will de-dimensionalise the ``zgrid``,
``xgrid``, and ``ygrid`` data you provide in meters, hence the name of the file is
often ``dimlessGBxboundaries.dat``.


II. Creating Cleo's Initial Superdroplet Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cleo's superdroplets are initialised by reading a "SDs init" ``.dat`` binary file.
To create this ``.dat`` binary file, you goal is to run:

.. code-block:: python

    from cleopy import geninitconds

    [...]

    geninitconds.generate_initial_superdroplet_conditions(
            initattrsgen,
            /path/to/dimlessSDsinit.dat,
            /path/to/your/config.yaml,
            /path/to/libs/cleoconstants.hpp,
            /path/to/dimlessGBxboundaries.dat,
            nsupers,
            numconc,
            numconc_tolerance=X,
            isprintinfo=True/False,
            isfigures=[True/False, True/False],
            savefigpath=/path/to/folder/to/save/figures/,
            gbxs2plt=[gbxindexes_of_superdroplets_to_plot],
            savelabel="optional_label_for_saved_figures",
        )

Note the call to ``generate_initial_superdroplet_conditions`` will de-dimensionalise the
superdroplet attributes you provide, hence the name of the file is
often ``dimlessSDsinit.dat``.

This function requires 3 key inputs ``initattrsgen``, ``nsupers`` and ``numconc``:

* ``initattrsgen`` is the attribute generator, it consists of 5 inidividual generators:
    ``radiigen``, ``dryradiigen``, ``coord3gen``, ``coord1gen``, ``coord2gen``, and a probability
    distirbution for droplet sizes, ``xiprobdist``, which is used togetehr with ``radiigen`` to
    generate each superdroplet's radius and multiplicity. See Cleo's examples, or the script
    ``create_initsuperdropsbinary_script.py`` for what generators are available and how they
    can be used (or create one yourself!).
* ``nsupers`` is the number of superdroplets per gridbox, it can either be a single number, or a
    dictionary of ``gbxindex: nsupers`` items for the number of superdroplets in each gridbox. You
    can also use the function ``nsupers_at_domain_top`` or ``nsupers_at_domain_base``from
    ``cleopy.crdgens`` to generate superdroplets only at the top/bottom of the
    domain (see e.g. the ``rainshaft1d`` example).
* ``numconc`` is the number concentration of droplets (m^-3) in each gridbox. currently it can only
    be a single number (i.e. the same for all gridboxes), but there is also the option to pass
    ``xi_by_pressure=True`` into the creation of ``initattrsgen`` to allow the number concentration
    in each gridbox to be consistent with the change in air density with height in a domain
    (see e.g. the ``constthermo2d`` example). The optional setting ``numconc_tolerance`` can be used
    to allow the multiplicities of the superdroplets in each gridbox to not exactly add to the
    desired numconc*gridbox_volume; this allows you a bit of rounding error if necesary.


III. Creating Thermodynamics Conditions for ``fromfile`` Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are using ``fromfile`` coupled dynamics, you need to create a ``.dat`` binary file for each
of thermodynamic field you are simulating. In a 0-D setup, these are pressure, temperature, qvap,
and qcond (``dimlessthermo_press.dat``, ``_temp.dat``, ``_qvap.dat``, and ``_qcond.dat``
respectively``). In a 1-D setup, you additionally need wvel (``_wvel.dat``), in a 2-D setup,
wvel and uvel (``_uvel.dat``), and in a 3-D setup, wvel, uvel, and vvel (``_vvel.dat``).

To create these ``.dat`` binary files, you goal is to run:

.. code-block:: python

    from cleopy import geninitconds

    [...]

    geninitconds.generate_thermodynamics_conditions_fromfile(
            /path/to/dimlessthermo.dat,
            thermodyngen,
            /path/to/your/config.yaml,
            /path/to/libs/cleoconstants.hpp,
            /path/to/dimlessGBxboundaries.dat,
            isfigures=[True/False, True/False],
            savefigpath=/path/to/folder/to/save/figures/,
            savelabel="optional_label_for_saved_figures",
            press_ref=None, # optional setting used only for plotting potential temperature
        )

Note for the function ``generate_thermodynamics_conditions_fromfile`` you pass a single filename,
``/path/to/dimlessthermo.dat``, to name all the thermodynamics files. The internals of this
function will split the path at ``.dat`` and append the relevant ``_press.dat``, ``_temp.dat`` etc.
to write out each thermodynamics file which your ``thermodyngen`` generator creates.

This function has one key argument, ``thermodyngen``, which is analogous to ``initattrsgen``.
``thermodyngen`` is a generator of thermodynamics and consists of two generators, a
``thermog`` (or "thermo-generator") for pressure, temperature, qvap, and qcond fields, and a
``windsg`` (or "winds-generator") for wvel, uvel and vvel fields. See Cleo's examples which
use ``fromfile`` dynamics or browse the script ``create_thermobinaries_script.py`` for the
different generators that are available and how to use them (or create one yourself!).

That's it, that covers the creation of all of Cleo's input files!

.. tip::
   To avoid buggy code, it's a good idea to read settings/arguments for creating the binary files
   (for your grid, initial superdroplets, and potentially thermodynamics) from a Cleo config YAML
   file rather than hard-coding variables in a Python script (and therefore having to ensure
   yourself that the settings/files which Cleo reads from the YAML config stay consistent with
   what your Python script is doing). Indeed, this is how all of the ``XXX_inputfiles.py`` scripts
   in Cleo's examples work, for example to ensure the ``nsupers`` passed to
   ``generate_initial_superdroplet_conditions`` is consisent with the ``maxnsupers`` setting in the
   YAML config which the Cleo executable reads in.


.. _building_cleo:

Step 5: Build, Compile and Run Cleo
-------------------------------------

Use Cleo's scripts
~~~~~~~~~~~~~~~~~~~

You are now ready to build, compile and sucessfully run your Cleo experiment!

To build and compile your executable, it's easiest to add your experiment to Cleo's bash
scripts in ``CLEO/scripts``, like what's done for Cleo's examples. Then the handling of the
compilers, modules, build flags, etc. for your computer architecture is sorted out by the
scripts automatically. Very roughly speaking, those scripts eventually call something like:

.. code-block:: console

  $ cmake -S . -B ./build -DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian [YAC_flags, Kokkos_flags, ...]
  $ cmake --build ./build --target your_executable

To run your experiment, you then need to execute

.. code-block:: console

  $ ./build/experiments/your_new_experiment/your_executable /path/to/your/config.yaml

Where ``your_executable`` matches the name you used for your exectuable in your ``CMakeLists.txt``.

Rather than actually run Cleo from a bash script (``./[...]/your_executable /path/to/your/config.yaml``),
and therefore have to ensure yourself that you load the relevant run-time modules/settings for
your computer architecture, you might want to write a Python driver script like in Cleo's examples
(e.g. ``rainshaft1d.py``), which includes not only running the executable (this step) but also
generating the input data (step 4) and post-processing/plotting it (step 6) (with control over
each steup via the ``--do_inputfiles`` / ``--do_run_executable`` / ``--do_plot_results`` flags).
You would then use the bash scripts in exactly the same way as the examples, namely to build and
compile Cleo and then to launch your python driver script.


.. _cleo_cmake_flags:

Cleo's CMake flags
~~~~~~~~~~~~~~~~~~~

During the build configuration step, ``cmake -S . -B ./build [...]``, you can set various CMake
flags, e.g. with ``-D<VARIABLE>=<value>``, to have more control over which of Cleo's optional
libraries, examples, and roughpaper targets get configured. The following flags are available:

- ``CLEO_COUPLED_DYNAMICS``: Selects which :doc:`coupled dynamics <../deep_dive/coupled_dynamics>`
  library/libraries to build. Options are ``null``, ``fromfile``, ``cvode``, ``yac``, ``numpy``, or
  ``all``. The default (leaving the flag unset), builds both ``cvode`` and ``fromfile``.
- ``CLEO_DOMAIN``: Selects domain geometry. Options are ``cartesian``. The default builds all
  options (currently only ``cartesian``).
- ``CLEO_NO_PYBINDINGS``: if true, skips configuring Cleo's Python bindings. Default is false.
- ``CLEO_NO_EXAMPLES``: if true, skips configuring Cleo's ``examples`` targets. Default is false.
- ``CLEO_NO_ROUGHPAPER``: if true, skips configuring Cleo's ``roughpaper`` targets. Default is false.

Additionally...

- You will also want to pass ``CLEO_CMAKE_BUILD_TYPE="Release"`` to benefit from performance
  optimisations rather than debugging. The default is otherwise ``CLEO_CMAKE_BUILD_TYPE="Debug"``.
- You must also pass various flags for Kokkos and YAC. See Cleo's bash scripts for the flags
  Cleo currently uses.
- Note: Cleo's configuration also requires ``MPI`` (``find_package(MPI REQUIRED COMPONENTS C)``),
  ``YAC`` (``find_package(YAC REQUIRED)``), and ``CMAKE_CXX_STANDARD`` to ``20``.

.. note::
   Several of Cleo's subdirectories are added with ``EXCLUDE_FROM_ALL`` (e.g. the examples, the
   domain and coupled-dynamics libraries, the Python bindings, ...) — this means even if your
   configuration includes them, CMake will only build those targets if something explicitly depends
   on them (e.g. in your executable's ``target_link_libraries()``). They will not be built
   automatically as part of a plain ``make``/``cmake --build``. You don't need to do anything about
   this other than make sure every target you link against is one you've actually named.



Step 6: Post-Process / Plot the Output
--------------------------------------------

Cleo writes output as to a Zarr store, readable with the ``xarray`` Python package
(``xr.open_dataset(/path/to/your/dataset.zarr, engine="zarr")``) :cite:`bayley2026a`.
However, the raw superdroplet data is stored in the
`CF convention's ragged array format <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html#_contiguous_ragged_array_representation>`_
which is pretty difficult to handle by yourself. Instead of tackling superdroplet data directly with
``xarray`` or ``numpy``, it's best to get familiar with the ``awkward`` Python package
`Awkward Array <https://awkward-array.org/doc/main/index.html>`_
which is designed exactly for handling ragged arrays.

To make things significantly easier, albeit sometimes slower, ``cleopy`` also has a module,
``sdmout_src``, which wraps various Python classes around, and sometimes re-shapes, your raw
dataset. In addition, there is the additional module  ``plotcleo`` which is used by Cleo's examples
to plot Cleo's output after it's been post-processed using ``sdmout_src``. ``plotcleo`` is not
available through ``PyPI``, but it can easily be installed locally (see
:doc:`../usage/examples/examples`). Exploring how the ``_plotting.py`` scripts handle and plot
Cleo's output may help you work with your dataset yourself.


Some Common Errors
--------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Step
     - Likely cause and fix
   * - Config file
     - ``std::invalid_argument`` or YAML parse errors usually mean a required field
       (``required_config_params.hpp``) is missing, or a path is wrong (this is when absolute paths
       are generally easier to debug that relative ones) — check every path in the YAML resolves
       from wherever you run the executable, not from the config file's own directory.
   * - Input files
     - A "wrong file size" or similar exception when Cleo reads a binary file almost always means
       a mismatch between what ``cleopy`` generated and what ```ngbxs`` or ``maxnsupers`` in
       the YAML say to expect — regenerate both from the same source of truth (a single Python
       config dict, as Cleo's examples do) rather than editing your YAML and python input files
       script independently - see the tip at the end of :ref:`create_input` for how this can be
       done.
   * - Wiring in ``main_impl.hpp``
     - A very long compiler error, usually naming a ``concept`` and ``requires`` clause, means
       some type you passed doesn't satisfy the concept expected at that point (e.g. a
       ``CouplingComms`` that doesn't match your chosen ``CoupledDynamics`` type — see the
       pitfall in :doc:`coupled_dynamics`). Read the *first* concept-related error in the output;
       later ones are often knock-on effects of the first.
   * - Build/CMake
     - "undefined reference" at link time usually means a needed ``coupldyn_*``/``cartesiandomain``
       library isn't listed in your executable's ``target_link_libraries`` (in your
       ``CMakeLists.txt``), even though the header compiled fine (see the pitfall in
       :doc:`main_cpp`). You may also have forgotten to ``#include`` at the top of your ``main.cpp``
       all of the ``.hpp`` files for the structures / classes / functions etc. that you then use.
   * - Runtime
     - A silent hang or a Kokkos-related abort often means the coupling timestep of your
       coupled-dynamics object and your SDM disagree — ``RunCLEO``'s constructor checks this
       (``check_coupling()`` in ``libs/runcleo/runcleo.hpp``) and throws
       ``"coupling timestep of dynamics solver and CLEO SDM are not equal"`` if so; make sure
       ``COUPLTSTEP`` in your YAML and however your ``create_coupldyn`` computes its coupling
       step agree.
   * - Runtime
     - A simulation fails due to the Newton Raphson method for condensation / evaporation, with
       an error such as ``No root converged upon within max number of iterations``, see the tip in
       :ref:`create_config` for a suggested solution.
   * - Output
     - Missing or empty output variables usually mean the corresponding observer wasn't included
       in your ``create_observer`` combination (:doc:`microphysics_and_observers`) — check the
       ``>>``-chain includes it.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
