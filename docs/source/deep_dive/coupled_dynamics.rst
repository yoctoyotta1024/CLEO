.. _coupled_dynamics_page:

Coupled Dynamics
==================

What "Coupled Dynamics" Means
--------------------------------

Cleo models super-droplet microphysics and super-droplet transport, but it does **not** itself
solve for how pressure, temperature, vapour/liquid mixing ratio and wind evolve in time, i.e. it is
not responsible for advecting thermodynamics, nor parameterising radiation nor turbulence — that is
the job of whatever Cleo is *coupled* to. Cleo calls this "external" source of thermodynamic and
wind fields its ``CoupledDynamics``, and the bridge between the coupled dynamics and Cleo's SDM,
is called the ``CouplingComms``. The coupling comms communicates fields either to Cleo
(in a one-way coupling) or both to and from Cleo (in a two-way coupling) . See the programming
guide's :ref:`dynamics_coupling` page, or B2026a :cite:`bayley2026a`, Sect. 2.4, for more detail
on this coupling framework.

The Two C++ Concepts Involved
-------------------------------

Two C++ concepts, both in ``libs/runcleo/``, define what's required for the ``CoupledDynamics``
and ``CouplingComms`` respectively.

The ``CoupledDynamics`` concept looks like this:

.. code-block:: cpp

    template <typename CD>
    concept CoupledDynamics = requires(CD cd, unsigned int t) {
      { cd.prepare_to_timestep() } -> std::same_as<void>;
      { cd.get_couplstep() } -> std::convertible_to<unsigned int>;
      { cd.on_step(t) } -> std::same_as<bool>;
      { cd.run_step(t, t) } -> std::same_as<void>;
    };

So, any type which satisfies these requirements is a valid ``CoupledDynamics`` — regardless of what
physics it actually implements. ``prepare_to_timestep`` allows you to configure your dynamics
model before timestepping begins but also after its construction; ``get_couplstep`` returns the
time interval for the coupling to/from Cleo; ``on_step`` says whether the current
timestep will run the coupled dynamics; and ``run_step`` allows you to run your dynamics model from
Cleo's time-stepping loop (although you may not need that). The requirements that make up the
concept of a ``CoupledDynamics`` ensure that any type that acts as coupled dynamics for Cleo
will fit into the timestepping routines in ``RunCleo``, as depicted in :doc:`../intro/timestepping`.

The ``CouplingComms`` concept looks slightly more involved because it is templated on *both* the
``GridboxMaps`` type (``GbxMaps``) and the ``CoupledDynamics`` type (``CD``) of your particular
setup (because, in general, to exchange a field between a specific dynamics solver and Cleo, the
coupler needs to know both the geometry of Cleo's domain and what the solver is). Overlooking the
templating though, you can see the ``CouplingComms`` concept actually only consists of two simple
functions:

.. code-block:: cpp

    template <typename Comms, typename GbxMaps, typename CD>
    concept CouplingComms = requires(Comms s, GbxMaps& gbxmaps, CD& coupldyn, viewh_gbx h_gbxs,
                                    dualview_constgbx gbxs, SupersInDomain& allsupers) {
      {
        s.template send_dynamics<GbxMaps, CD>(gbxmaps, gbxs, allsupers, coupldyn)
      } -> std::same_as<void>;
      { s.template receive_dynamics<GbxMaps, CD>(gbxmaps, coupldyn, h_gbxs) } -> std::same_as<void>;
    };

``receive_dynamics`` puts the dynamics solver's current fields into Cleo's gridboxes;
``send_dynamics`` sends Cleo's gridbox state, or change in state, to the dynamics solver. One-way
coupling is simply the case where ``send_dynamics`` does nothing. (Alternatively the dynamics
solver could ignore the data incoming from ``send_dynamics``, but this is not Cleo's responsibiltiy
to control nor know about.)


The Time-Ordering of the Coupling
----------------------------------

On each *coupling* timestep, the communication of thermodynamics *from* the coupled dynamics *to*
Cleo occurs at the start of the step. This is followed by the advance of the timestep of the
coupled dynamics, and then that of SDM (although this can in principle be done synchronously). At
the end of a *coupling* time-step, the communication *from* Cleo *to* the coupled dynamics occurs.
Note that in-between coupling timesteps, Cleo and the coupled dynamics may advance sub-timesteps.
To understand this time-stepping sequence in more detail, refer to
:ref:`main_cpp_timestepping_runcleo`.


.. _coupled_dynamics_available_options:

The Available Options
------------------------

Each option for the actual ``CoupledDynamics`` type you may want to use has its
own ``libs/coupldyn_*`` directory, alongside its associated ``CouplingComms`` type
(see :doc:`code_structure`), and is selected during the configuration of your Cleo build via
the CMake flag: ``-DCLEO_COUPLED_DYNAMICS=<option>``.

Here is a summary of the available options:

.. list-table::
   :header-rows: 1
   :widths: 16 14 70

   * - Directory → ``CoupledDynamics``, ``CouplingComms``
     - CMake flag
     - What it does
   * - ``coupldyn_null`` → ``NullDynamics``, ``NullDynComms``
     - ``null``
     - No dynamics and no coupling at all. Gridboxes keep
       whatever thermodynamic state they were initialised with. Useful for box-model microphysics
       tests where advection/dynamics are irrelevant (e.g. the collision-box-model examples).
   * - ``coupldyn_fromfile`` → ``FromFileDynamics``, ``FromFileComms``
     - ``fromfile``
     - One-way coupling: pressure, temperature, vapour/liquid mixing ratio and wind are read from
       pre-computed binary files (see :ref:`create_input`), one value per gridbox
       per coupling step. Cleo's own effect on these fields is never sent back — this is
       "piggybacking" a pre-existing dynamical solution. Used by e.g. the ``rainshaft1d`` and
       ``constthermo2d`` examples.
   * - ``coupldyn_cvode`` → ``CvodeDynamics``, ``CouplingComms``
     - ``cvode``
     - Two-way coupling to a 0-D adiabatic-parcel ODE model integrated with SUNDIALS' CVODE
       solver (see :doc:`../usage/installation/extern`). Used by the ``adiabaticparcel``
       examples, where condensation's latent heating/moistening genuinely feeds back into the
       parcel's temperature and vapour content.
   * - ``coupldyn_yac`` → ``YacCartesianDynamics``, ``YacComms``
     - ``yac``
     - Coupling via YAC (see :doc:`../usage/installation/extern`) to an external dynamical
       core (e.g. ICON), for large-eddy-simulations one-way or two-way coupled — the
       motivation for this setup is explained in B2026a :cite:`bayley2026a`, Sect. 1. Used by Cleo's
       ``bubble`` example, it is currently only compatible with Cleo's ``cartesiandomain`` grid
       and ICON's 100km, 5km resolution toroidal grid.
   * - (for python bindings) ``coupldyn_numpy`` → ``NumpyDynamics`` ``NumpyComms``
     - ``numpy``
     - Nested in ``cleo_python_bindings``, ``coupldyn_numpy`` contains a coupled dynamics type and
       coupler such that when Cleo is exposed to Python (via pybind11 python bindings;
       see :doc:`../usage/installation/extern`), a numpy-based coupled model in Python can drive it
       directly -- sending/receiving fields as numpy arrays to/from Cleo's C++ code.
       See ``examples/python_bindings``.

Note also that multiple options can be built
into the same binary (``-DCLEO_COUPLED_DYNAMICS=all``) — the choice of *which one* is instantiated
is still made in your ``main.cpp`` at compile time (see below), not at runtime.

Worked Example: Switching Coupled Dynamics
----------------------------------------------

As the top-level
~~~~~~~~~~~~~~~~

The *only* place a switch between coupled-dynamics options needs to happen is in your ``main.cpp``.

E.g. Consider you are using ``Null`` dynamics (like the ``boxmodelcollisions`` examples),
you could have something like:

.. code-block:: cpp

    #include "coupldyn_null/nulldynamics.hpp"
    #include "coupldyn_null/nulldyncomms.hpp"

    /* [...] */

    inline CoupledDynamics auto create_coupldyn(const Config& config,
                                                const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep,
                                                const unsigned int t_end) {
      return NullDynamics(couplstep);
    }

    /* [...] */

    int main(int argc, char* argv[]) {
    /* [...] */

    /* Solver of dynamics coupled to CLEO SDM */
    CoupledDynamics auto coupldyn =
        create_coupldyn(config, sdm.gbxmaps, tsteps.get_couplstep(), tsteps.get_t_end());

    /* coupling between coupldyn and SDM */
    const CouplingComms<CartesianMaps, NullDynamics> auto comms = NullDynComms{};

    /* [...] */
    }

To change to using a different dynamics, you need only change the libraries you include, the
return of ``create_coupldyn`` and your ``const CouplingComms<...>`` object. Everywhere else in
your ``main.cpp`` and in the rest of Cleo, the ``coupldyn`` and ``comms`` which you have made are
refered to only through the ``CoupledDynamics`` and ``CouplingComms`` concepts
(e.g. as shown in the above code snippet: ``CoupledDynamics auto coupldyn = create_coupldyn(...)``).
So no other code needs to change when you swap your coupled dynamics -- the compiler will figure out
the rest!

E.g. If you instead used ``FromFile`` (like the ``constthermo2d`` example), you would write:

.. code-block:: cpp

    #include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"
    #include "coupldyn_fromfile/fromfilecomms.hpp"

    /* [...] */

    inline CoupledDynamics auto create_coupldyn(const Config& config,
                                                const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep,
                                                const unsigned int t_end) {
      const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
      const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});
      const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

      return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
    }

    /* [...] */

    int main(int argc, char* argv[]) {
    /* [...] */

    /* Solver of dynamics coupled to CLEO SDM */
    CoupledDynamics auto coupldyn =
        create_coupldyn(config, sdm.gbxmaps, tsteps.get_couplstep(), tsteps.get_t_end());

    /* coupling between coupldyn and SDM */
    const CouplingComms<CartesianMaps, FromFileDynamics> auto comms = FromFileComms{};

    /* [...] */
    }

.. admonition:: Common pitfall
   :class: warning

   Note also the concept ``CouplingComms<...>`` type must match *both* your chosen ``GridboxMaps``
   type and your chosen ``CoupledDynamics`` type — for example
   ``CouplingComms<...> auto comms = FromFileComms{}``
   only satisfies ``CouplingComms<CartesianMaps, FromFileDynamics>``, not
   ``CouplingComms<CartesianMaps, NullDynamics>``. The compiler will refuse to build if the type
   of ``comms`` is inconsistent with the concept, which is of course the concept doing its job, but
   the resulting error message can be long -- look for which concept's ``requires`` clause failed.

*Note:* Whilst it's true that nothing else in your C++ code needs to change, you will also have to
edit which ``coupldyn_*`` library you link against in your executable's ``CMakeLists.txt``,
and which CMake ``CLEO_COUPLED_DYNAMICS`` value you configure your build with. You may also need to
change your YAML model configuration and generate different input files to run Cleo with a different
coupled dynamics. See :doc:`main_cpp` for where ``create_coupldyn`` and the ``CouplingComms``
instance fit into the rest of ``main.cpp``, and :doc:`tutorial_new_example` for how to build and run
an executable yourself from scratch.

**To summarise...**

When you switch your coupled dynamics, ``coupldyn``, remember to also switch the
coupled comms, ``comms``, and which libraries you ``#include`` at the top of your ``main.cpp``.
You also need to remember to switch which ``coupldyn_*`` library you link against in your
executable's ``CMakeLists.txt``, and which CMake ``CLEO_COUPLED_DYNAMICS`` value you configure
your build with. To run your new setup, you will also need to change your configuration and
input files accordingly.


In a bit more detail
~~~~~~~~~~~~~~~~~~~~~~

Compare the ``NullDynamics`` type with the ``FromFileDynamics`` type.
With ``coupldyn_null``, everything is trivial:

.. code-block:: cpp

    /* libs/coupldyn_null/nulldynamics.hpp */
    struct NullDynamics { // satisfies CoupledDynamics, does nothing
     private:
      const unsigned int interval;

     public:
      explicit NullDynamics(const unsigned int couplstep) : interval(couplstep) {}
      void prepare_to_timestep() const {}
      unsigned int get_couplstep() const { return interval; }
      bool on_step(const unsigned int t_mdl) const { return false; }
      void run_step(const unsigned int t_mdl, const unsigned int t_next) const {}
    };

    /* libs/coupldyn_null/nulldyncomms.hpp */
    struct NullDynComms {  // satisfies CouplingComms, does nothing
      template <typename GbxMaps, typename CD = NullDynComms>
      KOKKOS_INLINE_FUNCTION void receive_dynamics(const GbxMaps&, const NullDynamics&,
                                                   const viewh_gbx) const {}
      template <typename GbxMaps, typename CD = NullDynComms>
      KOKKOS_INLINE_FUNCTION void send_dynamics(const GbxMaps&, const viewh_constgbx,
                                                NullDynamics&) const {}
    };

With ``coupldyn_fromfile``, ``FromFileDynamics`` instead owns a ``CartesianDynamics`` object
holding vectors of pressure/temperature/vapour/liquid/wind read from binary files, and
``run_step`` advances its read-position by one coupling step each time it's called:

.. code-block:: cpp

    /* libs/coupldyn_fromfile/fromfile_cartesian_dynamics.hpp */
    struct FromFileDynamics { // satisfies CoupledDynamics, reads vectors initialised from files
     private:
      const unsigned int interval;
      std::shared_ptr<CartesianDynamics> dynvars;

      void run_dynamics(const unsigned int t_mdl) const { dynvars->increment_position(); }

     public:
      FromFileDynamics(const OptionalConfigParams::FromFileDynamicsParams& config,
                       const unsigned int couplstep, const std::array<size_t, 3> ndims,
                       const unsigned int nsteps)
          : interval(couplstep), dynvars(std::make_shared<CartesianDynamics>(config, ndims, nsteps)) {}

      auto get_couplstep() const { return interval; }
      void prepare_to_timestep() const {}
      bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }
      void run_step(const unsigned int t_mdl, const unsigned int t_next) const {
        if (on_step(t_mdl)) { run_dynamics(t_mdl); }
      }

      double get_press(const size_t ii) const { return dynvars->get_press(ii); }
      // ... get_temp, get_qvap, get_qcond, get_wvel, get_uvel, get_vvel similarly
    };

    /* libs/coupldyn_fromfile/fromfilecomms.hpp */
    struct FromFileComms { // satisfies CouplingComms, one-way coupling to Cleo's gridboxes
    private:
      void update_gridbox_state(const FromFileDynamics& ffdyn, const size_t ii, Gridbox& gbx) const;

    public:
      template <typename GbxMaps, typename CD = FromFileDynamics>
      void send_dynamics(const GbxMaps& gbxmaps, const dualview_constgbx gbxs,
                        const SupersInDomain& allsupers, FromFileDynamics& ffdyn) const {}

      template <typename GbxMaps, typename CD = FromFileDynamics>
      void receive_dynamics(const GbxMaps& gbxmaps, const FromFileDynamics& ffdyn,
                            const viewh_gbx h_gbxs) const  {
        const size_t ngbxs(h_gbxs.extent(0));
        Kokkos::parallel_for(
            "receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs), [=, *this](const size_t ii) {
              update_gridbox_state(ffdyn, gbxmaps.local_to_global_gridbox_index(ii), h_gbxs(ii));
            });
      }
    };

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
