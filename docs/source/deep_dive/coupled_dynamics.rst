.. _coupled_dynamics_page:

Coupled Dynamics
==================

What "Coupled Dynamics" Means
--------------------------------

Cleo models super-droplet microphysics and super-droplet transport, but it does **not** itself
solve for how pressure, temperature, vapour/liquid mixing ratio and wind evolve in time — that is
the job of whatever it is *coupled* to. Cleo calls this external (or internal, in the two-way
case) source of thermodynamic and wind fields its **coupled dynamics**, and the mechanism for
exchanging fields between it and Cleo's SDM its **coupling comms**.

Physically, this matters because a super-droplet's growth by condensation depends on the local
supersaturation, which depends on temperature and vapour content, which in a real cloud evolve
because of advection, radiation, turbulence, and the feedback of condensation itself onto the
thermodynamic state. Computationally, it matters because it is what decides whether Cleo is
running "one-way coupled" (dynamics affects the super-droplets, but the super-droplets' effect on
the dynamics is discarded — sometimes called *piggybacking*), "two-way coupled" (the
super-droplets' effect feeds back into the dynamics solver), or fully standalone (no dynamics at
all, useful for box-model microphysics tests). See the design paper, Sect. 2.4, and the numerics
paper, Sect. 2, for the physical framing of these modes.

The Two Concepts Involved
----------------------------

Two C++ concepts, both in ``libs/runcleo/``, define what's required:

``CoupledDynamics`` (``coupleddynamics.hpp``) — the dynamics solver itself:

.. code-block:: cpp

    template <typename CD>
    concept CoupledDynamics = requires(CD cd, unsigned int t) {
      { cd.prepare_to_timestep() } -> std::same_as<void>;
      { cd.get_couplstep() } -> std::convertible_to<unsigned int>;
      { cd.on_step(t) } -> std::same_as<bool>;
      { cd.run_step(t, t) } -> std::same_as<void>;
    };

``CouplingComms`` (``couplingcomms.hpp``) — how gridbox state and dynamics-solver state are
exchanged, templated on *both* the ``GridboxMaps`` type and the ``CoupledDynamics`` type it
connects (because, in general, extracting a field from a specific dynamics solver needs to know
both the solver's own representation and the domain's geometry):

.. code-block:: cpp

    template <typename Comms, typename GbxMaps, typename CD>
    concept CouplingComms = requires(Comms s, GbxMaps& gbxmaps, CD& coupldyn, viewh_gbx h_gbxs) {
      { s.template send_dynamics<GbxMaps, CD>(gbxmaps, h_gbxs, coupldyn) } -> std::same_as<void>;
      { s.template receive_dynamics<GbxMaps, CD>(gbxmaps, coupldyn, h_gbxs) } -> std::same_as<void>;
    };

``receive_dynamics`` pulls the dynamics solver's current fields into Cleo's gridboxes;
``send_dynamics`` pushes Cleo's gridbox state (after microphysics has acted on it) back out to
the dynamics solver. One-way coupling is simply the case where ``send_dynamics`` (and/or the
dynamics solver's own response to it) does nothing.

The Available Options
------------------------

Each option lives in its own ``libs/coupldyn_*`` directory (see :doc:`code_structure`) and is
selected at CMake configuration time via ``-DCLEO_COUPLED_DYNAMICS=<option>``
(``libs/CMakeLists.txt``):

.. list-table::
   :header-rows: 1
   :widths: 16 14 70

   * - Directory
     - CMake value
     - What it does
   * - ``coupldyn_null``
     - ``null``
     - No dynamics and no coupling at all (``NullDynamics``/``NullDynComms``). Gridboxes keep
       whatever thermodynamic state they were initialised with. Useful for box-model microphysics
       tests where advection/dynamics are irrelevant (e.g. the collision-box-model examples).
   * - ``coupldyn_fromfile``
     - ``fromfile``
     - One-way coupling: pressure, temperature, vapour/liquid mixing ratio and wind are read from
       pre-computed binary files (``FromFileDynamics``/``FromFileComms``), one value per gridbox
       per coupling step. Cleo's own effect on these fields is never sent back — this is
       "piggybacking" a pre-existing dynamical solution. Used by e.g. the ``rainshaft1d`` and
       ``constthermo2d`` examples.
   * - ``coupldyn_cvode``
     - ``cvode``
     - Two-way coupling to a 0-D adiabatic-parcel ODE model integrated with SUNDIALS' CVODE
       solver (see :doc:`../usage/installation/extern`). Used by the ``adiabaticparcel``
       examples, where condensation's latent heating/moistening genuinely feeds back into the
       parcel's temperature and vapour content.
   * - ``coupldyn_yac``
     - ``yac``
     - Coupling via `YAC <https://dkrz-sw.gitlab-pages.dkrz.de/yac>`_ to an external dynamical
       core (e.g. ICON), for large-domain, large-eddy-simulation-scale two-way coupled runs — the
       regime the design paper's Sect. 1 motivates Cleo for. Requires YAC and YAXT to be
       installed (:doc:`../usage/installation/extern`).
   * - (``cleo_python_bindings``)
     - ``numpy``
     - Not a separate ``coupldyn_*`` directory, but a related option: exposes Cleo via pybind11
       so a numpy/Python-based coupled model can drive it directly, receiving/sending fields as
       numpy arrays rather than through a C++ ``CoupledDynamics`` type. See
       ``examples/python_bindings``.

Choosing between them is a trade-off between physical realism and cost/complexity: ``null`` and
``fromfile`` are cheap and simple but cannot represent feedback from microphysics onto dynamics;
``cvode`` gives full two-way feedback but only for a single 0-D column; ``yac`` gives full
two-way, spatially-resolved feedback but requires an external dynamical core and its installation
overhead. Note also that multiple options can be built into the same binary
(``-DCLEO_COUPLED_DYNAMICS=all``) — the choice of *which one* is instantiated is still made in
your ``main.cpp`` at compile time (see below), not at runtime.

Worked Example: Switching Coupled Dynamics
----------------------------------------------

Compare the ``NullDynamics`` case with the ``FromFileDynamics`` case used by ``rainshaft1d``.
With ``coupldyn_null``, everything is trivial (``libs/coupldyn_null/nulldynamics.hpp`` and
``nulldyncomms.hpp``):

.. code-block:: cpp

    struct NullDynamics {
     private:
      const unsigned int interval;

     public:
      explicit NullDynamics(const unsigned int couplstep) : interval(couplstep) {}
      void prepare_to_timestep() const {}
      unsigned int get_couplstep() const { return interval; }
      bool on_step(const unsigned int t_mdl) const { return false; }
      void run_step(const unsigned int t_mdl, const unsigned int t_next) const {}
    };

    struct NullDynComms {  // satisfies CouplingComms, does nothing
      template <typename GbxMaps, typename CD = NullDynComms>
      KOKKOS_INLINE_FUNCTION void receive_dynamics(const GbxMaps&, const NullDynamics&,
                                                   const viewh_gbx) const {}
      template <typename GbxMaps, typename CD = NullDynComms>
      KOKKOS_INLINE_FUNCTION void send_dynamics(const GbxMaps&, const viewh_constgbx,
                                                NullDynamics&) const {}
    };

With ``coupldyn_fromfile``, ``FromFileDynamics`` (``libs/coupldyn_fromfile/
fromfile_cartesian_dynamics.hpp``) instead owns a ``CartesianDynamics`` object holding vectors of
pressure/temperature/vapour/liquid/wind read from binary files, and ``run_step`` advances its
read-position by one coupling step each time it's called:

.. code-block:: cpp

    struct FromFileDynamics {
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

And in ``main_rshaft1d.cpp``, this is how it's constructed and selected — this is the *only*
place a switch between coupled-dynamics options needs to happen in your ``main.cpp``:

.. code-block:: cpp

    inline CoupledDynamics auto create_coupldyn(const Config& config, const CartesianMaps& gbxmaps,
                                                const unsigned int couplstep,
                                                const unsigned int t_end) {
      const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
      const std::array<size_t, 3> ndims({h_ndims(0), h_ndims(1), h_ndims(2)});
      const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

      return FromFileDynamics(config.get_fromfiledynamics(), couplstep, ndims, nsteps);
      // to switch to no coupled dynamics instead, this whole function would instead just be:
      //   return NullDynamics(couplstep);
    }

Everywhere else in ``main.cpp`` that uses the resulting ``coupldyn`` (constructing ``RunCLEO``,
see :doc:`main_cpp`) refers to it only through the ``CoupledDynamics`` concept
(``CoupledDynamics auto coupldyn = create_coupldyn(...)``), so it does not need to change when
you swap which concrete type ``create_coupldyn`` returns — except that:

.. admonition:: Common pitfall
   :class: warning

   The ``CouplingComms`` type must match *both* your chosen ``GridboxMaps`` type and your chosen
   ``CoupledDynamics`` type — for example ``FromFileComms`` only satisfies
   ``CouplingComms<CartesianMaps, FromFileDynamics>``, not
   ``CouplingComms<CartesianMaps, NullDynamics>``. When you switch coupled dynamics you must also
   switch the comms type (and, in ``libs/CMakeLists.txt``/your executable's ``CMakeLists.txt``,
   which ``coupldyn_*`` library you link against and which CMake ``CLEO_COUPLED_DYNAMICS`` value
   you configure with) — the compiler will refuse to build if these are inconsistent, which is
   the concept system doing its job, but the resulting error message can be long; look for which
   concept's ``requires`` clause failed.

See :doc:`main_cpp` for where ``create_coupldyn`` and the ``CouplingComms`` instance fit into the
rest of ``main.cpp``, and :doc:`tutorial_new_example` for how the ``FromFileDynamicsParams`` seen
above are populated from your YAML configuration file.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
