.. _microphysics_and_observers:

Combining Microphysics and Observers
=======================================

This page assumes you've read the conceptual introduction to monoids in
:doc:`../intro/monoids` and the note on C++ concepts in :doc:`code_structure`. Here we go one
level more concrete: what ``MicrophysicalProcess`` and ``Observer`` actually require, how the
``>>`` operator combines them, and a fully worked, annotated example of each.

Microphysics and Observers as Software Abstractions
-----------------------------------------------------

Physically, a "microphysical process" is something that changes droplets — condensation,
collision-coalescence, breakup. Physically, an "observer" is something that reports on the
state of the simulation — writing super-droplet radii to disk, say. As **software**
abstractions in Cleo, both are C++20 concepts: named sets of required member functions, checked
by the compiler, with no shared base class or inheritance involved.

The ``MicrophysicalProcess`` concept (``libs/superdrops/microphysicalprocess.hpp``):

.. code-block:: cpp

    template <typename P>
    concept MicrophysicalProcess =
        requires(P p, const TeamMember& tm, const unsigned int t, subviewd_supers supers,
                 State& state, const NullSDMMonitor mo) {
          { p.next_step(t) } -> std::convertible_to<unsigned int>;
          { p.on_step(t) } -> std::same_as<bool>;
          { p.run_step(tm, t, supers, state, mo) } -> std::same_as<subviewd_supers>;
        };

Any type with these three functions *is* a microphysical process — regardless of what physics it
implements. ``next_step`` says when it next wants to run; ``on_step`` says whether "now" is one
of those times; ``run_step`` actually mutates the super-droplets in one gridbox's ``State`` and
returns the (possibly resized, e.g. after coalescence removes a droplet) view of super-droplets.
This is the same "action" pattern used for the outer timestepping loop, described in
:doc:`../intro/timestepping`.

The ``Observer`` concept (``libs/observers/observers.hpp``) is structurally similar but larger,
because an observer also needs one-time setup/teardown hooks around the whole run:

.. code-block:: cpp

    template <typename Obs>
    concept Observer = requires(Obs obs, unsigned int t, const viewd_constgbx d_gbxs,
                                const subviewd_constsupers d_supers) {
      { obs.next_obs(t) } -> std::convertible_to<unsigned int>;
      { obs.on_step(t) } -> std::same_as<bool>;
      { obs.before_timestepping(d_gbxs, d_supers) } -> std::same_as<void>;
      { obs.after_timestepping() } -> std::same_as<void>;
      { obs.at_start_step(t, d_gbxs, d_supers) } -> std::same_as<void>;
      { obs.get_sdmmonitor() };
    };

``before_timestepping``/``after_timestepping`` run once, at the very start/end of a run (e.g. to
write fixed metadata). ``at_start_step`` is where an observer actually does its
observing/writing, called once per outer-loop timestep (see :doc:`main_cpp`).
``get_sdmmonitor()`` returns an ``SDMMonitor`` — a separate, smaller concept
(``libs/superdrops/sdmmonitor.hpp``) for observers that need to record something *during*
microphysics itself (e.g. how much mass condensed), rather than only at the start of a step; if
an observer doesn't need this it returns ``NullSDMMonitor{}``.

How Combination Works: the Monoid Machinery
-----------------------------------------------

Both concepts have an associated "combined" struct and an ``operator>>`` overload. For
microphysics (``microphysicalprocess.hpp``):

.. code-block:: cpp

    template <MicrophysicalProcess Microphys1, MicrophysicalProcess Microphys2>
    struct CombinedMicrophysicalProcess {
     private:
      Microphys1 a;
      Microphys2 b;

     public:
      CombinedMicrophysicalProcess(const Microphys1 a, const Microphys2 b) : a(a), b(b) {}

      unsigned int next_step(const unsigned int subt) const {
        return Kokkos::min(a.next_step(subt), b.next_step(subt));   // whichever wants to run sooner
      }

      bool on_step(const unsigned int subt) const { return a.on_step(subt) || b.on_step(subt); }

      subviewd_supers run_step(const TeamMember& team_member, const unsigned int subt,
                                subviewd_supers supers, State& state, const SDMMonitor auto mo) const {
        supers = a.run_step(team_member, subt, supers, state, mo);   // a acts first, on a's output
        supers = b.run_step(team_member, subt, supers, state, mo);   // then b acts on a's result
        return supers;
      }
    };

    auto operator>>(const MicrophysicalProcess auto a, const MicrophysicalProcess auto b) {
      return CombinedMicrophysicalProcess(a, b);
    }

``CombinedMicrophysicalProcess`` itself satisfies the ``MicrophysicalProcess`` concept (it has
``next_step``/``on_step``/``run_step`` with the right signatures) — which is precisely what makes
it a **monoid**: combining two members of the set produces another member of the same set, so you
can chain arbitrarily many together (``a >> b >> c`` combines left-to-right). Note the order
matters for ``run_step``: ``a >> b`` means *a runs, then b runs on a's result* — not that they
run simultaneously or in some unspecified order.

``NullMicrophysicalProcess`` is the monoid's null element — it satisfies the concept but its
``run_step`` is a no-op — completing the three monoid properties (create, combine, null) from
:doc:`../intro/monoids`.

``Observer`` follows exactly the same pattern (``CombinedObserver``/``CombineObservers``/
``operator>>`` in ``observers.hpp``), with one extra piece: combining two observers also combines
their ``SDMMonitor``s via ``CombinedSDMMonitor``, so that ``get_sdmmonitor()`` on the combination
returns a monitor that runs *both* constituent monitors.

Worked Example: Combining Two Microphysical Processes
--------------------------------------------------------

This is (lightly trimmed) real code from ``examples/rainshaft1d/src/main_rshaft1d.cpp``:

.. code-block:: cpp

    inline MicrophysicalProcess auto create_microphysics(const Config& config,
                                                         const Timesteps& tsteps) {
      // (1) condensation/evaporation, with its own constant sub-timestep
      const auto c = config.get_condensation();
      const MicrophysicalProcess auto cond =
          Condensation(tsteps.get_condstep(), &step2dimlesstime, c.do_alter_thermo, c.maxniters,
                       c.rtol, c.atol, c.MINSUBTSTEP, &realtime2dimless);

      // (2) collision-coalescence, with Long's hydrodynamic collision efficiency, and its own
      //     (in general different) constant sub-timestep
      const PairProbability auto coalprob =
          LongHydroProb(config.get_coalescence().constcoaleff.coaleff);
      const MicrophysicalProcess auto colls =
          CollCoal(tsteps.get_collstep(), &step2realtime, coalprob);

      // (3) combine: condensation acts on the super-droplets first, then collision-coalescence
      //     acts on condensation's result, each on its own adaptive sub-timestep
      return cond >> colls;
    }

Line by line:

1. ``Condensation`` is constructed with its own sub-timestep (``tsteps.get_condstep()``) and
   solver tolerances read from the YAML config (see :doc:`tutorial_new_example`). It satisfies
   ``MicrophysicalProcess`` by construction (its class implements the three required functions).
2. ``LongHydroProb`` satisfies a different, narrower concept, ``PairProbability`` (the
   probability that a randomly-sampled pair of super-droplets collides), which ``CollCoal`` takes
   as one of its constructor arguments; ``CollCoal`` itself is what satisfies
   ``MicrophysicalProcess``.
3. ``cond >> colls`` produces a single ``CombinedMicrophysicalProcess<Condensation, CollCoal>``
   value. Its exact C++ type is long and you never need to write it out — this is exactly what
   the ``MicrophysicalProcess auto`` return-type syntax buys you (see :doc:`code_structure`'s
   note on concepts). At runtime, on every sub-timestep of the combined process, condensation
   runs on whatever sub-timestep multiple it's due, and collision-coalescence runs on whatever
   (in general different) sub-timestep multiple *it's* due — each governed independently by its
   own ``next_step``/``on_step``, combined adaptively as shown above.

Worked Example: Combining Observers
--------------------------------------

Also from ``main_rshaft1d.cpp``, combining several observers, including a helper function that
itself builds one observer out of several smaller "data-collector" pieces:

.. code-block:: cpp

    template <typename Dataset, typename Store>
    inline Observer auto create_superdrops_observer(const unsigned int interval, Dataset& dataset,
                                                    Store& store, const int maxchunk) {
      // each of these collects one super-droplet attribute into the output dataset
      CollectDataForDataset<Dataset> auto sdid       = CollectSdId(dataset, maxchunk);
      CollectDataForDataset<Dataset> auto sdgbxindex = CollectSdgbxindex(dataset, maxchunk);
      CollectDataForDataset<Dataset> auto xi         = CollectXi(dataset, maxchunk);
      CollectDataForDataset<Dataset> auto radius     = CollectRadius(dataset, maxchunk);
      CollectDataForDataset<Dataset> auto msol       = CollectMsol(dataset, maxchunk);
      CollectDataForDataset<Dataset> auto coord3     = CollectCoord3(dataset, maxchunk);

      // combined into one "collect everything" object, then wrapped in an Observer
      const auto collect_sddata = coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid;
      return SuperdropsObserver(interval, dataset, store, maxchunk, collect_sddata);
    }

    template <typename Dataset, typename Store>
    inline Observer auto create_observer(const Config& config, const Timesteps& tsteps,
                                         Dataset& dataset, Store& store) {
      const auto obsstep  = tsteps.get_obsstep();
      const auto maxchunk = config.get_maxchunk();
      const auto ngbxs    = config.get_ngbxs();

      const Observer auto obs0 = StreamOutObserver(obsstep * 10, &step2realtime);       // progress to stdout
      const Observer auto obs1 = TimeObserver(obsstep, dataset, store, maxchunk, &step2dimlesstime);
      const Observer auto obs2 = GbxindexObserver(dataset, store, maxchunk, ngbxs);
      const Observer auto obs3 = NsupersObserver(obsstep, dataset, maxchunk, ngbxs);
      const Observer auto obs4 = MassMomentsObserver(obsstep, dataset, store, maxchunk, ngbxs);
      const Observer auto obssd = create_superdrops_observer(obsstep, dataset, store, maxchunk);

      return obssd >> obs4 >> obs3 >> obs2 >> obs1 >> obs0;
    }

There are two levels of combination happening here: ``CollectDataForDataset`` is a *third*, more
specialised monoid (satisfying a concept of the same name, defined in
``observers/collect_data_for_dataset.hpp``) used specifically to build up which columns get
written into one Zarr dataset. ``coord3 >> msol >> radius >> xi >> sdgbxindex >> sdid`` combines
six single-attribute collectors into one, which ``SuperdropsObserver`` then wraps to actually
satisfy the full ``Observer`` concept. That resulting ``Observer`` is then itself combined with
five other independent observers (progress printing, time, gridbox indices, super-droplet counts,
mass moments) via the *same* ``>>`` you saw for microphysics — because ``Observer`` is a monoid
in exactly the same structural sense.

.. tip::
   Order matters for readability even when it doesn't change results: it's common in Cleo's
   examples to build the "innermost"/most specific observer or process first and combine
   outward, as shown above, so the final ``return`` statement reads as a short summary of
   everything being combined.

Compatibility Rules and Constraints
--------------------------------------

* You can only combine two values that satisfy the *same* concept — you cannot ``>>`` a
  ``MicrophysicalProcess`` with an ``Observer``; there is a separate ``operator>>`` overload
  (and a separate ``Combined...`` struct) defined for each monoid, and the C++ compiler will
  reject a call whose argument types don't match any overload's constraints.
* Order is significant for the *effect* of combination (as noted above for microphysics), even
  though it doesn't matter for which sub-timestep the combination runs on (``next_step`` and
  ``on_step`` use ``min``/logical-``or``, which are commutative).
* An empty/no-op branch of a combination should use the monoid's null member
  (``NullMicrophysicalProcess{}``, ``NullObserver{}``) rather than, say, an ``if`` statement
  guarding whether to call ``create_microphysics`` at all — this keeps the *type* of your SDM
  fixed at compile time regardless of runtime configuration choices, which is the whole point of
  the monoid pattern (design paper, Sect. 2.3).
* Because combination happens via templates, the *type* of your combined microphysics/observer
  is fixed once you write your ``main.cpp`` and compile — you cannot pick a different combination
  at runtime from the YAML config file alone (though you can, and generally should, use the YAML
  config to set the *parameters* of each already-selected process — see
  :doc:`tutorial_new_example`).

See also :doc:`code_structure` for what a C++ concept/template is if any of this is unfamiliar,
and :doc:`main_cpp` for how ``create_microphysics``/``create_observer`` fit into the rest of
``main.cpp``.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
