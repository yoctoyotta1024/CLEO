.. _external_project:

Using Cleo as a Dependency in an External Project
=====================================================

.. admonition:: Ground-truth gap
   :class: warning

   As of the code inspected for this documentation, Cleo's top-level ``CMakeLists.txt`` and
   ``libs/CMakeLists.txt`` do not define any ``install()`` rules or export a
   ``CLEOConfig.cmake``/``find_package(CLEO)`` package. This means Cleo cannot currently be
   consumed the way you'd consume, say, a system-installed library with
   ``find_package(SomeLib REQUIRED)``. What follows describes the integration route that *is*
   supported today — building Cleo's source tree as a CMake subdirectory of your own project,
   the same way Cleo's own top-level ``CMakeLists.txt`` builds ``examples/`` and ``roughpaper/``
   as subdirectories of itself — and flags where you would need to extend Cleo's build system
   (adding ``install()``/export rules) if you want true "linked external library" style
   packaging instead. If you need that and it doesn't exist yet, please
   :ref:`contact us <contact>` or open an issue.

Fetching Cleo
----------------

Add Cleo as a Git submodule of your project (this also pins a specific version/commit, which is
recommended given Cleo is under active development):

.. code-block:: console

  $ git submodule add https://github.com/yoctoyotta1024/CLEO.git extern/CLEO
  $ git -C extern/CLEO checkout <tag-or-commit>
  $ git submodule update --init --recursive

Alternatively, use CMake's ``FetchContent`` to pull Cleo directly in your ``CMakeLists.txt``
(shown in the next section) instead of using a submodule — either works, since both end up
giving CMake a path to Cleo's source tree.

Building Against Cleo from an External CMake Project
---------------------------------------------------------

Because Cleo has no installable package, the working pattern is: bring Cleo's source tree into
your build with ``add_subdirectory`` (or ``FetchContent_MakeAvailable``, which does the same
thing after downloading), exactly as Cleo's own top-level ``CMakeLists.txt`` does for its
``libs/`` directory. A minimal external project's ``CMakeLists.txt``:

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.18.0)
    project(myclimateapp LANGUAGES CXX C)

    find_package(MPI REQUIRED COMPONENTS C)   # Cleo's own CMakeLists.txt also requires this

    # Option A: submodule already present in extern/CLEO
    set(CLEO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/CLEO)
    set(CLEO_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/CLEO)
    set(CLEO_NO_EXAMPLES true)      # don't build Cleo's own examples as part of your build
    set(CLEO_NO_ROUGHPAPER true)    # don't build Cleo's roughpaper scratch executable either
    add_subdirectory(${CLEO_SOURCE_DIR} ${CLEO_BINARY_DIR})

    # Option B (alternative to A): fetch Cleo instead of vendoring it as a submodule
    # include(FetchContent)
    # FetchContent_Declare(cleo
    #   GIT_REPOSITORY https://github.com/yoctoyotta1024/CLEO.git
    #   GIT_TAG <tag-or-commit>)
    # set(CLEO_NO_EXAMPLES true)
    # set(CLEO_NO_ROUGHPAPER true)
    # FetchContent_MakeAvailable(cleo)
    # # FetchContent_MakeAvailable sets cleo_SOURCE_DIR / cleo_BINARY_DIR for you;
    # # Cleo's own CMakeLists.txt still expects CLEO_SOURCE_DIR / CLEO_BINARY_DIR to be set:
    # set(CLEO_SOURCE_DIR ${cleo_SOURCE_DIR})
    # set(CLEO_BINARY_DIR ${cleo_binary_dir})

    add_executable(myapp src/main.cpp)

    # link exactly as Cleo's own examples do (see examples/rainshaft1d/src/CMakeLists.txt) -
    # pick the CLEOLIBS + coupled-dynamics/domain libraries your code actually uses
    set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)
    target_link_libraries(myapp PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
    target_link_libraries(myapp PUBLIC Kokkos::kokkos)
    target_include_directories(myapp PRIVATE "${CLEO_SOURCE_DIR}/libs")

    set_target_properties(myapp PROPERTIES
      CMAKE_CXX_STANDARD_REQUIRED ON CMAKE_CXX_EXTENSIONS ON CXX_STANDARD 20)

Then configure and build your project the same way you would configure Cleo itself (see
:doc:`../usage/installation/installation` and :doc:`../usage/installation/extern` for compiler,
MPI, and optional-dependency requirements — YAC/YAXT, CVODE, or pybind11, only if your chosen
Cleo libraries need them):

.. code-block:: console

  $ cmake -S . -B build -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc \
      -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DCLEO_COUPLED_DYNAMICS=fromfile
  $ cmake --build build --target myapp

Toolchain and Dependency Minimums
-------------------------------------

Carried over from :doc:`../usage/installation/requirements`:

* CMake >= 3.18.0
* A C++ compiler with C++20 support, plus MPI (``find_package(MPI REQUIRED COMPONENTS C)``)
* Kokkos and yaml-cpp are built automatically by Cleo's own ``CMakeLists.txt`` from
  ``extern/kokkos``/``extern/yaml-cpp`` when you ``add_subdirectory`` it — you do not need to
  install these yourself.
* YAC/YAXT (only if you link ``coupldyn_yac``), CVODE (only if you link ``coupldyn_cvode``), and
  pybind11 (only if you build ``cleo_python_bindings``) must be installed separately — see
  :doc:`../usage/installation/extern`.

Common Pitfalls
-------------------

.. admonition:: Version mismatches
   :class: warning

   Because there is no installed/versioned package, "which version of Cleo am I building against"
   is entirely determined by which commit your submodule/``FetchContent`` points at. Pin it
   explicitly (a tag or commit hash, not a moving branch) and record it in your own project's
   version control, the same way you would a vendored dependency.

.. admonition:: Missing ``CLEO_SOURCE_DIR``/``CLEO_BINARY_DIR``
   :class: warning

   Cleo's ``CMakeLists.txt`` files reference ``${CLEO_SOURCE_DIR}``/``${CLEO_BINARY_DIR}``
   directly (rather than CMake's automatic ``PROJECT_SOURCE_DIR``), so if you use
   ``FetchContent`` you must set these two variables yourself from the ``..._SOURCE_DIR``/
   ``..._BINARY_DIR`` variables ``FetchContent_MakeAvailable`` defines, exactly as shown as
   Option B above — Cleo's build will otherwise fail with unset-variable or path errors.

.. admonition:: Build-type / Kokkos backend mismatches
   :class: warning

   If your outer project and Cleo's ``add_subdirectory`` end up configured with different Kokkos
   backends (e.g. your project forces ``Kokkos_ENABLE_CUDA=ON`` globally while Cleo's inherited
   cache still defaults to serial, or vice versa) you will typically see link-time errors deep in
   Kokkos rather than an obvious top-level message. Set Kokkos-related ``-D`` flags once, at your
   project's top-level ``cmake`` invocation, so they propagate consistently into Cleo's
   ``add_subdirectory`` call.

.. admonition:: Excluding Cleo's own examples/roughpaper
   :class: warning

   Forgetting ``CLEO_NO_EXAMPLES``/``CLEO_NO_ROUGHPAPER`` will pull Cleo's own example and
   roughpaper targets into your build graph too (though, because Cleo adds them with
   ``EXCLUDE_FROM_ALL``, they won't be built by a plain ``cmake --build build`` unless
   explicitly named) — harmless, but worth knowing about if your build feels slower to configure
   than expected.

----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
