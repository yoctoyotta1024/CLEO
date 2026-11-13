.. _external_project:

Using Cleo in an External Project
===================================

This page covers two ways to use Cleo as a dependency in an external project.

.. note::
   This page assumes you've already read the installation and dependencies documentation for
   building Cleo itself (compiler, CMake, Kokkos, and MPI requirements) — the same requirements
   apply when building Cleo as a dependency of your own project, since you are compiling exactly
   the same Cleo source code either way.

Currently, Cleo's top-level ``CMakeLists.txt`` does not define any ``install()`` rules or export a
``CLEOConfig.cmake`` / ``find_package(CLEO)`` package. This means Cleo cannot currently be used in
the way you might recognise, e.g. for a system-installed library, with ``find_package(SomeLib REQUIRED)``.
Instead, you build Cleo's libraries *as part of* your own project — the same way Cleo does for
``Kokkos`` and ``yaml-cpp``, by pointing ``add_subdirectory()`` to a fetched, or cloned, checkout
of their source code (see the ``CLEO/extern`` directory).

There are two equivalent ways to get a copy of Cleo's source into your project so you can
``add_subdirectory()`` it: a local clone/submodule, or CMake's ``FetchContent`` module to clone it
automatically. Both end up doing the same thing — making Cleo's library targets (defined in
``libs/CMakeLists.txt``) available for your executable to link against.

Option A: Local Clone or Git Submodule
---------------------------------------

Clone Cleo (or add it as a ``git submodule``) somewhere inside your own project's repository, then
point ``add_subdirectory()`` at it from your project's top-level ``CMakeLists.txt``, exactly as
Cleo's own outermost ``CMakeLists.txt`` does for ``extern/kokkos`` and ``extern/yaml-cpp``:

.. code-block:: console

   ### from the root of your own project's repository
   git submodule add https://github.com/yoctoyotta1024/CLEO.git extern/CLEO
   git submodule update --init --recursive


.. code-block:: cmake

   # in your top-level CMakeLists.txt
   add_subdirectory(extern/CLEO)

.. tip::
   Pin the version by checking out a specific tag or commit inside the submodule (``cd extern/CLEO &&
   git checkout <tag_or_commit>``) and committing that submodule reference, so anyone who clones your
   project gets the same Cleo version you built and tested against.

Option B: FetchContent
------------------------

Alternatively, do exactly what Cleo does to install Kokkos and yaml-cpp. Let CMake fetch Cleo for
you at configure time, with no manual clone or submodule bookkeeping, using the ``FetchContent``
module. Inside your ``extern/CLEO`` directory, rather than keeping a full-clone of CLEO, have a
single ``CMakeLists.txt`` with the following contents:

.. code-block:: cmake

   # set cmake version
   if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
      cmake_minimum_required(VERSION 3.18.0)
   endif()

   include(FetchContent)
   FetchContent_Declare(
     CLEO
     DOWNLOAD_EXTRACT_TIMESTAMP TRUE
     GIT_REPOSITORY https://github.com/yoctoyotta1024/CLEO.git
     GIT_TAG        <tag_or_commit>   # pin a specific released version or commit
   )
   FetchContent_MakeAvailable(CLEO)

   message(STATUS "CLEO installation in: ${CMAKE_CURRENT_BINARY_DIR}")

.. tip::
   ``GIT_TAG`` lets you pin the Cleo version you want to use — set it to a specific tag or
   commit, never to a moving branch name like ``main``, so your build is reproducible.


What's Made Available
----------------------

Both of the above options gives you the same thing: during your project's build configuration,
Cleo's outermost ``CMakeLists.txt`` will be run. In turn this will run the  ``CMakeLists.txt``
files in ``libs``, ``extern``, ``examples``, etc. to define all of Cleo's library and
executable targets for you to link against.

If you want to exclude some of Cleo's targets from your configuration (e.g. those in Cleo's examples,
roughpaper, or certain libraries) you can pass Cleo's CMake flags just like you would when you
build Cleo standalone. Most likely, you will want to pass ``-DCLEO_NO_ROUGHPAPER=true``
and ``-DCLEO_NO_EXAMPLES=true``, since you probably don't need Cleo's own executables, only its
libraries. See :ref:`cleo_cmake_flags` for more information.


Linking Cleo's Libraries For your Executable
---------------------------------------------

The following library targets are always built for you to link against:
- ``configuration``
- ``initialise``
- ``superdrops``
- ``gridboxes``
- ``runcleo``
- ``observers``
- ``zarr``

Depending on the CMake flags you pass for ``CLEO_DOMAIN``, ``CLEO_COUPLED_DYNAMICS`` and
``CLEO_NO_PYBINDINGS`` (see :ref:`cleo_cmake_flags`), the following targets may also be built:
- ``cartesiandomain``
- ``coupldyn_cvode``
- ``coupldyn_fromfile``
- ``coupldyn_yac``
- ``cleo_python_bindings``

A minimal external executable links against all of Cleo's library targets which it uses in its own
source code. E.g. via `target_link_libraries(your_executable PRIVATE configuration initialise [...])``.
The worked example below shows how this can be done for one domain and one coupled-dynamics choice.

.. caution::
   **Common Pitfall:** forgetting to list one of these targets in your executable's
   ``target_link_libraries()`` produces undefined-reference linker errors for the symbols that
   library provides, not a configure-time error — since ``EXCLUDE_FROM_ALL`` targets are silently
   skipped until something references them. If you get undefined-reference errors at link time,
   check whether the missing symbol belongs to one of the domain/coupled-dynamics libraries above.

Worked Example: an External Project's ``CMakeLists.txt``
-------------------------------------------------------------

A minimal external project with one executable, ``my_cleo_app.cpp``, that includes Cleo headers
and links against Cleo's libraries, using the local-clone/submodule approach from Option A:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.18.0)

   project(your_executable LANGUAGES CXX C)

   set(CMAKE_CXX_STANDARD 20)
   set(CMAKE_CXX_STANDARD_REQUIRED ON)

   # skip Cleo's own examples/roughpaper/py-bindings; we only want its libraries
   set(CLEO_NO_EXAMPLES ON)
   set(CLEO_NO_ROUGHPAPER ON)
   set(CLEO_NO_PYBINDINGS ON)

   # choose which of Cleo's optional libraries to build
   set(CLEO_DOMAIN "cartesian")
   set(CLEO_COUPLED_DYNAMICS "fromfile")

   # pull in Cleo's own CMakeLists.txt (either through clone or FetchContent; Option A/B above)
   add_subdirectory(extern/CLEO)

   add_executable(your_executable your_executable.cpp)
   target_link_libraries(your_executable PRIVATE
     configuration
     initialise
     superdrops
     gridboxes
     runcleo
     zarr
     observers
     cartesiandomain
     coupldyn_fromfile
   )

Annotated, block by block:

1. ``cmake_minimum_required``/``project``: your own project's usual boilerplate; the CMake version
   must still satisfy Cleo's own ``cmake_minimum_required(VERSION 3.18.0)``.
2. ``CMAKE_CXX_STANDARD 20``: set *before* ``add_subdirectory(extern/CLEO)``, so your own sources
   and Cleo's are compiled with the same C++ standard — Cleo sets this itself for its own targets,
   but setting it in your own project too keeps your sources consistent with it.
3. The three ``CLEO_NO_*`` variables: skip configuring parts of Cleo you don't need as a dependency
   (Cleo's own examples, roughpaper, and Python bindings).
4. ``CLEO_DOMAIN``/``CLEO_COUPLED_DYNAMICS``: choose which optional libraries get configured — here,
   the Cartesian domain and the from-file coupled dynamics (see :doc:`code_structure`
   for what these mean and how to pick between them).
5. ``add_subdirectory(extern/CLEO)``: through either option A or option B above, runs Cleo's
   top-level ``CMakeLists.txt`` (and therefore ``libs/CMakeLists.txt``) as part of your own
   project's configure step — this is what actually defines all Cleo's library targets listed below.
6. ``add_executable``/``target_link_libraries``: your own executable, linked against every Cleo
   library target it needs — the seven "always built" targets, plus ``cartesiandomain`` and
   ``coupldyn_fromfile``, which match the options set in step 4.


Common Pitfalls
-----------------

Often people struggle to pass the correct flags (for Kokkos, MPI, YAC etc.) to configure their
project with Cleo. To avoid this, we reccomend following similar scripting to build your project
as we do for Cleo stand-alone (see :ref:`building_cleo`). Some common issues include:

- **Build type.** Set ``CMAKE_BUILD_TYPE`` (e.g. ``Release``) consistently for your whole
  configure — if you don't set it at all, Cleo and your own project both fall back to CMake's
  (unoptimized) default, which is a correctness-neutral but easy-to-miss performance pitfall for
  Kokkos-based code in particular.
- **Kokkos and YAC flags.** Any ``-DKokkos_ARCH_...``/``-DKokkos_ENABLE_...`` or ``YAC`` flags you
  pass to ``cmake`` apply globally to the whole configure, including Cleo's fetched Kokkos — set
  them once, when you configure your own project's build.
- **MPI compiler wrappers.** Since Cleo's top-level ``CMakeLists.txt`` unconditionally requires MPI,
  you'll generally want to configure your whole project with MPI compiler wrappers as you would
  for building Cleo standalone (e.g. ``-DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc``).
- **Missing target causing linker error.** See the caution box above: an omitted
  domain/coupled-dynamics target fails at link time with undefined references, not at
  ``cmake``'s configure time.


Using Cleo through Python
---------------------------

To use Cleo through Python you will no longer need to create an executable or ``CMakeLists.txt``
for your project, simply building and compiling Cleo alone is sufficient. You need to compile
the ``cleo_python_bindings`` target.

Once you've done that, you can import the bindings in Python just like a normal module, e.g.

.. code-block:: python

   import cleo_python_bindings as cleo
   from cleo_python_bindings import coupldyn_numpy

From there you can use the parts of Cleo which the bindings expose to drive Cleo through Python
numpy arrays.

Have a look at Cleo's ``python_bindings`` example to see how this can work in practice:
``CLEO/examples/python_bindings/`` and :doc:`the instructions to run it<../usage/examples/examples>`.
In this example (run twice, once with and once without microphysics), we also created a Python
class, ``CleoSDM``, to control the time-stepping sequence in a similar way to how the class,
``RunCleo``, works in Cleo's C++ code.


----

.. note::
   This page was generated by Claude (Anthropic), at the request of Clara Bayley, based on
   Cleo's design and numerics papers (Bayley et al. 2026, GMD 19, 6099-6120 and 6121-6142) and
   the Cleo source code. Please check it for accuracy and open an issue or get in touch if you
   spot anything wrong or out of date.
