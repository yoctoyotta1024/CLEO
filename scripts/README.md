# Cleo's bash scripts

Bash scripts to build, compile, run and plot Cleo's examples on a "vanilla" machine, on DKRZ's
Levante and on JSC's JUPITER. `common/` holds everything that is the same on every machine, and
each machine has its own directory with its job scripts (`cpu.sh`, `gpu.sh`), compiler and Kokkos
flags (`build_flags.sh`), runtime environment (`runtime_settings.sh`) and packages (`helpers/`).

- To run the examples, see the docs page for your machine under *Usage → Examples*
  (`docs/source/usage/examples/`).
- For how the scripts work and where to change things (compiler flags, modules, adding an example
  or a machine), see *The Bash Scripts in Detail* (`docs/source/usage/examples/bashscripts.rst`).
