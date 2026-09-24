# Cleo's bash scripts

Bash scripts to build, compile, run and plot Cleo's examples on a "vanilla" computer, on DKRZ's
Levante and on JSC's JUPITER. `common/` holds everything that is the same on every computer, and
each computer has its own directory with its job scripts (`cpu.sh`, `gpu.sh`), compiler and Kokkos
flags (`build_flags.sh`), runtime environment (`runtime_settings.sh`) and packages (`helpers/`).

- To run the examples, see the docs page for your computer under *Usage → Examples*
  (`docs/source/usage/examples/`).
- For how the scripts work and where to change things (compiler flags, modules, adding an example
  or a computer), see *The Bash Scripts in Detail* (`docs/source/usage/examples/bashscripts.rst`).
