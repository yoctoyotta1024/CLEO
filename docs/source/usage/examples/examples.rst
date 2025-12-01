.. _examples:

Examples
========

There are various examples of CLEO, with different build configurations, domains, microphysics,
coupling, and super-droplet motion etc. They can be found in the ``CLEO/examples`` directory. If you
would like to a copy of the reference solutions please :ref:`contact us <contact>`.

Before being able to run the examples you will need to locally install the ``plotcleo`` python
package from the ``examples/exampleplotting/`` directory. E.g.

.. code-block:: console

  $ uv build examples/exampleplotting/plotcleo
  $ uv pip install examples/exampleplotting/plotcleo/dist/plotcleo-[version].tar.gz
  $ uv run python -c "import plotcleo"


Each example can be run by building CLEO, compiling the relevant executable, and then running the
example's Python script. There are bash scripts to help you to do all this relatively smoothly on
DKRZ's Levante HPC or on a generic/arbitrary, so-called "vanilla", computer:

.. toctree::
   :maxdepth: 1
   :caption: Running Examples on Different Computers:

   examples_vanilla
   examples_levante
