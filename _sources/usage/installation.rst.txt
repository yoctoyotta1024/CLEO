.. _installation:

Installation
============

First clone CLEO's GitHub repository. Everything will be much easier for you if you clone CLEO in
your home directory. It’s not essential, but if you choose to do otherwise you may have to change
some extra paths in the bash and Python scripts.

.. code-block:: console

  $ git clone https://github.com/yoctoyotta1024/CLEO.git


Everything will also be much easier for you if you set up and use an environment with the packages
we suggest in our environment.yml file, e.g. using Conda/Micromamba:

.. code-block:: console

  $ micromamba create --file=environment.yml
  $ micromamba activate cleoenv

*Note*: on Levante you need to first load an openmpi package to install mpi4py, e.g. via
```module load openmpi/4.1.2-gcc-11.2.0```

Finally we suggest you use pre-commit. You can install our hooks via:

.. code-block:: console

  $ pre-commit install


That's it, you're done!

Now maybe you want to run :doc:`some examples <examples>` or try out
the :doc:`quickstart<quickstart>`...
