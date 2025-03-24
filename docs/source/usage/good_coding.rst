Good Coding Practices
=====================

.. _environment:

Packages and Environment
------------------------
We provide an environment yaml file for creating an envirnoment with all the packages you need to
use CLEO (to use PySD, to build the documentation, or to run pre-commit etc.).

The easiest way to use it is to create a new Conda or Micromamba environment, e.g.

.. code-block:: console

  $ conda env create --file=environment.yml

or

.. code-block:: console

  $ micromamba env create --file=environment.yml

And then activate this environment everytime you use CLEO, e.g.

.. code-block:: console

  $ micromamba activate cleoenv

We kindly ask that you also :ref:`contact us <contact>` or `open a new
issue <https://github.com/yoctoyotta1024/CLEO/issues/new>`_ on our GitHub repository if you discover
this environment is missing any packages.

Pre-Commit
----------
We use pre-commit to check our code for simple issues before submission to code review, such as
before pushing to a GitHub repository. We reccomend you use it too.

Pre-Commit has "hooks" which can be used for a vast range of measures to ensure you follow good
coding practices. For example, we use pre-commit to find and correct mistakes such as missing
semicolons, trailing whitespaces, bad indentations etc., and to ensure we conform to
our chosen code style guide.

You can learn more by checking out `pre-commit <https://pre-commit.com/>`_'s comprehensive
documentation.

In case you want to use our pre-commit hooks, first load your environment with pre-commit installed,
then install the pre-commit hooks, e.g.

.. code-block:: console

  $ pre-commit install

Pre-commit should then run automatically upon ```git commit``.


Commit Specification
--------------------
We use `conventional commits <https://www.conventionalcommits.org/>`_ as a a lightweight convention
on top of commit messages to ensure they are meaningful, provide a useful history and can be used
by automated tools such as `Cocogitto-bot <https://github.com/apps/cocogitto-bot>`_.

Code Style / Formatting
-----------------------
For Python, we use `ruff <https://docs.astral.sh/ruff/>`_` for formatting and linting. Ruff checks
are something like the combination of several Python linters (Flake8, isort, pydocstyle etc.) and
the black formatter. We obey the default settings of ruff except we ignore E402 errors.

For C++ we obey the Google C++ Style Guide with:
  | IndentWidth: 2
  | TabWidth: 2
  | UseTab: Never
  | ColumnLimit: 100
