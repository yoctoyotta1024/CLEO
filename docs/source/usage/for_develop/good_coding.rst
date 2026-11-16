Good Coding Practices
=====================

.. _environment:

Packages and Environment
------------------------
We provide a ``pyproject.toml`` file for version control over the python and pip-installable
packages you need to use Cleo (e.g to use our ``cleopy`` python package, to run the examples,
to build the documentation, or to run pre-commit etc.).

The recommended way to setup your python environment with all the correct packages is with
`uv <https://docs.astral.sh/uv/>`_. To install ``uv``, if you want, you can use our
Conda/Micromamba environment e.g. via

.. code-block:: console

  $ mamba create --file=environment.yml

which addionally installs Cleo's non-python dependencies available on conda-forge, e.g.
doxygen (which is a non-python requirement for locally building Cleo's documentation).

Once you have installed ``uv`` you can install Cleo's python version and required dependencies.
For only the dependencies required by ``cleopy`` you can do

.. code-block:: console

  $ uv sync --no-dev

Or, if you would like to have the dependencies for Cleo's examples and for Cleo's development
(e.g. for ``pre-commit`` and documentation), you can do

.. code-block:: console

  $ uv sync --extra examples --extra yac

Sometimes it is useful to know where uv has installed your python
(e.g. for :ref:`Cleo's examples <examples>`), in which case you can run

.. code-block:: console

  $ uv python find

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

  $ uv run pre-commit install

Pre-commit should then run automatically upon ``git commit``.


Commit Specification
--------------------
We use `conventional commits <https://www.conventionalcommits.org/>`_ as a a lightweight convention
on top of commit messages to ensure they are meaningful, provide a useful history and can be used
by automated tools such as `Cocogitto <https://github.com/cocogitto/cocogitto>`_.


Code Style / Formatting
-----------------------
For Python, we use `ruff <https://docs.astral.sh/ruff/>`_ for formatting and linting. Ruff checks
are something like the combination of several Python linters (Flake8, isort, pydocstyle etc.) and
the black formatter. We obey the default settings of ruff.

For C++ we obey the Google C++ Style Guide with:
  | IndentWidth: 2
  | TabWidth: 2
  | UseTab: Never
  | ColumnLimit: 100


Auto-Updating Packages
-----------------------

We use dependabot to try to keep our our CI, python (through ``uv``) and pre-commit
dependencies up-to-date. Nevertheless, things can become out of date.

To automatically update some of Cleo's dependencies to their latest versions, run the
following commands:

* for python package dependencies (controlled by Cleo's ``uv.lock``), run ``uv lock --upgrade``

* for pre-commit hook repositories (controlled by ``.pre-commit-config.yaml``),
  run ``pre-commit autoupdate``

*Note:* these commands do not update any manually-installed dependencies, e.g. Kokkos or YAC, nor
the packages installed in a mamba environment (controlled by ``environment.yml``). These commands
also do not update any of the CI ``uses`` dependencies, e.g. for GitHub ``actions`` or ``astral-sh``.
