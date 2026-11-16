Contributing Guidelines
==========================

How to Contribute?
-------------------

To contribute to the Cleo repository: fork the repository, create a branch for your
contribution, and align with the following practices:

1. :ref:`Be consistent with our general coding rules<general_code_rules>`

2. :ref:`Format and lint your code correctly<format_lint>`

3. :ref:`Use pre-commit to check your commits<pre_commit>`

4. :ref:`Use conventional commit messages<conventional_commits>`

5. :ref:`Create unit tests for your code<unit_tests>`

6. :ref:`Update Cleo's documentation with your contributions<contrib_our_docs>`

You can then contribute via pull request (PR) back into Cleo's main code repository.

The Rules on Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Your PRs *must* obey the practices listed above** (see more details below).

* **Keep your PRs *up-to-date* by rebasing against main**. Until a PR is merged into Cleo's main
  branch, is the responsibilty of the owner of the PR to keep it up-to-date and resolve any
  code-conflicts.

* ***Note the 3x1 Rule:* One Branch; One Contribution; One PR**
    * One branch should contain one "contribution".
    * One "contribution" may be large, but should be as small as reasonably possible.
    * One PR should contain one "contribution". A PR will not be accepted if it is deemed
      to contain many contributions that could instead logically be split into smaller PRs from
      seperate branches.

* **Contributions from a ``main`` branch will never be accepted.**


.. _general_code_rules:

General Coding Rules
----------------------

* Blend-in with whatever you change or add to.

* A new class/struct deserves a new file.

* As much code should be in ``.cpp`` rather than ``.hpp`` files as is reasonably possible.

* Improve code readability wherever possible, but take caution not to be overly verbose.

* Favour more-involved structral changes which will make the overall code better over more-elongated
  code and exceptional case handling.

* Name variables/functions/classes appropriately. Favour readable and descriptive names over shorter
  code. At the same time, take caution not to be overly verbose.

* All formulas, constants, physical values, numerical methods, etc. should be referenced (to
  papers, people or other code-bases, etc.) in clearly visible doc-strings near the relevant code.

* Correct and improve doc-strings wherever possible. Keep comments and doc-strings short.

* Avoid adding comments about future actions.

* Never add commented-out code to the code base. You need a *very* good reason to do this.


.. _format_lint:

Formatting and Linting your Code
----------------------------------

Python
~~~~~~~

For Python, we use `ruff <https://docs.astral.sh/ruff/>`_ for formatting and linting. Ruff checks
are something like the combination of several Python linters (Flake8, isort, pydocstyle etc.) and
the black formatter. We obey the default settings of ruff and these are automatically enforced for
you if you use our pre-commit hooks (see below).

C++
~~~~

For C++, we obey the Google C++ Style Guide with:
  | IndentWidth: 2
  | TabWidth: 2
  | UseTab: Never
  | ColumnLimit: 100

We recommend to enforce these rules with the `ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`_
tool, which can be integrated in several text/code editors (e.g. Vim, Emacs, VSCode, ...).


.. _pre_commit:

Pre-Commit
----------
We use pre-commit to check our code for simple issues before submission to code review, such as
before pushing to a GitHub repository. Pull requests that fail our pre-commit checks will fail
in our continuous integration, therefore, you must use pre-commit too if you want to contribute
to Cleo's code base.

Pre-Commit has "hooks" which can be used for a vast range of measures to ensure you follow good
coding practices. For example, we use pre-commit to find and correct mistakes such as missing
semicolons, trailing whitespaces, bad indentations etc., and to ensure we conform to
our chosen code style guide.

To use our pre-commit hooks, first load your environment with pre-commit installed
(see e.g. :ref:`environment`), then install the pre-commit hooks, with ``pre-commit install``, e.g.

.. code-block:: console

  $ uv run pre-commit install

Pre-commit should then run automatically upon ``git commit``.

You can learn more by checking out `pre-commit <https://pre-commit.com/>`_'s comprehensive
documentation.


.. _conventional_commits:

Commit Specification
--------------------
We use `conventional commits <https://www.conventionalcommits.org/>`_ as a a lightweight convention
on top of commit messages to ensure they are meaningful, provide a useful history and can be used
by automated versioning tools such as `Cocogitto <https://github.com/cocogitto/cocogitto>`_.

In short, all your commits must be structured as follows:

.. code-block:: console

  <type>[optional scope]: <description>

  [optional body]

  [optional footer(s)]


* The ``<type>`` may be ``fix:``, ``feat:``, ``refactor:``, ``chore:``,  ``ci:``, ``test:``,  ``perf:``, ...
  (see `conventional commits <https://www.conventionalcommits.org/>`_ for more options).

* A breaking change must be marked with a footer ``BREAKING CHANGE: <description>`` and/or a ``!`` after
  the type/scope.

* AI generated commits must be marked with a footer: ``generated with <AI model>``.


.. _unit_tests:

Unit Tests
----------

More on this comming soon! (Last updated: Nov. 2026.)


.. _contrib_our_docs:

Contributing to Cleo's Documentation
--------------------------------------

Your PR *must* also update Cleo's docs if your code contribution requires new documentation or
corrections to the existing documentation. If you do not believe your PR should update Cleo's
docs, you will need to explain why compellingly in your PR description.

You can see :doc:`the documentation on the docs<ourdocs>` for instructions how to locally
build the documentation. This enables you to visualise, as if it were published, the documentation
you edit.
