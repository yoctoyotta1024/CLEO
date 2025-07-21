"""
----- CLEO -----
File: conf.py
Project: source
Created Date: Monday 20th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
configuration for CLEO documentation made using Sphinx
"""

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import pathlib
import sys

# sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "CLEO"
license = "BSD 3-Clause"
copyright = "(2023) MPI-M, Clara Bayley"
author = "Clara Bayley & Other Developers"
release = "0.38.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.duration",
    "sphinx_copybutton",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "breathe",
    "sphinx.ext.viewcode",
]

# configuration of citations using bibtex file(s)
bibtex_bibfiles = ["./references.bib"]
bibtex_reference_style = "label"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Integrate doxygen with sphinx via breathe
breathe_projects = {
    "libs": "../build/doxygen/xml/",
    "runcleo": "../build/doxygen/xml/",
    "superdrops": "../build/doxygen/xml/",
    "zarr": "../build/doxygen/xml/",
    "observers": "../build/doxygen/xml/",
    "initialise": "../build/doxygen/xml/",
}

breathe_default_project = "proj"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# To include the date of the last visit in Sphinx documentation, use the last
# updated feature of Sphinx. This feature automatically adds the last modification
# date of the source file to the rendered HTML output.
html_last_updated_fmt = "%d %B %Y"
