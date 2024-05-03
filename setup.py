'''
----- CLEO -----
File: setup.py
Project: CLEO
Created Date: Thursday 12th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 1st May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
'''


from setuptools import setup, find_packages

setup(
    name='CLEO',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # 'python>=3.12',
        'pytest',
        'sphinx',
        'furo',
        'sphinx_copybutton',
        'sphinxcontrib-bibtex',
        'breathe',
        'numpy',
        'scipy',
        'matplotlib',
        'netcdf4',
        'xarray',
        'awkward',
        'zarr',
    ],
)
