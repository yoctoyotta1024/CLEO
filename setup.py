'''
----- CLEO -----
File: setup.py
Project: CLEOfire
Created Date: Thursday 12th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Thursday 23rd November 2023
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
    version='4.0',
    packages=find_packages(),
    install_requires=[
        'pytest',
        'sphinx',
        'matplotlib',
    ],
)
