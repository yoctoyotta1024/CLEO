"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: thermodyngen.py
Project: thermobinary_src
Created Date: Wednesday 26th March 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Thermodyanmics generator to generate winds, temperature, pressure, qvap and qcond fields
from thermo and winds generators in order to create files which CLEO can read for fromfile dynamics.
"""


class ThermodynamicsGenerator:
    """create thermodynamics from winds generator (for wvel, vvel and uvel) and thermo generator
    (for temperature, pressure, qvap and qcond) flow field"""

    def __init__(
        self,
        thermogen,
        windsgen,
    ):
        self.thermogen = thermogen
        self.windsgen = windsgen

    """
    create thermodynamics from winds generator (for wvel, vvel and uvel) and thermo generator
    (for temperature, pressure, qvap and qcond) flow field

    THERMODATA is dictionary with keys: "temp", "press" "qvap" "qcond"
    WINDSDATA is dictionary with keys: "wwvel", "vvel" and "uvel"
    THERMODYNDATA is combination of two dictionaries. Note if keys are not distinct error thrown,
    or data from key in THERMODATA is overwritten by data from matching key in WINDDATA
    """

    def generate_thermodyn(self, gbxbounds, ndims, ntime):
        THERMODATA = self.thermogen.generate_thermo(gbxbounds, ndims, ntime)  # dict
        WINDDATA = self.windsgen.generate_winds(
            gbxbounds, ndims, ntime, THERMODATA
        )  # dict

        matching_keys = set(THERMODATA.keys()) & set(WINDDATA.keys())
        assert (
            not matching_keys
        ) and "THERMODATA and WINDDATA cannot have matching keys, unsafe overwriting"

        THERMODYNDATA = {**THERMODATA, **WINDDATA}

        return THERMODYNDATA
