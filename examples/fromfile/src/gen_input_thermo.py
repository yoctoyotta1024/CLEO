"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: gen_input_thermo.py
Project: src
Created Date: Monday 25th March 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Python functions used by fromfile example to make thermo and wind fields
for CLEO to run example with 3-D time-varying thermodynamics.
"""

import numpy as np

from cleopy.gbxboundariesbinary_src import read_gbxboundaries as rgrid


class TimeVarying3DThermodyn:
    """create some sinusoidal thermodynamics that varies in time and is
    hetergenous throughout 3D domain"""

    def __init__(
        self, PRESSz0, TEMPz0, qvapz0, qcondz0, WMAX, Zlength, Xlength, VMAX, Ylength
    ):
        ### parameters of profile ###
        self.PRESSz0 = PRESSz0  # pressure at z=0m [Pa]
        self.TEMPz0 = TEMPz0  # temperature at z=0m [K]
        self.qvapz0 = qvapz0  # vapour mass mixing ratio at z=0m [Kg/Kg]
        self.qcondz0 = qcondz0  # liquid mass mixing ratio at z=0m [Kg/Kg]
        self.dimless_omega = np.pi / 4.0  # ~ frequency of time modulation []

        self.WMAX = WMAX  # max velocities constant [m/s]
        self.Zlength = Zlength  # wavelength of velocity modulation in z direction [m]
        self.Xlength = Xlength  # wavelength of velocity modulation in x direction [m]
        self.VMAX = VMAX  # max horizontal (y) velocity
        self.Ylength = Ylength  # wavelength of velocity modulation in y direction [m]

    def idealised_flowfield2D(self, gbxbounds, ndims):
        zfaces, xcens_z, ycens_z = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")
        zcens_x, xfaces, ycens_x = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "x")

        ztilda = self.Zlength / np.pi
        xtilda = self.Xlength / (2 * np.pi)
        wamp = 2 * self.WMAX

        WVEL = wamp * np.sin(zfaces / ztilda) * np.sin(xcens_z / xtilda)
        UVEL = (
            wamp * xtilda / ztilda * np.cos(zcens_x / ztilda) * np.cos(xfaces / xtilda)
        )

        # modulation in y direction
        WVEL *= 1.0 + 0.5 * np.cos(self.Ylength / np.pi * ycens_z)
        UVEL *= 1.0 + 0.5 * np.cos(self.Ylength / np.pi * ycens_x)

        return WVEL, UVEL

    def gen_3dvvelocity(self, gbxbounds, ndims):
        zcens_y, xcens_y, yfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "y")
        zxmod = zcens_y / self.Zlength + xcens_y / self.Xlength
        VVEL = self.VMAX * (zxmod + np.cos(self.Ylength / np.pi * yfaces))

        return VVEL

    def generate_timevarying_3dwinds(self, gbxbounds, ndims, ntime, THERMODATA):
        # time modulation factor for variables at each timestep
        tmod = np.full(ntime, -0.5)
        tmod = np.power(tmod, np.array(range(0, ntime, 1)))

        WVEL, UVEL = self.idealised_flowfield2D(gbxbounds, ndims)
        VVEL = self.gen_3dvvelocity(gbxbounds, ndims)

        THERMODATA["WVEL"] = np.outer(tmod, WVEL).flatten()
        THERMODATA["UVEL"] = np.outer(tmod, UVEL).flatten()
        THERMODATA["VVEL"] = np.outer(tmod, VVEL).flatten()

        return THERMODATA

    def generate_3dsinusoidal_variable(self, gbxbounds, ndims, amp):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

        ztilda = self.Zlength / np.pi / 3.0
        xtilda = self.Xlength / np.pi / 4.0
        ytilda = self.Ylength / np.pi / 1.33

        return amp + 0.25 * amp * (
            np.sin(zfulls / ztilda) * np.sin(xfulls / xtilda) + np.sin(yfulls / ytilda)
        )

    def generate_thermodyn(self, gbxbounds, ndims, ntime):
        PRESS = self.generate_3dsinusoidal_variable(gbxbounds, ndims, self.PRESSz0)
        TEMP = self.generate_3dsinusoidal_variable(gbxbounds, ndims, self.TEMPz0)
        qvap = self.generate_3dsinusoidal_variable(gbxbounds, ndims, self.qvapz0)
        qcond = self.generate_3dsinusoidal_variable(gbxbounds, ndims, self.qcondz0)

        tmod = np.cos(self.dimless_omega * np.arange(0.0, ntime, 1.0))
        tmod = 1 + 0.5 * tmod

        THERMODATA = {
            "PRESS": np.outer(tmod, PRESS).flatten(),
            "TEMP": np.outer(tmod, TEMP).flatten(),
            "qvap": np.outer(tmod, qvap).flatten(),
            "qcond": np.outer(tmod, qcond).flatten(),
        }

        THERMODATA = self.generate_timevarying_3dwinds(
            gbxbounds, ndims, ntime, THERMODATA
        )

        return THERMODATA
