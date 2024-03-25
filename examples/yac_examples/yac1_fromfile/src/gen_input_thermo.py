'''
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: gen_input_thermo.py
Project: src
Created Date: Monday 25th March 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 25th March 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Python functions used by yac1_fromfile example to make thermo and wind fields
for CLEO to run example with 3-D time-varying thermodynamics.
'''


import sys
import numpy as np

sys.path.append("../../../..") # for imports from pySD package
from pySD import cxx2py
from pySD.thermobinary_src.create_thermodynamics import thermoinputsdict
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid

class TimeVarying3DThermo:
  ''' create some sinusoidal thermodynamics that varies in time and is
  hetergenous throughout 3D domain '''

  def __init__(self, PRESSz0, TEMPz0, qvapz0, qcondz0, Noscs,
               WMAX, Zlength, Xlength, VMAX, Ylength):

    ### parameters of profile ###
    self.PRESSz0 = PRESSz0 # pressure at z=0m [Pa]
    self.TEMPz0 = TEMPz0 # temperature at z=0m [K]
    self.qvapz0 = qvapz0 # vapour mass mixing ratio at z=0m [Kg/Kg]
    self.qcondz0 = qcondz0 # liquid mass mixing ratio at z=0m [Kg/Kg]
    self.Noscs = Noscs # number of oscilations to have in simulations (= T_END / oscilation_period)

    self.WMAX = WMAX  # max velocities constant [m/s]
    self.Zlength = Zlength # wavelength of velocity modulation in z direction [m]
    self.Xlength = Xlength # wavelength of velocity modulation in x direction [m]
    self.VMAX = VMAX  # max horizontal (y) velocity
    self.Ylength = Ylength # wavelength of velocity modulation in y direction [m]

  def idealised_flowfield2D(self, gbxbounds, ndims):

    zfaces, xcens_z = rgrid.coords_forgridboxfaces(gbxbounds, ndims, 'z')[0:2]
    zcens_x, xfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, 'x')[0:2]

    ztilda = self.Zlength / np.pi
    xtilda = self.Xlength / (2*np.pi)
    wamp = 2 * self.WMAX

    WVEL = wamp * np.sin(zfaces / ztilda) * np.sin(xcens_z / xtilda)
    UVEL = wamp * xtilda / ztilda * np.cos(zcens_x/ztilda) * np.cos(xfaces/xtilda)

    return WVEL, UVEL

  def generate_timevarying_3dwinds(self, gbxbounds, ndims, ntime, THERMODATA):

    # time modulation factor for variables at each timestep
    tmod = np.full(ntime, -0.5)
    tmod = np.power(tmod, np.array(range(0, ntime, 1)))

    WVEL, UVEL = self.idealised_flowfield2D(gbxbounds, ndims)

    shape_yface = int((ndims[2]+1)*ndims[0]*ndims[1])
    VVEL = np.full(shape_yface, self.VMAX)

    THERMODATA["WVEL"] = np.outer(tmod, WVEL).flatten()
    THERMODATA["UVEL"] = np.outer(tmod, UVEL).flatten()
    THERMODATA["VVEL"] = np.outer(tmod, VVEL).flatten()

    return THERMODATA

  def generate_thermo(self, gbxbounds, ndims, ntime):

    zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

    shape_cen = int(np.prod(ndims))
    PRESS = np.full(shape_cen, self.PRESSz0)
    TEMP = np.full(shape_cen, self.TEMPz0)
    qvap = np.full(shape_cen, self.qvapz0)
    qcond = np.full(shape_cen, self.qcondz0)

    dimless_omega = 2.0 * np.pi / self.Noscs
    tmod = np.cos(dimless_omega * np.arange(0.0, ntime, 1.0))

    THERMODATA = {
      "PRESS": np.outer(tmod, PRESS).flatten(),
      "TEMP": np.tile(TEMP, ntime),
      "qvap": np.tile(qvap, ntime),
      "qcond": np.tile(qcond, ntime),
    }

    THERMODATA = self.generate_timevarying_3dwinds(gbxbounds, ndims, ntime, THERMODATA)

    return THERMODATA
