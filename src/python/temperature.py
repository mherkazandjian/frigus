# -*- coding: utf-8 -*-
"""
radiation temperature and gas temperature as a function of redshift/cosmic time.

The data of the gas temperature is read from a pre-computed file
 (Read/temp_gas.txt) computed by C. Coppola.
"""

import pylab
import numpy
from numpy import *
from scipy.constants import *

from scipy import interpolate


def gen_T_g_interpolation_function():
    """contrunct the interpolation function for the gas temperature"""
    temperature = loadtxt("Read/temp_gas.txt")
    z_read = temperature[:, 0]
    tg_read = temperature[:, 1]
    
    # itg = interpolate.splrep(z_read, tg_read, s=0)
    
    log10_T_g = interpolate.interp1d(z_read, numpy.log10(tg_read))

    return log10_T_g

log10_T_g = gen_T_g_interpolation_function()


def T_r(z):
    """ given the redshift, it returns the radiation temperature"""
    tr=2.726*(1.0+z)
    return tr


def T_g(z):
    """ given the redshift, it returns the gas temperature"""
    # tg  = interpolate.splev(z,itg,der=0)

    # .. todo:: 
    # there shuould be a check that the interpolated temperature is not negative

    return 10.0**log10_T_g(z)


if __name__ == "__main__":

    print 'testing the temperature interpolation'

    temperature = loadtxt("Read/temp_gas.txt")
    z_read = temperature[:,0]
    tg_read = temperature[:,1]

    pylab.xlim([-1, 4000])
    #pylab.xscale('log')
    pylab.yscale('log')
    pylab.figure(1)
    pylab.plot(z_read, tg_read,'r+')
    zgrid = numpy.linspace(z_read.min(), z_read.max(), 1000.0)
    pylab.plot(zgrid, T_g(zgrid), 'b--')

    #tg  = interpolate.splev(z,itg,der=0)

    pylab.show()
