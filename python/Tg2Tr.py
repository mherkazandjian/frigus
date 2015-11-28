# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 18:03:31 2015

@author: carla
"""

import pylab
import numpy
from   numpy import *
from   scipy.constants import *
from Newton import Newton

from scitools.std import *
from scipy import interpolate

def gen_T_g_interpolation_function():
    '''contrunct the interpolation function for the gas temperature'''
    temperature  = loadtxt("Read/temp_gas.txt")  
    z_read       = temperature[:,0]
    tg_read      = temperature[:,1]
    
    #itg   = interpolate.splrep(z_read,tg_read,s=0) 
    
    log10_T_g = interpolate.interp1d(z_read, numpy.log10(tg_read))

    return log10_T_g

log10_T_g = gen_T_g_interpolation_function()

def T_r(z):
    ''' given the redshift, it returns the radiation temperature'''
    tr=2.726*(1.0+z)
    return tr

def T_g(z):
    ''' given the redshift, it returns the gas temperature'''
    #tg  = interpolate.splev(z,itg,der=0)

    # .. todo:: 
    # there shuould be a check that the interpolated temperature is not negative

    return 10.0**log10_T_g(z)

def f(x):
    return x**2 - 1
def F(gamma):
    return f(gamma) - xi
def dFdx(gamma):
    return (F(gamma+h) - F(gamma-h))/(2*h)

h = 1E-6
x = linspace(0.01, 3, 21)
g = zeros(len(x))
for i in range(len(x)):
xi = x[i]
# Compute start value (use last g[i-1] if possible)
if i == 0:
gamma0 = x[0]
else:
gamma0 = g[i-1]
gamma, n, F_value = Newton(F, gamma0, dFdx)
g[i] = gamma
plot(x, f(x), ’r-’, x, g, ’b-’,
title=’f1’, legend=(’original’, ’inverse’))




#if __name__ == "__main__":

#    print 'testing the temperature interpolation'
#
#    temperature  = loadtxt("Read/temp_gas.txt")  
#    z_read       = temperature[:,0]
#    tg_read      = temperature[:,1]

#    pylab.xlim([-1, 4000])
#    #pylab.xscale('log')
#    pylab.yscale('log')
#    pylab.figure(1)
#    pylab.plot(z_read, tg_read,'r+')
#    zgrid = numpy.linspace(z_read.min(), z_read.max(), 1000.0)
#    pylab.plot(zgrid, T_g(zgrid), 'b--')

#    #tg  = interpolate.splev(z,itg,der=0)

#    pylab.show()
