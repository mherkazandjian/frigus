# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 23:45:39 2015

@author: carla
"""

import read_ei
import read_cr
from read_cr import label
import read_levels
import population
from population import reduce_vj_repr
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer

# read the energy levels of H2
en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")

# read the collisional rates for H2 with H
cr, T, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")

# read the einstein coefficients for the H2 transitions
A = read_ei.read_einstein()

# reduce the level representation from 2D indexing to 1D indexing
data = reduce_vj_repr(en_H2, A, cr, T, ini, fin, vj_unique)

Tracer()()
#
matrix = population.computeRateMatrix()
#nvj = population.solveEquilibrium()

#print A,c,h,kB,en_H2

print 'done reading!'

#plt.yscale('log')
#plt.xscale('log')
#plt.plot(T, cr[ini[:,1][0],ini[:,1][1],fin[:,1][0],fin[:,1][1]])
#plt.show()
