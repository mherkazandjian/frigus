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
from numpy import log10

nc = 2.e-6

# read the energy levels of H2
en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")

# read the collisional rates for H2 with H
cr, T, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")

# read the einstein coefficients for the H2 transitions
A = read_ei.read_einstein()

# reduce the level representation from 2D indexing to 1D indexing
lin_data = reduce_vj_repr(en_H2, A, cr, T, ini, fin, vj_unique)
en_l, a_eins_l, cr_l, ini_l, fin_l, vj_unique_l, g = lin_data

matrix = population.computeRateMatrix(en_l,
                                      a_eins_l,
                                      cr_l,
                                      ini_l,
                                      fin_l,
                                      vj_unique_l,
                                      g,
                                      T,
                                      nc)

nvj = population.solveEquilibrium(matrix)

#print A,c,h,kB,en_H2

# plotting the population densities

def plot_vj_populations(v, j, nvj):
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(v, j, nvj, 'o')
    ax.set_xlabel('v')
    ax.set_ylabel('j')
    ax.set_zlabel('n(v,j)')
    pylab.show()

nvj[0] = 0
plot_vj_populations(vj_unique[0], vj_unique[1], log10(nvj.flatten()))

cf = population.coolingFunction(en_l,a_eins_l,nvj)
print 'done reading!'

#plt.yscale('log')
#plt.xscale('log')
#plt.plot(T, cr[ini[:,1][0],ini[:,1][1],fin[:,1][0],fin[:,1][1]])
#plt.show()
