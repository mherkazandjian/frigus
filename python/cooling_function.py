# -*- coding: utf-8 -*-
"""
Calculates the cooling function of H2 colliding with H using collisional data
by F. Lique.

In calculating the cooling function, the required data are:

  - energy levels of H2 (vibrational and rotational)
  - collisional coefficients of H2 with H (K_ij)
  - radiative coefficients (A_ij, B_ij, B_ji)

Limitations

  - The smallest data set of (A, B, K) determines the number of states to be
    inserted in the model.
"""

import read_ei
import read_cr
import read_levels
import population
from population import reduce_vj_repr, coolingFunction, fit_glover
import matplotlib.pyplot as plt

from numpy import zeros
import pdb
from numpy import log10, unique

# read the energy levels (v, j, energy)
levels = read_levels.read_levels_lique('Read/H2Xvjlevels_francois_mod.cs')
# print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
# for level in levels:
#     print('{:<3}{:<3}{:<10}'.format(level['v'], level['j'], level['E']))

# read the einstein coefficients for the H2 transitions
A = read_ei.read_einstein()
pdb.set_trace()

# read the Einstein coefficients
nc = 1.e9
'''density of the colliding species, in units of 1.e3 cm-3 as in Lipovka'''

# read the energy levels of H2
en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")

# read the collisional rates for H2 with H
cr, T, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")


# creating the array for the radiation temperatures;
# at the moment, same size of tkin but constant value @ 2.73
tcmb = zeros(T.size)
tcmb = 2*T
#tcmb[:] = 2.*T[:]


# reduce the level representation from 2D indexing to 1D indexing
lin_data = reduce_vj_repr(en_H2, A, cr, T, ini, fin, vj_unique)
en_l, a_eins_l, cr_l, ini_l, fin_l, vj_unique_l, g = lin_data

cf = zeros(T.size)
fit = zeros(T.size)
for itemp in range(T.size):
    matrix = population.computeRateMatrix(en_l,
                                      a_eins_l,
                                      cr_l,
                                      ini_l,
                                      fin_l,
                                      vj_unique_l,
                                      g,
                                      T,
                                      tcmb,
                                      nc,
                                      itemp,
                                      itemp)

    nvj = population.solveEquilibrium(matrix)
    cooling_function = population.coolingFunction(nvj,
                                              en_l,
                                              a_eins_l,
                                              ini_l,
                                              fin_l,
                                              vj_unique_l)
    cf[itemp] = cooling_function
    fit[itemp] = fit_glover(T[itemp])

    print cooling_function*1e13 # to have the output in erg cm-3 s-1


plt.plot(T, cf*1e13, 'x', T, fit, 'o')
plt.xscale('log')
plt.yscale('log')
#plt.axes(xrange())
#yaxes: 1e-29:1e-22
#plt.legend('reaction rate [cm^{-3}s^{-1}]')
plt.show()



## plotting the population densities

def plot_vj_populations(v, j, nvj):
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(v, j, nvj, 'o')
    ax.axes.xaxis.set_ticks(unique(v))
    ax.axes.yaxis.set_ticks(unique(j)[::2])
    ax.set_xlabel('v')
    ax.set_ylabel('j')
    ax.set_zlabel('n(v,j)')
    pylab.show()

plot_vj_populations(vj_unique[0], vj_unique[1], log10(nvj.flatten()))


print 'done reading!'