# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""

import pylab

from astropy import units as u

import matplotlib.pyplot as plt

import population
from population import (fit_glover,
                        cooling_rate_at_steady_state,
                        population_density_at_steady_state)
import utils

import pdb

from readers import DataLoader

# density of the colliding species, in m^3
nc_H = 1e14 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 3000.0 * u.Kelvin

species_data = DataLoader().load('H2_lique')

pop_dens_equ = population_density_at_steady_state(species_data,T_kin, nc_H)
cooling_rate = cooling_rate_at_steady_state(species_data, T_kin, nc_H)

utils.load_ascii_matrix_data()

if False:
    print('this')
    lambda_vs_T_kin = []
    T_rng = species_data.raw_data.collision_rates_T_range
    for T_kin in T_rng:
        print(T_kin)
        lambda_vs_T_kin += [cooling_rate_at_steady_state(species_data,
                                                         T_kin, nc_H)]

    lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)
    lambda_vs_T_kin_glover = u.Quantity([fit_glover(T_kin) for T_kin in
                                         T_rng.value])

    pylab.loglog(T_rng.value, lambda_vs_T_kin.si.value, '-o', label='cooling H2')
    # pylab.loglog(T_rng.value, lambda_vs_T_kin_glover.si.value,
    #            'r--', label='cooling H2 glover')
    pylab.legend()
    pylab.show()

pdb.set_trace()


# creating the array for the radiation temperatures;
# at the moment, same size of tkin but constant value @ 2.73
# tcmb = zeros(T.size)
# tcmb = 2*T
#tcmb[:] = 2.*T[:]

# cf = zeros(T.size)
# fit = zeros(T.size)
# for itemp in range(T.size):
#     cf[itemp] = cooling_function
#     fit[itemp] = fit_glover(T[itemp])
#
#     print cooling_function*1e13 # to have the output in erg cm-3 s-1


# plt.plot(T, cf*1e13, 'x', T, fit, 'o')
# plt.xscale('log')
# plt.yscale('log')
#plt.axes(xrange())
#yaxes: 1e-29:1e-22
#plt.legend('reaction rate [cm^{-3}s^{-1}]')
# plt.show()

## plotting the population densities

# def plot_vj_populations(v, j, nvj):
#     import pylab
#     from mpl_toolkits.mplot3d import Axes3D
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(v, j, nvj, 'o')
#     ax.axes.xaxis.set_ticks(unique(v))
#     ax.axes.yaxis.set_ticks(unique(j)[::2])
#     ax.set_xlabel('v')
#     ax.set_ylabel('j')
#     ax.set_zlabel('n(v,j)')
#     pylab.show()

# plot_vj_populations(vj_unique[0], vj_unique[1], log10(nvj.flatten()))


print 'done reading!'