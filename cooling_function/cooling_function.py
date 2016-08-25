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

import pylab

import numpy
from numpy import zeros, log10, unique

from astropy import units as u

import matplotlib.pyplot as plt

# .. todo:: organize these imports
import read_einstien_coefficient
from read_collision_coefficients import read_collision_coefficients
import read_energy_levels
import population
from population import cooling_rate, fit_glover

import pdb

#
# .. todo:: use proper units and standardize them or use astropy units..etc..
# .. todo:: to make sure that there are no errors done in unit conversion
#

# density of the colliding species, in m^3
nc_H = 1e14 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 3000.0 * u.Kelvin

# read the energy levels (v, j, energy)
#
energy_levels = read_energy_levels.read_levels_lique(
                       'Read/H2Xvjlevels_francois_mod.cs')
# en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")
# print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
# for level in levels:
#     print('{:<3}{:<3}{:<10}'.format(level['v'], level['j'], level['E']))

#
# read the einstein coefficients for the H2 transitions
#
A, A_info_nnz = read_einstien_coefficient.read_einstein()

# read the collisional rates for H2 with H
collision_rates, T_rng, collision_rates_info_nnz = read_collision_coefficients(
                                                      "Read/Rates_H_H2.dat")

# find the maximum v and j from the Einstein and collisional rates data sets
# and adjust the labels of the energy levels according to that
v_max_data, j_max_data = population.find_v_max_j_max_from_data(
                                          A_info_nnz,
                                          collision_rates_info_nnz)
energy_levels.set_labels(v_max=v_max_data+1)

#
# reduce the Einstein coefficients to a 2D matrix (construct the A matrix)
# [n_levels, n_levels]
# A_reduced_slow = population.reduce_einstein_coefficients_slow(A_info_nnz,
#                                                               energy_levels)
A_matrix = population.reduce_einstein_coefficients(A, energy_levels)

# getting the collisional de-excitation matrix (K_dex) (for all tabulated
# values)  [n_level, n_level, n_T_kin_values]
K_dex_matrix = population.reduce_collisional_coefficients_slow(
                                          collision_rates_info_nnz,
                                          energy_levels)

# compute the interpolator that produces K_dex at a certain temperature
K_dex_matrix_interpolator = population.compute_K_dex_matrix_interpolator(
                                                         K_dex_matrix, T_rng)

def cooling_rate_at_steady_state_T_kin_nc(T_kin, nc):
    return population.cooling_rate_at_steady_state(A_matrix,
                                                   energy_levels,
                                                   K_dex_matrix_interpolator,
                                                   T_kin,
                                                   nc)

cooling_rate = cooling_rate_at_steady_state_T_kin_nc(T_kin, nc_H)

print('this')
lambda_vs_T_kin = []
for T_kin in T_rng:
    print(T_kin)
    lambda_vs_T_kin += [cooling_rate_at_steady_state_T_kin_nc(T_kin, nc_H)]

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