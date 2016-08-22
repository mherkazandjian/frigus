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
import read_einstien_coefficient
from read_collision_coefficients import read_collision_coefficients
import read_levels
import population
from population import cooling_rate, fit_glover
import matplotlib.pyplot as plt

import numpy
from numpy import zeros
import pdb
from numpy import log10, unique

#
# .. todo:: use proper units and standardize them or use astropy units..etc..
# .. todo:: to make sure that there are no errors done in unit conversion
#

# density of the colliding species, in m^3
nc = 1.e14

# the kinetic temperature of the gas
T_kin = 3000.0

# read the energy levels (v, j, energy)
#
energy_levels = read_levels.read_levels_lique(
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
collision_rates, T, collision_rates_info_nnz = read_collision_coefficients(
                                                      "Read/Rates_H_H2.dat")

# find the maximum v and j from the Einstein and collisional rates data sets
# and adjust the labels of the energy levels according to that
v_max_data, j_max_data = population.find_v_max_j_max_from_data(
                                          A_info_nnz,
                                          collision_rates_info_nnz)
energy_levels.set_labels(v_max=v_max_data+1)

#
# reduce the Einstein coefficients to a 2D matrix (construct the A matrix)
#
# A_reduced_slow = population.reduce_einstein_coefficients_slow(A_info_nnz,
#                                                               energy_levels)
A_matrix = population.reduce_einstein_coefficients(A, energy_levels)

# compute the detla E matrix
delta_e_matrix = population.compute_delta_energy_matrix(energy_levels)

# compute the stimulated emission and absorption coefficients matrix
B_J_nu_matrix = population.compute_B_J_nu_matrix_from_A_matrix(energy_levels,
                                                               A_matrix,
                                                               T_kin)


# getting the collisional de-excitation matrix (K_dex) (for all tabulated
# values
K_dex_matrix = population.reduce_collisional_coefficients_slow(
                                                 collision_rates_info_nnz,
                                                 energy_levels)

# get the K matrix for a certain temperature in the tabulated range
K_matrix = population.compute_K_matrix_from_K_dex_matrix(energy_levels,
                                                         K_dex_matrix,
                                                         T,
                                                         T_kin)

# compute the M matrix that can be used to compute the equilibrium state of the
# levels (see notebook)
O_matrix = (A_matrix + B_J_nu_matrix + K_matrix*nc).T

D_matrix = numpy.zeros(O_matrix.shape, 'f8')
D_matrix[numpy.diag_indices(D_matrix.shape[0])] = -O_matrix.sum(axis=0)

M_matrix = O_matrix + D_matrix

# solve the equilibrium population densities
n_vj_equilibrium = population.solveEquilibrium(M_matrix)

# compute the cooling rate
c_rate = cooling_rate(n_vj_equilibrium, energy_levels, A_matrix)

pdb.set_trace()


# creating the array for the radiation temperatures;
# at the moment, same size of tkin but constant value @ 2.73
tcmb = zeros(T.size)
tcmb = 2*T
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

# plot_vj_populations(vj_unique[0], vj_unique[1], log10(nvj.flatten()))


print 'done reading!'