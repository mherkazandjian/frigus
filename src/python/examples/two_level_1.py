# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""

import pylab
from numpy import exp, fabs, logspace, array, vstack

from astropy import units as u
from astropy.constants import k_B

from frigus.population import population_density_at_steady_state
from frigus.readers import DataLoader

import pdb

# density of the colliding species, in m^3
nc = 1e14 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 3000.0 * u.Kelvin
T_rad = 30.0 * u.Kelvin

species_data = DataLoader().load('two_levels_system')

def analytic_solution_no_radiation_field(species_data, T_kin, n_collider):
    """

    :param T_rad: The radiation temperature
    :param T_kin: The kinetic temperature
    :param n_collider: The density of the colliding species
    :return: fractional population densities of the levels
    """
    A_10 = species_data.A_matrix[1, 0]
    K_10 = species_data.K_dex_matrix_interpolator(T_kin)[1, 0]
    n_c = n_collider
    g_0, g_1 = species_data.energy_levels.data['g']
    E_0, E_1 = species_data.energy_levels.data['E']

    K_01 = (g_1/g_0)*K_10*exp(-fabs(E_1 - E_0)/(k_B*T_kin))

    r = n_c * K_01 / (n_c * K_10 + A_10)
    x_0 = 1.0 / (1.0 + r)
    x_1 = r / (1.0 + r)

    return x_0, x_1

T_kin_range = logspace(2.0, 6.0, 100)*u.K

x_exact = array(
    [analytic_solution_no_radiation_field(species_data, T_kin, nc)
     for T_kin in T_kin_range], 'f8')

x_numerical = array(
    [population_density_at_steady_state(species_data, T_kin, 0.0, nc)
    for T_kin in T_kin_range], 'f8')[:, :, 0]

if True:
    # plot the relative population densities as a function of temperatures
    pylab.figure()
    pylab.loglog(T_kin_range, x_exact[:, 0], 'r-', label='x_0 exact')
    pylab.loglog(T_kin_range, x_numerical[:, 0], 'ro', label='x_0 numerical')
    pylab.loglog(T_kin_range, x_exact[:, 1], 'b-', label='x_1 exact')
    pylab.loglog(T_kin_range, x_numerical[:, 1], 'bo', label='x_1 numerical')
    pylab.legend()
    pylab.title('fractional population densities of levels 0 and level 1')
    pylab.show()

rel_errors = fabs(1.0 - x_numerical / x_exact)

assert rel_errors.max() < 1e-13

if True:
    # plot the relative errors in the relative population densities as a
    # function of temperatures
    pylab.figure()
    pylab.loglog(T_kin_range, rel_errors[:, 0], 'r-', label='x_0')
    pylab.loglog(T_kin_range, rel_errors[:, 1], 'k-', label='x_1')
    pylab.legend()
    pylab.title('relative errors between numerical and exact solution')
    pylab.show()


print 'done'