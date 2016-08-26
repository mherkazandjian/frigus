# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""

import pylab
from numpy import exp, fabs, logspace, array

from astropy import units as u
from astropy.constants import k_B

from population import population_density_at_steady_state
from readers import DataLoader

import pdb

# density of the colliding species, in m^3
nc = 1e14 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 3000.0 * u.Kelvin

species_data = DataLoader().load('two_level_1')

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

x_0_exact = []
x_0_numerical = []
x_1_exact = []
x_1_numerical = []
for T_kin in T_kin_range:
    x_exact = analytic_solution_no_radiation_field(species_data, T_kin, nc)
    x_numerical = population_density_at_steady_state(species_data, T_kin, nc)
    x_0_exact.append(x_exact[0])
    x_1_exact.append(x_exact[1])
    x_0_numerical.append(x_numerical[0])
    x_1_numerical.append(x_numerical[1])

x_0_exact = array(x_0_exact).flatten()
x_1_exact = array(x_1_exact).flatten()
x_0_numerical = array(x_0_numerical).flatten()
x_1_numerical = array(x_1_numerical).flatten()

if False:
    pylab.figure()
    pylab.loglog(T_kin_range, x_0_exact, 'r-')
    pylab.loglog(T_kin_range, x_0_numerical, 'ro')
    pylab.loglog(T_kin_range, x_1_exact, 'b-')
    pylab.loglog(T_kin_range, x_1_numerical, 'bo')
    pylab.show()
rel_error_x0 = fabs(1.0 - x_0_numerical / x_0_exact)
rel_error_x1 = fabs(1.0 - x_1_numerical / x_1_exact)

assert rel_error_x0.max() < 1e-13
assert rel_error_x1.max() < 1e-13

if False:
    pylab.figure()
    pylab.loglog(T_kin_range, rel_error_x0, 'r-')
    pylab.loglog(T_kin_range, rel_error_x0, 'k-')
    pylab.show()


print 'done'