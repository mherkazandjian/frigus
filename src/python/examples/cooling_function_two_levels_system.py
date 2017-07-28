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
nc_H = 1e6 * u.meter ** -3
T_kin = 100.0 * u.Kelvin
T_rad = 0.0 * u.Kelvin

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

nc_H = 1e6 * u.meter ** -3
T_rad = 0.0 * u.Kelvin

lambda_vs_T_kin = []
pop_dens_vs_T_kin = []

T_rng = species_data.raw_data.collision_rates_T_range
for T_kin in T_rng:

    print(T_kin)

    pop_dens_vs_T_kin += [population_density_at_steady_state(species_data,
                                                                 T_kin,
                                                                 T_rad,
                                                                 nc_H)]

    lambda_vs_T_kin += [cooling_rate_at_steady_state(species_data,
                                                         T_kin,
                                                         T_rad,
                                                         nc_H)]


lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)
lambda_vs_T_kin_analytic = \
    u.Quantity([analytic_solution_no_radiation_field(species_data, T_kin, nc_H)\
                for T_kin in T_rng.value])


plt.ion()
fig, axs = plt.subplots(2)

axs[0].loglog(
        T_rng.value, lambda_vs_T_kin.si.value,
        '-o', label='cooling H2 (computed)')
axs[0].loglog(
        T_rng.value, lambda_vs_T_kin_analytic.si.value,
        'r--', label='cooling H2 (analytic)')

pop_dens_vs_T_kin = numpy.array(pop_dens_vs_T_kin)[:, :, 0]
axs[1].plot(
        T_rng.value,
        lambda_vs_T_kin.si.value / lambda_vs_T_kin_analytic.si.value,
        '-', label='lambda / lambda_analytic')

axs[1].plot(
        T_rng.value,
        pop_dens_vs_T_kin[:, 0] / pop_dens_vs_T_kin[:, 1],
        '--', label='x_v_0_j_0 / x_v_0_j_1')

plt.legend()
plt.show()

print('done')
