# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling function using the
Lique data.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy

from astropy import units as u

from frigus.population import (fit_glover,
                               cooling_rate_at_steady_state,
                               population_density_at_steady_state)

from frigus.readers import DataLoader

species_data = DataLoader().load('H2_lique')

# Calculate the population density and the cooling rate per particle
# for one value of kinetic temperature and gas density of H
if False:

    # density of the colliding species, in m^3
    nc_H = 1e6 * u.meter ** -3
    T_kin = 100.0 * u.Kelvin
    T_rad = 0.0 * u.Kelvin

    pop_dens_equ = population_density_at_steady_state(
        species_data,
        T_kin,
        T_rad,
        nc_H)

    cooling_rate = cooling_rate_at_steady_state(
        species_data,
        T_kin,
        T_rad,
        nc_H)

if True:

    # density of the colliding species, in m^3
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
    lambda_vs_T_kin_glover = u.Quantity([fit_glover(T_kin) for T_kin in
                                         T_rng.value]) * nc_H

    plt.ion()
    fig, axs = plt.subplots(2)

    axs[0].loglog(
        T_rng.value, lambda_vs_T_kin.si.value,
        '-o', label='cooling H2')
    axs[0].loglog(
        T_rng.value, lambda_vs_T_kin_glover.si.value,
        'r--', label='cooling H2 glover')

    pop_dens_vs_T_kin = numpy.array(pop_dens_vs_T_kin)[:, :, 0]
    axs[1].plot(
        T_rng.value,
        lambda_vs_T_kin.si.value / lambda_vs_T_kin_glover.si.value,
        '-', label='lambda / lambda_glover')

    axs[1].plot(
        T_rng.value,
        pop_dens_vs_T_kin[:, 0] / pop_dens_vs_T_kin[:, 1],
        '--', label='x_v_0_j_0 / x_v_0_j_1')

    plt.legend()
    plt.show()

print('done')
