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

x_fit = []
y_fit = []
lambda_fit = []

plt.ion()
fig, axs = plt.subplots(2)

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
    # nc_H = 1e6 * u.meter ** -3
    nc_H_rng = [1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e24]
    T_rad = 0.0 * u.Kelvin

    T_rng = species_data.raw_data.collision_rates_T_range
        for nc_H in nc_H_rng:
            pop_dens_vs_T_kin = []

            d = []
            nc_H = nc_H * u.meter ** -3
            for T_kin in T_rng:
                print(T_kin, nc_H)

            pop_dens_vs_T_kin += [population_density_at_steady_state(species_data,
                                                   T_kin,
                                                   T_rad,
                                                   nc_H)]

            lambda_grid += [cooling_rate_at_steady_state(species_data,
                                                         T_kin,
                                                         T_rad,
                                                         nc_H)]

            x_fit.append(T_kin.value)
            y_fit.append(nc_H.cgs.value)


        lambda_grid = u.Quantity(lambda_grid)
        lambda_fit.append(lambda_grid.value)

        lambda_vs_T_kin_glover = u.Quantity([fit_glover(T_kin) for T_kin in
                                             T_rng.value]) * nc_H


        axs[0].loglog(
            T_rng.value, lambda_grid.si.value,
            '-o', label='cooling H2')
        axs[0].loglog(
            T_rng.value, lambda_vs_T_kin_glover.si.value,
            'r--', label='cooling H2 glover')

        pop_dens_vs_T_kin = numpy.array(pop_dens_vs_T_kin)[:, :, 0]
        axs[1].plot(
            T_rng.value,
            lambda_grid.si.value / lambda_vs_T_kin_glover.si.value,
            '-', label='lambda / lambda_glover')

        axs[1].plot(
            T_rng.value,
            pop_dens_vs_T_kin[:, 0] / pop_dens_vs_T_kin[:, 1],
            '--', label='x_v_0_j_0 / x_v_0_j_1')

        plt.legend()
        plt.show()

print('done')
