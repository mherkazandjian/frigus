# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""
from __future__ import print_function
import pylab

from astropy import units as u

from population import (fit_lipovka,
                        cooling_rate_at_steady_state,
                        population_density_at_steady_state)

from readers import DataLoader

species_data = DataLoader().load('HD_lipovka')

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
        nc_H,
        debug=True)

    cooling_rate = cooling_rate_at_steady_state(species_data,
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

    lambda_vs_T_kin_lipovka = fit_lipovka(T_rng, nc_H)

    pylab.ion()
    fig, axs = pylab.subplots()

    axs.loglog(
        T_rng.value, lambda_vs_T_kin.si.value,
        '-o', label='cooling H2')

    axs.loglog(
        T_rng.value, lambda_vs_T_kin_lipovka.si.value,
        'r--', label='cooling H2 lipovka')

    pylab.legend()
    pylab.show()

print('done')