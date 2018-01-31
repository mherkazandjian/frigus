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

from frigus.readers.dataset import DataLoader

species_data = DataLoader().load('H2_lique')

plt.ion()
fig, axs = plt.subplots(figsize=(8, 8))

# density of the colliding species, in m^3
nc_H_rng = [
               1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14
           ] * u.meter ** -3
T_rad = 0.0 * u.Kelvin

T_rng = numpy.logspace(2, 3.2, 10) * u.Kelvin
# T_rng = species_data.raw_data.collision_rates_T_range

plot_markers = [(2 + i//2, 1 + i % 2, 0) for i in range(16)]

for nc_index, nc_H in enumerate(nc_H_rng):

    lambda_vs_T_kin = []
    pop_dens_vs_T_kin = []

    for T_kin in T_rng:
        print(T_kin, nc_H)

        lambda_vs_T_kin += [
            cooling_rate_at_steady_state(
                species_data,
                T_kin,
                T_rad,
                nc_H)
        ]

    lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)

    axs.loglog(
        T_rng.value, lambda_vs_T_kin.to(u.Joule / u.second).value,
        '-x', color='black', marker=plot_markers[nc_index], label=nc_H
    )

axs.set_xlabel('T$_\mathrm{kin}$ [K]')
axs.set_ylabel('cooling function [erg s$^{-1}$]')

lambda_vs_T_kin_glover = u.Quantity(
    [fit_glover(T_kin) for T_kin in T_rng.value]
) * (1.0e6 * u.meter**-3)

axs.loglog(
    T_rng.value, lambda_vs_T_kin_glover.to(u.Joule / u.second).value,
    'r--', color='black', label='cooling H2 glover'
)


plt.legend()
plt.show()

print('done')
