# -*- coding: utf-8 -*-
"""
Dalculate equilibrium population of species and the cooling function using
the Lipovka data
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u

from frigus.cooling_function.fits import fit_lipovka
from frigus.population import cooling_rate_at_steady_state
from frigus.readers.dataset import DataLoader

species_data = DataLoader().load('HD_lipovka')

plt.ion()
fig, axs = plt.subplots(figsize=(8, 8))


# density of the colliding species, in m^3
nc_h_rng = [
               1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14
           ] * u.meter ** -3
t_rad = 0.0 * u.Kelvin

t_rng = np.logspace(2, 3.2, 10) * u.Kelvin
# T_rng = species_data.raw_data.collision_rates_T_range

plot_markers = [(2 + i//2, 1 + i % 2, 0) for i in range(16)]

for nc_index, nc_h in enumerate(nc_h_rng):

    lambda_vs_t_kin = []
    pop_dens_vs_t_kin = []

    for t_kin in t_rng:
        print(t_kin, nc_h)

        lambda_vs_t_kin += [
            cooling_rate_at_steady_state(
                species_data,
                t_kin,
                t_rad,
                nc_h)
        ]

    lambda_vs_t_kin = u.Quantity(lambda_vs_t_kin)

    lambda_vs_t_kin_lipovka = fit_lipovka(t_rng, nc_h)

    axs.loglog(
        t_rng.value, lambda_vs_t_kin.cgs.value,
        '-x', color='black', marker=plot_markers[nc_index], label=nc_h
    )

    axs.loglog(
        t_rng.value, lambda_vs_t_kin_lipovka.cgs.value,
        'r--', color='black', label=''
    )

    axs.set_xlabel('T$_\mathrm{kin}$ [K]')
    axs.set_ylabel('cooling function [erg s$^{-1}$]')

plt.legend()
plt.show()

print('done')
