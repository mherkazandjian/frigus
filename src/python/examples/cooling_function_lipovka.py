# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling function using
the Lipovka data
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u

from frigus.population import (
    fit_lipovka,
    cooling_rate_at_steady_state
)

from frigus.readers import DataLoader

species_data = DataLoader().load('HD_lipovka')

plt.ioff()
fig, axs = plt.subplots()

# density of the colliding species, in m^3
# nc_H = 1e6 * u.meter ** -3
nc_H_rng = [1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14]
T_rad = 0.0 * u.Kelvin

T_rng = np.logspace(2, 3.2, 10) * u.Kelvin
# T_rng = species_data.raw_data.collision_rates_T_range
for nc_H in nc_H_rng:
    lambda_vs_T_kin = []
    pop_dens_vs_T_kin = []
    nc_H = nc_H * u.meter ** -3
    for T_kin in T_rng:
        print(T_kin, nc_H)


        lambda_vs_T_kin += [
            cooling_rate_at_steady_state(
                species_data,
                T_kin,
                T_rad,
                nc_H)]

    lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)

    lambda_vs_T_kin_lipovka = fit_lipovka(T_rng, nc_H)

    axs.loglog(
        T_rng.value, lambda_vs_T_kin.si.value,
        '-x', color = 'black', label='')

    axs.loglog(
        T_rng.value, lambda_vs_T_kin_lipovka.si.value,
        'r--', color = 'black', label='')

    axs.set_xlabel('T$_\mathrm{kin}$ [K]')
    axs.set_ylabel('cooling function [erg s$^{-1}$]')

    plt.legend()

plt.show()

print('done')
