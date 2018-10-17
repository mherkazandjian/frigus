# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling function using
the data by Benjamin Desrorreaux et al 2018
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u

from frigus.cooling_function.fits import fit_lipovka, fit_flower
from frigus.population import cooling_rate_at_steady_state
from frigus.readers.dataset import DataLoader

species_data = DataLoader().load('HD_galileo_project')

# density of the colliding species, in m^3
#nc_h_rng = [
#               1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14
#           ] * u.meter ** -3
nc_h_rng = [
               1.e6
           ] * u.meter ** -3
t_rad = 0.0 * u.Kelvin

t_rng = np.logspace(1, 3.1, 30) * u.Kelvin
#T_rng = species_data.raw_data.collision_rates_T_range
#T_rng_old = species_data_old.raw_data.collision_rates_T_range



plot_markers = [(2 + i//2, 1 + i % 2, 0) for i in range(16)]

fig, axs = plt.subplots(figsize=(10, 10))

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

    lambda_vs_t_kin = lambda_vs_t_kin.to(u.erg / u.second)


    lambda_vs_t_kin_lipovka = fit_lipovka(t_rng[t_rng.value >=100], nc_h)

    lambda_vs_t_kin_flower = fit_flower(nc_h, t_rng[t_rng.value >=30])

    #lambda_vs_t_kin_flower = lambda_vs_t_kin_flower.to(u.eV / u.second)

    axs.set_yscale('linear')
    axs.set_xscale('linear')

    axs.loglog(
        t_rng.value, lambda_vs_t_kin.value,
        'o', color='black', label='present calculation'
    )
    axs.loglog(
        t_rng[t_rng.value >= 100.], lambda_vs_t_kin_lipovka.value,
        '-.x', color='red', marker=plot_markers[nc_index], label='fit by Lipovka et al. 2005'
    )
    axs.loglog(
        t_rng[t_rng.value >= 30.], lambda_vs_t_kin_flower.value,
        'r--', color='green', marker=plot_markers[nc_index], label='fit on data by Flower et al. 2000'
    )



#    axs.loglog(
#        t_rng.value, lambda_vs_t_kin_lipovka.cgs.value,
#        'r--', color='black', label=''
#    )

    axs.set_xlabel('T$_\mathrm{kin}$ [K]')
    axs.set_ylabel('cooling function [erg s$^{-1}$]')

plt.legend()
plt.show()

print('done')
