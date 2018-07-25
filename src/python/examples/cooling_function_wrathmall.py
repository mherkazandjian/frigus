# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling function using the
Lique data.

.. warning:: this results in negative population densities. fix it.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy

from astropy import units as u

from frigus.cooling_function.fits import fit_glover
from frigus.population import (cooling_rate_at_steady_state,
                               population_density_at_steady_state)

from frigus.readers.dataset import DataLoader

species_data = DataLoader().load('H2_wrathmall')

# density of the colliding species, in m^3
nc_h = 1e14 * u.meter ** -3
t_rad = 0.0 * u.Kelvin

lambda_vs_t_kin = []
pop_dens_vs_t_kin = []

t_rng = species_data.raw_data.collision_rates_t_range
for i, t_kin in enumerate(t_rng):

    #print(t_kin)

    pop_dens_vs_t_kin += [
        population_density_at_steady_state(species_data, t_kin, t_rad, nc_h)
    ]

    lambda_vs_t_kin += [
        cooling_rate_at_steady_state(species_data, t_kin, t_rad, nc_h)
    ]

    print(t_kin,
           pop_dens_vs_t_kin[i][1],
           pop_dens_vs_t_kin[i][3],
           pop_dens_vs_t_kin[i][5],
           pop_dens_vs_t_kin[i][7],
           pop_dens_vs_t_kin[i][9],
           pop_dens_vs_t_kin[i][11],
           pop_dens_vs_t_kin[i][13],
           pop_dens_vs_t_kin[i][15],
           pop_dens_vs_t_kin[i][17]
          )

lambda_vs_t_kin = u.Quantity(lambda_vs_t_kin)
lambda_vs_t_kin_glover = u.Quantity(
    [
        fit_glover(T_kin) for T_kin in t_rng.value
    ]) * nc_h

plt.ion()
fig, axs = plt.subplots(2)

axs[0].loglog(
    t_rng.value, lambda_vs_t_kin.cgs.value,
    '-o', label='cooling H2')
axs[0].loglog(
    t_rng.value, lambda_vs_t_kin_glover.cgs.value,
    'r--', label='cooling H2 glover')

pop_dens_vs_t_kin = numpy.array(pop_dens_vs_t_kin)[:, :, 0]
axs[1].plot(
    t_rng.value,
    lambda_vs_t_kin.cgs.value / lambda_vs_t_kin_glover.cgs.value,
    '-', label='lambda / lambda_glover')

axs[1].plot(
    t_rng.value,
    pop_dens_vs_t_kin[:, 0] / pop_dens_vs_t_kin[:, 1],
    '--', label='x_v_0_j_0 / x_v_0_j_1')





plt.legend()
plt.show()

print('done')
