# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling function using the
Lique data.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy

from astropy import units as u

from frigus.cooling_function.fits import fit_glover
from frigus.population import (cooling_rate_at_steady_state,
                               population_density_at_steady_state,
                               compute_transition_rate_matrix, cooling_rate)
from frigus.solvers.linear import solve_equilibrium
from frigus.readers.dataset import DataLoader
from frigus.utils import display_matrix
from frigus.cooling_function.fits import fit_coppola

import matplotlib
from matplotlib.ticker import ScalarFormatter

from multivariate_fit import fit_lambda, func

from decimal import Decimal

species_data_h2_h = DataLoader().load('H2_lique')

# load the data of He + H2 w/wo  radiative data
species_data_h2_he = DataLoader().load('HeH2')
species_data_h2_he_zero_a = DataLoader().load('HeH2')
species_data_h2_he_zero_a.a_matrix *= 0.0
species_data_h2_he_zero_a.raw_data.a *= 0.0
species_data_h2_he_zero_a.raw_data.a_info_nnz = list(
    species_data_h2_he_zero_a.raw_data.a_info_nnz
)
species_data_h2_he_zero_a.raw_data.a_info_nnz[4] *= 0.0

# load the data of H+ + H2 w/wo  radiative data
species_data_h2_hp = DataLoader().load('HpH2')
species_data_h2_hp_zero_a = DataLoader().load('HpH2')
species_data_h2_hp_zero_a.a_matrix *= 0.0
species_data_h2_hp_zero_a.raw_data.a *= 0.0
species_data_h2_hp_zero_a.raw_data.a_info_nnz = list(
    species_data_h2_hp_zero_a.raw_data.a_info_nnz
)
species_data_h2_hp_zero_a.raw_data.a_info_nnz[4] *= 0.0

species_data_wrathmall = DataLoader().load('H2_wrathmall')

plt.ion()
fig, axs = plt.subplots(figsize=(8, 8))

# density of the colliding species, in m^3
nc_h_rng = [
               1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14
           ] * u.meter ** -3

# nc_he_rng = nc_h_rng * 1.e-1

t_rad = 0.0 * u.Kelvin

t_rng = numpy.logspace(2, 3.6, 10) * u.Kelvin
# T_rng = species_data.raw_data.collision_rates_T_range

plot_markers = [(2 + i // 2, 1 + i % 2, 0) for i in range(16)]

x_fit = []
y_fit = []
lambda_fit = []

for nc_index, nc_h in enumerate(nc_h_rng):

    lambda_vs_t_kin = []
    lambda_vs_t_kin_wrathmall = []
    pop_dens_vs_t_kin = []

    for i, t_kin in enumerate(t_rng):
        print(t_kin, nc_h)

        m_h2_h = compute_transition_rate_matrix(species_data_h2_h,
                                                t_kin,
                                                t_rad,
                                                nc_h)

        m_h2_he_without_a = compute_transition_rate_matrix(
            species_data_h2_he_zero_a,
            t_kin,
            t_rad,
            nc_h * 1.e-1)

        m_h2_hp_without_a = compute_transition_rate_matrix(
            species_data_h2_hp_zero_a,
            t_kin,
            t_rad,
            nc_h * 2.e-4)

        m_h2_wrathmall = compute_transition_rate_matrix(species_data_wrathmall,
                                                t_kin,
                                                t_rad,
                                                nc_h)


        m = m_h2_h + m_h2_he_without_a + m_h2_hp_without_a

        m_wrathmall = m_h2_wrathmall

        x_equilibrium = solve_equilibrium(m.si.value)
        x_equilibrium_wrathmall = solve_equilibrium(m_wrathmall.si.value)


        sca = cooling_rate(
            x_equilibrium,
            species_data_h2_h.energy_levels,
            species_data_h2_h.a_matrix
        )

        lambda_vs_t_kin += [cooling_rate(
            x_equilibrium,
            species_data_h2_h.energy_levels,
            species_data_h2_h.a_matrix
            )
        ]
        lambda_vs_t_kin_wrathmall += [cooling_rate(
            x_equilibrium_wrathmall,
            species_data_wrathmall.energy_levels,
            species_data_wrathmall.a_matrix
            )
        ]
        x_fit.append(t_kin.value)
        y_fit.append(nc_h.cgs.value)
        lambda_fit.append(sca.to(u.erg/u.s).value)

    lambda_vs_t_kin = u.Quantity(lambda_vs_t_kin)
    lambda_vs_t_kin_wrathmall = u.Quantity(lambda_vs_t_kin_wrathmall)


    lambda_vs_t_kin_lique = fit_coppola(t_rng, nc_h)

    axs.loglog(
        t_rng.value, lambda_vs_t_kin.to(u.erg / u.second).value,
        '-x', color='black', marker=plot_markers[nc_index], label=''
    )

    axs.loglog(
         t_rng.value, lambda_vs_t_kin_lique.value,
         'r-.', color='green', label=''
    )

    axs.set_xlabel('T$_\mathrm{kin}$ [K]')
    axs.set_ylabel('cooling function [erg $\cdot$ s$^{-1}$]')

    # axs.loglog(
    #     t_rng.value, lambda_vs_t_kin_wrathmall.to(u.erg / u.second).value,
    #      '-.', color='green', marker=plot_markers[nc_index]
    # )

popt, pcov = fit_lambda(x_fit, y_fit, lambda_fit)


lambda_vs_t_kin_glover = u.Quantity(
    [fit_glover(_t_kin) for _t_kin in t_rng.value]
) * (1.0e6 * u.meter ** -3)

axs.loglog(
    t_rng.value, lambda_vs_t_kin_glover.to(u.erg / u.second).value,
    'r--', color='blue'
)

#import pdb; pdb.set_trace()
axs.set_yticks([
                1e-29, 2e-29, 3e-29, 4e-29, 5e-29, 6e-29, 7e-29, 8e-29, 9e-29,
                1e-28, 2e-28, 3e-28, 4e-28, 5e-28, 6e-28, 7e-28, 8e-28, 9e-28,
                1e-27, 2e-27, 3e-27, 4e-27, 5e-27, 6e-27, 7e-27, 8e-27, 9e-27,
                1e-26, 2e-26, 3e-26, 4e-26, 5e-26, 6e-26, 7e-26, 8e-26, 9e-26,
                1e-25, 2e-25, 3e-25, 4e-25, 5e-25, 6e-25, 7e-25, 8e-25, 9e-25,
                1e-24, 2e-24, 3e-24, 4e-24, 5e-24, 6e-24, 7e-24, 8e-24, 9e-24,
                1e-23, 2e-23, 3e-23, 4e-23, 5e-23, 6e-23, 7e-23, 8e-23, 9e-23,
                1e-22, 2e-22, 3e-22, 4e-22, 5e-22, 6e-22, 7e-22, 8e-22, 9e-22,
                1e-21, 2e-21, 3e-21, 4e-21, 5e-21, 6e-21, 7e-21, 8e-21, 9e-21,
                1e-20, 2e-20, 3e-20, 4e-20, 5e-20, 6e-20, 7e-20, 8e-20, 9e-20,
                1e-19, 2e-19, 3e-19, 4e-19, 5e-19, 6e-19, 7e-19, 8e-19, 9e-19,
                1e-18
            ])
# matplotlib.pyplot.ylim(ymin=3e-36, ymax=3e-26)

axs.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())

plt.legend()
plt.show()

print('done')
