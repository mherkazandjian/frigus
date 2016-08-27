# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""

import pylab
import numpy

from astropy import units as u

import matplotlib.pyplot as plt

import population
from population import (fit_glover,
                        cooling_rate_at_steady_state,
                        population_density_at_steady_state)
import utils

import pdb

from readers import DataLoader

# density of the colliding species, in m^3
nc_H = 1e2 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 100.0 * u.Kelvin
T_rad = 30.0 * u.Kelvin

species_data = DataLoader().load('H2_lique')

pop_dens_equ = population_density_at_steady_state(species_data,
                                                  T_kin,
                                                  T_rad,
                                                  nc_H)
cooling_rate = cooling_rate_at_steady_state(species_data,
                                            T_kin,
                                            T_rad,
                                            nc_H)

# utils.load_ascii_matrix_data()

if True:
    print('this')
    lambda_vs_T_kin = []
    T_rng = species_data.raw_data.collision_rates_T_range
    for T_kin in T_rng:
        print(T_kin)
        lambda_vs_T_kin += [cooling_rate_at_steady_state(species_data,
                                                         T_kin, T_rad, nc_H)]

    lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)
    lambda_vs_T_kin_glover = u.Quantity([fit_glover(T_kin) for T_kin in
                                         T_rng.value])

    pylab.loglog(T_rng.value, lambda_vs_T_kin.si.value, '-o', label='cooling H2')
    pylab.loglog(T_rng.value, lambda_vs_T_kin_glover.si.value,
                 'r--', label='cooling H2 glover')
    pylab.legend()
    pylab.show()

print 'done'