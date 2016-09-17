# -*- coding: utf-8 -*-
"""
dump the cooling rate for grids of gas density, kinetic temperature and
radiation temperatures.
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

# utils.load_ascii_matrix_data()

# nc    = 1e6, 1e7, 1e8 -> 1e12
# T_rad = 0.0
# T_kin = francios list


#####
# T_rad   T_kin    lambda    fit_glover
#####
with open('/home/mher/dropbox/Dropbox/us/cooling_function/mher/cooling_rate.out', 'w') as fobj:

    # the header of the file
    fobj.write('#{:<12} {:<12} {:<12} {:<12}\n'.format(
        'n(H)', 'T_rad', 'T_kin', 'rate'
    ))

    fobj.write('#{:<12} {:<12} {:<12} {:<12}\n'.format(
        u.m**-3, u.Kelvin, u.Kelvin, u.Joule / u.second
    ))

    # the data
    for nc_H in numpy.logspace(6, 12, 7)*u.meter**-3:
        for T_rad in numpy.array([0.0], 'f8')*u.Kelvin:
            for T_kin in species_data.raw_data.collision_rates_T_range:
                cooling_rate = cooling_rate_at_steady_state(species_data,
                                             T_kin, T_rad, nc_H)
                fobj.write('{:e} {:e} {:e} {:e}\n'.format(
                    nc_H.value,
                    T_rad.value,
                    T_kin.value,
                    cooling_rate.to(u.Joule / u.second).value))

print 'done'