# -*- coding: utf-8 -*-
"""
dump the cooling rate for a grid of gas density, kinetic temperature and
radiation temperatures.
"""
from __future__ import print_function
import pylab
import numpy

from astropy import units as u

from frigus.population import fit_glover, cooling_rate_at_steady_state
from frigus.readers.dataset import DataLoader

# density of the colliding species, in m^3
nc_H = 1e2 * u.meter**-3

# the kinetic temperature of the gas
T_kin = 100.0 * u.Kelvin
T_rad = 0.0 * u.Kelvin

species_data = DataLoader().load('H2_lique')

#####
# T_rad   T_kin    lambda    fit_glover
#####
with open('cooling_rate.out', 'w') as fobj:

    # the header of the file
    fobj.write('#{:<12} {:<12} {:<12} {:<12}\n'.format(
        'n(H)', 'T_rad', 'T_kin', 'rate'
    ))

    fobj.write('#{:<12} {:<12} {:<12} {:<12}\n'.format(
        u.m**-3, u.Kelvin, u.Kelvin, u.Joule / u.second
    ))

    fig, axs = pylab.subplots(1)

    # the data
    for nc_H in numpy.logspace(6, 12, 7)*u.meter**-3:
        for T_rad in numpy.array([0.0], 'f8')*u.Kelvin:
            cooling_rate_vs_T_kin = []
            for T_kin in species_data.raw_data.collision_rates_T_range:

                cooling_rate = cooling_rate_at_steady_state(
                    species_data,
                    T_kin,
                    T_rad,
                    nc_H
                )

                fobj.write('{:e} {:e} {:e} {:e}\n'.format(
                    nc_H.value,
                    T_rad.value,
                    T_kin.value,
                    cooling_rate.to(u.Joule / u.second).value)
                )

                cooling_rate_vs_T_kin.append(
                    cooling_rate.to(u.Joule / u.second).value
                )

            pylab.loglog(
                species_data.raw_data.collision_rates_T_range.value,
                cooling_rate_vs_T_kin, '-', label='{}'.format(nc_H)
            )


nc_H = 1e6 * u.meter**-3
lambda_vs_T_kin = []
for T_kin in species_data.raw_data.collision_rates_T_range:
    lambda_vs_T_kin += [
        cooling_rate_at_steady_state(
            species_data,
            T_kin,
            T_rad,
            nc_H
        )
    ]

pylab.loglog(
    species_data.raw_data.collision_rates_T_range.value,
    u.Quantity(lambda_vs_T_kin).cgs.value,
    'ko',
    label='lambda compared with glover'
)

lambda_vs_T_kin_glover = u.Quantity(
    [
        fit_glover(T_kin) for T_kin in
        species_data.raw_data.collision_rates_T_range.value
    ]
)*nc_H

pylab.loglog(
    species_data.raw_data.collision_rates_T_range.value,
    lambda_vs_T_kin_glover.cgs.value,
    'kx',
    label='glover'
)

pylab.legend()
pylab.show()
print('done')
